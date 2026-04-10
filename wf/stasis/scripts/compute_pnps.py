#!/usr/bin/env python3
"""
Compute pN/pS (Ka/Ks) from per-gene codon alignments using KaKs_Calculator 2.0.

Method: YN (Yang & Nielsen 2000, Mol. Biol. Evol. 17:32-43)
  - Maximum likelihood estimation of Ka and Ks
  - Accounts for transition/transversion rate bias
  - Accounts for unequal base composition
  - Standard for pairwise dN/dS in microbial evolution

Genetic code: 11 (Bacterial, Archaeal, and Plant Plastid Code)

Input:  40 per-gene codon FASTA files from 06_trimmed/{marker}.codon.fna
Output: aggregate_tsv  — per-pair aggregate Ka/Ks across all genes
        per_gene_tsv   — per-gene Ka/Ks for each pair

Key pairs:
  ancient_consensus vs s17_ancient_relative  (primary stasis comparison)
  ancient_consensus vs modern_consensus      (intra-Kap K)
  ancient_consensus vs NAY3300025461b7       (closest permafrost relative)
  modern_consensus  vs s17_ancient_relative
"""

import sys
import csv
import subprocess
import tempfile
import shutil
import re
from pathlib import Path

import numpy as np
from Bio import SeqIO


KAKS_BIN = "KaKs_Calculator"

PAIRS = [
    ("ancient_consensus",   "s17_ancient_relative"),
    ("ancient_consensus",   "modern_consensus"),
    ("ancient_consensus",   "NAY3300025461b7"),
    ("modern_consensus",    "s17_ancient_relative"),
]


def strip_marker_suffix(name: str) -> str:
    for sep in ["_COG", "_DNGNGWU", "_FAM"]:
        idx = name.find(sep)
        if idx >= 0:
            return name[:idx]
    return name


def load_alignment(fasta_path: Path) -> dict:
    seqs = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        taxon = strip_marker_suffix(rec.id)
        seqs[taxon] = str(rec.seq).upper()
    return seqs


def clean_codon_pair(seq1: str, seq2: str) -> tuple[str, str]:
    """Remove codon positions with gaps or ambiguous bases from either sequence."""
    assert len(seq1) == len(seq2)
    trim = len(seq1) - (len(seq1) % 3)
    seq1, seq2 = seq1[:trim], seq2[:trim]

    keep1, keep2 = [], []
    for i in range(0, len(seq1), 3):
        c1 = seq1[i:i+3]
        c2 = seq2[i:i+3]
        if '-' in c1 or '-' in c2:
            continue
        if 'N' in c1 or 'N' in c2:
            continue
        keep1.append(c1)
        keep2.append(c2)
    return "".join(keep1), "".join(keep2)


def write_axt(pairs_seqs: list[tuple[str, str, str]], axt_path: Path):
    """Write AXT format: >name\\nseq1\\nseq2 for each pair."""
    with open(axt_path, "w") as fh:
        for name, s1, s2 in pairs_seqs:
            fh.write(f">{name}\n{s1}\n{s2}\n")


def run_kaks(axt_path: Path, out_path: Path, method: str = "YN") -> bool:
    """Run KaKs_Calculator; return True on success."""
    cmd = [KAKS_BIN, "-i", str(axt_path), "-o", str(out_path),
           "-m", method, "-c", "11"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0


def parse_kaks_output(out_path: Path) -> dict:
    """
    Parse KaKs_Calculator output TSV.
    Returns {pair_name: {Ka, Ks, Ka_Ks, n_sites}} or {} on failure.
    """
    results = {}
    if not out_path.exists():
        return results
    with open(out_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            name = row.get("Sequence", "").strip().lstrip(">")
            if not name:
                continue
            try:
                ka   = float(row["Ka"])
                ks   = float(row["Ks"])
                # Ka/Ks may be "NA" or numeric
                ratio_raw = row.get("Ka/Ks", "NA").strip()
                ratio = float(ratio_raw) if ratio_raw not in ("NA", "") else float('nan')
                # N-Sites is the number of non-synonymous sites (proxy for alignment length)
                n_sites_raw = row.get("N-Sites", "0").strip()
                n_sites = float(n_sites_raw) if n_sites_raw else 0.0
            except (ValueError, KeyError):
                continue
            results[name] = {
                "Ka": ka, "Ks": ks, "Ka_Ks": ratio,
                "n_sites": n_sites
            }
    return results


def collect_concat_seqs(gene_files: list[Path]) -> dict:
    """
    Collect and concatenate cleaned sequences for all pairs across all genes.
    Returns {pair_label: {"s1": str, "s2": str, "n_genes": int, "n_codons": int}}.
    """
    concat = {f"{ta}_vs_{tb}": {"s1": [], "s2": [], "n_genes": 0, "n_codons": 0}
              for ta, tb in PAIRS}

    for gf in gene_files:
        seqs = load_alignment(gf)
        for taxon_a, taxon_b in PAIRS:
            label = f"{taxon_a}_vs_{taxon_b}"
            if taxon_a not in seqs or taxon_b not in seqs:
                continue
            s1, s2 = clean_codon_pair(seqs[taxon_a], seqs[taxon_b])
            if len(s1) < 3:
                continue
            concat[label]["s1"].append(s1)
            concat[label]["s2"].append(s2)
            concat[label]["n_genes"] += 1
            concat[label]["n_codons"] += len(s1) // 3

    return {
        label: {
            "s1": "".join(v["s1"]),
            "s2": "".join(v["s2"]),
            "n_genes": v["n_genes"],
            "n_codons": v["n_codons"],
        }
        for label, v in concat.items()
        if v["n_codons"] > 0
    }


def run_concat_kaks(concat_seqs: dict, tmpdir: Path) -> dict:
    """
    Run KaKs_Calculator on concatenated super-sequences, one AXT file per pair.
    KaKs_Calculator 2.0 has a buffer overflow when multiple long sequences share
    one AXT file — running one pair per file avoids this.
    Returns {pair_label: {Ka, Ks, Ka_Ks, n_genes, n_codons}}.
    """
    results = {}
    for label, v in concat_seqs.items():
        if len(v["s1"]) < 9:
            continue
        axt_path = tmpdir / f"{label}.axt"
        out_path = tmpdir / f"{label}.kaks"
        write_axt([(label, v["s1"], v["s2"])], axt_path)
        run_kaks(axt_path, out_path)
        parsed = parse_kaks_output(out_path)
        r = parsed.get(label, {})
        results[label] = {
            "Ka":       r.get("Ka", float("nan")),
            "Ks":       r.get("Ks", float("nan")),
            "Ka_Ks":    r.get("Ka_Ks", float("nan")),
            "n_genes":  v["n_genes"],
            "n_codons": v["n_codons"],
        }
    return results


def process_gene_level(gene_files: list[Path], tmpdir: Path) -> list[dict]:
    """
    Per-gene KaKs_Calculator for genes that have ≥1 substitution.
    Returns list of {pair_label: result} dicts (one per gene).
    """
    gene_results = []
    for gf in gene_files:
        gene_name = re.sub(r"\.codon$", "", gf.stem)
        seqs = load_alignment(gf)

        axt_entries = []
        for taxon_a, taxon_b in PAIRS:
            if taxon_a not in seqs or taxon_b not in seqs:
                continue
            s1, s2 = clean_codon_pair(seqs[taxon_a], seqs[taxon_b])
            if len(s1) < 9:
                continue
            # Only include if there's at least 1 substitution
            if s1 == s2:
                continue
            label = f"{taxon_a}_vs_{taxon_b}"
            axt_entries.append((label, s1, s2))

        if not axt_entries:
            gene_results.append({})
            continue

        axt_path = tmpdir / f"{gene_name}.axt"
        out_path = tmpdir / f"{gene_name}.kaks"
        write_axt(axt_entries, axt_path)
        run_kaks(axt_path, out_path)

        parsed = parse_kaks_output(out_path)
        gr = {}
        for label, s1, _ in axt_entries:
            if label in parsed:
                parsed[label]["n_codons"] = len(s1) // 3
                gr[label] = parsed[label]
                gr[label]["gene"] = gene_name
        gene_results.append(gr)

    return gene_results


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <codon_dir> <aggregate_tsv> [per_gene_tsv]",
              file=sys.stderr)
        sys.exit(1)

    codon_dir    = Path(sys.argv[1])
    output_tsv   = Path(sys.argv[2])
    per_gene_tsv = Path(sys.argv[3]) if len(sys.argv) > 3 else \
                   output_tsv.parent / (output_tsv.stem + "_per_gene.tsv")

    gene_files = sorted(codon_dir.glob("*.codon.fna"))
    if not gene_files:
        print(f"No *.codon.fna files in {codon_dir}", file=sys.stderr)
        sys.exit(1)
    print(f"Found {len(gene_files)} gene alignments", file=sys.stderr)

    tmpdir = Path(tempfile.mkdtemp(prefix="kaks_"))
    try:
        # Primary: concatenated super-sequence per pair (handles low divergence)
        print("Concatenating sequences across genes...", file=sys.stderr)
        concat_seqs = collect_concat_seqs(gene_files)
        concat_results = run_concat_kaks(concat_seqs, tmpdir)

        # Secondary: per-gene for genes with ≥1 substitution (supplementary)
        print("Running per-gene KaKs (genes with substitutions only)...", file=sys.stderr)
        gene_results = process_gene_level(gene_files, tmpdir)

        # Write per-gene TSV
        per_gene_rows = []
        for gr in gene_results:
            for pair_label, r in gr.items():
                ta, tb = pair_label.split("_vs_", 1)
                per_gene_rows.append({
                    "gene":     r.get("gene", ""),
                    "pair":     pair_label,
                    "taxon_a":  ta,
                    "taxon_b":  tb,
                    "Ka":       r["Ka"],
                    "Ks":       r["Ks"],
                    "omega":    r["Ka_Ks"],
                    "n_codons": r.get("n_codons", 0),
                })

        per_gene_tsv.parent.mkdir(parents=True, exist_ok=True)
        with open(per_gene_tsv, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=["gene", "pair", "taxon_a", "taxon_b",
                            "Ka", "Ks", "omega", "n_codons"],
                delimiter="\t"
            )
            writer.writeheader()
            writer.writerows(per_gene_rows)

        # Write aggregate (concatenated) results
        print("\n=== pN/pS results (YN method, code 11, KaKs_Calculator 2.0) ===",
              file=sys.stderr)
        print("=== Concatenated super-sequence across all 40 marker genes ===\n",
              file=sys.stderr)
        agg_rows = []
        for taxon_a, taxon_b in PAIRS:
            pair_label = f"{taxon_a}_vs_{taxon_b}"
            r = concat_results.get(pair_label, {})
            ka = r.get("Ka", float("nan"))
            ks = r.get("Ks", float("nan"))
            om = r.get("Ka_Ks", float("nan"))
            ka_str = f"{ka:.6e}" if not np.isnan(ka) else "NA"
            ks_str = f"{ks:.6e}" if not np.isnan(ks) else "NA"
            om_str = f"{om:.4f}" if not np.isnan(om) else "NA"
            print(f"  {pair_label}:", file=sys.stderr)
            print(f"    Ka={ka_str}  Ks={ks_str}  pN/pS={om_str}  "
                  f"({r.get('n_genes',0)} genes, {r.get('n_codons',0)} codons)",
                  file=sys.stderr)
            agg_rows.append({
                "pair":         pair_label,
                "taxon_a":      taxon_a,
                "taxon_b":      taxon_b,
                "Ka":           ka_str,
                "Ks":           ks_str,
                "omega":        om_str,
                "method":       "YN",
                "genetic_code": 11,
                "n_genes":      r.get("n_genes", 0),
                "total_codons": r.get("n_codons", 0),
            })

        output_tsv.parent.mkdir(parents=True, exist_ok=True)
        with open(output_tsv, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=["pair", "taxon_a", "taxon_b", "Ka", "Ks", "omega",
                            "method", "genetic_code", "n_genes", "total_codons"],
                delimiter="\t"
            )
            writer.writeheader()
            writer.writerows(agg_rows)

        print(f"\nAggregate: {output_tsv}", file=sys.stderr)
        print(f"Per-gene:  {per_gene_tsv}", file=sys.stderr)

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()
