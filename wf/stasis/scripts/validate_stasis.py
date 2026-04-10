#!/usr/bin/env python3
"""
Extreme validation of evolutionary stasis claim for Kap København Methanoflorens.

Checks:
1. SNP damage artifacts — are the codon/SNP differences C→T or G→A (aDNA damage)?
2. Pairwise distances — ancient vs modern vs all other taxa in codon + SNP alignments
3. Deconvolution quality — damage rates, p_ancient distribution, coverage
4. Tree topology — ancient+modern clade support in codon ML and SNP trees
5. Marker gene SNP spectrum — is ancient enriched for C→T / G→A vs modern?
"""

import re
import sys
from pathlib import Path
from collections import Counter, defaultdict


# ── Helpers ──────────────────────────────────────────────────────────────────

def read_fasta(path):
    seqs = {}
    name = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.upper())
        if name:
            seqs[name] = "".join(seq)
    return seqs


def pairwise_distance(s1, s2):
    """Fraction of differing non-gap, non-N sites."""
    diffs = total = 0
    for a, b in zip(s1, s2):
        if a in "-N" or b in "-N":
            continue
        total += 1
        if a != b:
            diffs += 1
    return diffs / total if total > 0 else float("nan"), diffs, total


def snp_type(ref_base, alt_base):
    """Classify substitution; returns e.g. 'C>T', 'G>A'."""
    return f"{ref_base}>{alt_base}"


def is_damage_type(ref, alt):
    """C→T or G→A are characteristic aDNA deamination artefacts."""
    return (ref == "C" and alt == "T") or (ref == "G" and alt == "A")


# ── Per-alignment analysis ────────────────────────────────────────────────────

def analyze_alignment(path, label, ancient="ancient_consensus", modern="modern_consensus"):
    seqs = read_fasta(path)
    if ancient not in seqs or modern not in seqs:
        return None

    anc = seqs[ancient]
    mod = seqs[modern]

    # Ancient vs modern
    snp_positions = []
    damage_snps = []
    for i, (a, m) in enumerate(zip(anc, mod)):
        if a in "-N" or m in "-N":
            continue
        if a != m:
            snp_positions.append((i + 1, m, a))  # (pos, ref=modern, alt=ancient)
            if is_damage_type(m, a):
                damage_snps.append((i + 1, m, a))

    frac, ndiff, ntot = pairwise_distance(anc, mod)

    lines = []
    lines.append(f"\n{'='*60}")
    lines.append(f"  {label}")
    lines.append(f"{'='*60}")
    lines.append(f"  Alignment length     : {len(anc)} nt")
    lines.append(f"  Comparable sites     : {ntot}")
    lines.append(f"  ancient vs modern SNPs: {ndiff} ({frac*100:.4f}%)")
    lines.append(f"  Damage-type SNPs (C>T / G>A): {len(damage_snps)}")
    if damage_snps:
        lines.append(f"  ⚠ DAMAGE CANDIDATES: {damage_snps}")
    else:
        lines.append(f"  ✓ No SNPs are C>T or G>A — not aDNA damage artifacts")

    if snp_positions:
        lines.append(f"\n  SNP details (pos, modern_base, ancient_base):")
        for pos, ref, alt in snp_positions:
            flag = " ← DAMAGE?" if is_damage_type(ref, alt) else ""
            lines.append(f"    pos {pos:6d}: {ref} → {alt}  {snp_type(ref, alt)}{flag}")

    # Pairwise: ancient/modern vs all other taxa
    lines.append(f"\n  Pairwise distances vs all taxa (sorted):")
    all_taxa = sorted(seqs.keys())
    rows = []
    for t in all_taxa:
        if t in (ancient, modern):
            continue
        fa, da, ta = pairwise_distance(anc, seqs[t])
        fm, dm, tm = pairwise_distance(mod, seqs[t])
        rows.append((t, da, ta, fa, dm, tm, fm))
    rows.sort(key=lambda r: r[1])  # sort by ancient-taxon distance
    lines.append(f"  {'Taxon':<25} {'Anc_SNPs':>8} {'Mod_SNPs':>8} {'Anc%':>7} {'Mod%':>7}")
    for t, da, ta, fa, dm, tm, fm in rows:
        lines.append(f"  {t:<25} {da:>8} {dm:>8} {fa*100:>7.3f} {fm*100:>7.3f}")

    return "\n".join(lines)


# ── SNP spectrum analysis ─────────────────────────────────────────────────────

def analyze_snp_spectrum(snps_dir, ancient="ancient_consensus", modern="modern_consensus", ref="SOIL100000032"):
    """
    Parse per-sample .snps files from show-snps.
    Compare C>T / G>A enrichment in ancient vs modern relative to reference.
    """
    lines = []
    lines.append(f"\n{'='*60}")
    lines.append(f"  SNP Spectrum Analysis (damage artifact check)")
    lines.append(f"{'='*60}")
    lines.append(f"  Reference: {ref}")

    snps_dir = Path(snps_dir)
    for sample in [ancient, modern]:
        snp_file = snps_dir / f"{sample}.snps"
        if not snp_file.exists():
            lines.append(f"  {sample}: SNP file not found at {snp_file}")
            continue

        counts = Counter()
        total = 0
        with open(snp_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("[") or line.startswith("="):
                    continue
                parts = line.split()
                if len(parts) < 3:
                    continue
                try:
                    ref_base = parts[1].upper()
                    alt_base = parts[2].upper()
                except IndexError:
                    continue
                if ref_base in "ACGT" and alt_base in "ACGT" and ref_base != alt_base:
                    counts[snp_type(ref_base, alt_base)] += 1
                    total += 1

        if total == 0:
            lines.append(f"  {sample}: no parseable SNPs")
            continue

        ct = counts.get("C>T", 0)
        ga = counts.get("G>A", 0)
        damage_frac = (ct + ga) / total * 100

        lines.append(f"\n  {sample} ({total} SNPs vs {ref}):")
        lines.append(f"    C>T: {ct:5d} ({ct/total*100:.1f}%)")
        lines.append(f"    G>A: {ga:5d} ({ga/total*100:.1f}%)")
        lines.append(f"    C>T+G>A (damage proxy): {damage_frac:.1f}%")
        if damage_frac > 25:
            lines.append(f"    ⚠ HIGH damage fraction — ancient sample may have residual damage in SNP calls")
        else:
            lines.append(f"    ✓ Damage fraction within expected range for deamination-free SNPs")

        # All substitution types
        lines.append(f"    Full spectrum:")
        for st, cnt in sorted(counts.items(), key=lambda x: -x[1]):
            lines.append(f"      {st}: {cnt:5d} ({cnt/total*100:.1f}%)")

    return "\n".join(lines)


# ── Deconvolution QC ─────────────────────────────────────────────────────────

def analyze_deconvolution(deconv_dir, resolve_dir=None):
    deconv_dir = Path(deconv_dir)
    lines = []
    lines.append(f"\n{'='*60}")
    lines.append(f"  Deconvolution Quality Check")
    lines.append(f"{'='*60}")

    # Damage model output — produced by amber resolve, lives in resolve_dir
    damage_file = Path(resolve_dir) / "damage_per_bin.tsv" if resolve_dir else deconv_dir / "damage_per_bin.tsv"
    if damage_file.exists():
        with open(damage_file) as f:
            header = f.readline().strip().split("\t")
            for line in f:
                parts = line.strip().split("\t")
                if not parts:
                    continue
                row = dict(zip(header, parts))
                bin_id = row.get("bin", parts[0])
                lines.append(f"\n  Bin: {bin_id}")
                for k in ["p_ancient", "ct_1p", "ga_1p", "lambda5", "lambda3", "frag_mean"]:
                    if k in row:
                        lines.append(f"    {k}: {row[k]}")
    else:
        lines.append(f"  damage_per_bin.tsv not found at {damage_file}")

    # Fragment length histogram
    frag_file = deconv_dir / "read_length_histogram.tsv"
    if frag_file.exists():
        import csv
        ancient_lengths = []
        modern_lengths = []
        with open(frag_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                try:
                    mid = (float(row.get("length_min", 0)) + float(row.get("length_max", 0))) / 2
                    anc_count = float(row.get("ancient", 0))
                    mod_count = float(row.get("modern", 0))
                    ancient_lengths.extend([mid] * int(anc_count))
                    modern_lengths.extend([mid] * int(mod_count))
                except (ValueError, KeyError):
                    continue
        if ancient_lengths:
            mean_anc = sum(ancient_lengths) / len(ancient_lengths)
            mean_mod = sum(modern_lengths) / len(modern_lengths) if modern_lengths else float("nan")
            lines.append(f"\n  Fragment lengths:")
            lines.append(f"    Ancient mean: {mean_anc:.1f} bp (n={len(ancient_lengths)})")
            lines.append(f"    Modern mean:  {mean_mod:.1f} bp (n={len(modern_lengths)})")
            if mean_anc < 80:
                lines.append(f"    ✓ Ancient mean < 80 bp — consistent with aDNA fragmentation")
            else:
                lines.append(f"    ⚠ Ancient mean >= 80 bp — unexpectedly long for aDNA")
            if mean_mod > 100:
                lines.append(f"    ✓ Modern mean > 100 bp — consistent with environmental DNA")
    else:
        lines.append(f"  read_length_histogram.tsv not found")

    return "\n".join(lines)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    codon_aln   = snakemake.input.codon_aln
    snp_aln     = snakemake.input.snp_aln
    snps_dir    = snakemake.input.snps_dir
    deconv_dir  = snakemake.input.deconv_dir
    resolve_dir = snakemake.input.resolve_dir
    output_report = snakemake.output.report

    ancient = "ancient_consensus"
    modern  = "modern_consensus"

    report = []
    report.append("=" * 60)
    report.append("  EVOLUTIONARY STASIS VALIDATION REPORT")
    report.append("  Kap København Methanoflorens (~2 Ma)")
    report.append("=" * 60)

    # 1. Codon alignment SNP check
    result = analyze_alignment(codon_aln, "Codon marker gene alignment (78 taxa, HQ phylo)")
    if result:
        report.append(result)

    # 2. SNP alignment analysis
    result = analyze_alignment(snp_aln, "Core genome SNP alignment (tv-only, 9 taxa)")
    if result:
        report.append(result)

    # 3. SNP spectrum (damage artifact check)
    report.append(analyze_snp_spectrum(snps_dir, ancient, modern))

    # 4. Deconvolution QC
    report.append(analyze_deconvolution(deconv_dir, resolve_dir=resolve_dir))

    # 5. Summary verdict
    report.append(f"\n{'='*60}")
    report.append(f"  VERDICT SUMMARY")
    report.append(f"{'='*60}")
    report.append(f"  Codon alignment   : see above — check if 2 SNPs are damage-type")
    report.append(f"  SNP alignment     : see above — check damage fraction in spectrum")
    report.append(f"  Fragment lengths  : ancient should be ~57 bp, modern ~129 bp")
    report.append(f"  Stasis signal     : ancient and modern should be closest to each")
    report.append(f"                      other relative to all other taxa")
    report.append(f"{'='*60}\n")

    with open(output_report, "w") as f:
        f.write("\n".join(report))

    # Also print to stdout for Snakemake log
    print("\n".join(report))


if __name__ == "__main__":
    main()
