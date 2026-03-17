#!/usr/bin/env python3
"""
Genome-wide pN/pS from MUMmer SNPs + Prodigal annotations.

Strategy: reference-based reconstruction.
All samples are aligned to SOIL100000032 via nucmer. For each pair, reconstruct
gene sequences for both samples by applying their respective SNP sets to the
shared reference, then compare codon-by-codon.

Method: Nei-Gojobori 1986 with Jukes-Cantor correction (standard NG86).
  Ka = -3/4 × ln(1 - 4/3 × pN)
  Ks = -3/4 × ln(1 - 4/3 × pS)

Genetic code 11 (Bacterial/Archaeal).

Usage:
    python genome_pnps.py <snps_dir> <prodigal_dir> <ref_name> <output_tsv>
                          <pair1_a,pair1_b> [pair2_a,pair2_b ...]

Example:
    python genome_pnps.py 10_snp_tree/snps 02_prodigal SOIL100000032 pnps.tsv \
        ancient_consensus,s17_ancient_relative \
        ancient_consensus,modern_consensus
"""

import sys
import csv
import math
import re
from pathlib import Path
from collections import defaultdict

from Bio import SeqIO


# Genetic code 11 (Bacterial/Archaeal/Plant Plastid)
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    # Code 11 start codons: TTG, CTG, ATT, ATC, ATA all encode M as start
    # For internal codons, standard table applies
}
BASES = list('ACGT')


def revcomp(seq: str) -> str:
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq.upper()))


def translate_codon(codon: str) -> str:
    return CODON_TABLE.get(codon.upper(), 'X')


def synonymous_sites_per_codon(codon: str) -> float:
    """
    NG86: for each position in the codon, compute the fraction of synonymous
    single-nucleotide changes. Sum over 3 positions = synonymous sites in codon.
    """
    codon = codon.upper()
    if len(codon) != 3 or any(b not in BASES for b in codon):
        return 0.75  # default for ambiguous codons
    aa_ref = translate_codon(codon)
    if aa_ref == '*':
        return 0.0

    s = 0.0
    for pos in range(3):
        ref_base = codon[pos]
        synonymous = 0
        total = 0
        for alt_base in BASES:
            if alt_base == ref_base:
                continue
            mut_codon = codon[:pos] + alt_base + codon[pos+1:]
            aa_alt = translate_codon(mut_codon)
            if aa_alt == '*':
                continue  # stop codons excluded from site counting
            total += 1
            if aa_alt == aa_ref:
                synonymous += 1
        if total > 0:
            s += synonymous / total
    return s


def count_sites_in_gene(gene_seq: str) -> tuple[float, float]:
    """Return (S_sites, N_sites) for a full CDS sequence (NG86 method)."""
    s_total = 0.0
    n_total = 0.0
    for i in range(0, len(gene_seq) - len(gene_seq) % 3, 3):
        codon = gene_seq[i:i+3]
        if 'N' in codon or '-' in codon:
            continue
        s = synonymous_sites_per_codon(codon)
        s_total += s
        n_total += (3 - s)
    return s_total, n_total


def is_synonymous(codon_ref: str, codon_alt: str) -> bool | None:
    """
    True = synonymous, False = non-synonymous, None = skip (stop/ambiguous).
    """
    codon_ref = codon_ref.upper()
    codon_alt = codon_alt.upper()
    if any(b not in BASES for b in codon_ref + codon_alt):
        return None
    aa_ref = translate_codon(codon_ref)
    aa_alt = translate_codon(codon_alt)
    if aa_ref == '*' or aa_alt == '*':
        return None
    return aa_ref == aa_alt


def jukes_cantor(p: float) -> float:
    """JC correction: d = -3/4 × ln(1 - 4p/3). Returns nan if saturated."""
    arg = 1.0 - 4.0 * p / 3.0
    if arg <= 0:
        return float('nan')
    return -0.75 * math.log(arg)


def load_genome(fasta_path: Path) -> dict:
    """Return {contig_id: sequence}."""
    return {rec.id: str(rec.seq).upper()
            for rec in SeqIO.parse(fasta_path, "fasta")}


def load_gff(gff_path: Path) -> list:
    """
    Parse Prodigal GFF. Returns list of dicts with keys:
    contig, start (1-based), end (1-based, inclusive), strand, gene_id.
    """
    genes = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue
            contig   = fields[0]
            start    = int(fields[3])   # 1-based
            end      = int(fields[4])   # 1-based inclusive
            strand   = fields[6]
            attr     = fields[8]
            gene_id  = re.search(r'ID=([^;]+)', attr)
            gene_id  = gene_id.group(1) if gene_id else f"{contig}_{start}"
            # Skip partial genes at contig edges
            partial  = re.search(r'partial=([01]+)', attr)
            if partial and partial.group(1) not in ('00',):
                continue
            genes.append({
                'contig': contig, 'start': start, 'end': end,
                'strand': strand, 'gene_id': gene_id
            })
    return genes


def load_snps(snps_path: Path) -> dict:
    """
    Parse show-snps -CT output.
    Returns {contig: {pos_1based: (ref_base, qry_base)}}.
    Skips indels (positions where either base is '.').
    """
    snps = defaultdict(dict)
    with open(snps_path) as fh:
        for line in fh:
            if line.startswith('[') or line.startswith('NUCMER') or \
               line.startswith('/') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            try:
                pos    = int(fields[0])
                ref_b  = fields[1].upper()
                qry_b  = fields[2].upper()
                contig = fields[8]     # reference contig (TAGS col 1)
            except (ValueError, IndexError):
                continue
            if ref_b == '.' or qry_b == '.':
                continue   # skip indels
            if ref_b not in BASES or qry_b not in BASES:
                continue
            snps[contig][pos] = (ref_b, qry_b)
    return dict(snps)


def apply_snps_to_gene(ref_seq: str, gene_start: int, snps_on_contig: dict) -> str:
    """
    Apply sample SNPs to reference gene sequence.
    gene_start is 1-based. snps_on_contig is {pos_1based: (ref_b, qry_b)}.
    Returns modified gene sequence (sample's version of the gene).
    """
    seq = list(ref_seq)
    gene_len = len(ref_seq)
    for pos, (ref_b, qry_b) in snps_on_contig.items():
        # Convert global 1-based pos to local 0-based index
        idx = pos - gene_start
        if 0 <= idx < gene_len:
            if seq[idx] == ref_b:   # sanity check
                seq[idx] = qry_b
    return ''.join(seq)


def compare_genes_pairwise(gene_a: str, gene_b: str) -> tuple[int, int, int]:
    """
    Compare two gene sequences codon-by-codon.
    Returns (n_syn, n_nonsyn, n_codons_compared).
    Skips codons with ambiguous bases or stop codons.
    """
    n_syn = n_nonsyn = n_compared = 0
    length = min(len(gene_a), len(gene_b))
    length -= length % 3

    for i in range(0, length, 3):
        c_a = gene_a[i:i+3].upper()
        c_b = gene_b[i:i+3].upper()
        if c_a == c_b:
            n_compared += 1
            continue
        result = is_synonymous(c_a, c_b)
        if result is None:
            continue
        n_compared += 1
        if result:
            n_syn += 1
        else:
            n_nonsyn += 1
    return n_syn, n_nonsyn, n_compared


def compute_pnps_pair(
    sample_a: str, sample_b: str,
    snps_dir: Path, prodigal_dir: Path,
    ref_name: str, ref_genome: dict, ref_genes: list
) -> dict:
    """
    Compute pN/pS between sample_a and sample_b using SOIL-aligned SNPs.
    """
    snps_a_path = snps_dir / f"{sample_a}.snps"
    snps_b_path = snps_dir / f"{sample_b}.snps"

    if not snps_a_path.exists():
        print(f"  Missing: {snps_a_path}", file=sys.stderr)
        return {}
    if not snps_b_path.exists():
        print(f"  Missing: {snps_b_path}", file=sys.stderr)
        return {}

    print(f"  Loading SNPs: {sample_a} ({snps_a_path.stat().st_size//1024}kb), "
          f"{sample_b} ({snps_b_path.stat().st_size//1024}kb)", file=sys.stderr)

    snps_a = load_snps(snps_a_path)
    snps_b = load_snps(snps_b_path)

    s_sites = n_sites = 0.0
    s_subs = n_subs = 0
    n_genes_used = 0

    for gene in ref_genes:
        contig  = gene['contig']
        start   = gene['start']
        end     = gene['end']
        strand  = gene['strand']

        if contig not in ref_genome:
            continue

        # Extract reference gene sequence (0-based slicing, 1-based coords)
        ref_gene_seq = ref_genome[contig][start-1:end]
        if len(ref_gene_seq) < 3:
            continue

        # Apply SNPs from each sample to get their gene sequences
        snps_a_contig = snps_a.get(contig, {})
        snps_b_contig = snps_b.get(contig, {})

        gene_a = apply_snps_to_gene(ref_gene_seq, start, snps_a_contig)
        gene_b = apply_snps_to_gene(ref_gene_seq, start, snps_b_contig)

        # Reverse complement for minus-strand genes
        if strand == '-':
            gene_a = revcomp(gene_a)
            gene_b = revcomp(gene_b)

        # Count synonymous/non-synonymous sites from sample_a sequence
        s, n = count_sites_in_gene(gene_a)
        if s + n < 3:
            continue

        # Compare pairwise
        syn, nonsyn, n_compared = compare_genes_pairwise(gene_a, gene_b)
        if n_compared == 0:
            continue

        s_sites  += s
        n_sites  += n
        s_subs   += syn
        n_subs   += nonsyn
        n_genes_used += 1

    if s_sites == 0 or n_sites == 0:
        return {}

    # Nei-Gojobori proportions
    pS = s_subs / s_sites
    pN = n_subs / n_sites

    # Jukes-Cantor correction
    Ks = jukes_cantor(pS)
    Ka = jukes_cantor(pN)
    omega = Ka / Ks if (not math.isnan(Ka) and not math.isnan(Ks) and Ks > 0) \
            else float('nan')

    return {
        'sample_a':     sample_a,
        'sample_b':     sample_b,
        'Ka':           Ka,
        'Ks':           Ks,
        'omega':        omega,
        'pN':           pN,
        'pS':           pS,
        'N_subs':       n_subs,
        'S_subs':       s_subs,
        'N_sites':      round(n_sites, 1),
        'S_sites':      round(s_sites, 1),
        'n_genes':      n_genes_used,
        'method':       'NG86+JC',
        'ref':          ref_name,
    }


def main():
    if len(sys.argv) < 6:
        print(__doc__, file=sys.stderr)
        sys.exit(1)

    snps_dir    = Path(sys.argv[1])
    prodigal_dir = Path(sys.argv[2])
    ref_name    = sys.argv[3]
    output_tsv  = Path(sys.argv[4])
    pairs_args  = sys.argv[5:]

    # Parse pairs
    pairs = []
    for p in pairs_args:
        parts = p.split(',')
        if len(parts) != 2:
            print(f"Invalid pair: {p} (expected a,b)", file=sys.stderr)
            sys.exit(1)
        pairs.append((parts[0].strip(), parts[1].strip()))

    # Load reference
    ref_fasta = prodigal_dir / f"{ref_name}.fna"    # CDS fna from prodigal? No — need genome
    # Try genome fasta in 01_genomes (standard workflow location)
    ref_genome_path = prodigal_dir.parent / "01_genomes" / f"{ref_name}.fna"
    if not ref_genome_path.exists():
        # Fallback: look for it anywhere nearby
        for candidate in [
            prodigal_dir.parent / f"{ref_name}.fna",
            prodigal_dir / f"{ref_name}.fa",
        ]:
            if candidate.exists():
                ref_genome_path = candidate
                break

    if not ref_genome_path.exists():
        print(f"Reference genome not found: {ref_genome_path}", file=sys.stderr)
        sys.exit(1)

    ref_gff_path = prodigal_dir / f"{ref_name}.gff"
    if not ref_gff_path.exists():
        print(f"Reference GFF not found: {ref_gff_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Reference: {ref_name}", file=sys.stderr)
    print(f"  Genome: {ref_genome_path}", file=sys.stderr)
    print(f"  GFF: {ref_gff_path}", file=sys.stderr)

    ref_genome = load_genome(ref_genome_path)
    ref_genes  = load_gff(ref_gff_path)
    print(f"  Contigs: {len(ref_genome)}, Genes: {len(ref_genes)}", file=sys.stderr)

    # Compute pN/pS for each pair
    results = []
    for sample_a, sample_b in pairs:
        pair_label = f"{sample_a}_vs_{sample_b}"
        print(f"\nProcessing: {pair_label}", file=sys.stderr)
        r = compute_pnps_pair(
            sample_a, sample_b,
            snps_dir, prodigal_dir,
            ref_name, ref_genome, ref_genes
        )
        if r:
            results.append(r)
            print(f"  Ka={r['Ka']:.4e}  Ks={r['Ks']:.4e}  pN/pS={r['omega']:.4f}  "
                  f"(N_subs={r['N_subs']}, S_subs={r['S_subs']}, "
                  f"N_sites={r['N_sites']:.0f}, S_sites={r['S_sites']:.0f}, "
                  f"genes={r['n_genes']})", file=sys.stderr)
        else:
            print(f"  No result", file=sys.stderr)

    # Print summary
    print("\n=== pN/pS SUMMARY (NG86 + JC, genome-wide) ===", file=sys.stderr)
    for r in results:
        print(f"  {r['sample_a']} vs {r['sample_b']}: "
              f"pN/pS={r['omega']:.3f}  Ka={r['Ka']:.3e}  Ks={r['Ks']:.3e}  "
              f"({r['N_subs']}N/{r['S_subs']}S subs, {r['n_genes']} genes)",
              file=sys.stderr)

    # Write TSV
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ['sample_a', 'sample_b', 'Ka', 'Ks', 'omega', 'pN', 'pS',
                  'N_subs', 'S_subs', 'N_sites', 'S_sites', 'n_genes',
                  'method', 'ref']
    with open(output_tsv, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(results)
    print(f"\nResults: {output_tsv}", file=sys.stderr)


if __name__ == '__main__':
    main()
