#!/usr/bin/env python3
"""
Count non-damage SNPs between focal genome pairs from direct pairwise nucmer alignments.

Produces: focal_pairwise_snp_counts.tsv
  pair                   total_snps  damage_snps  nondamage_snps  tv_snps
  ancient_vs_modern      ...         ...          ...             ...
  ancient_vs_s17         ...         ...          ...             ...

Damage-type SNPs (C→T and G→A) are excluded from the reported count because
they are characteristic aDNA deamination artefacts and do not reflect genuine
evolutionary divergence.
"""

from pathlib import Path

TRANSITIONS = {('A','G'), ('G','A'), ('C','T'), ('T','C')}
DAMAGE      = {('C','T'), ('G','A')}


def count_snps(snps_file):
    total = dmg = tv = 0
    with open(snps_file) as f:
        for i, line in enumerate(f):
            if i < 3:          # skip 2-line header + column header
                continue
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            ref = fields[1].upper()
            alt = fields[2].upper()
            if ref not in 'ACGT' or alt not in 'ACGT':
                continue
            total += 1
            if (ref, alt) in DAMAGE:
                dmg += 1
            if (ref, alt) not in TRANSITIONS:
                tv += 1
    return total, dmg, tv


pairs      = snakemake.params.pairs          # list of (name, snps_file) tuples
output_tsv = snakemake.output.tsv

rows = []
for name, snps_file in pairs:
    total, dmg, tv = count_snps(snps_file)
    nondamage = total - dmg
    rows.append((name, total, dmg, nondamage, tv))

with open(output_tsv, 'w') as f:
    f.write("pair\ttotal_snps\tdamage_snps\tnondamage_snps\ttv_snps\n")
    for row in rows:
        f.write("\t".join(str(x) for x in row) + "\n")

# Print summary to log
print("=== Focal pairwise SNP counts (non-damage = reported values) ===")
for name, total, dmg, nondamage, tv in rows:
    pct_dmg = 100 * dmg / total if total else 0
    print(f"  {name}:")
    print(f"    total={total}  damage={dmg} ({pct_dmg:.1f}%)  non-damage={nondamage}  Tv-only={tv}")
