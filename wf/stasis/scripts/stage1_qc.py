#!/usr/bin/env python3
"""Generate Stage 1 QC summary."""

from Bio import SeqIO
import re

# Snakemake variables
tree_file = snakemake.input.tree
aln_file = snakemake.input.aln
output_file = snakemake.output[0]

# Count taxa and alignment stats
seqs = list(SeqIO.parse(aln_file, "fasta"))
aln_len = len(seqs[0].seq) if seqs else 0

# Check for ancient/modern
has_ancient = any("ancient" in s.id for s in seqs)
has_modern = any("modern" in s.id for s in seqs)

# Read tree and extract branch lengths for consensus
with open(tree_file) as f:
    tree_str = f.read()

# Simple regex to find branch lengths
anc_match = re.search(r'ancient_consensus:([0-9.e-]+)', tree_str)
mod_match = re.search(r'modern_consensus:([0-9.e-]+)', tree_str)

anc_branch = float(anc_match.group(1)) if anc_match else None
mod_branch = float(mod_match.group(1)) if mod_match else None

with open(output_file, "w") as out:
    out.write("=== Stage 1 QC Summary ===\n\n")
    out.write(f"Taxa: {len(seqs)}\n")
    out.write(f"Alignment length: {aln_len} bp\n")
    out.write(f"Ancient consensus present: {has_ancient}\n")
    out.write(f"Modern consensus present: {has_modern}\n")
    out.write(f"\nBranch lengths (IQ-TREE):\n")
    out.write(f"  ancient_consensus: {anc_branch if anc_branch else 'N/A'}\n")
    out.write(f"  modern_consensus: {mod_branch if mod_branch else 'N/A'}\n")
    if anc_branch and mod_branch and anc_branch > 0:
        out.write(f"  Ratio (modern/ancient): {mod_branch/anc_branch:.2f}x\n")
