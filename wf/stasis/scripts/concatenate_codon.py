#!/usr/bin/env python3
"""Concatenate all marker codon alignments."""

from collections import defaultdict
from pathlib import Path

# Snakemake variables
input_files = snakemake.input
output_aln = snakemake.output.aln
output_partitions = snakemake.output.partitions
log_file = snakemake.log[0]
samples = snakemake.params.samples

# Collect sequences per sample
sample_seqs = defaultdict(list)
partition_info = []
current_pos = 1

for marker_file in sorted(input_files):
    marker = Path(marker_file).stem.replace(".codon", "")

    marker_seqs = {}
    with open(marker_file) as f:
        name = None
        seq = []
        for line in f:
            if line.startswith(">"):
                if name:
                    marker_seqs[name.rsplit("_", 1)[0]] = "".join(seq)
                name = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if name:
            marker_seqs[name.rsplit("_", 1)[0]] = "".join(seq)

    if not marker_seqs:
        continue

    # Get alignment length
    aln_len = len(next(iter(marker_seqs.values())))

    # Add to concatenation
    for sample in samples:
        if sample in marker_seqs:
            sample_seqs[sample].append(marker_seqs[sample])
        else:
            sample_seqs[sample].append("N" * aln_len)

    # Record partition
    partition_info.append(f"DNA, {marker} = {current_pos}-{current_pos + aln_len - 1}")
    current_pos += aln_len

# Write concatenated alignment
with open(output_aln, "w") as out:
    for sample in sorted(sample_seqs.keys()):
        concat = "".join(sample_seqs[sample])
        out.write(f">{sample}\n{concat}\n")

# Write partition file
with open(output_partitions, "w") as out:
    out.write("\n".join(partition_info))

with open(log_file, "w") as logf:
    logf.write(f"Concatenated {len(partition_info)} markers\n")
    logf.write(f"Total alignment length: {current_pos - 1} bp\n")
    logf.write(f"Samples: {len(sample_seqs)}\n")
