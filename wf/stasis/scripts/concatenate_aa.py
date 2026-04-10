#!/usr/bin/env python3
"""Concatenate trimmed AA marker alignments into a protein supermatrix."""

from collections import defaultdict
import os

try:
    input_files = list(snakemake.input)
    output_file = snakemake.output[0]
    log_file    = snakemake.log[0]
    samples     = snakemake.params.samples
except NameError:
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--inputs",  nargs="+", required=True)
    p.add_argument("--output",  required=True)
    p.add_argument("--log",     required=True)
    a = p.parse_args()
    input_files = a.inputs
    output_file = a.output
    log_file    = a.log
    samples     = None

marker_data    = {}
marker_lengths = {}
all_taxa       = set()

for f in sorted(input_files):
    marker = os.path.basename(f).replace(".aln.faa", "").replace(".trimmed.faa", "")
    suffix = f"_{marker}"
    seqs   = {}
    cur_id = None
    cur_seq = []
    with open(f) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_id:
                    seqs[cur_id] = "".join(cur_seq)
                raw_id = line[1:].split()[0]
                # Strip _{marker} suffix to recover sample name
                cur_id = raw_id[:-len(suffix)] if raw_id.endswith(suffix) else raw_id
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id:
            seqs[cur_id] = "".join(cur_seq)
    if seqs:
        marker_lengths[marker] = len(next(iter(seqs.values())))
        marker_data[marker]    = seqs
        all_taxa.update(seqs.keys())

sequences = {t: "" for t in all_taxa}
for marker in sorted(marker_data):
    seqs    = marker_data[marker]
    gap_seq = "-" * marker_lengths[marker]
    for taxon in all_taxa:
        sequences[taxon] += seqs.get(taxon, gap_seq)

os.makedirs(os.path.dirname(output_file), exist_ok=True)
with open(output_file, "w") as out:
    for taxon in sorted(sequences):
        out.write(f">{taxon}\n{sequences[taxon]}\n")

total_len = sum(marker_lengths.values())
with open(log_file, "w") as logf:
    logf.write(f"Concatenated {len(marker_data)} markers ({total_len} AA positions) "
               f"for {len(all_taxa)} taxa\n")
