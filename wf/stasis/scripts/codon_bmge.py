#!/usr/bin/env python3
"""Codon-aware alignment trimmer for concatenated codon supermatrices.

Removes codon triplets (3-nt atomic units) that exceed a gap fraction
threshold, preserving reading frame. Outputs a trimmed alignment and an
updated partition file with corrected marker boundaries.

Rationale: BMGE on the concatenated supermatrix (not per-marker) gives more
stable entropy estimates and retains more signal. This script mirrors that
approach but operates on codon triplets so IQ-TREE's partition model remains
valid on the output.

Gap definition: a triplet is gapped for a taxon if it is exactly "---"
(back_translate.py produces "---" for AA alignment gaps). Ambiguous "NNN"
padding for absent markers is treated as present (missing data, not a gap).
"""

import re
import os

try:
    input_aln    = snakemake.input.aln
    input_parts  = snakemake.input.partitions
    output_aln   = snakemake.output.aln
    output_parts = snakemake.output.partitions
    log_file     = snakemake.log[0]
    max_gap      = snakemake.params.get("max_gap", 0.2)
    strip_third  = snakemake.params.get("strip_third", False)
except NameError:
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--aln",               required=True)
    p.add_argument("--partitions",        required=True)
    p.add_argument("--output-aln",        required=True, dest="output_aln")
    p.add_argument("--output-partitions", required=True, dest="output_parts")
    p.add_argument("--log",               required=True)
    p.add_argument("--max-gap",           type=float, default=0.2)
    p.add_argument("--strip-third",       action="store_true", default=False)
    a = p.parse_args()
    input_aln    = a.aln
    input_parts  = a.partitions
    output_aln   = a.output_aln
    output_parts = a.output_parts
    log_file     = a.log
    max_gap      = a.max_gap
    strip_third  = a.strip_third

# --- Read alignment (order matters for output) ---
order = []
seqs  = {}
cur_id  = None
cur_seq = []
with open(input_aln) as fh:
    for line in fh:
        line = line.rstrip()
        if line.startswith(">"):
            if cur_id:
                seqs[cur_id] = "".join(cur_seq)
            cur_id = line[1:].split()[0]
            order.append(cur_id)
            cur_seq = []
        else:
            cur_seq.append(line)
    if cur_id:
        seqs[cur_id] = "".join(cur_seq)

aln_len = len(next(iter(seqs.values())))
n_taxa  = len(seqs)
assert aln_len % 3 == 0, f"Alignment length {aln_len} is not divisible by 3"

# --- Read partitions (1-based inclusive: "DNA, NAME = START-END") ---
markers = []  # (name, start0, end0) — 0-based, end exclusive
with open(input_parts) as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        m = re.match(r"DNA,\s*(\S+)\s*=\s*(\d+)-(\d+)", line)
        if m:
            markers.append((m.group(1), int(m.group(2)) - 1, int(m.group(3))))

# --- Filter: keep codon triplets with gap fraction <= max_gap ---
kept_triplets   = []   # list of (aln_start0, aln_end0) for each kept triplet
new_part_bounds = []   # (name, new_start1, new_end1) for output partitions

new_pos = 1  # 1-based running position in trimmed alignment

for marker_name, marker_start0, marker_end0 in markers:
    marker_len = marker_end0 - marker_start0
    assert marker_len % 3 == 0, (
        f"Marker {marker_name} span {marker_len} nt is not divisible by 3"
    )
    marker_kept = 0
    for codon_start in range(marker_start0, marker_end0, 3):
        codon_end = codon_start + 3
        gap_count = sum(
            1 for seq in seqs.values()
            if seq[codon_start:codon_end] == "---"
        )
        if gap_count / n_taxa <= max_gap:
            if strip_third:
                # Keep only 1st and 2nd positions (drop saturated 3rd position)
                kept_triplets.append((codon_start,     codon_start + 1))
                kept_triplets.append((codon_start + 1, codon_start + 2))
                marker_kept += 1  # counts codons (2 nt each in output)
            else:
                kept_triplets.append((codon_start, codon_end))
                marker_kept += 1

    nt_per_codon = 2 if strip_third else 3
    if marker_kept > 0:
        new_part_bounds.append(
            (marker_name, new_pos, new_pos + marker_kept * nt_per_codon - 1)
        )
        new_pos += marker_kept * nt_per_codon

# --- Build and write trimmed alignment ---
os.makedirs(os.path.dirname(output_aln) or ".", exist_ok=True)
with open(output_aln, "w") as out:
    for taxon in order:
        seq = seqs[taxon]
        trimmed = "".join(seq[s:e] for s, e in kept_triplets)
        out.write(f">{taxon}\n{trimmed}\n")

# --- Write updated partition file ---
with open(output_parts, "w") as out:
    out.write("\n".join(
        f"DNA, {name} = {s}-{e}" for name, s, e in new_part_bounds
    ))

n_before   = aln_len // 3
n_codons   = len(kept_triplets) if not strip_third else len(kept_triplets) // 2
nt_out     = sum(e - s for s, e in kept_triplets)
with open(log_file, "w") as logf:
    mode = "1st+2nd positions only" if strip_third else "all positions"
    logf.write(
        f"Codon BMGE ({mode}): {n_before} codons in → {n_codons} codons out "
        f"({100 * n_codons / n_before:.1f}% retained, max_gap={max_gap})\n"
    )
    logf.write(f"Alignment: {n_taxa} taxa × {nt_out} nt\n")
    logf.write(f"Markers: {len(new_part_bounds)} of {len(markers)} retained boundaries\n")
