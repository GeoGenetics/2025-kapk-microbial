#!/usr/bin/env python3
"""Extract best-hit marker sequences for each sample from a combined hmmsearch tblout."""

from Bio import SeqIO

# Snakemake variables (try snakemake object, fall back to CLI args)
try:
    faa_file    = snakemake.input.faa
    tbl_file    = snakemake.input.tbl
    output_file = snakemake.output[0]
    markers     = snakemake.params.markers
    sample      = snakemake.wildcards.sample
    log_file    = snakemake.log[0]
except NameError:
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--faa",     required=True)
    p.add_argument("--tbl",     required=True)
    p.add_argument("--output",  required=True)
    p.add_argument("--sample",  required=True)
    p.add_argument("--log",     required=True)
    p.add_argument("--markers", nargs="+", required=True)
    a = p.parse_args()
    faa_file    = a.faa
    tbl_file    = a.tbl
    output_file = a.output
    sample      = a.sample
    log_file    = a.log
    markers     = a.markers

# Load all proteins
proteins = {r.id: r for r in SeqIO.parse(faa_file, "fasta")}

# Parse combined tblout: collect best hit per marker (query name = col 2)
# tblout columns: target_name target_acc query_name query_acc evalue score ...
best_hits = {}  # marker -> (prot_id, score)
with open(tbl_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 6:
            continue
        prot_id = fields[0]
        marker  = fields[2]
        score   = float(fields[5])
        if prot_id not in proteins:
            continue
        if marker not in best_hits or score > best_hits[marker][1]:
            best_hits[marker] = (prot_id, score)

# Build output sequences
extracted = []
for marker in markers:
    if marker not in best_hits:
        continue
    prot_id, _ = best_hits[marker]
    seq = proteins[prot_id]
    seq.id = f"{sample}_{marker}"
    seq.description = prot_id  # keep original gene ID for back-translation
    extracted.append(seq)

SeqIO.write(extracted, output_file, "fasta")

with open(log_file, "w") as logf:
    logf.write(f"Extracted {len(extracted)}/{len(markers)} markers\n")
