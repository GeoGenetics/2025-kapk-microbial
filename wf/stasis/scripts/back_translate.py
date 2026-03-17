#!/usr/bin/env python3
"""Back-translate AA alignment to codon alignment."""

from Bio import SeqIO
from pathlib import Path

# Snakemake variables (try snakemake object, fall back to CLI args)
try:
    aa_aln_file   = snakemake.input.aa_aln
    nt_seqs_files = snakemake.input.nt_seqs
    output_file   = snakemake.output[0]
    log_file      = snakemake.log[0]
except NameError:
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--aa-aln",  required=True)
    p.add_argument("--output",  required=True)
    p.add_argument("--log",     required=True)
    p.add_argument("--nt-seqs", nargs="+", required=True)
    a = p.parse_args()
    aa_aln_file   = a.aa_aln
    output_file   = a.output
    log_file      = a.log
    nt_seqs_files = a.nt_seqs

# Build sample -> nt_file mapping
sample_to_nt_file = {}
for f in nt_seqs_files:
    sample = Path(f).stem  # e.g., "ancient_consensus" from "ancient_consensus.fna"
    sample_to_nt_file[sample] = f

# Back-translate
count = 0
missing = []

with open(output_file, "w") as out:
    for record in SeqIO.parse(aa_aln_file, "fasta"):
        # Extract sample and marker from ID (format: sample_marker)
        # Handle markers with underscores (e.g., COG0085-FAM000453_finalBA)
        parts = record.id.split("_")
        if len(parts) < 2:
            continue

        # Try to find which part is the sample by checking sample_to_nt_file
        sample = None
        marker = None
        for i in range(1, len(parts)):
            candidate_sample = "_".join(parts[:i])
            if candidate_sample in sample_to_nt_file:
                sample = candidate_sample
                marker = "_".join(parts[i:])
                break

        if not sample:
            missing.append(record.id)
            continue

        # Get original gene ID from description (set by extract_markers)
        original_gene_id = record.description.split()[1] if len(record.description.split()) > 1 else None

        # Load nucleotide sequences for this sample
        nt_file = sample_to_nt_file[sample]
        nt_seqs = {r.id.split()[0]: str(r.seq) for r in SeqIO.parse(nt_file, "fasta")}

        # Find the matching nucleotide sequence
        aa_ungapped = str(record.seq).replace("-", "").replace("*", "")
        aa_len = len(aa_ungapped)
        expected_nt_len = aa_len * 3

        nt_seq = None

        # First try exact gene ID match
        if original_gene_id and original_gene_id in nt_seqs:
            nt_seq = nt_seqs[original_gene_id]
        else:
            # Fallback: find by length match
            for gene_id, seq in nt_seqs.items():
                # Allow for stop codon difference
                if abs(len(seq) - expected_nt_len) <= 6:
                    nt_seq = seq
                    break

        if not nt_seq:
            missing.append(record.id)
            continue

        # Back-translate with gaps
        codon_aln = []
        nt_idx = 0
        for aa in str(record.seq):
            if aa == "-":
                codon_aln.append("---")
            elif aa == "*":
                # Stop codon
                if nt_idx + 3 <= len(nt_seq):
                    codon_aln.append(nt_seq[nt_idx:nt_idx+3])
                    nt_idx += 3
                else:
                    codon_aln.append("NNN")
            else:
                if nt_idx + 3 <= len(nt_seq):
                    codon_aln.append(nt_seq[nt_idx:nt_idx+3])
                    nt_idx += 3
                else:
                    codon_aln.append("NNN")

        out.write(f">{record.id}\n{''.join(codon_aln)}\n")
        count += 1

with open(log_file, "w") as logf:
    logf.write(f"Back-translated {count} sequences\n")
    if missing:
        logf.write(f"Missing: {len(missing)}\n")
        for m in missing[:10]:
            logf.write(f"  {m}\n")
