#!/usr/bin/env python3
"""
Select the ancient bin for deconvolution.

Criteria (from config):
  - GTDB-Tk classification contains gtdbtk_clade (e.g. "Bog-38")
  - Completeness >= min_completeness
  - Contamination <= max_contamination
  - GUNC pass (if require_gunc_pass: true)

Tiebreak: highest completeness (or highest p_ancient if tiebreak: "ancient")
"""

import csv
import sys

output_file      = snakemake.output.selected_bin
gtdbtk_clade     = snakemake.params.gtdbtk_clade
min_completeness = float(snakemake.params.min_completeness)
max_contamination = float(snakemake.params.max_contamination)
require_gunc_pass = bool(snakemake.params.require_gunc_pass)
tiebreak         = snakemake.params.tiebreak

# --- Load CheckM2 ---
checkm = {}
with open(snakemake.input.checkm) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        checkm[row['Name']] = {
            'completeness': float(row['Completeness']),
            'contamination': float(row['Contamination']),
        }

# --- Load GUNC ---
gunc = {}
with open(snakemake.input.gunc) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        gunc[row['genome']] = row['pass.GUNC'].strip().lower() in ('true', 'pass', '1')

# --- Load GTDB-Tk (ar53 + bac120) ---
gtdbtk = {}
for path in [snakemake.input.gtdbtk_ar53, snakemake.input.gtdbtk_bac120]:
    try:
        with open(path) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                name = row.get('user_genome', '')
                if name:
                    gtdbtk[name] = row.get('classification', '')
    except Exception:
        pass

# --- Apply filters ---
candidates = []
for bin_name, qc in checkm.items():
    if qc['completeness'] < min_completeness:
        continue
    if qc['contamination'] > max_contamination:
        continue
    if require_gunc_pass and not gunc.get(bin_name, False):
        continue
    classification = gtdbtk.get(bin_name, '')
    if gtdbtk_clade not in classification:
        continue
    candidates.append({
        'bin': bin_name,
        'completeness': qc['completeness'],
        'contamination': qc['contamination'],
        'classification': classification,
    })

if not candidates:
    raise RuntimeError(
        f"No bin found matching: clade={gtdbtk_clade}, "
        f"comp>={min_completeness}, cont<={max_contamination}, "
        f"gunc_pass={require_gunc_pass}\n"
        f"Available GTDB-Tk classifications:\n" +
        "\n".join(f"  {k}: {v}" for k, v in sorted(gtdbtk.items()))
    )

# Tiebreak
candidates.sort(key=lambda x: x[tiebreak if tiebreak != 'completeness' else 'completeness'],
                reverse=True)
selected = candidates[0]

print(f"[select_ancient_bin] Selected: {selected['bin']}", file=sys.stderr)
print(f"  Completeness: {selected['completeness']:.1f}%", file=sys.stderr)
print(f"  Contamination: {selected['contamination']:.2f}%", file=sys.stderr)
print(f"  GTDB-Tk: {selected['classification']}", file=sys.stderr)
if len(candidates) > 1:
    print(f"  ({len(candidates)} candidates — picked by {tiebreak})", file=sys.stderr)

with open(output_file, 'w') as f:
    f.write(selected['bin'] + '\n')
