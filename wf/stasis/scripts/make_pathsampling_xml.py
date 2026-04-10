#!/usr/bin/env python3
"""
Create BEAST2 path-sampling XMLs for marginal likelihood estimation.

PRIMARY comparison (model selection at 2 Ma — the scientifically correct test):
  ps_strict_2Ma   vs  ps_2clock_stasis  vs  ps_ucld_2Ma
  2×ln BF(2clock vs strict) is the formal Bayesian test of the stasis hypothesis.

SECONDARY comparison (fixed-age sensitivity under strict clock — shows rate anomaly):
  ps_0Ma  ps_0.5Ma  ps_1Ma  ps_2Ma
  These confirm the ancient branch is anomalous under a strict clock, but should NOT
  be interpreted as the molecular data "preferring" a younger age — they just reflect
  that a rate-homogeneous model fits better when the outlier branch is shorter.
"""
import re
import sys

XMLDIR = "/maps/projects/caeg/people/kbd606/scratch/kapk-assm/amber/paper_run/phylogenomics_79taxa_bak/09_beast2/stasis_analysis/no_dec"
OUTDIR = "/maps/projects/caeg/people/kbd606/scratch/kapk-assm/amber/paper_run/beast2_nodec"
PSDIR  = "/maps/projects/caeg/people/kbd606/scratch/kapk-assm/amber/paper_run/beast2_nodec"

# (source_xml, ps_label, run_dir_name)
# PRIMARY: model comparison at 2 Ma fixed tip date
MODELS_PRIMARY = [
    ("codon_tipdating_nodec.xml",        "ps_strict_2Ma",      "ps_strict_2Ma"),
    ("codon_2clock_stasis_nodec.xml",    "ps_2clock_stasis",   "ps_2clock_stasis"),
    ("codon_ucld_tipdating_nodec.xml",   "ps_ucld_2Ma",        "ps_ucld_2Ma"),
]

# SECONDARY: fixed-age sensitivity under strict clock (already run; kept for reference)
MODELS_SECONDARY = [
    ("codon_tipdating_0Ma_nodec.xml",    "ps_0Ma",   "ps_0Ma"),
    ("codon_tipdating_0.5Ma_nodec.xml",  "ps_0.5Ma", "ps_0.5Ma"),
    ("codon_tipdating_1Ma_nodec.xml",    "ps_1Ma",   "ps_1Ma"),
    ("codon_tipdating_nodec.xml",        "ps_2Ma",   "ps_2Ma"),
]

MODELS = MODELS_PRIMARY + MODELS_SECONDARY

# Path sampling settings
NR_STEPS      = 64
CHAIN_LENGTH  = 2_000_000   # per step
ALPHA         = 0.3         # power posterior schedule
PREBURNIN     = 1_000_000   # pre-burnin samples
BURNIN_PCT    = 50          # % of chain used as burnin per step

PS_CMD = "cd $(dir)\n$(java) -cp $(java.class.path) beast.pkgmgmt.launcher.BeastLauncher $(resume/overwrite) -seed $(seed) beast.xml"

def make_ps_xml(src_xml, label, run_dir):
    with open(f"{XMLDIR}/{src_xml}") as f:
        content = f.read()

    # Find the run element opening tag line
    # Pattern: <run id="mcmc" spec="MCMC" chainLength="...">
    run_open_pat = re.compile(
        r'<run\s+id="mcmc"\s+spec="MCMC"\s+chainLength="\d+">'
    )
    m = run_open_pat.search(content)
    if not m:
        print(f"ERROR: could not find run element in {src_xml}")
        sys.exit(1)

    ps_rootdir = f"{PSDIR}/{run_dir}"

    ps_open = (
        f'<run spec="modelselection.inference.PathSampler"\n'
        f'         nrOfSteps="{NR_STEPS}"\n'
        f'         chainLength="{CHAIN_LENGTH}"\n'
        f'         alpha="{ALPHA}"\n'
        f'         rootdir="{ps_rootdir}"\n'
        f'         burnInPercentage="{BURNIN_PCT}"\n'
        f'         preBurnin="{PREBURNIN}"\n'
        f'         deleteOldLogs="true">\n'
        f'    {PS_CMD}\n'
        f'    <mcmc id="mcmc" spec="MCMC" chainLength="{CHAIN_LENGTH}">'
    )

    content = run_open_pat.sub(ps_open, content, count=1)

    # Replace closing </run> with </mcmc>\n    </run>
    # The closing tag is the last </run> in the file
    last_run_close = content.rfind('</run>')
    if last_run_close == -1:
        print(f"ERROR: could not find closing </run> in {src_xml}")
        sys.exit(1)
    content = content[:last_run_close] + '    </mcmc>\n    </run>' + content[last_run_close+6:]

    # Write output
    outfile = f"{OUTDIR}/{label}.xml"
    with open(outfile, 'w') as f:
        f.write(content)
    print(f"Written: {outfile}")
    return outfile


for src_xml, label, run_dir in MODELS:
    make_ps_xml(src_xml, label, run_dir)

print("Done. All path-sampling XMLs created.")
