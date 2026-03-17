#!/usr/bin/env python3
"""
Generate BEAST2 XML for tip-dating with BDSS tree prior.
"""

import re
from pathlib import Path

def read_fasta(path):
    """Read FASTA file."""
    seqs = {}
    name = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name:
            seqs[name] = "".join(seq)
    return seqs

def read_newick(path):
    """Read Newick tree file."""
    with open(path) as f:
        tree = f.read().strip()
    # Remove bootstrap values
    tree = re.sub(r'\)[0-9.]+:', '):', tree)
    return tree

def main():
    # Get parameters from snakemake
    alignment_path = snakemake.input.aln
    tree_path = snakemake.input.tree
    output_xml = snakemake.output.xml
    output_script = snakemake.output.script

    chain_length = snakemake.params.chain_length
    log_every = snakemake.params.log_every
    tip_dates = snakemake.params.tip_dates
    clock_rate_mean = snakemake.params.clock_rate_mean

    # Read alignment
    seqs = read_fasta(alignment_path)
    taxa = sorted(seqs.keys())
    aln_length = len(next(iter(seqs.values())))

    # Read starting tree
    starting_tree = read_newick(tree_path)

    # Build tip dates string
    tip_date_parts = []
    for taxon in taxa:
        date = tip_dates.get(taxon, 0)
        tip_date_parts.append(f"{taxon}={date}")
    tip_date_str = ",".join(tip_date_parts)

    # Generate XML
    xml = f'''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<beast namespace="beast.base.core:beast.base.inference:beast.base.evolution.alignment:
beast.base.evolution.tree:beast.base.evolution.tree.coalescent:
beast.pkgmgmt:beast.base.evolution.operator:beast.base.inference.operator:
beast.base.evolution.sitemodel:beast.base.evolution.substitutionmodel:
beast.base.evolution.likelihood:beast.base.evolution.branchratemodel:
beast.base.evolution.speciation:beast.base.inference.distribution:beast.base.inference.parameter" required="BEAST.base v2.7.0" version="2.7">

    <!-- Alignment -->
    <data id="alignment" spec="Alignment" dataType="nucleotide">
'''

    for taxon in taxa:
        seq = seqs[taxon]
        xml += f'        <sequence id="seq_{taxon}" spec="Sequence" taxon="{taxon}" value="{seq}"/>\n'

    xml += f'''    </data>

    <!-- Codon position partitions -->
    <data id="pos12" spec="FilteredAlignment" filter="1-{aln_length}\\3,2-{aln_length}\\3" data="@alignment"/>
    <data id="pos3" spec="FilteredAlignment" filter="3-{aln_length}\\3" data="@alignment"/>

    <!-- Tip dates -->
    <trait id="dateTrait" spec="beast.base.evolution.tree.TraitSet"
           traitname="date-backward" units="year"
           value="{tip_date_str}">
        <taxa id="taxonSet" spec="TaxonSet" alignment="@alignment"/>
    </trait>

    <run id="mcmc" spec="MCMC" chainLength="{chain_length}">
        <state id="state" spec="State" storeEvery="{log_every}">
            <tree id="Tree" spec="beast.base.evolution.tree.Tree" name="stateNode">
                <trait idref="dateTrait"/>
                <taxonset idref="taxonSet"/>
            </tree>
            <parameter id="clockRate" spec="parameter.RealParameter" name="stateNode">{clock_rate_mean}</parameter>
            <parameter id="birthRate" spec="parameter.RealParameter" name="stateNode">1.0</parameter>
            <parameter id="deathRate" spec="parameter.RealParameter" name="stateNode">0.5</parameter>
            <parameter id="samplingRate" spec="parameter.RealParameter" name="stateNode">0.01</parameter>
            <parameter id="kappa12" spec="parameter.RealParameter" lower="0.0" name="stateNode">2.0</parameter>
            <parameter id="kappa3" spec="parameter.RealParameter" lower="0.0" name="stateNode">2.0</parameter>
            <parameter id="gammaShape12" spec="parameter.RealParameter" lower="0.0" name="stateNode">1.0</parameter>
            <parameter id="gammaShape3" spec="parameter.RealParameter" lower="0.0" name="stateNode">1.0</parameter>
            <parameter id="freqs12" spec="parameter.RealParameter" dimension="4" lower="0.0" upper="1.0" name="stateNode">0.25 0.25 0.25 0.25</parameter>
            <parameter id="freqs3" spec="parameter.RealParameter" dimension="4" lower="0.0" upper="1.0" name="stateNode">0.25 0.25 0.25 0.25</parameter>
        </state>

        <init id="startingTree" spec="beast.base.evolution.tree.TreeParser"
              newick="{starting_tree}"
              initial="@Tree" taxa="@alignment"
              IsLabelledNewick="true" adjustTipHeights="false"/>

        <distribution id="posterior" spec="CompoundDistribution">
            <distribution id="prior" spec="CompoundDistribution">
                <!-- Birth-Death Serial Sampling tree prior (BDSKY 1.5.1) -->
                <distribution id="BDSS" spec="bdsky.evolution.speciation.BirthDeathSkylineModel"
                              tree="@Tree"
                              birthRate="@birthRate"
                              deathRate="@deathRate"
                              samplingRate="@samplingRate"
                              origin="120000000"/>

                <!-- Clock rate prior -->
                <distribution id="clockRatePrior" spec="beast.base.inference.distribution.Prior" x="@clockRate">
                    <distr id="clockRatePriorDist" spec="beast.base.inference.distribution.LogNormalDistributionModel" meanInRealSpace="true" M="{clock_rate_mean}" S="1.0"/>
                </distribution>

                <!-- Other priors -->
                <distribution id="birthRatePrior" spec="beast.base.inference.distribution.Prior" x="@birthRate">
                    <distr id="birthRatePriorDist" spec="beast.base.inference.distribution.Uniform" lower="0.0" upper="100.0"/>
                </distribution>
                <distribution id="kappa12Prior" spec="beast.base.inference.distribution.Prior" x="@kappa12">
                    <distr id="kappa12PriorDist" spec="beast.base.inference.distribution.LogNormalDistributionModel" M="1.0" S="1.25" meanInRealSpace="false"/>
                </distribution>
                <distribution id="kappa3Prior" spec="beast.base.inference.distribution.Prior" x="@kappa3">
                    <distr id="kappa3PriorDist" spec="beast.base.inference.distribution.LogNormalDistributionModel" M="1.0" S="1.25" meanInRealSpace="false"/>
                </distribution>
                <distribution id="gammaShape12Prior" spec="beast.base.inference.distribution.Prior" x="@gammaShape12">
                    <distr id="gammaShape12PriorDist" spec="beast.base.inference.distribution.Exponential" mean="1.0"/>
                </distribution>
                <distribution id="gammaShape3Prior" spec="beast.base.inference.distribution.Prior" x="@gammaShape3">
                    <distr id="gammaShape3PriorDist" spec="beast.base.inference.distribution.Exponential" mean="1.0"/>
                </distribution>
            </distribution>

            <distribution id="likelihood" spec="CompoundDistribution" useThreads="true">
                <distribution id="treeLikelihood12" spec="ThreadedTreeLikelihood" data="@pos12" tree="@Tree">
                    <siteModel id="siteModel12" spec="SiteModel" gammaCategoryCount="4" shape="@gammaShape12">
                        <substModel id="hky12" spec="HKY" kappa="@kappa12">
                            <frequencies id="freqs12Dist" spec="Frequencies" frequencies="@freqs12"/>
                        </substModel>
                    </siteModel>
                    <branchRateModel id="StrictClock" spec="StrictClockModel" clock.rate="@clockRate"/>
                </distribution>
                <distribution id="treeLikelihood3" spec="ThreadedTreeLikelihood" data="@pos3" tree="@Tree">
                    <siteModel id="siteModel3" spec="SiteModel" gammaCategoryCount="4" shape="@gammaShape3">
                        <substModel id="hky3" spec="HKY" kappa="@kappa3">
                            <frequencies id="freqs3Dist" spec="Frequencies" frequencies="@freqs3"/>
                        </substModel>
                    </siteModel>
                    <branchRateModel idref="StrictClock"/>
                </distribution>
            </distribution>
        </distribution>

        <!-- Operators -->
        <operator id="clockRateScaler" spec="ScaleOperator" parameter="@clockRate" scaleFactor="0.5" weight="3.0"/>
        <operator id="birthRateScaler" spec="ScaleOperator" parameter="@birthRate" scaleFactor="0.75" weight="3.0"/>
        <operator id="deathRateScaler" spec="ScaleOperator" parameter="@deathRate" scaleFactor="0.75" weight="1.0"/>
        <operator id="samplingRateScaler" spec="ScaleOperator" parameter="@samplingRate" scaleFactor="0.75" weight="1.0"/>
        <operator id="kappa12Scaler" spec="ScaleOperator" parameter="@kappa12" scaleFactor="0.75" weight="1.0"/>
        <operator id="kappa3Scaler" spec="ScaleOperator" parameter="@kappa3" scaleFactor="0.75" weight="1.0"/>
        <operator id="gammaShape12Scaler" spec="ScaleOperator" parameter="@gammaShape12" scaleFactor="0.75" weight="1.0"/>
        <operator id="gammaShape3Scaler" spec="ScaleOperator" parameter="@gammaShape3" scaleFactor="0.75" weight="1.0"/>
        <operator id="freqs12Delta" spec="DeltaExchangeOperator" parameter="@freqs12" delta="0.01" weight="1.0"/>
        <operator id="freqs3Delta" spec="DeltaExchangeOperator" parameter="@freqs3" delta="0.01" weight="1.0"/>

        <operator id="treeScaler" spec="ScaleOperator" tree="@Tree" scaleFactor="0.95" weight="3.0"/>
        <operator id="treeRootScaler" spec="ScaleOperator" tree="@Tree" rootOnly="true" scaleFactor="0.95" weight="3.0"/>
        <operator id="UniformOperator" spec="beast.base.evolution.operator.Uniform" tree="@Tree" weight="30.0"/>
        <operator id="SubtreeSlide" spec="SubtreeSlide" tree="@Tree" weight="15.0"/>
        <operator id="Narrow" spec="Exchange" tree="@Tree" weight="15.0"/>
        <operator id="Wide" spec="Exchange" tree="@Tree" isNarrow="false" weight="3.0"/>
        <operator id="WilsonBalding" spec="WilsonBalding" tree="@Tree" weight="3.0"/>

        <operator id="updownClock" spec="UpDownOperator" scaleFactor="0.75" weight="3.0">
            <up idref="clockRate"/>
            <down idref="Tree"/>
        </operator>

        <!-- Loggers -->
        <logger id="tracelog" spec="Logger" fileName="codon_tipdating.log" logEvery="{log_every}" model="@posterior" sort="smart">
            <log idref="posterior"/>
            <log idref="likelihood"/>
            <log idref="prior"/>
            <log idref="treeLikelihood12"/>
            <log idref="treeLikelihood3"/>
            <log id="TreeHeight" spec="beast.base.evolution.tree.TreeHeightLogger" tree="@Tree"/>
            <log idref="clockRate"/>
            <log idref="birthRate"/>
            <log idref="deathRate"/>
            <log idref="samplingRate"/>
            <log idref="kappa12"/>
            <log idref="kappa3"/>
            <log idref="gammaShape12"/>
            <log idref="gammaShape3"/>
            <log idref="BDSS"/>
        </logger>

        <logger id="screenlog" spec="Logger" logEvery="100000">
            <log idref="posterior"/>
            <log idref="clockRate"/>
            <log idref="TreeHeight"/>
        </logger>

        <logger id="treelog" spec="Logger" fileName="codon_tipdating.trees" logEvery="{log_every}" mode="tree">
            <log id="TreeWithMetaDataLogger" spec="beast.base.evolution.TreeWithMetaDataLogger" tree="@Tree"/>
        </logger>
    </run>
</beast>
'''

    with open(output_xml, 'w') as f:
        f.write(xml)

    # Create submission script
    script = f'''#!/bin/bash
#SBATCH --job-name=beast2_phylo
#SBATCH --partition=compregular
#SBATCH --gres=gpu:a100:2
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=beast2_%j.out
#SBATCH --error=beast2_%j.err

source /maps/projects/fernandezguerra/apps/opt/conda/etc/profile.d/conda.sh
conda activate beast2

# BEAGLE GPU library path
export LD_LIBRARY_PATH="/maps/projects/caeg/people/kbd606/scratch/kapk-assm/119_50/phylogenomics/analysis/beast2_correct/beagle_build/beagle-lib/build/libhmsbeagle:$LD_LIBRARY_PATH"

cd $(dirname $0)

echo "=== BEAST2 Tip-Dating ==="
echo "Node: $(hostname)"
echo "Date: $(date)"
nvidia-smi

# Use both GPUs: -beagle_order 0,1 assigns partitions to each GPU
beast -beagle_GPU -beagle_order 0,1 -threads 8 -overwrite codon_tipdating.xml

echo ""
echo "=== Done ==="
echo "Date: $(date)"
'''

    with open(output_script, 'w') as f:
        f.write(script)

    # Also create CPU-only version for clusters without GPU
    script_cpu = script.replace('--gres=gpu:a100:1', '')
    script_cpu = script_cpu.replace('-beagle_GPU', '-beagle_CPU -beagle_SSE')
    script_cpu = script_cpu.replace('nvidia-smi', '# No GPU')
    script_cpu = script_cpu.replace('beast2_phylo', 'beast2_cpu')

    cpu_script_path = str(output_script).replace('run_beast.sh', 'run_beast_cpu.sh')
    with open(cpu_script_path, 'w') as f:
        f.write(script_cpu)

if __name__ == "__main__":
    main()
