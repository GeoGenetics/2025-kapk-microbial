#!/usr/bin/env python3
"""
Generate BEAST2 XML for SNP-based tip-dating with ASC correction.

Uses tv-only SNP alignment (9 taxa) with:
  - Ascertainment bias correction via constantSiteWeights (invariant site counts
    from reference genome minus SNP positions)
  - Clock rate prior recentered to ~1e-9/site/year (appropriate for ancient bacteria)
  - Coalescent constant-size tree prior (simpler than BDSS for 9 taxa / 1 ancient tip)
  - Single HKY+G4 partition (no codon partitioning for SNP data)
  - adjustTipHeights=true for date-consistent starting tree
"""

import re
from pathlib import Path
from collections import Counter


def read_fasta(path):
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
                seq.append(line.upper())
        if name:
            seqs[name] = "".join(seq)
    return seqs


def count_invariant_sites(ref_fasta_path, n_variable_sites):
    """
    Count A/C/G/T in the reference genome.
    Subtract variable site counts proportionally to get invariant site weights.
    Returns (countA, countC, countG, countT) for constant sites.
    """
    seqs = read_fasta(ref_fasta_path)
    genome = "".join(seqs.values()).upper()
    total = Counter(c for c in genome if c in "ACGT")
    total_acgt = sum(total.values())
    # Fraction of genome that is invariant
    invariant = max(total_acgt - n_variable_sites, 0)
    frac = invariant / total_acgt if total_acgt > 0 else 1.0
    return (
        int(total["A"] * frac),
        int(total["C"] * frac),
        int(total["G"] * frac),
        int(total["T"] * frac),
    )


def read_newick(path):
    with open(path) as f:
        tree = f.read().strip()
    # Take first tree if .contree / multi-tree
    if "\n" in tree:
        tree = tree.split("\n")[0]
    # Remove bootstrap values
    tree = re.sub(r'\)[0-9.]+:', '):', tree)
    return tree


def main():
    alignment_path  = snakemake.input.aln
    ref_genome_path = snakemake.input.ref
    tree_path       = snakemake.input.tree
    output_xml      = snakemake.output.xml

    chain_length     = snakemake.params.chain_length
    log_every        = snakemake.params.log_every
    tip_dates        = snakemake.params.tip_dates
    clock_rate_mean  = snakemake.params.clock_rate_mean  # expected ~1e-9

    seqs = read_fasta(alignment_path)
    taxa = sorted(seqs.keys())
    aln_length = len(next(iter(seqs.values())))

    # Ascertainment bias correction: count invariant sites from reference
    cA, cC, cG, cT = count_invariant_sites(ref_genome_path, aln_length)

    starting_tree = read_newick(tree_path)

    # Build tip dates trait
    tip_date_parts = [f"{t}={tip_dates.get(t, 0)}" for t in taxa]
    tip_date_str = ",".join(tip_date_parts)

    xml = f'''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<beast namespace="beast.base.core:beast.base.inference:beast.base.evolution.alignment:
beast.base.evolution.tree:beast.base.evolution.tree.coalescent:
beast.pkgmgmt:beast.base.evolution.operator:beast.base.inference.operator:
beast.base.evolution.sitemodel:beast.base.evolution.substitutionmodel:
beast.base.evolution.likelihood:beast.base.evolution.branchratemodel:
beast.base.inference.distribution:beast.base.inference.parameter" required="BEAST.base v2.7.0" version="2.7">

    <!-- SNP alignment (tv-only, {aln_length} sites) -->
    <data id="raw_alignment" spec="Alignment" dataType="nucleotide">
'''
    for taxon in taxa:
        xml += f'        <sequence id="seq_{taxon}" spec="Sequence" taxon="{taxon}" value="{seqs[taxon]}"/>\n'

    xml += f'''    </data>

    <!-- ASC correction: add invariant site weights from reference genome -->
    <data id="alignment" spec="FilteredAlignment" data="@raw_alignment"
          filter="1-{aln_length}"
          constantSiteWeights="{cA} {cC} {cG} {cT}"/>

    <!-- Tip dates -->
    <trait id="dateTrait" spec="beast.base.evolution.tree.TraitSet"
           traitname="date-backward" units="year"
           value="{tip_date_str}">
        <taxa id="taxonSet" spec="TaxonSet" alignment="@alignment"/>
    </trait>

    <!-- TreeIntervals defined at beast level so coalescent can reference it -->
    <treeIntervals id="TreeIntervals" spec="beast.base.evolution.tree.TreeIntervals" tree="@Tree"/>

    <run id="mcmc" spec="MCMC" chainLength="{chain_length}">
        <state id="state" spec="State" storeEvery="{log_every}">
            <tree id="Tree" spec="beast.base.evolution.tree.Tree" name="stateNode">
                <trait idref="dateTrait"/>
                <taxonset idref="taxonSet"/>
            </tree>
            <parameter id="clockRate"   spec="parameter.RealParameter" name="stateNode">{clock_rate_mean}</parameter>
            <parameter id="popSize"     spec="parameter.RealParameter" lower="0.0" name="stateNode">1000000.0</parameter>
            <parameter id="kappa"       spec="parameter.RealParameter" lower="0.0" name="stateNode">2.0</parameter>
            <parameter id="gammaShape"  spec="parameter.RealParameter" lower="0.0" name="stateNode">1.0</parameter>
            <parameter id="freqs"       spec="parameter.RealParameter" dimension="4" lower="0.0" upper="1.0" name="stateNode">0.25 0.25 0.25 0.25</parameter>
        </state>

        <!-- Date-consistent starting tree -->
        <init id="startingTree" spec="beast.base.evolution.tree.TreeParser"
              newick="{starting_tree}"
              initial="@Tree" taxa="@alignment"
              IsLabelledNewick="true" adjustTipHeights="true"/>

        <distribution id="posterior" spec="CompoundDistribution">
            <distribution id="prior" spec="CompoundDistribution">

                <!-- Coalescent constant-size tree prior -->
                <distribution id="coalescent" spec="beast.base.evolution.tree.coalescent.Coalescent" treeIntervals="@TreeIntervals">
                    <populationModel id="constPop" spec="beast.base.evolution.tree.coalescent.ConstantPopulation" popSize="@popSize"/>
                </distribution>

                <!-- Population size: Jeffreys (1/x) approximated as LogNormal -->
                <distribution id="popSizePrior" spec="beast.base.inference.distribution.Prior" x="@popSize">
                    <distr spec="beast.base.inference.distribution.OneOnX"/>
                </distribution>

                <!-- Clock rate prior: LogNormal centred at {clock_rate_mean} /site/year
                     S=1.5 gives ~3 orders of magnitude uncertainty -->
                <distribution id="clockRatePrior" spec="beast.base.inference.distribution.Prior" x="@clockRate">
                    <distr spec="beast.base.inference.distribution.LogNormalDistributionModel"
                           meanInRealSpace="true" M="{clock_rate_mean}" S="1.5"/>
                </distribution>

                <distribution id="kappaPrior" spec="beast.base.inference.distribution.Prior" x="@kappa">
                    <distr spec="beast.base.inference.distribution.LogNormalDistributionModel" M="1.0" S="1.25" meanInRealSpace="false"/>
                </distribution>
                <distribution id="gammaShapePrior" spec="beast.base.inference.distribution.Prior" x="@gammaShape">
                    <distr spec="beast.base.inference.distribution.Exponential" mean="1.0"/>
                </distribution>
            </distribution>

            <distribution id="likelihood" spec="CompoundDistribution">
                <distribution id="treeLikelihood" spec="ThreadedTreeLikelihood" data="@alignment" tree="@Tree">
                    <siteModel id="siteModel" spec="SiteModel" gammaCategoryCount="4" shape="@gammaShape">
                        <substModel id="hky" spec="HKY" kappa="@kappa">
                            <frequencies id="freqsDist" spec="Frequencies" frequencies="@freqs"/>
                        </substModel>
                    </siteModel>
                    <branchRateModel id="StrictClock" spec="StrictClockModel" clock.rate="@clockRate"/>
                </distribution>
            </distribution>
        </distribution>

        <!-- Operators -->
        <operator id="clockRateScaler"  spec="ScaleOperator"         parameter="@clockRate"  scaleFactor="0.5"  weight="3.0"/>
        <operator id="popSizeScaler"    spec="ScaleOperator"         parameter="@popSize"    scaleFactor="0.75" weight="3.0"/>
        <operator id="kappaScaler"      spec="ScaleOperator"         parameter="@kappa"      scaleFactor="0.75" weight="1.0"/>
        <operator id="gammaScaler"      spec="ScaleOperator"         parameter="@gammaShape" scaleFactor="0.75" weight="1.0"/>
        <operator id="freqsDelta"       spec="DeltaExchangeOperator" parameter="@freqs"      delta="0.01"       weight="1.0"/>
        <operator id="treeScaler"       spec="ScaleOperator"         tree="@Tree"            scaleFactor="0.95" weight="3.0"/>
        <operator id="treeRootScaler"   spec="ScaleOperator"         tree="@Tree"            rootOnly="true"    scaleFactor="0.95" weight="3.0"/>
        <operator id="UniformOperator"  spec="beast.base.evolution.operator.Uniform" tree="@Tree" weight="30.0"/>
        <operator id="SubtreeSlide"     spec="SubtreeSlide"          tree="@Tree"            weight="15.0"/>
        <operator id="Narrow"           spec="Exchange"              tree="@Tree"            weight="15.0"/>
        <operator id="Wide"             spec="Exchange"              tree="@Tree"            isNarrow="false"   weight="3.0"/>
        <operator id="WilsonBalding"    spec="WilsonBalding"         tree="@Tree"            weight="3.0"/>
        <operator id="updownClock"      spec="UpDownOperator"        scaleFactor="0.75"      weight="3.0">
            <up idref="clockRate"/>
            <down idref="Tree"/>
        </operator>

        <!-- Loggers -->
        <logger id="tracelog" spec="Logger" fileName="snp_tipdating.log" logEvery="{log_every}" model="@posterior" sort="smart">
            <log idref="posterior"/>
            <log idref="likelihood"/>
            <log idref="prior"/>
            <log idref="treeLikelihood"/>
            <log id="TreeHeight" spec="beast.base.evolution.tree.TreeHeightLogger" tree="@Tree"/>
            <log idref="clockRate"/>
            <log idref="popSize"/>
            <log idref="kappa"/>
            <log idref="gammaShape"/>
        </logger>
        <logger id="screenlog" spec="Logger" logEvery="100000">
            <log idref="posterior"/>
            <log idref="clockRate"/>
            <log idref="TreeHeight"/>
        </logger>
        <logger id="treelog" spec="Logger" fileName="snp_tipdating.trees" logEvery="{log_every}" mode="tree">
            <log id="TreeWithMetaDataLogger" spec="beast.base.evolution.TreeWithMetaDataLogger" tree="@Tree"/>
        </logger>
    </run>
</beast>
'''

    with open(output_xml, 'w') as f:
        f.write(xml)


if __name__ == "__main__":
    main()
