rule reads_dmg_distribution_metaDMG:
    input:
        bam=(
            config["rdir"]
            + "/reads-dmg-distribution-mags-filtering/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{mag_db}"
            + "/{dr_fraction}"
            + "/{smp}.dmg-reads.dedup.filtered.bam"
        ),
    output:
        mdmg_outdir=(
            config["rdir"]
            + "/reads-dmg-distribution-mags-metaDMG/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/local"
            + "/{mag_db}"
            + "/{dr_fraction}"
            + "/{smp}.tp-mdmg.weight-{weight}.tar.gz"
        ),
        mdmg_csv=(
            config["rdir"]
            + "/reads-dmg-distribution-mags-metaDMG/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/local"
            + "/{mag_db}"
            + "/{dr_fraction}"
            + "/{smp}.tp-mdmg.weight-{weight}.csv.gz"
        ),
        mdmg_misincorporation=(
            config["rdir"]
            + "/reads-dmg-distribution-mags-metaDMG/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/local"
            + "/{mag_db}"
            + "/{dr_fraction}"
            + "/{smp}.tp-mdmg.misincorporation.weight-{weight}.txt.gz"
        ),
        # mdmg_pmd=(
        #     config["rdir"]
        #     + "/taxonomic-profiling-dmg/{smp}.tp-mdmg.{dmg_mode}.pmd.weight-{weight}.txt.gz"
        # ),
    threads: config["mdmg_threads"]
    params:
        mdmg_bin=config["mdmg_bin"],
        mdmg_cpp=config["mdmg_cpp"],
        mdmg_config="{smp}-config.weight-{weight}.yaml",
        mdmg_fixncbi=config["mdmg_fixncbi"],
        mdmg_weight="{weight}",
        mdmg_names=config["mdmg_names"],
        mdmg_nodes=config["mdmg_nodes"],
        mdmg_acc2tax=config["mdmg_acc2tax"],
        mdmg_cores=config["mdmg_cores"],
        mdmg_cores_pr_fit=config["mdmg_cores_pr_fit"],
        mdmg_outdir="{smp}.tp-mdmg.weight-{weight}",
        mdmg_positions=config["mdmg_positions"],
        mdmg_lcarank=config["mdmg_lcarank"],
        mdmg_csv=(
            config["rdir"]
            + "/reads-dmg-distribution-mags-metaDMG/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/local"
            + "/{mag_db}"
            + "/{dr_fraction}"
            + "/{smp}.tp-mdmg.weight-{weight}.csv"
        ),
        mdmg_simscorelow=config["mdmg_simscorelow"],
        mdmg_extra=config["mdmg_extra"],
        mdmg_libprep=lambda wildcards: sample_table_read.libprep[wildcards.smp],
        mdmg_misincorporation=(
            config["rdir"]
            + "/reads-dmg-distribution-mags-metaDMG/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/local"
            + "/{mag_db}"
            + "/{dr_fraction}"
            + "/{smp}.tp-mdmg.misincorporation.weight-{weight}.txt"
        ),
        # mdmg_pmd="{smp}.tp-mdmg.{dmg_mode}.weight-{weight}/pmd/{smp}.pmd.txt.gz",
        rdir=config["rdir"]
        + "/reads-dmg-distribution-mags-metaDMG/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/local"
        + "/{mag_db}"
        + "/{dr_fraction}",
        wdir=config["wdir"],
        label="{smp}",
    log:
        config["rdir"]
        + "/logs/reads-dmg-distribution-mags-metaDMG/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/local"
        + "/{mag_db}"
        + "/{dr_fraction}"
        + "/{smp}.tp-mdmg.weight-{weight}.log",
    benchmark:
        (
            config["rdir"]
            + "/benchmarks/reads-dmg-distribution-mags-metaDMG/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/local"
            + "/{mag_db}"
            + "/{dr_fraction}"
            + "/{smp}.tp-mdmg.weight-{weight}.bmk"
        )
    conda:
        "../envs/metaDMG.yaml"
    message:
        """--- Damage estimation with metaDMG."""
    shell:
        """
        set -x

        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        # check if BAM file has size > 0
        if [ ! -s {input.bam} ]; then
            echo "BAM file is empty"
            touch {output.mdmg_outdir}
            touch {output.mdmg_csv}
            touch {output.mdmg_misincorporation}
            cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
            exit 0
        fi

        if [ -d {params.mdmg_outdir} ]; then
            rm -rf {params.mdmg_outdir}.tmp
            mkdir -p {params.mdmg_outdir}.tmp
        else
            mkdir -p {params.mdmg_outdir}.tmp
        fi

        cd {params.mdmg_outdir}.tmp || {{ echo "Cannot change dir"; exit 1; }}

        # Run metaDMG
        if [[ {params.mdmg_libprep} == "double" ]]; then
            {params.mdmg_bin} config  \
                --overwrite \
                --config-file {params.mdmg_config} \
                --custom-database \
                --weight-type {params.mdmg_weight} \
                --names {params.mdmg_names} \
                --nodes {params.mdmg_nodes} \
                --acc2tax {params.mdmg_acc2tax} \
                --metaDMG-cpp {params.mdmg_cpp} \
                --parallel-samples {params.mdmg_cores} \
                --cores-per-sample {params.mdmg_cores_pr_fit} \
                --output-dir {params.mdmg_outdir} \
                --max-position {params.mdmg_positions}  \
                --lca-rank {params.mdmg_lcarank} \ \
                --min-similarity-score {params.mdmg_simscorelow} \
                --damage-mode local {params.mdmg_extra} \
                {input.bam}
        else
            {params.mdmg_bin} config  \
                --overwrite \
                --config-file {params.mdmg_config} \
                --custom-database \
                --weight-type {params.mdmg_weight} \
                --names {params.mdmg_names} \
                --nodes {params.mdmg_nodes} \
                --acc2tax {params.mdmg_acc2tax} \
                --metaDMG-cpp {params.mdmg_cpp} \
                --parallel-samples {params.mdmg_cores} \
                --cores-per-sample {params.mdmg_cores_pr_fit} \
                --output-dir {params.mdmg_outdir} \
                --max-position {params.mdmg_positions}  \
                --lca-rank {params.mdmg_lcarank} \
                --bayesian \
                --min-similarity-score {params.mdmg_simscorelow} \
                --damage-mode local \
                --forward-only {params.mdmg_extra} \
                {input.bam}
        fi

        {params.mdmg_bin} compute {params.mdmg_config}

        # check if folder exists, if not, exit
        if [ ! -d {params.mdmg_outdir}/results ]; then
            echo "Cannot find output folder"
            touch {output.mdmg_outdir}
            touch {output.mdmg_csv}
            touch {output.mdmg_misincorporation}
            exit 0
        fi

        {params.mdmg_bin} convert --add-fit-predictions --output {params.mdmg_csv} --results {params.mdmg_outdir}/results

        # Get misincorporation file
        {params.mdmg_bin} mismatch-to-mapDamage {params.mdmg_outdir}/mismatches/{params.label}.mismatches.parquet
        mv misincorporation.txt {params.mdmg_misincorporation}

        tar cvfz {params.mdmg_outdir}.tar.gz {params.mdmg_outdir}
        mv {params.mdmg_outdir}.tar.gz {params.rdir}/
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        gzip {params.mdmg_csv}
        gzip {params.mdmg_misincorporation}

        rm -rf {params.mdmg_outdir}.tmp

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
