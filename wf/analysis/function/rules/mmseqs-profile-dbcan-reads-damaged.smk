rule mmseqs_profile_dbcan_reads_damaged:
    input:
        reads_dr=config["rdir"]
        + "/read-dmg-rank/standard/k1000/94/weight-1/{smp}.damaged.derep.fasta.gz",
    output:
        out_mm=config["rdir"]
        + "/dbcan-reads-damaged/{smp}.profile.damaged.dbcan_no-multimap.tsv.gz",
        out_cov=config["rdir"]
        + "/dbcan-reads-damaged/{smp}.profile.damaged.dbcan_cov-stats.tsv.gz",
        out_abun=(
            config["rdir"]
            + "/dbcan-reads-damaged/{smp}.profile.damaged.dbcan_group-abundances-anvio.tsv.gz"
        ),
        out_abun_agg=(
            config["rdir"]
            + "/dbcan-reads-damaged/{smp}.profile.damaged.dbcan_group-abundances-agg.tsv.gz"
        ),
        out=config["rdir"] + "/dbcan-reads-damaged/{smp}.profile.damaged.dbcan.tsv.gz",
    threads: config["mmseqs2_threads"]
    params:
        mmseqs_bin=config["mmseqs_bin"],
        xfilter_bin=config["xfilter_bin"],
        out=config["rdir"] + "/dbcan-reads-damaged/{smp}.profile.damaged.dbcan.tsv",
        dbcan_db=config["dbcan_db"],
        dbcan_gene_list=config["dbcan_gene_list"],
        mmseqs_dbcan_parms=config["mmseqs_dbcan_parms"],
        xfilter_parms=config["xfilter_parms"],
        annotation_source=config["dbcan_annotation_source"],
        rdir=config["rdir"] + "/dbcan-reads-damaged",
        wdir=config["wdir"],
        label="{smp}.profile.damaged.dbcan",
        tmp="{smp}-tmp",
    conda:
        "../envs/kegg-xfilter.yaml"
    log:
        config["rdir"] + "/logs/dbcan-reads-damaged/{smp}.dbcan-reads-damaged.log",
    benchmark:
        config["rdir"] + "/benchmarks/dbcan-reads-damaged/{smp}.dbcan-reads-damaged.bmk"
    message:
        """--- Annotate dbcans on reads."""
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        rm -rf {params.tmp}

        if [[ -s {input.reads_dr} ]]; then

            {params.mmseqs_bin} easy-search {input.reads_dr} \
                {params.dbcan_db} \
                {params.out} \
                {params.tmp} \
                --threads {threads} \
                {params.mmseqs_dbcan_parms}

            {params.xfilter_bin} \
                -i {params.out} \
                -p {params.label} \
                -m {params.dbcan_gene_list} \
                -t {threads} \
                --anvio \
                --annotation-source {params.annotation_source} \
                {params.xfilter_parms}

            if [ ! -f {output.out_abun} ]; then
                touch {output.out_abun}
                touch {output.out_mm}
                touch {output.out_cov}
                touch {output.out_abun_agg}
                touch {output.out}
            fi

            rm -rf {params.tmp}
            pigz -f {params.out}
        else
            touch {output.out_mm}
            touch {output.out_cov}
            touch {output.out_abun}
            touch {output.out_abun_agg}
            touch {output.out}
        fi

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
