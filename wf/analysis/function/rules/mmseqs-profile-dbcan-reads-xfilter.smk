rule mmseqs_profile_dbcan_reads_xfilter:
    input:
        dbcan_results=config["rdir"] + "/dbcan-reads/{smp}.profile.dbcan.tsv.gz",
    output:
        out_mm=config["rdir"] + "/dbcan-reads/{smp}.profile.dbcan_no-multimap.tsv.gz",
        out_cov=config["rdir"] + "/dbcan-reads/{smp}.profile.dbcan_cov-stats.tsv.gz",
        out_abun=(
            config["rdir"]
            + "/dbcan-reads/{smp}.profile.dbcan_group-abundances-anvio.tsv.gz"
        ),
        out_abun_agg=(
            config["rdir"]
            + "/dbcan-reads/{smp}.profile.dbcan_group-abundances-agg.tsv.gz"
        ),
    threads: config["mmseqs2_threads"]
    params:
        xfilter_bin=config["xfilter_bin"],
        dbcan_gene_list=config["dbcan_gene_list"],
        mmseqs_dbcan_parms=config["mmseqs_dbcan_parms"],
        xfilter_parms=config["xfilter_parms"],
        annotation_source=config["dbcan_annotation_source"],
        rdir=config["rdir"] + "/dbcan-reads",
        wdir=config["wdir"],
        label="{smp}.profile.dbcan",
        tmp="{smp}-tmp",
    conda:
        "../envs/kegg-xfilter.yaml"
    log:
        config["rdir"] + "/logs/dbcan-reads/{smp}.dbcan-reads-filter.log",
    benchmark:
        config["rdir"] + "/benchmarks/dbcan-reads/{smp}.dbcan-reads-filter.bmk"
    message:
        """--- Annotate dbcans on reads."""
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        {params.xfilter_bin} \
            -i {input.dbcan_results} \
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
        fi

        rm -rf {params.tmp}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
