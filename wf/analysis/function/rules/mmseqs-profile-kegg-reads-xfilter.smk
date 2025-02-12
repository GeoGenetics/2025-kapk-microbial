rule mmseqs_profile_kegg_reads_xfilter:
    input:
        kegg_results=config["rdir"] + "/kegg-reads/{smp}.profile.kegg.tsv.gz",
    output:
        out_mm=config["rdir"] + "/kegg-reads/{smp}.profile.kegg_no-multimap.tsv.gz",
        out_cov=config["rdir"] + "/kegg-reads/{smp}.profile.kegg_cov-stats.tsv.gz",
        out_abun=(
            config["rdir"]
            + "/kegg-reads/{smp}.profile.kegg_group-abundances-anvio.tsv.gz"
        ),
        out_abun_agg=(
            config["rdir"]
            + "/kegg-reads/{smp}.profile.kegg_group-abundances-agg.tsv.gz"
        ),
    threads: config["mmseqs2_threads"]
    params:
        xfilter_bin=config["xfilter_bin"],
        kegg_gene_list=config["kegg_gene_list"],
        mmseqs_kegg_parms=config["mmseqs_kegg_parms"],
        xfilter_parms=config["xfilter_parms"],
        annotation_source=config["annotation_source"],
        rdir=config["rdir"] + "/kegg-reads",
        wdir=config["wdir"],
        label="{smp}.profile.kegg",
        tmp="{smp}-tmp",
    conda:
        "../envs/kegg-xfilter.yaml"
    log:
        config["rdir"] + "/logs/kegg-reads/{smp}.kegg-reads-filter.log",
    benchmark:
        config["rdir"] + "/benchmarks/kegg-reads/{smp}.kegg-reads-filter.bmk"
    message:
        """--- Annotate KEGGs on reads."""
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        {params.xfilter_bin} \
            -i {input.kegg_results} \
            -p {params.label} \
            -m {params.kegg_gene_list} \
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
