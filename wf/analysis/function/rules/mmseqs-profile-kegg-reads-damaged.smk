rule mmseqs_profile_kegg_reads_damaged:
    input:
        reads=config["rdir"]
        + "/read-dmg-rank/standard/k1000/94/weight-1/{smp}.damaged.derep.fasta.gz",
    output:
        out_mm=config["rdir"]
        + "/kegg-reads-damaged/{smp}.profile.damaged.kegg_no-multimap.tsv.gz",
        out_cov=config["rdir"]
        + "/kegg-reads-damaged/{smp}.profile.damaged.kegg_cov-stats.tsv.gz",
        out_abun=(
            config["rdir"]
            + "/kegg-reads-damaged/{smp}.profile.damaged.kegg_group-abundances-anvio.tsv.gz"
        ),
        out_abun_agg=(
            config["rdir"]
            + "/kegg-reads-damaged/{smp}.profile.damaged.kegg_group-abundances-agg.tsv.gz"
        ),
        out=config["rdir"] + "/kegg-reads-damaged/{smp}.profile.damaged.kegg.tsv.gz",
    threads: config["mmseqs2_threads"]
    params:
        mmseqs_bin=config["mmseqs_bin"],
        xfilter_bin=config["xfilter_bin"],
        out=config["rdir"] + "/kegg-reads-damaged/{smp}.profile.damaged.kegg.tsv",
        kegg_db=config["kegg_db"],
        kegg_gene_list=config["kegg_gene_list"],
        mmseqs_kegg_parms=config["mmseqs_kegg_parms"],
        xfilter_parms=config["xfilter_parms"],
        annotation_source=config["annotation_source"],
        rdir=config["rdir"] + "/kegg-reads-damaged",
        wdir=config["wdir"],
        label="{smp}.profile.damaged.kegg",
        tmp="{smp}-tmp",
    conda:
        "../envs/kegg-xfilter.yaml"
    log:
        config["rdir"] + "/logs/kegg-reads-damaged/{smp}.kegg-reads-damaged.log",
    benchmark:
        config["rdir"] + "/benchmarks/kegg-reads-damaged/{smp}.kegg-reads-damaged.bmk"
    message:
        """--- Annotate KEGGs on reads."""
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        rm -rf {params.tmp}

        if [[ -s {input.reads} ]]; then

            {params.mmseqs_bin} easy-search {input.reads} \
                {params.kegg_db} \
                {params.out} \
                {params.tmp} \
                --threads {threads} \
                {params.mmseqs_kegg_parms}

            {params.xfilter_bin} \
                -i {params.out} \
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
