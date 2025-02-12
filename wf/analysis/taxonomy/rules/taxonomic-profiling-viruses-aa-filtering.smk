def get_viruses_db_aa_map(db):
    for entry in config["viruses_db_aa"]:
        for key, value in entry.items():
            if key == db:
                return value["mapping_list"]


rule mmseqs_profile_viruses_aa_filter:
    input:
        viral_results=config["rdir"]
        + "/taxonomic-profiling-mapping-aa/"
        + "{indata}"
        + "/{viruses_db_aa}"
        + "/{smp}.tsv.gz",
    output:
        out_mm=config["rdir"]
        + "/taxonomic-profiling-filtering-aa/"
        + "{indata}"
        + "/{viruses_db_aa}"
        + "/{smp}.profile.viral_no-multimap.tsv.gz",
        out_cov=config["rdir"]
        + "/taxonomic-profiling-filtering-aa/"
        + "{indata}"
        + "/{viruses_db_aa}"
        + "/{smp}.profile.viral_cov-stats.tsv.gz",
        out_abun_agg=(
            config["rdir"]
            + "/taxonomic-profiling-filtering-aa/"
            + "{indata}"
            + "/{viruses_db_aa}"
            + "/{smp}.profile.viral_group-abundances-agg.tsv.gz"
        ),
        out_abun=(
            config["rdir"]
            + "/taxonomic-profiling-filtering-aa/"
            + "{indata}"
            + "/{viruses_db_aa}"
            + "/{smp}.profile.viral_group-abundances.tsv.gz"
        ),
    threads: config["mmseqs2_threads"]
    params:
        xfilter_bin=config["xfilter_bin"],
        viral_gene_list=lambda wc: get_viruses_db_aa_map(db=wc.viruses_db_aa),
        xfilter_parms=config["xfilter_parms"],
        rdir=config["rdir"]
        + "/taxonomic-profiling-filtering-aa/"
        + "{indata}"
        + "/{viruses_db_aa}",
        wdir=config["wdir"],
        label="{smp}.profile.viral",
        tmp="{smp}-tmp",
    conda:
        "../envs/xfilter.yaml"
    log:
        config["rdir"]
        + "/taxonomic-profiling-filtering-aa/"
        + "{indata}"
        + "/{viruses_db_aa}"
        + "/{smp}.viral-reads-filter.log",
    benchmark:
        config["rdir"]
        +"/taxonomic-profiling-filtering-aa/"
        +"{indata}"
        +"/{viruses_db_aa}"
        +"/{smp}.viral-reads-filter.bmk"
    message:
        """--- Annotate viruses on reads."""
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        {params.xfilter_bin} \
            -i {input.viral_results} \
            -p {params.label} \
            -m {params.viral_gene_list} \
            -t {threads} \
            {params.xfilter_parms}

        if [ ! -f {output.out_abun_agg} ]; then
            touch {output.out_mm}
            touch {output.out_cov}
            touch {output.out_abun_agg}
        fi

        rm -rf {params.tmp}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
