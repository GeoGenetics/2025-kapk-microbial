def get_viruses_db_aa(db):
    for entry in config["viruses_db_aa"]:
        for key, value in entry.items():
            if key == db:
                return value["mmseqs_db"]


rule mmseqs_profile_viruses_aa_map:
    input:
        reads=config["rdir"] + "/read-derep/{smp}.fa.gz",
    output:
        out=config["rdir"]
        + "/taxonomic-profiling-mapping-aa/"
        + "{indata}"
        + "/{viruses_db_aa}"
        + "/{smp}.tsv.gz",
    threads: config["mmseqs2_threads"]
    params:
        mmseqs_bin=config["mmseqs_bin"],
        viral_db=lambda wc: get_viruses_db_aa(db=wc.viruses_db_aa),
        mmseqs_viral_parms=config["mmseqs_viral_parms"],
        out=config["rdir"]
        + "/taxonomic-profiling-mapping-aa/"
        + "{indata}"
        + "/{viruses_db_aa}"
        + "/{smp}.tsv",
        rdir=config["rdir"]
        + "/taxonomic-profiling-mapping-aa/"
        + "{indata}"
        + "/{viruses_db_aa}",
        wdir=config["wdir"],
        tmp="{smp}-tmp",
    conda:
        "../envs/xfilter.yaml"
    log:
        config["rdir"]
        + "/taxonomic-profiling-mapping-aa/"
        + "{indata}"
        + "/{viruses_db_aa}/{smp}.log",
    benchmark:
        config["rdir"]
        +"/taxonomic-profiling-mapping-aa/"
        +"{indata}"
        +"/{viruses_db_aa}/{smp}.bmk"
    message:
        """--- Annotate viruses on reads."""
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        rm -rf {params.tmp}

        {params.mmseqs_bin} easy-search {input.reads} \
            {params.viral_db} \
            {params.out} \
            {params.tmp} \
            --threads {threads} \
            {params.mmseqs_viral_parms}

        pigz -p {threads} {params.out}

        rm -rf {params.tmp}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
