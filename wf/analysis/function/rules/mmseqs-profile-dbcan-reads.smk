rule mmseqs_profile_dbcan_reads:
    input:
        reads=config["rdir"] + "/read-derep/{smp}.fa.gz",
    output:
        out=config["rdir"] + "/dbcan-reads/{smp}.profile.dbcan.tsv.gz",
    threads: config["mmseqs2_threads"]
    params:
        mmseqs_bin=config["mmseqs_bin"],
        dbcan_db=config["dbcan_db"],
        mmseqs_dbcan_parms=config["mmseqs_dbcan_parms"],
        out=config["rdir"] + "/dbcan-reads/{smp}.profile.dbcan.tsv",
        rdir=config["rdir"] + "/dbcan-reads",
        wdir=config["wdir"],
        label="{smp}.profile.dbcan",
        tmp="{smp}-tmp",
    conda:
        "../envs/kegg-xfilter.yaml"
    log:
        config["rdir"] + "/logs/dbcan-reads/{smp}.dbcan-reads.log",
    benchmark:
        config["rdir"] + "/benchmarks/dbcan-reads/{smp}.dbcan-reads.bmk"
    message:
        """--- Annotate dbcans on reads."""
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        rm -rf {params.tmp}

        {params.mmseqs_bin} easy-search {input.reads} \
            {params.dbcan_db} \
            {params.out} \
            {params.tmp} \
            --threads {threads} \
            {params.mmseqs_dbcan_parms}

        pigz -p {threads} {params.out}

        rm -rf {params.tmp}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
