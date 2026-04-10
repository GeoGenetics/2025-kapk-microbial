rule mmseqs_profile_kegg_reads:
    input:
        reads=config["rdir"] + "/read-derep/{smp}.fa.gz",
    output:
        out=config["rdir"] + "/kegg-reads/{smp}.profile.kegg.tsv.gz",
    threads: config["mmseqs2_threads"]
    params:
        mmseqs_bin=config["mmseqs_bin"],
        kegg_db=config["kegg_db"],
        mmseqs_kegg_parms=config["mmseqs_kegg_parms"],
        out=config["rdir"] + "/kegg-reads/{smp}.profile.kegg.tsv",
        rdir=config["rdir"] + "/kegg-reads",
        wdir=config["wdir"],
        label="{smp}.profile.kegg",
        tmp="{smp}-tmp",
    conda:
        "../envs/kegg-xfilter.yaml"
    log:
        config["rdir"] + "/logs/kegg-reads/{smp}.kegg-reads.log",
    benchmark:
        config["rdir"] + "/benchmarks/kegg-reads/{smp}.kegg-reads.bmk"
    message:
        """--- Annotate KEGGs on reads."""
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        rm -rf {params.tmp}

        {params.mmseqs_bin} easy-search {input.reads} \
            {params.kegg_db} \
            {params.out} \
            {params.tmp} \
            --threads {threads} \
            {params.mmseqs_kegg_parms}

        pigz -p {threads} {params.out}

        rm -rf {params.tmp}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
