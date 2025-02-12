rule rename_reads:
    input:
        reads=(
            lambda wildcards: config["sdir"]
            + "/"
            + sample_table_read.file[wildcards.smp]
        ),
    output:
        renamed_reads=config["rdir"] + "/read-renamed/{smp}.fq.gz",
        renamed_reads_mapping=config["rdir"] + "/read-renamed/{smp}.mapping.tsv.gz",
    threads: config["read_rename_threads"]
    params:
        rename_reads=config["rename_reads"],
        seqkit_bin=config["seqkit_bin"],
        rename_bin=config["rename_bin"],
        rename_parms=config["rename_parms"],
        minlen=config["read_minlen"],
        renamed_reads=config["rdir"] + "/read-renamed/{smp}.tmp.fq.gz",
        raw_ids=config["rdir"] + "/read-renamed/{smp}.raw.ids",
        pattern=".+",
        subst="{smp}_\{nr\}",
    conda:
        "../envs/qc.yaml"
    log:
        config["rdir"] + "/logs/read-renamed/{smp}.read-renamed.log",
    benchmark:
        config["rdir"] + "/benchmarks/read-renamed/{smp}.read-renamed.bmk"
    message:
        """--- Rename reads"""
    shell:
        """
        # Derep and rename reads
        if [[ {params.rename_reads} == "True" ]]; then
        {params.seqkit_bin} seq \
            -j {threads} -m {params.minlen} \
            {input.reads} \
            | {params.seqkit_bin} \
                replace \
                -j {threads} \
                -p {params.pattern} \
                -r {params.subst} \
                --nr-width 12 \
                -o {output.renamed_reads} 

        # Get mapping between ids
        {params.seqkit_bin} seq \
            -j {threads} -m {params.minlen} \
            {input.reads} \
            | {params.seqkit_bin} fx2tab -n -o {params.raw_ids}

        paste {params.raw_ids} <({params.seqkit_bin} fx2tab -n {output.renamed_reads}) \
            | gzip > {output.renamed_reads_mapping}

        rm -rf {params.renamed_reads} {params.raw_ids}

        else
            {params.seqkit_bin} seq \
                -j {threads} -m {params.minlen} \
                {input.reads} -o {output.renamed_reads}
            touch {output.renamed_reads_mapping}
        fi
        """
