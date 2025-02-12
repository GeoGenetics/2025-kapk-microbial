rule derep_reads:
    input:
        reads=config["rdir"] + "/read-extension/{smp}.extended.fastq.gz",
        stats=config["rdir"] + "/stats/{smp}.stats-extension.txt",
    output:
        derep_reads=config["rdir"] + "/read-derep/{smp}.fa.gz",
        derep_stats=config["rdir"] + "/stats/{smp}.stats-derep.txt",
    threads: config["read_derep_threads"]
    params:
        derep_reads=lambda wildcards: sample_table_read.derep_reads[wildcards.smp],
        use_vsearch=config["use_vsearch"],
        seqkit_bin=config["seqkit_bin"],
        vsearch_bin=config["vsearch_bin"],
        derep_minlen=config["read_minlen"],
        derep_parms=config["derep_parms"],
        derep_tsv=config["rdir"] + "/read-derep/{smp}.derep.tsv.gz",
    conda:
        "../envs/qc.yaml"
    log:
        config["rdir"] + "/logs/read-derep/{smp}.read-derep.log",
    benchmark:
        config["rdir"] + "/benchmarks/read-derep/{smp}.read-derep.bmk"
    message:
        """--- Dereplicate reads"""
    shell:
        """
        if [[ {params.derep_reads} == "TRUE" ]]; then
            if [[ {params.use_vsearch} == "TRUE" ]]; then
                # Derep and rename reads
                {params.vsearch_bin} \
                    --fastx_uniques \
                    {input.reads} \
                    --minseqlength {params.derep_minlen} \
                    {params.derep_parms} \
                | sed -e 's|;size=|__|' \
                | gzip > {output.derep_reads} 
            else
                {params.seqkit_bin} rmdup -j {threads} -s -i \
                    -D {params.derep_tsv} {input.reads} \
                | {params.seqkit_bin} replace -p "\s.+" -o {output.derep_reads} 
            fi

            {params.seqkit_bin} stats -j {threads} -T -a --quiet \
                {output.derep_reads} -o {output.derep_stats}
        else
            ln -sf {input.reads} {output.derep_reads}
            ln -sf {input.stats} {output.derep_stats}
        fi
        """
