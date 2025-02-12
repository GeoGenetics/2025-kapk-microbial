rule extension_stats:
    input:
        reads=config["rdir"] + "/read-extension/{smp}.extended.fastq.gz",
    output:
        stats=config["rdir"] + "/stats/{smp}.stats-extension.txt",
    threads: config["seqkit_threads"]
    params:
        seqkit_bin=config["seqkit_bin"],
    conda:
        "../envs/qc.yaml"
    log:
        config["rdir"] + "/logs/stats/{smp}.stats-extension.log",
    benchmark:
        config["rdir"] + "/benchmarks/stats/{smp}.stats-extension.bmk"
    message:
        """--- Get initial stats"""
    shell:
        """
        {params.seqkit_bin} stats -j {threads} -T -a --quiet \
        {input.reads} -o {output.stats}
        """
