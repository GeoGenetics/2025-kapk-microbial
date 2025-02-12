rule initial_stats:
    input:
        reads=config["rdir"] + "/read-renamed/{smp}.fq.gz",
    output:
        stats=config["rdir"] + "/stats/{smp}.stats-initial.txt",
    threads: config["seqkit_threads"]
    params:
        seqkit_bin=config["seqkit_bin"],
    conda:
        "../envs/qc.yaml"
    log:
        config["rdir"] + "/logs/stats/{smp}.stats.log",
    benchmark:
        config["rdir"] + "/benchmarks/stats/{smp}.stats-initial.bmk"
    message:
        """--- Get initial stats"""
    shell:
        """
        {params.seqkit_bin} stats -j {threads} -T -a --quiet \
        {input.reads} -o {output.stats}
        """
