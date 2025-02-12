rule extend_reads:
    input:
        reads=config["rdir"] + "/read-renamed/{smp}.fq.gz",
        stats=config["rdir"] + "/stats/{smp}.stats-initial.txt",
    output:
        extended=config["rdir"] + "/read-extension/{smp}.extended.fastq.gz",
        stats=config["rdir"] + "/stats/{smp}.stats-extension.txt",
    threads: config["read_extension_threads"]
    params:
        extend_reads=lambda wildcards: sample_table_read.extend_reads[wildcards.smp],
        seqkit_bin=config["seqkit_bin"],
        tadpole_bin=config["tadpole_bin"],
        loglog_bin=config["loglog_bin"],
        mode=config["tadpole_mode"],
        length=config["tadpole_length"],
        memory=config["tadpole_memory"],
        k=config["tadpole_k"],
        ibb=config["ibb"],
        ecc=config["ecc"],
        ecco=config["ecco"],
        conservative=config["conservative"],
        seed=config["loglog_seed"],
        bits=config["loglog_bits"],
        maxmem=config["maxmem"],
        hashes=config["loglog_hashes"],
        table_cap=config["loglog_table_cap"],
        trimends=config["trimends"],
        prefilter=config["prefilter"],
    conda:
        "../envs/qc.yaml"
    log:
        config["rdir"] + "/logs/read-extension/{smp}.read-extension.log",
    benchmark:
        config["rdir"] + "/benchmarks/read-extension/{smp}.read-extension.bmk"
    message:
        """--- READ EXTENSION"""
    shell:
        """
        if [[ {params.extend_reads} == "TRUE" ]]; then
            set -x
            MEM=$({params.loglog_bin} seed={params.seed} k={params.k} in={input.reads} ignorebadquality 2> >(grep Cardinality) \
                | awk -vP={params.table_cap} -vB={params.bits} -vH={params.hashes} '{{print int( (((B*H)/8)*$2)/P )}}')
            echo ${{MEM}}
            if [[ ${{MEM}} -gt {params.maxmem} ]]; then
                MEM={params.maxmem}
            fi

            {params.tadpole_bin} \
                -Xmx{params.memory} \
                k={params.k} \
                in={input.reads} out={output.extended} mode={params.mode} ibb={params.ibb} prefilter={params.prefilter} \
                el={params.length} er={params.length} threads={threads} overwrite=true trimends={params.trimends} \
                ecc={params.ecc} ecco={params.ecco} filtermem=${{MEM}} conservative={params.conservative} ignorebadquality \
                >> {log} 2>&1

            {params.seqkit_bin} stats -j {threads} -T -a --quiet \
                {output.extended} -o {output.stats}
        else
            ln -s {input.reads} {output.extended}
            ln -s {input.stats} {output.stats}
        fi
        """
