rule agp_predict:
    """Predict genes from reads using AGP (DART)"""
    input:
        reads=config["rdir"] + "/read-derep/{smp}.fa.gz",
    output:
        faa=config["rdir"] + "/agp-predict/{smp}.faa",
        gff=config["rdir"] + "/agp-predict/{smp}.gff",
        summary=config["rdir"] + "/agp-predict/{smp}.json",
    threads: config.get("agp_threads", 4)
    params:
        agp_bin=config["agp_bin"],
        rdir=config["rdir"] + "/agp-predict",
        wdir=config["wdir"],
        domain=config.get("agp_domain", "gtdb"),
    conda:
        "../envs/agp.yaml"
    log:
        config["rdir"] + "/logs/agp-predict/{smp}.agp-predict.log",
    benchmark:
        config["rdir"] + "/benchmarks/agp-predict/{smp}.agp-predict.bmk"
    message:
        """--- AGP predict: {wildcards.smp}"""
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        {params.agp_bin} predict \
            -i {input.reads} \
            -o {output.gff} \
            --fasta-aa {output.faa} \
            --summary {output.summary} \
            --domain {params.domain} \
            --threads {threads} \
            -v 2>&1 | tee {log}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
