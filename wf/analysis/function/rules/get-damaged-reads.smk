rule get_damaged_reads:
    input:
        reads=config["rdir"] + "/read-derep/{smp}.fa.gz",
        bam=(config["rdir"] + "/taxonomic-profiling/{smp}.dedup.filtered.bam"),
        mdmg_csv=(
            config["rdir"]
            + "/taxonomic-profiling-dmg/{smp}.tp-mdmg.local.weight-1.csv.gz"
        ),
    output:
        dmg_reads=config["rdir"] + "/read-dmg/{smp}.damaged.fastq.gz",
        nondmg_reads=config["rdir"] + "/read-dmg/{smp}.nondamaged.fastq.gz",
        dmg_reads_derep=config["rdir"] + "/read-dmg/{smp}.damaged.derep.fasta.gz",
    threads: config["seqkit_threads"]
    params:
        seqkit_bin=config["seqkit_bin"],
        dread_bin=config["dread_bin"],
        dread_filter=config["dread_filter"],
        label="{smp}",
        rdir=config["rdir"] + "/read-dmg",
        wdir=config["wdir"],
    conda:
        "../envs/get-dmg-reads.yaml"
    log:
        config["rdir"] + "/logs/read-dmg/{smp}.read-dmg.log",
    benchmark:
        config["rdir"] + "/benchmarks/read-dmg/{smp}.read-dmg.bmk"
    message:
        """--- Get damaged reads reads"""
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}
        {params.dread_bin} -r {input.mdmg_csv} -b {input.bam} -f '{params.dread_filter}' --prefix {params.label}

        {params.seqkit_bin} grep  -j {threads} -f <(seqkit fx2tab -n {output.dmg_reads}) -o {output.dmg_reads_derep}  <({params.seqkit_bin} replace -p '(\S+)__(\d+)' -r '${{1}}' {input.reads})
        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}

        """
