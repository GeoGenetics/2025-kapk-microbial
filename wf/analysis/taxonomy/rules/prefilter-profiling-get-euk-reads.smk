rule prefilter_profiling_get_euk_reads:
    input:
        reads=config["rdir"] + "/read-renamed/{smp}.fq.gz",
        reads_derep=config["rdir"] + "/read-derep/{smp}.fa.gz",
        bac_ids=config["rdir"] + "/prefilter-profiling-coarse/{smp}_Bacteria.txt.gz",
        arc_ids=config["rdir"] + "/prefilter-profiling-coarse/{smp}_Archaea.txt.gz",
        vir_ids=config["rdir"] + "/prefilter-profiling-coarse/{smp}_Viruses.txt.gz",
    output:
        fastq=config["rdir"] + "/read-euk/{smp}.read-euk.fq.gz",
    threads: config["prefilter_get_reads_threads"]
    params:
        awk_bin=config["awk_bin"],
        seqkit_bin=config["seqkit_bin"],
        taxonkit_bin=config["taxonkit_bin"],
        read_ids=config["rdir"] + "/read-euk/read-euk-ids.txt",
        label="{smp}",
        rdir=config["rdir"] + "/read-euk",
        wdir=config["wdir"],
    log:
        config["rdir"] + "/logs/read-euk/{smp}.read-euk.log",
    benchmark:
        config["rdir"] + "/benchmarks/read-euk/{smp}.read-euk.bmk"
    conda:
        "../envs/bam-filter.yaml"
    message:
        """--- Get eukaryotic reads."""
    shell:
        """
        set -x
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        rm -rf {params.read_ids}

        if [[ -s {input.bac_ids} ]]; then
            zcat {input.bac_ids} >> {params.read_ids}
        fi

        if [[ -s {input.arc_ids} ]]; then
            zcat {input.arc_ids} >> {params.read_ids}
        fi

        if [[ -s {input.vir_ids} ]]; then
            zcat {input.vir_ids} >> {params.read_ids}
        fi

        # We find those reads that are not classified as arc/bact/vir
        if [[ {params.use_vsearch} == "TRUE" ]]; then
            {params.seqkit_bin} grep \
                -v \
                -f {params.read_ids} \
                -j {threads} {input.reads_derep} \
            | {params.seqkit_bin} fx2tab \
            | sed 's/__[0-9]\+$//' \
            | cut -f1 > {params.read_ids}
        else
            {params.seqkit_bin} grep \
                -v \
                -f {params.read_ids} \
                -j {threads} {input.reads_derep} \
            | {params.seqkit_bin} fx2tab -n > {params.read_ids}
        fi

        {params.seqkit_bin} grep \
            -f {params.read_ids} \
            -j {threads} \
            -o {output.fastq} {input.reads}    

        rm -rf {params.read_ids}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
