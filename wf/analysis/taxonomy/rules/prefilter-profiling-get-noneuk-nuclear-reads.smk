rule prefilter_profiling_get_noneuk_nuclear_reads:
    input:
        reads=config["rdir"] + "/read-renamed/{smp}.fq.gz",
        reads_derep=config["rdir"] + "/read-derep/{smp}.fa.gz",
        bac_ids=config["rdir"] + "/prefilter-profiling-coarse/{smp}_Bacteria.txt.gz",
        arc_ids=config["rdir"] + "/prefilter-profiling-coarse/{smp}_Archaea.txt.gz",
        vir_ids=config["rdir"] + "/prefilter-profiling-coarse/{smp}_Viruses.txt.gz",
        root_ids=config["rdir"] + "/prefilter-profiling-coarse/{smp}_root.txt.gz",
        euk_ids=config["rdir"]
        + "/prefilter-profiling-organelles/{smp}_Eukaryota.txt.gz",
        unclassified_ids=config["rdir"]
        + "/prefilter-profiling-coarse/{smp}_unclassified.txt.gz",
    output:
        fastq=config["rdir"] + "/read-noneuk-nuclear/{smp}.read-noneuk-nuclear.fq.gz",
    threads: config["prefilter_get_reads_threads"]
    params:
        awk_bin=config["awk_bin"],
        use_vsearch=config["use_vsearch"],
        seqkit_bin=config["seqkit_bin"],
        taxonkit_bin=config["taxonkit_bin"],
        read_ids=config["rdir"]
        + "/read-noneuk-nuclear/{smp}-read-noneuk-nuclear-ids.txt",
        include_noneuk_nuclear_unclassified=config[
            "include_noneuk_nuclear_unclassified"
        ],
        fastq_derep=config["rdir"]
        + "/read-noneuk-nuclear/{smp}.read-noneuk-nuclear-derep.fq.gz",
        label="{smp}",
        rdir=config["rdir"] + "/read-noneuk-nuclear",
        wdir=config["wdir"],
    log:
        config["rdir"] + "/logs/read-noneuk-nuclear/{smp}.read-euk-noneuk-nuclear.log",
    benchmark:
        config["rdir"] + "/benchmarks/read-noneuk-nuclear/{smp}.read-noneuk-nuclear.bmk"
    conda:
        "../envs/bam-filter.yaml"
    message:
        """--- Get all non-eukaryotic nuclear reads."""
    shell:
        """
        set -x
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        # Extract reads
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

        if [[ -s {input.euk_ids} ]]; then
            zcat {input.euk_ids} >> {params.read_ids}
        fi

        if [[ {params.include_noneuk_nuclear_unclassified} == "TRUE" ]]; then
            if [[ -s {input.unclassified_ids} ]]; then
                zcat {input.unclassified_ids} >> {params.read_ids}
            elif [[ -s {input.root_ids} ]]; then
                zcat {input.root_ids} >> {params.read_ids}
            fi
        fi

        # get original representative reads from the derep_reads
        if [[ {params.use_vsearch} == "TRUE" ]]; then
            {params.seqkit_bin} grep \
                -f <(sed 's/__[0-9]\+$//' {params.read_ids} | sort -u --parallel {threads} -S20%) \
                -j {threads} \
                -o {output.fastq} {input.reads}
        else
            {params.seqkit_bin} grep \
                -f <(sort -u --parallel {threads} -S20% {params.read_ids}) \
                -j {threads} \
                -o {output.fastq} {input.reads}
        fi

        # Clean up
        #rm -rf {params.read_ids} {params.fastq_derep}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
