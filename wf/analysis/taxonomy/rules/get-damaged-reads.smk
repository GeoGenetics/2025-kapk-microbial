rule get_damaged_reads:
    wildcard_constraints:
        indata="prefilter|standard|raw",
    input:
        reads=config["rdir"] + "/read-derep/{smp}.fa.gz",
        bowtie2_bam_dedup_filt=(
            config["rdir"]
            + "/taxonomic-profiling-filtering/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/{smp}.dedup.filtered.bam"
        ),
        mdmg_csv=(
            config["rdir"]
            + "/taxonomic-profiling-dmg/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/local"
            + "/{smp}.tp-mdmg.weight-{weight}.csv.gz"
        ),
    output:
        dmg_reads=(
            config["rdir"]
            + "/read-dmg/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{smp}.damaged.fastq.gz"
        ),
        nondmg_reads=(
            config["rdir"]
            + "/read-dmg/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{smp}.nondamaged.fastq.gz"
        ),
        dmg_reads_derep=(
            config["rdir"]
            + "/read-dmg/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{smp}.damaged.derep.fasta.gz"
        ),
    threads: config["seqkit_threads"]
    params:
        seqkit_bin=config["seqkit_bin"],
        dread_bin=config["dread_bin"],
        dread_filter=config["dread_filter"],
        label="{smp}",
        bowtie2_bam_dedup_filt=config["rdir"]
        + "/read-dmg/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/weight-{weight}"
        + "/{smp}.dedup.filtered.bam",
        rdir=config["rdir"]
        + "/read-dmg/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/weight-{weight}",
        wdir=config["wdir"],
        use_vsearch=config["use_vsearch"],
        dread_sort_memory=config["dread_sort_memory"],
    conda:
        "../envs/get-dmg-reads.yaml"
    log:
        config["rdir"]
        + "/logs"
        + "/read-dmg/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/weight-{weight}"
        + "/{smp}.read-dmg.log",
    benchmark:
        config["rdir"]
        +"/logs"
        +"/read-dmg/"
        +"{indata}"
        +"/k{tp_bowtie2_k}"
        +"/{read_ani}"
        +"/weight-{weight}"
        +"/{smp}.read-dmg.bmk"
    message:
        """--- Get damaged reads reads"""
    shell:
        """
        set -x 
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        # Get damaged and non-damaged reads
        ln -sf {input.bowtie2_bam_dedup_filt} {params.bowtie2_bam_dedup_filt}        
        # Check that both needed files exist if any is missing touch output files
        if [[ -s {input.mdmg_csv} && -s {input.bowtie2_bam_dedup_filt} ]]; then
            {params.dread_bin} \
                -m {input.mdmg_csv} \
                -b {params.bowtie2_bam_dedup_filt} \
                -f '{params.dread_filter}' \
                --prefix {params.label} \
                --threads {threads} \
                --sort-memory {params.dread_sort_memory}

                # If we get reads back
                # Get the reads that are damaged from the derep set
                if [[ -s {output.dmg_reads} ]]; then
                    # Did we use vsearch?
                    if [[ {params.use_vsearch} == "TRUE" ]]; then
                        {params.seqkit_bin} grep \
                            -j {threads} \
                            -f <(seqkit fx2tab -n {output.dmg_reads}) \
                            -o {output.dmg_reads_derep}  \
                            <({params.seqkit_bin} replace -p '(\S+)__(\d+)' -r '${{1}}' {input.reads})
                    else
                        {params.seqkit_bin} grep \
                            -j {threads} \
                            -f <(seqkit fx2tab -n {output.dmg_reads}) \
                            -o {output.dmg_reads_derep} {input.reads}
                    fi
                else
                    touch {output.dmg_reads}
                    touch {output.dmg_reads_derep}
                fi

                # Check if we didn't get any non-damaged
                if [[ ! -s {output.nondmg_reads} ]]; then
                    touch {output.nondmg_reads}
                fi
        else
            touch {output.dmg_reads}
            touch {output.nondmg_reads}
            touch {output.dmg_reads_derep}
        fi

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
