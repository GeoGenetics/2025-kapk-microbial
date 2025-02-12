rule get_damaged_reads_rank:
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
        fb_data=(
            config["rdir"]
            + "/taxonomic-profiling-filtering/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/{smp}.dedup_stats-filtered.tsv.gz"
        ),
    output:
        dmg_reads_bact=(
            config["rdir"]
            + "/read-dmg-rank/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{smp}.d__Bacteria.damaged.fastq.gz"
        ),
        nondmg_reads_bact=(
            config["rdir"]
            + "/read-dmg-rank/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{smp}.d__Bacteria.non-damaged.fastq.gz"
        ),
        multi_reads_bact=(
            config["rdir"]
            + "/read-dmg-rank/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{smp}.d__Bacteria.multi.fastq.gz"
        ),
        dmg_reads_arc=(
            config["rdir"]
            + "/read-dmg-rank/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{smp}.d__Archaea.damaged.fastq.gz"
        ),
        nondmg_reads_arc=(
            config["rdir"]
            + "/read-dmg-rank/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{smp}.d__Archaea.non-damaged.fastq.gz"
        ),
        multi_reads_arc=(
            config["rdir"]
            + "/read-dmg-rank/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{smp}.d__Archaea.multi.fastq.gz"
        ),
        dmg_reads_vir=(
            config["rdir"]
            + "/read-dmg-rank/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{smp}.d__Viruses.damaged.fastq.gz"
        ),
        nondmg_reads_vir=(
            config["rdir"]
            + "/read-dmg-rank/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{smp}.d__Viruses.non-damaged.fastq.gz"
        ),
        multi_reads_vir=(
            config["rdir"]
            + "/read-dmg-rank/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{smp}.d__Viruses.multi.fastq.gz"
        ),
        dmg_reads_euk=(
            config["rdir"]
            + "/read-dmg-rank/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{smp}.d__Eukaryota.damaged.fastq.gz"
        ),
        nondmg_reads_euk=(
            config["rdir"]
            + "/read-dmg-rank/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{smp}.d__Eukaryota.non-damaged.fastq.gz"
        ),
        multi_reads_euk=(
            config["rdir"]
            + "/read-dmg-rank/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{smp}.d__Eukaryota.multi.fastq.gz"
        ),
        discarded=(
            config["rdir"]
            + "/read-dmg-rank/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{smp}.discarded.fastq.gz"
        ),
        dmg_reads=(
            config["rdir"]
            + "/read-dmg-rank/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{smp}.damaged.fastq.gz"
        ),
        dmg_reads_derep=(
            config["rdir"]
            + "/read-dmg-rank/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{smp}.damaged.derep.fasta.gz"
        ),
        nondmg_reads=(
            config["rdir"]
            + "/read-dmg-rank/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{smp}.non-damaged.fastq.gz"
        ),
        nondmg_reads_derep=(
            config["rdir"]
            + "/read-dmg-rank/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{smp}.non-damaged.derep.fasta.gz"
        ),
    threads: config["seqkit_threads"]
    params:
        seqkit_bin=config["seqkit_bin"],
        dread_bin=config["dread_bin"],
        dread_filter_mdmg=config["dread_filter_mdmg"],
        dread_rank=config["dread_rank"],
        dread_tax=config["dread_tax"],
        dread_filter_fb=config["dread_filter_fb"],
        label="{smp}",
        bowtie2_bam_dedup_filt=config["rdir"]
        + "/read-dmg-rank/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/weight-{weight}"
        + "/{smp}.dedup.filtered.bam",
        rdir=config["rdir"]
        + "/read-dmg-rank/"
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
        + "/read-dmg-rank/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/weight-{weight}"
        + "/{smp}.read-dmg-rank.log",
    benchmark:
        config["rdir"]
        +"/logs"
        +"/read-dmg-rank/"
        +"{indata}"
        +"/k{tp_bowtie2_k}"
        +"/{read_ani}"
        +"/weight-{weight}"
        +"/{smp}.read-dmg-rank.bmk"
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
                -f '{params.dread_filter_mdmg}' \
                --fb-data {input.fb_data} \
                --fb-filter '{params.dread_filter_fb}' \
                --prefix {params.label} \
                --threads {threads} \
                --sort-memory {params.dread_sort_memory} \
                --rank '{params.dread_rank}' \
                --taxonomy-file {params.dread_tax}


                # Check we get all ouput files
                if [[ ! -s {output.dmg_reads_bact} ]]; then
                    touch {output.dmg_reads_bact}
                fi
                if [[ ! -s {output.nondmg_reads_bact} ]]; then
                    touch {output.nondmg_reads_bact}
                fi
                if [[ ! -s {output.multi_reads_bact} ]]; then
                    touch {output.multi_reads_bact}
                fi

                if [[ ! -s {output.dmg_reads_vir} ]]; then
                    touch {output.dmg_reads_vir}
                fi
                if [[ ! -s {output.nondmg_reads_vir} ]]; then
                    touch {output.nondmg_reads_vir}
                fi
                if [[ ! -s {output.multi_reads_vir} ]]; then
                    touch {output.multi_reads_vir}
                fi

                if [[ ! -s {output.dmg_reads_arc} ]]; then
                    touch {output.dmg_reads_arc}
                fi
                if [[ ! -s {output.nondmg_reads_arc} ]]; then
                    touch {output.nondmg_reads_arc}
                fi
                if [[ ! -s {output.multi_reads_arc} ]]; then
                    touch {output.multi_reads_arc}
                fi

                if [[ ! -s {output.dmg_reads_euk} ]]; then
                    touch {output.dmg_reads_euk}
                fi
                if [[ ! -s {output.nondmg_reads_euk} ]]; then
                    touch {output.nondmg_reads_euk}
                fi
                if [[ ! -s {output.multi_reads_euk} ]]; then
                    touch {output.multi_reads_euk}
                fi

                if [[ ! -s {output.discarded} ]]; then
                    touch {output.discarded}
                fi

                # combine all reads
                cat {output.dmg_reads_bact} {output.dmg_reads_vir} {output.dmg_reads_arc} > {output.dmg_reads}
                cat {output.nondmg_reads_bact} {output.nondmg_reads_vir} {output.nondmg_reads_arc} > {output.nondmg_reads}

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

                if [[ -s {output.nondmg_reads} ]]; then
                    # Did we use vsearch?
                    if [[ {params.use_vsearch} == "TRUE" ]]; then
                        {params.seqkit_bin} grep \
                            -j {threads} \
                            -f <(seqkit fx2tab -n {output.nondmg_reads}) \
                            -o {output.nondmg_reads_derep}  \
                            <({params.seqkit_bin} replace -p '(\S+)__(\d+)' -r '${{1}}' {input.reads})
                    else
                        {params.seqkit_bin} grep \
                            -j {threads} \
                            -f <(seqkit fx2tab -n {output.nondmg_reads}) \
                            -o {output.nondmg_reads_derep} {input.reads}
                    fi
                else
                    touch {output.nondmg_reads}
                    touch {output.nondmg_reads_derep}
                fi
        else
            touch {output.dmg_reads_bact}
            touch {output.nondmg_reads_bact}
            touch {output.multi_reads_bact}

            touch {output.dmg_reads_vir}
            touch {output.nondmg_reads_vir}
            touch {output.multi_reads_vir}

            touch {output.dmg_reads_arc}
            touch {output.nondmg_reads_arc}
            touch {output.multi_reads_arc}

            touch {output.dmg_reads_euk}
            touch {output.nondmg_reads_euk}
            touch {output.multi_reads_euk}

            touch {output.dmg_reads}
            touch {output.dmg_reads_derep}
            touch {output.nondmg_reads}
            touch {output.nondmg_reads_derep}

            touch {output.discarded}
        fi

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
