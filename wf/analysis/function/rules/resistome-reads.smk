rule resistome_reads:
    input:
        reads_derep=config["rdir"] + "/read-derep/{smp}.fa.gz",
        reads=config["rdir"] + "/read-renamed/{smp}.fq.gz",
        stats=config["rdir"] + "/stats/{smp}.stats-initial.txt",
    output:
        bam=config["rdir"] + "/resistome-reads/{smp}.megares.bam",
        stats=config["rdir"] + "/resistome-reads/{smp}.megares_stats.tsv.gz",
        stats_filtered=config["rdir"]
        + "/resistome-reads/{smp}.megares_stats-filtered.tsv.gz",
    threads: config["resistome_reads_threads"]
    params:
        bowtie2_bin=config["bowtie2_bin"],
        seqkit_bin=config["seqkit_bin"],
        samtools_bin=config["samtools_bin"],
        filter_bam_bin=config["filter_bam_bin"],
        megares_db_bw=config["megares_db_bw"],
        bowtie2_params=config["resistome_bowtie2_params"],
        filter_sort_mem=config["filter_sort_mem"],
        resistome_thres=config["resistome_thres"],
        resistome_filter_nreads=config["resistome_filter_nreads"],
        fastq=config["rdir"] + "/resistome-reads/{smp}.read-reps.fq.gz",
    conda:
        "../envs/resistome.yaml"
    log:
        config["rdir"] + "/logs/resistome-reads/{smp}.resistome.log",
    benchmark:
        config["rdir"] + "/benchmarks/resistome-reads/{smp}.resistome.bmk"
    message:
        """--- Resistome analysis"""
    shell:
        """
        # Extract reads
        {params.seqkit_bin} grep \
            -f <({params.seqkit_bin} fx2tab -j {threads} -n {input.reads_derep} | awk '{{split($1, a,"__"); print a[1]}}') \
            -j {threads} \
            -o {params.fastq} {input.reads} 

        IS_FASTQ=$(awk 'NR==2{{print $2}}' {input.stats})

        if [ ${{IS_FASTQ}} = "FASTQ" ]; then
            # Align reads
            {params.bowtie2_bin} -x {params.megares_db_bw} \
                -p {threads} \
                -q -U {params.fastq} \
                | {params.samtools_bin} view -bS \
                | {params.samtools_bin} sort -@ {threads} > {output.bam}
        else
            # Align reads
            {params.bowtie2_bin} -x {params.megares_db_bw}  \
                -p {threads} \
                -f -U {params.fastq} \
                | {params.samtools_bin} view -bS \
                | {params.samtools_bin} sort -@ {threads} >{output.bam}
        fi

        {params.filter_bam_bin} -m {params.filter_sort_mem} \
            -t {threads} \
            -n {params.resistome_filter_nreads} \
            -b 0 \
            -c 0 \
            -B {params.resistome_thres} \
            --only-stats-filtered \
            {output.bam}

        if [ ! -f {output.stats} ]; then
            touch {output.stats}
        fi
        if [ ! -f {output.stats_filtered} ]; then
            touch {output.stats_filtered}
        fi

        """
