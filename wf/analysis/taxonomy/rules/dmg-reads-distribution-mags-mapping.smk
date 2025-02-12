def get_dr_db(db):
    for entry in config["db_mag_distribution"]:
        for key, value in entry.items():
            if key == db:
                return value["bowtie2_db"]


rule reads_dmg_distribution_mapping:
    wildcard_constraints:
        indata="prefilter|standard|raw",
        weight="0|1",
        dr_fraction="damaged|non-damaged",
    input:
        reads_dmg=config["rdir"]
        + "/read-dmg-rank/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/weight-{weight}"
        + "/{smp}.{dr_fraction}.fastq.gz",
    output:
        bowtie2_bam_dedup=config["rdir"]
        + "/reads-dmg-distribution-mags-mapping/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/weight-{weight}"
        + "/{mag_db}"
        + "/{dr_fraction}"
        + "/{smp}.dmg-reads.dedup.bam",
        bowtie2_bam_dedup_metrics=(
            config["rdir"]
            + "/reads-dmg-distribution-mags-mapping/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{mag_db}"
            + "/{dr_fraction}"
            + "/{smp}.dmg-reads.dedup.metrics"
        ),
    threads: config["bowtie2_dr_threads"]
    params:
        awk_bin=config["awk_bin"],
        seqkit_bin=config["seqkit_bin"],
        bowtie2_bin=config["bowtie2_bin"],
        use_vsearch=config["use_vsearch"],
        fastq=config["tp_bowtie2_tmp"]
        + "/reads-dmg-distribution-mags-mapping/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/weight-{weight}"
        + "/{mag_db}"
        + "/{dr_fraction}"
        + "/{smp}.dmg-reads.fq.gz",
        bowtie2_rd_db=lambda wc: get_dr_db(db=wc.mag_db),
        tp_bowtie2_parms=config["tp_bowtie2_parms"],
        tp_bowtie2_k="{tp_bowtie2_k}",
        samtools_bin=config["samtools_bin"],
        samtools_view_parms=config["samtools_view_parms"],
        samtools_sort_parms=config["samtools_sort_parms"],
        picard_jar=config["picard_jar"],
        picard_max_records_in_ram=config["picard_max_records_in_ram"],
        java_opts=config["picard_java_opts"],
        rdir=config["rdir"]
        + "/reads-dmg-distribution-mags-mapping/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/weight-{weight}"
        + "/{mag_db}"
        + "/{dr_fraction}",
        wdir=config["wdir"],
        bowtie2_tmp_dir=config["tp_bowtie2_tmp"]
        + "/reads-dmg-distribution-mags-mapping/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/weight-{weight}"
        + "/{mag_db}"
        + "/{dr_fraction}"
        + "/{smp}",
        bowtie2_sam_tmp=config["tp_bowtie2_tmp"]
        + "/reads-dmg-distribution-mags-mapping/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/weight-{weight}"
        + "/{mag_db}"
        + "/{dr_fraction}"
        + "/{smp}"
        + "/{smp}.dmg-reads.tmp.sam",
        bowtie2_bam=config["rdir"]
        + "/reads-dmg-distribution-mags-mapping/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/weight-{weight}"
        + "/{mag_db}"
        + "/{dr_fraction}"
        + "/{smp}.dmg-reads.bam",
        bowtie2_bam_dedup_tmp=config["rdir"]
        + "/reads-dmg-distribution-mags-mapping/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/weight-{weight}"
        + "/{mag_db}"
        + "/{dr_fraction}"
        + "/{smp}.dmg-reads.dedup.tmp.bam",
        profiling_data_input="{indata}",
        remove_bam=config["remove_bam"],
        remove_fastq=config["remove_fastq"],
    log:
        config["rdir"]
        + "/logs/reads-dmg-distribution-mags-mapping/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/weight-{weight}"
        + "/{mag_db}"
        + "/{dr_fraction}"
        + "/{smp}.dmg-reads-mapping.log",
    benchmark:
        (
            config["rdir"]
            + "/benchmarks/reads-dmg-distribution-mags-mapping/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{mag_db}"
            + "/{dr_fraction}"
            + "/{smp}.dmg-reads-mapping.bmk"
        )
    conda:
        "../envs/bam-filter.yaml"
    message:
        """--- High-res taxonomic profiling."""
    shell:
        """
        set -x
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        # create

        if [[ ! -s {input.reads_dmg} ]]; then
            echo "No reads to map"
            touch {output.bowtie2_bam_dedup}
            touch {output.bowtie2_bam_dedup_metrics}
            exit 0
        fi


        if [[ -f {params.fastq} ]]; then
            rm -rf {params.fastq}
        fi

        # check if tmp dir exists, if not remove it
        if [[ -d {params.bowtie2_tmp_dir} ]]; then
            rm -rf {params.bowtie2_tmp_dir}
            mkdir -p {params.bowtie2_tmp_dir}
        else
            mkdir -p {params.bowtie2_tmp_dir}
        fi

        IS_FASTQ=$({params.seqkit_bin} head -n 1000 {input.reads_dmg} | {params.seqkit_bin} stats -T | awk 'NR==2{{print $2}}')

        # find out if bowtie2 is running
        # Very naive, if there's a process already running
        # with the DB needed use memory mapping instead of loading the DB

        ps -Af | grep bowtie2 | grep "{params.bowtie2_rd_db}" > /dev/null
        if [ $? -eq 0 ]; then
            echo "Bowtie2 is running."
            MM="--mm"
        else
         echo "Bowtie2 is not running."
            MM=""
        fi

        if [ ${{IS_FASTQ}} = "FASTQ" ]; then
            # Align reads
            {params.bowtie2_bin} -x {params.bowtie2_rd_db} \
                -t \
                -p {threads} \
                -q -U {input.reads_dmg} \
                -k {params.tp_bowtie2_k} \
                {params.tp_bowtie2_parms} ${{MM}} > {params.bowtie2_sam_tmp} 2> {log}
        else
            # Align reads
            {params.bowtie2_bin} -x {params.bowtie2_rd_db} \
                -t \
                -p {threads} \
                -f -U {input.reads_dmg} \
                -k {params.tp_bowtie2_k} \
                {params.tp_bowtie2_parms} ${{MM}} > {params.bowtie2_sam_tmp} 2> {log}
        fi

        # Sorting decoupled from mapping to save memory so we can run multiple samples in parallel
        {params.samtools_bin} sort -@ {threads} -T {params.bowtie2_tmp_dir} {params.samtools_sort_parms} -O BAM -o {params.bowtie2_bam} {params.bowtie2_sam_tmp} 

        rm -rf {params.bowtie2_sam_tmp}

        java {params.java_opts} -jar {params.picard_jar} MarkDuplicates \
            -INPUT {params.bowtie2_bam} \
            -OUTPUT {params.bowtie2_bam_dedup_tmp} \
            -METRICS_FILE {output.bowtie2_bam_dedup_metrics} \
            -ASO coordinate \
            -VALIDATION_STRINGENCY LENIENT \
            -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 \
            -MAX_RECORDS_IN_RAM {params.picard_max_records_in_ram} \
            -REMOVE_DUPLICATES TRUE  \
            >> {log} 2>&1

        {params.samtools_bin} sort -@ {threads} -T {params.bowtie2_tmp_dir} {params.samtools_sort_parms} {params.bowtie2_bam_dedup_tmp} > {output.bowtie2_bam_dedup}
        #{params.samtools_bin} index -@ {threads} {output.bowtie2_bam_dedup}

        if [[ {params.remove_bam} ==  TRUE ]]; then
            rm -rf {params.bowtie2_bam}
        fi
        if [[ {params.remove_fastq} == TRUE ]]; then
            rm -rf {params.fastq}
        else
            mv {params.fastq} {params.rdir}
        fi

        rm -rf {params.bowtie2_bam} {params.bowtie2_bam_dedup_tmp} {params.bowtie2_tmp_dir}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
