def get_mapping_input(wc):
    infile = wc.indata
    if infile == "prefilter":
        file = config["rdir"] + "/read-noneuk-nuclear/{smp}.read-noneuk-nuclear.fq.gz"
    elif infile == "standard":
        file = config["rdir"] + "/read-derep/{smp}.fa.gz"
    elif infile == "raw":
        file = config["rdir"] + "/read-renamed/{smp}.fq.gz"
    return file


def get_other_db(db):
    for entry in config["db_other"]:
        for key, value in entry.items():
            if key == db:
                return value["bowtie2_db"]


rule other_db_mapping:
    wildcard_constraints:
        indata="prefilter|standard|raw",
    input:
        reads_raw=config["rdir"] + "/read-renamed/{smp}.fq.gz",
        reads=lambda wc: get_mapping_input(wc),
    output:
        bowtie2_bam_dedup=config["rdir"]
        + "/other-mapping/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{other_db}"
        + "/{smp}.dedup.bam",
        bowtie2_bam_dedup_metrics=(
            config["rdir"]
            + "/other-mapping/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{other_db}"
            + "/{smp}.dedup.metrics"
        ),
    threads: config["bowtie2_tp_threads"]
    params:
        awk_bin=config["awk_bin"],
        seqkit_bin=config["seqkit_bin"],
        bowtie2_bin=config["bowtie2_bin"],
        use_vsearch=config["use_vsearch"],
        fastq=config["tp_bowtie2_tmp"]
        + "/other-mapping/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{other_db}"
        + "/{smp}"
        + "/{smp}.read-reps.fq.gz",
        bowtie2_other_db=lambda wc: get_other_db(db=wc.other_db),
        tp_bowtie2_parms=config["tp_bowtie2_parms"],
        tp_bowtie2_k="{tp_bowtie2_k}",
        samtools_bin=config["samtools_bin"],
        samtools_view_parms=config["samtools_view_parms"],
        samtools_sort_parms=config["samtools_sort_parms"],
        picard_jar=config["picard_jar"],
        picard_max_records_in_ram=config["picard_max_records_in_ram"],
        java_opts=config["picard_java_opts"],
        label="{smp}",
        rdir=config["rdir"]
        + "/other-mapping/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{other_db}",
        wdir=config["wdir"],
        bowtie2_tmp_dir=config["tp_bowtie2_tmp"]
        + "/other-mapping/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{other_db}"
        + "/{smp}",
        bowtie2_sam_tmp=config["tp_bowtie2_tmp"]
        + "/other-mapping/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{other_db}"
        + "/{smp}"
        + "/{smp}.tmp.sam",
        bowtie2_bam=config["rdir"]
        + "/other-mapping/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{other_db}"
        + "/{smp}.bam",
        bowtie2_bam_dedup_tmp=config["rdir"]
        + "/other-mapping/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{other_db}"
        + "/{smp}.dedup.tmp.bam",
        derep_reads=lambda wildcards: sample_table_read.derep_reads[wildcards.smp],
        profiling_data_input="{indata}",
        remove_bam=config["remove_bam"],
        remove_fastq=config["remove_fastq"],
    log:
        config["rdir"]
        + "/logs/other-mapping/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{other_db}"
        + "/{smp}.other-mapping.log",
    benchmark:
        (
            config["rdir"]
            + "/benchmarks/other-mapping/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{other_db}"
            + "/{smp}.other-mapping.bmk"
        )
    conda:
        "../envs/bam-filter.yaml"
    message:
        """--- Other DB mapping"""
    shell:
        """
        set -x
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        # create

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

        if [[ {params.profiling_data_input} == "raw" ]]; then
            ln -sf {input.reads_raw} {params.fastq}
        elif [[ {params.profiling_data_input} == "prefilter" ]]; then
            ln -sf {input.reads} {params.fastq}
        else
            if [[ {params.use_vsearch} == "TRUE" ]]; then
                {params.seqkit_bin} grep \
                    -f <({params.seqkit_bin} fx2tab -j {threads} -n {input.reads} | sed 's/__[0-9]\+$//') \
                    -j {threads} \
                    -o {params.fastq} {input.reads_raw}
            else
                {params.seqkit_bin} grep \
                    -f <({params.seqkit_bin} fx2tab -j {threads} -n {input.reads}) \
                    -j {threads} \
                    -o {params.fastq} {input.reads_raw}
            fi
        fi

        IS_FASTQ=$({params.seqkit_bin} head -n 1000 {params.fastq} | {params.seqkit_bin} stats -T | awk 'NR==2{{print $2}}')

        # find out if bowtie2 is running
        # Very naive, if there's a process already running
        # with the DB needed use memory mapping instead of loading the DB

        ps -Af | grep bowtie2 | grep "{params.bowtie2_other_db}" > /dev/null
        if [ $? -eq 0 ]; then
            echo "Bowtie2 is running."
            MM="--mm"
        else
         echo "Bowtie2 is not running."
            MM=""
        fi

        if [ ${{IS_FASTQ}} = "FASTQ" ]; then
            # Align reads
            {params.bowtie2_bin} -x {params.bowtie2_other_db} \
                -t \
                -p {threads} \
                -q -U {params.fastq} \
                -k {params.tp_bowtie2_k} \
                {params.tp_bowtie2_parms} ${{MM}} > {params.bowtie2_sam_tmp} 2> {log}
        else
            # Align reads
            {params.bowtie2_bin} -x {params.bowtie2_other_db} \
                -t \
                -p {threads} \
                -f -U {params.fastq} \
                -k {params.tp_bowtie2_k} \
                {params.tp_bowtie2_parms} ${{MM}} > samtools view -Sb > {params.bowtie2_sam_tmp} 2> {log}
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
