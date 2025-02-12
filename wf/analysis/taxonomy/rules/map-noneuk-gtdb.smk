rule map_noneuk_gtdb:
    input:
        read_ids=config["rdir"] + "/prefilter-profiling-hires/{smp}.prefilter-hires.ids",
        reads=config["rdir"] + "/read-renamed/{smp}.fq.gz",
        stats=config["rdir"] + "/stats/{smp}.stats-initial.txt",
    output:
        fastq=config["rdir"] + "/read-noneuk/{smp}.read-noneuk.fq.gz",
        out_none=config["rdir"]
        + "/taxonomic-profiling/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{smp}.woltka-none.tsv.gz",
        out_free=config["rdir"]
        + "/taxonomic-profiling/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{smp}.woltka-free.tsv.gz",
        out_species=(
            config["rdir"]
            + "/taxonomic-profiling/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{smp}.woltka-species.tsv.gz"
        ),
        out_genus=config["rdir"]
        + "/taxonomic-profiling/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{smp}.woltka-genus.tsv.gz",
        out_family=config["rdir"]
        + "/taxonomic-profiling/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{smp}.woltka-family.tsv.gz",
        out_order=config["rdir"]
        + "/taxonomic-profiling/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{smp}.woltka-order.tsv.gz",
        out_class=config["rdir"]
        + "/taxonomic-profiling/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{smp}.woltka-class.tsv.gz",
        out_phylum=config["rdir"]
        + "/taxonomic-profiling/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{smp}.woltka-phylum.tsv.gz",
        bowtie2_bam_dedup=(
            config["rdir"]
            + "/taxonomic-profiling/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{smp}.woltka.dedup.bam"
        ),
        bowtie2_bam_dedup_filt=(
            config["rdir"]
            + "/taxonomic-profiling/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{smp}.woltka.dedup.filtered.bam"
        ),
        bowtie2_bam_dedup_metrics=(
            config["rdir"]
            + "/taxonomic-profiling/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{smp}.woltka.dedup.metrics"
        ),
        woltka_stats=(
            config["rdir"]
            + "/taxonomic-profiling/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{smp}.woltka.dedup_stats.tsv.gz"
        ),
        woltka_stats_filt=(
            config["rdir"]
            + "/taxonomic-profiling/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{smp}.woltka.dedup_stats-filtered.tsv.gz"
        ),
    threads: config["woltka_threads"]
    params:
        woltka_bin=config["woltka_bin"],
        awk_bin=config["awk_bin"],
        seqkit_bin=config["seqkit_bin"],
        bowtie2_bin=config["bowtie2_bin"],
        filter_bam_bin=config["filter_bam_bin"],
        bowtie2_gtdb_db=config["bowtie2_gtdb_db"],
        woltka_bowtie2_parms=config["woltka_bowtie2_parms"],
        samtools_bin=config["samtools_bin"],
        samtools_view_parms=config["samtools_view_parms"],
        samtools_sort_parms=config["samtools_sort_parms"],
        picard_bin=config["picard_bin"],
        java_opts=config["picard_java_opts"],
        woltka_filter_sort_mem=config["woltka_filter_sort_mem"],
        woltka_filter_nreads=config["woltka_filter_nreads"],
        woltka_filter_exp_bratio=config["woltka_filter_exp_bratio"],
        woltka_filter_cov_evenness=config["woltka_filter_cov_evenness"],
        woltka_taxids=config["woltka_taxids"],
        woltka_nodes=config["woltka_nodes"],
        woltka_names=config["woltka_names"],
        woltka_sizes=config["woltka_sizes"],
        woltka_norm_scale=config["woltka_norm_scale"],
        label="{smp}",
        rdir=config["rdir"] + "/woltka-profiling",
        wdir=config["wdir"],
        bowtie2_bam=config["rdir"]
        + "/taxonomic-profiling/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{smp}.woltka.bam",
        bowtie2_bam_dedup_tmp=(
            config["rdir"]
            + "/taxonomic-profiling/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{smp}.woltka.tmp.bam"
        ),
        bowtie2_sam=(
            config["rdir"]
            + "/taxonomic-profiling/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{smp}.woltka.dedup.filtered.sam"
        ),
    log:
        config["rdir"] + "/logs/taxonomic-profiling/{smp}.woltka-profiling.log",
    benchmark:
        config["rdir"] + "/benchmarks/taxonomic-profiling/{smp}.woltka-profiling.bmk"
    conda:
        "../envs/woltka.yaml"
    message:
        """--- High-res taxonomic profiling with WOLTKA."""
    shell:
        """
        set -x
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        # Extract reads
        {params.seqkit_bin} grep \
            -f <(awk '{{split($1, a,"__"); print a[1]}}' {input.read_ids}) \
            -j {threads} \
            -o {output.fastq} {input.reads} 

        IS_FASTQ=$(awk 'NR==2{{print $2}}' {input.stats})

        if [ ${{IS_FASTQ}} = "FASTQ" ]; then
            # Align reads
            {params.bowtie2_bin} -x {params.bowtie2_gtdb_db} \
                -p {threads} \
                -q -U {output.fastq} \
                {params.woltka_bowtie2_parms} \
                | {params.samtools_bin} view {params.samtools_view_parms} \
                | {params.samtools_bin} sort -@ {threads} {params.samtools_sort_parms} > {params.bowtie2_bam}
        else
            # Align reads
            {params.bowtie2_bin} -x {params.bowtie2_gtdb_db} \
                -p {threads} \
                -f -U {output.fastq} \
                {params.woltka_bowtie2_parms} \
                | {params.samtools_bin} view {params.samtools_view_parms} \
                | {params.samtools_bin} sort -@ {threads} {params.samtools_sort_parms} > {params.bowtie2_bam}
        fi

        {params.picard_bin} {params.java_opts} MarkDuplicates \
            -INPUT {params.bowtie2_bam} \
            -OUTPUT {params.bowtie2_bam_dedup_tmp} \
            -METRICS_FILE {output.bowtie2_bam_dedup_metrics} \
            -AS TRUE \
            -VALIDATION_STRINGENCY LENIENT \
            -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 \
            -REMOVE_DUPLICATES TRUE  \
             >> {log} 2>&1

        {params.samtools_bin} sort -@ {threads} {params.samtools_sort_parms} {params.bowtie2_bam_dedup_tmp} > {output.bowtie2_bam_dedup}
        {params.samtools_bin} index {output.bowtie2_bam_dedup}

        touch -c {output.bowtie2_bam_dedup}.bai

        # Filter woltka mapping file
        {params.filter_bam_bin} -r {params.woltka_sizes} \
            -m {params.woltka_filter_sort_mem} \
            -t {threads} \
            -n {params.woltka_filter_nreads} \
            -b {params.woltka_filter_exp_bratio} \
            -c {params.woltka_filter_cov_evenness} \
            {output.bowtie2_bam_dedup}

        if [[ -f {output.bowtie2_bam_dedup_filt} ]]; then
            {params.samtools_bin} view -o {params.bowtie2_sam} {output.bowtie2_bam_dedup_filt}

            {params.woltka_bin} classify \
                -i {params.bowtie2_sam} \
                --map {params.woltka_taxids} \
                --nodes {params.woltka_nodes} \
                --names {params.woltka_names} \
                --rank none,free,species,genus,family,order,class,phylum \
                --sizes {params.woltka_sizes} \
                --scale {params.woltka_norm_scale} \
                --add-rank \
                --add-lineage \
                --output {params.label}-woltka \
                --outmap {params.label}-woltka_map \
                --to-tsv

            cat {params.label}-woltka/none.tsv | gzip > {output.out_none}
            cat {params.label}-woltka/free.tsv | gzip > {output.out_free}
            cat {params.label}-woltka/species.tsv | gzip > {output.out_species}
            cat {params.label}-woltka/genus.tsv | gzip > {output.out_genus}
            cat {params.label}-woltka/family.tsv | gzip > {output.out_family}
            cat {params.label}-woltka/order.tsv | gzip > {output.out_order}
            cat {params.label}-woltka/class.tsv | gzip > {output.out_class}
            cat {params.label}-woltka/phylum.tsv | gzip > {output.out_phylum}
        else
            touch {output.bowtie2_bam_dedup_filt}
            touch {output.out_none}
            touch {output.out_free}
            touch {output.out_species}
            touch {output.out_genus}
            touch {output.out_family}
            touch {output.out_order}
            touch {output.out_class}
            touch {output.out_phylum}
            touch {output.woltka_stats}
            touch {output.woltka_stats_filt}
        fi

        rm -rf {params.label}-woltka {params.bowtie2_sam} {params.bowtie2_bam_dedup_tmp} {params.bowtie2_bam}
        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
