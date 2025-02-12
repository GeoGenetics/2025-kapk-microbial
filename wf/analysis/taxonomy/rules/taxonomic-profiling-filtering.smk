rule taxonomic_profiling_filtering:
    wildcard_constraints:
        indata="prefilter|standard|raw",
    input:
        bowtie2_bam_dedup=config["rdir"]
        + "/taxonomic-profiling-mapping/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{smp}.dedup.bam",
    output:
        bowtie2_bam_dedup_filt=(
            config["rdir"]
            + "/taxonomic-profiling-filtering/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/{smp}.dedup.filtered.bam"
        ),
        tp_stats=config["rdir"]
        + "/taxonomic-profiling-filtering/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/{smp}.dedup_stats.tsv.gz",
        tp_stats_filt=(
            config["rdir"]
            + "/taxonomic-profiling-filtering/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/{smp}.dedup_stats-filtered.tsv.gz"
        ),
    threads: config["tp_threads"]
    params:
        filter_bam_bin=config["filter_bam_bin"],
        bowtie2_bam_dedup=config["rdir"]
        + "/taxonomic-profiling-filtering/{indata}/k{tp_bowtie2_k}/{read_ani}/{smp}.dedup.bam",
        tp_filter_sort_mem=config["tp_filter_sort_mem"],
        tp_filter_nreads=config["tp_filter_nreads"],
        tp_filter_exp_bratio=config["tp_filter_exp_bratio"],
        tp_filter_cov_evenness=config["tp_filter_cov_evenness"],
        tp_filter_read_ani="{read_ani}",
        tp_filter_entropy=config["tp_filter_entropy"],
        tp_filter_gini=config["tp_filter_gini"],
        tp_sizes=config["tp_sizes"],
        tp_norm_scale=config["tp_norm_scale"],
        tp_extra=config["tp_extra"],
        label="{smp}",
        profiling_data_input="{indata}",
        rdir=config["rdir"]
        + "/taxonomic-profiling-filtering/{indata}/k{tp_bowtie2_k}/{read_ani}",
        wdir=config["wdir"],
    log:
        config["rdir"]
        + "/logs/taxonomic-profiling-filtering/{indata}/k{tp_bowtie2_k}/{read_ani}/{smp}.taxonomic-profiling.log",
    benchmark:
        (
            config["rdir"]
            + "/benchmarks/taxonomic-profiling-filtering/{indata}/k{tp_bowtie2_k}/{read_ani}/{smp}.taxonomic-profiling.bmk"
        )
    conda:
        "../envs/bam-filter.yaml"
    message:
        """--- High-res taxonomic profiling."""
    shell:
        """
        set -x
        set +e
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        if [[ -f  {params.bowtie2_bam_dedup}.bai ]]; then
            rm -rf  {params.bowtie2_bam_dedup}.bai
        fi

        # Filter mapping file
        ln -sf {input.bowtie2_bam_dedup} {params.bowtie2_bam_dedup}

        {params.filter_bam_bin} -r {params.tp_sizes} \
            -N \
            -g {params.tp_filter_entropy} \
            -e {params.tp_filter_gini} \
            -m {params.tp_filter_sort_mem} \
            -t {threads} \
            -n {params.tp_filter_nreads} \
            -b {params.tp_filter_exp_bratio} \
            -c {params.tp_filter_cov_evenness} \
            -A {params.tp_filter_read_ani} \
            {params.tp_extra} {params.bowtie2_bam_dedup}
        errorcode=$?
        if [[ ! -s {output.bowtie2_bam_dedup_filt} ]] || [[ $errorcode -eq 1 ]]; then
            echo "No references passed the filters"
            touch {output.bowtie2_bam_dedup_filt}
            touch {output.tp_stats}
            touch {output.tp_stats_filt}
        fi

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
