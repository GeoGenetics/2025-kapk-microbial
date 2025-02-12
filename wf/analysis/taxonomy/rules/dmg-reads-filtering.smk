rule reads_dmg_filtering:
    wildcard_constraints:
        indata="prefilter|standard|raw",
        weight="0|1",
        dr_fraction="damaged|non-damaged",
    input:
        bowtie2_bam_dedup=config["rdir"]
        + "/reads-dmg-mapping/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/weight-{weight}"
        + "/{dr_fraction}"
        + "/{smp}.dmg-reads.dedup.bam",
    output:
        bowtie2_bam_dedup_filt=(
            config["rdir"]
            + "/reads-dmg-filtering/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{dr_fraction}"
            + "/{bm_read_ani}"
            + "/{smp}.dmg-reads.dedup.filtered.bam"
        ),
        tp_stats=config["rdir"]
        + "/reads-dmg-filtering/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/weight-{weight}"
        + "/{dr_fraction}"
        + "/{bm_read_ani}"
        + "/{smp}.dmg-reads.dedup_stats.tsv.gz",
        tp_stats_filt=(
            config["rdir"]
            + "/reads-dmg-filtering/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{dr_fraction}"
            + "/{bm_read_ani}"
            + "/{smp}.dmg-reads.dedup_stats-filtered.tsv.gz"
        ),
    threads: config["tp_threads"]
    params:
        filter_bam_bin=config["filter_bam_bin"],
        bowtie2_bam_dedup=config["rdir"]
        + "/reads-dmg-filtering/{indata}/k{tp_bowtie2_k}/{read_ani}/weight-{weight}/{dr_fraction}/{bm_read_ani}/{smp}.dmg-reads.dedup.bam",
        tp_filter_sort_mem=config["tp_filter_sort_mem"],
        tp_filter_nreads=config["tp_filter_nreads"],
        tp_filter_exp_bratio=config["tp_filter_exp_bratio_dmg"],
        tp_filter_cov_evenness=config["tp_filter_cov_evenness"],
        tp_filter_read_ani="{read_ani}",
        filter_read_ani="{bm_read_ani}",
        tp_filter_entropy=config["tp_filter_entropy"],
        tp_filter_gini=config["tp_filter_gini"],
        tp_sizes=config.get("tp_sizes"),
        tp_norm_scale=config["tp_norm_scale"],
        tp_extra=config["tp_extra"],
        profiling_data_input="{indata}",
        rdir=config["rdir"]
        + "/reads-dmg-filtering/{indata}/k{tp_bowtie2_k}/{read_ani}/weight-{weight}/{dr_fraction}/{bm_read_ani}",
        wdir=config["wdir"],
    log:
        config["rdir"]
        + "/logs/reads-dmg-filtering/{indata}/k{tp_bowtie2_k}/{read_ani}/weight-{weight}/{dr_fraction}/{bm_read_ani}/{smp}.dmg-reads.taxonomic-profiling.log",
    benchmark:
        (
            config["rdir"]
            + "/benchmarks/reads-dmg-filtering/{indata}/k{tp_bowtie2_k}/{read_ani}/weight-{weight}/{dr_fraction}/{bm_read_ani}/{smp}.dmg-reads.taxonomic-profiling.bmk"
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
            -A {params.filter_read_ani} \
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
