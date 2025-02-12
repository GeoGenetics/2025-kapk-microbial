rule prefilter_profiling_hires:
    input:
        read_ids=(
            config["rdir"]
            + "/prefilter-profiling-coarse/{smp}.prefilter-coarse.ids.gz"
        ),
        reads=config["rdir"] + "/read-derep/{smp}.fa.gz",
    output:
        k_out=config["rdir"] + "/prefilter-profiling-hires/{smp}.prefilter-hires.out.gz",
        k_report=(
            config["rdir"]
            + "/prefilter-profiling-hires/{smp}.prefilter-hires-filt.report.gz"
        ),
        read_ids=(
            config["rdir"] + "/prefilter-profiling-hires/{smp}.prefilter-hires.ids.gz"
        ),
        k_hires_out_tax_summary=(
            config["rdir"]
            + "/prefilter-profiling-hires/{smp}.prefilter-hires.tax-summary.tsv.gz"
        ),
        k_hires_out_tax=(
            config["rdir"]
            + "/prefilter-profiling-hires/{smp}.prefilter-hires.tax.tsv.gz"
        ),
        bac_ids=config["rdir"] + "/prefilter-profiling-hires/{smp}_Bacteria.txt.gz",
        arc_ids=config["rdir"] + "/prefilter-profiling-hires/{smp}_Archaea.txt.gz",
        vir_ids=config["rdir"] + "/prefilter-profiling-hires/{smp}_Viruses.txt.gz",
    threads: config["prefilter_k2_threads"]
    params:
        prefilter_bin=config["prefilter_bin"],
        awk_bin=config["awk_bin"],
        taxonkit_bin=config["taxonkit_bin"],
        seqkit_bin=config["seqkit_bin"],
        csvtk_bin=config["csvtk_bin"],
        datamash_bin=config["datamash_bin"],
        conifer_bin=config["conifer_bin"],
        prefilter_report_bin=config["prefilter_report_bin"],
        db=config["prefilter_k2_db_hires"],
        prefilter_k2_db_hires_taxonomy=config["prefilter_k2_db_hires_taxonomy"],
        prefilter_k2_db_hires_taxo=config["prefilter_k2_db_hires_taxo"],
        k_out_min=(
            config["rdir"]
            + "/prefilter-profiling-hires/{smp}.prefilter-hires-minimizer.out"
        ),
        k_out=config["rdir"] + "/prefilter-profiling-hires/{smp}.prefilter-hires.out",
        k_out_filt=(
            config["rdir"]
            + "/prefilter-profiling-hires/{smp}.prefilter-hires-filt.out"
        ),
        k_out_filt_tmp=(
            config["rdir"]
            + "/prefilter-profiling-hires/{smp}.prefilter-hires-filt.tmp"
        ),
        k_report_filt=(
            config["rdir"]
            + "/prefilter-profiling-hires/{smp}.prefilter-hires-filt.report"
        ),
        k_report=(
            config["rdir"] + "/prefilter-profiling-hires/{smp}.prefilter-hires.report"
        ),
        read_ids=config["rdir"] + "/prefilter-profiling-hires/{smp}.prefilter-hires.ids",
        k_confidence=config["k_confidence_hires"],
        k_rtl=config["k_rtl_hires"],
        k_min_hit_group=config["k_min_hit_group"],
        unique_kmers=config["unique_kmers"],
        label="{smp}",
        rdir=config["rdir"] + "/prefilter-profiling-hires",
        wdir=config["wdir"],
        k_out_fq=config["rdir"] + "/prefilter-profiling-hires/{smp}.prefilter-coarse.fq",
        header="label,name,taxonomy_id,taxonomy_lvl,kraken_assigned_reads,added_reads,new_est_reads,fraction_total_reads,tax_string",
    log:
        config["rdir"] + "/logs/prefilter-profiling-hires/{smp}.prefilter-hires.log",
    benchmark:
        (
            config["rdir"]
            + "/benchmarks/prefilter-profiling-hires/{smp}.prefilter-hires.bmk"
        )
    conda:
        "../envs/prefilter.yaml"
    message:
        """--- High-res taxonomic profiling with prefilter."""
    shell:
        """
        set -x
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        rm -rf {params.k_out_fq} {output.k_out} {output.k_report} {params.rdir}/{params.label}*gz

        # Extract reads
        {params.seqkit_bin} grep \
            -f <(zcat {input.read_ids}) \
            -j {threads} \
            {input.reads} -o {params.k_out_fq}

        # run prefilter
        {params.prefilter_bin} \
            --memory-mapping \
            --db {params.db} \
            --threads {threads} \
            --report {params.k_out_min} \
            --output {params.k_out} \
            --confidence {params.k_confidence} \
            --report-minimizer-data \
            --minimum-hit-groups {params.k_min_hit_group} \
            {params.k_out_fq} >> {log} 2>&1

        {params.conifer_bin} -i {params.k_out} \
            -d {params.prefilter_k2_db_hires_taxo} \
            -b -a > {params.k_out_filt}

        {params.awk_bin} -vFS="\t" -vR={params.k_rtl} -vA={params.k_out_filt_tmp} -vB={params.read_ids} '{{if ($1 == "C" && $7 > R){{print $2"\t"$3 > A}}else{{print $2 > B}}}}' {params.k_out_filt} 

        # get tax domain
        {params.taxonkit_bin} lineage -j {threads} -i 2 --data-dir {params.prefilter_k2_db_hires_taxonomy} {params.k_out_filt_tmp} \
            | {params.awk_bin} -vFS="\t" -vL={params.label} '{{print L"\t"$0}}' \
            | gzip > {output.k_hires_out_tax} 

        zcat {output.k_hires_out_tax} \
            | {params.awk_bin} -vFS="\t" -vL={params.label} '{{split($4, a, ";"); gsub("d__", "", a[1]); if(a[1]==""){{a[1]="Unclassified"}}; print $1 > L"_"a[1]".txt"}}'

        if [[ -s {params.label}_Bacteria.txt ]]; then
            cat {params.label}_Bacteria.txt >> {params.read_ids}
        else
            touch {params.label}_Bacteria.txt.gz
        fi

        if [[ -s {params.label}_Archaea.txt ]]; then
            cat {params.label}_Archaea.txt >> {params.read_ids}
        else
            touch {params.label}_Archaea.txt.gz
        fi

        if [[ -s {params.label}_root.txt ]]; then
            cat {params.label}_root.txt >> {params.read_ids}
        else
            touch {params.label}_root.txt.gz
        fi

        if [[ ! -s {params.label}_Viruses.txt ]]; then
            touch {params.label}_Viruses.txt.gz
        fi

        # # Create report
        # {params.prefilter_report_bin} {params.prefilter_k2_db_hires_taxo} {params.k_out_filt_tmp} {params.k_report_filt}

        # Use report with minimizers
        # filter out the ones with less N unique kmers
        awk -vF={params.unique_kmers} '$5 >= F' {params.k_out_min} | cut -f1-3,6-8 > {params.k_report_filt}

        # Get summary for hires
        zcat {output.k_hires_out_tax} \
            | {params.awk_bin} -vFS="\t" -vL={params.label} '{{split($4, a, ";"); gsub("d__", "", a[1]); if(a[1]==""){{a[1]="Unclassified"}};b[a[1]]++}}END{{for (i in b){{print L"\t"i"\t"b[i]}}}}' \
            | {params.csvtk_bin} add-header -t -n 'label,domain,counts' | gzip > {output.k_hires_out_tax_summary}

        pigz -p {threads} {params.k_out} {params.k_report_filt} {params.k_out_min} {params.label}_*.txt {params.read_ids}
        rm -rf {params.k_out_filt_tmp} {params.k_report} \
                {params.k_out_fq} {params.k_out_filt}
        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
