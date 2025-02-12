rule prefilter_profiling_coarse:
    input:
        reads=config["rdir"] + "/read-derep/{smp}.fa.gz",
    output:
        k_out=config["rdir"]
        + "/prefilter-profiling-coarse/{smp}.prefilter-coarse.out.gz",
        k_coarse_out_tax_sumamry=(
            config["rdir"]
            + "/prefilter-profiling-coarse/{smp}.prefilter-coarse.tax-summary.tsv.gz"
        ),
        k_coarse_out_report=(
            config["rdir"]
            + "/prefilter-profiling-coarse/{smp}.prefilter-coarse.report.txt.gz"
        ),
        k_coarse_out_tax=(
            config["rdir"]
            + "/prefilter-profiling-coarse/{smp}.prefilter-coarse.tax.tsv.gz"
        ),
        unclassified_ids=(
            config["rdir"] + "/prefilter-profiling-coarse/{smp}_unclassified.txt.gz"
        ),
        euk_ids=config["rdir"] + "/prefilter-profiling-coarse/{smp}_Eukaryota.txt.gz",
        root_ids=config["rdir"] + "/prefilter-profiling-coarse/{smp}_root.txt.gz",
        bac_ids=config["rdir"] + "/prefilter-profiling-coarse/{smp}_Bacteria.txt.gz",
        arc_ids=config["rdir"] + "/prefilter-profiling-coarse/{smp}_Archaea.txt.gz",
        vir_ids=config["rdir"] + "/prefilter-profiling-coarse/{smp}_Viruses.txt.gz",
    threads: config["prefilter_k2_threads"]
    params:
        prefilter_k2_bin=config["prefilter_k2_bin"],
        awk_bin=config["awk_bin"],
        taxonkit_bin=config["taxonkit_bin"],
        seqkit_bin=config["seqkit_bin"],
        conifer_bin=config["conifer_bin"],
        csvtk_bin=config["csvtk_bin"],
        datamash_bin=config["datamash_bin"],
        db=config["prefilter_k2_db_coarse"],
        prefilter_k2_db_coarse_taxonomy=config["prefilter_k2_db_coarse_taxonomy"],
        prefilter_k2_db_coarse_taxo=config["prefilter_k2_db_coarse_taxo"],
        k_out=config["rdir"] + "/prefilter-profiling-coarse/{smp}.prefilter-coarse.out",
        k_out_filt=(
            config["rdir"]
            + "/prefilter-profiling-coarse/{smp}.prefilter-coarse-filt.out"
        ),
        k_out_filt_tmp=(
            config["rdir"]
            + "/prefilter-profiling-coarse/{smp}.prefilter-coarse-filt.tmp"
        ),
        unclassified_ids=config["rdir"]
        + "/prefilter-profiling-coarse/{smp}_unclassified.txt",
        k_coarse_out_report=(
            config["rdir"]
            + "/prefilter-profiling-coarse/{smp}.prefilter-coarse.report.txt"
        ),
        k_confidence=config["k_confidence_coarse"],
        k_rtl=config["k_rtl_coarse"],
        k_min_hit_group=config["k_min_hit_group"],
        k_extra_parms_coarse=config["k_extra_parms_coarse"],
        label="{smp}",
        rdir=config["rdir"] + "/prefilter-profiling-coarse",
        wdir=config["wdir"],
    log:
        config["rdir"] + "/logs/prefilter-profiling-coarse/{smp}.prefilter-coarse.log",
    benchmark:
        (
            config["rdir"]
            + "/benchmarks/prefilter-profiling-coarse/{smp}.prefilter-coarse.bmk"
        )
    conda:
        "../envs/prefilter.yaml"
    message:
        """--- Coarse taxonomic profiling with prefilter."""
    shell:
        """
        set -x

        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}
        rm -rf {output.unclassified_ids} {output.k_out} {output.k_coarse_out_tax} {params.unclassified_ids}

        # run prefilter
        {params.prefilter_k2_bin} \
            --db {params.db} \
            --threads {threads} \
            --output {params.k_out} \
            --confidence {params.k_confidence} \
            --minimum-hit-groups {params.k_min_hit_group} \
            --report {params.k_coarse_out_report} \
            {params.k_extra_parms_coarse} \
            {input.reads} >> {log} 2>&1

        # filter prefilter output
        {params.conifer_bin} -i {params.k_out} \
            -d {params.prefilter_k2_db_coarse_taxo} \
            -b -a > {params.k_out_filt} 

        # split reads between classified (k_out_filt_tmp), non-classified (unclassified_ids)
        {params.awk_bin} -vFS=$'\t' -vR={params.k_rtl} -vA={params.k_out_filt_tmp} -vB={params.unclassified_ids} '{{if ($1 == "C" && $7 > R){{print $2"\t"$3 > A}}else{{print $2 > B}}}}' {params.k_out_filt} 

        # get tax domain of the classified ones
        {params.taxonkit_bin} lineage -j {threads} -i 2 --data-dir {params.prefilter_k2_db_coarse_taxonomy} {params.k_out_filt_tmp} \
            | {params.awk_bin} -vFS=$'\t' -vL={params.label} '{{print L"\t"$0}}' \
            | gzip > {output.k_coarse_out_tax} 

        # separate reads by domain
        zcat {output.k_coarse_out_tax} \
            | {params.awk_bin} -vFS=$'\t' -vL={params.label} '{{split($4, a, ";"); gsub("d__", "", a[1]); if(a[1]==""){{a[1]="Unclassified"}}; print $2 > L"_"a[1]".txt"}}'

        if [[ ! -s {params.label}_Bacteria.txt ]]; then
            echo "No Bacteria reads found"
            touch {output.bac_ids}
        fi

        if [[ ! -s {params.label}_Archaea.txt ]]; then
            echo "No Archaea reads found"
            touch {output.arc_ids}
        fi

        if [[ ! -s {params.label}_Viruses.txt ]]; then
            echo "No Viruses reads found"
            touch {output.vir_ids}
        fi

        if [[ ! -s {params.label}_root.txt ]]; then
            echo "No root reads found"
            touch {output.root_ids}
        fi

        if [[ ! -s {params.label}_Eukaryota.txt ]]; then
            echo "No Eukaryota reads found"
            touch {output.euk_ids}
        fi

        if [[ ! -s {params.unclassified_ids} ]]; then
            echo "No unclassified reads found"
            touch {output.unclassified_ids}
        fi

        # Get summary for coarse
        zcat {output.k_coarse_out_tax} \
            | {params.awk_bin} -vFS=$'\t' -vL={params.label} '{{split($4, a, ";"); gsub("d__", "", a[1]); if(a[1]==""){{a[1]="Unclassified"}};b[a[1]]++}}END{{for (i in b){{print L"\t"i"\t"b[i]}}}}' \
            | {params.csvtk_bin} add-header -t -n 'label,domain,counts' | gzip > {output.k_coarse_out_tax_sumamry}

        rm -rf {params.k_out_filt_tmp} {params.k_out_filt}
        pigz {params.k_out} {params.label}_*.txt {params.k_coarse_out_report}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
