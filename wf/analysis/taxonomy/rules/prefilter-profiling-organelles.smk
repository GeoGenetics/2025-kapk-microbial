rule prefilter_profiling_organelles:
    input:
        unclassified_ids=(
            config["rdir"] + "/prefilter-profiling-coarse/{smp}_unclassified.txt.gz"
        ),
        euk_ids=config["rdir"] + "/prefilter-profiling-coarse/{smp}_Eukaryota.txt.gz",
        root_ids=config["rdir"] + "/prefilter-profiling-coarse/{smp}_root.txt.gz",
        reads=config["rdir"] + "/read-derep/{smp}.fa.gz",
    output:
        k_out=config["rdir"]
        + "/prefilter-profiling-organelles/{smp}.prefilter-organelles.out.gz",
        k_organelles_out_tax=(
            config["rdir"]
            + "/prefilter-profiling-organelles/{smp}.prefilter-organelles.tax.tsv.gz"
        ),
        k_organelles_out_report=(
            config["rdir"]
            + "/prefilter-profiling-organelles/{smp}.prefilter-organelles.report.txt.gz"
        ),
        euk_ids=(
            config["rdir"] + "/prefilter-profiling-organelles/{smp}_Eukaryota.txt.gz",
        ),
        root_ids=(config["rdir"] + "/prefilter-profiling-organelles/{smp}_root.txt.gz",),
        unclassified_ids=(
            config["rdir"]
            + "/prefilter-profiling-organelles/{smp}_unclassified.txt.gz",
        ),
    threads: config["prefilter_k2_threads"]
    params:
        prefilter_k2_bin=config["prefilter_k2_bin"],
        awk_bin=config["awk_bin"],
        taxonkit_bin=config["taxonkit_bin"],
        seqkit_bin=config["seqkit_bin"],
        csvtk_bin=config["csvtk_bin"],
        datamash_bin=config["datamash_bin"],
        conifer_bin=config["conifer_bin"],
        db=config["prefilter_k2_db_organelles"],
        prefilter_k2_db_organelles_taxonomy=config[
            "prefilter_k2_db_organelles_taxonomy"
        ],
        prefilter_k2_db_organelles_taxo=config["prefilter_k2_db_organelles_taxo"],
        k_out=config["rdir"]
        + "/prefilter-profiling-organelles/{smp}.prefilter-organelles.out",
        k_out_filt=(
            config["rdir"]
            + "/prefilter-profiling-organelles/{smp}.prefilter-organelles-filt.out"
        ),
        k_out_filt_tmp=(
            config["rdir"]
            + "/prefilter-profiling-organelles/{smp}.prefilter-organelles-filt.tmp"
        ),
        read_ids_tmp=config["rdir"]
        + "/prefilter-profiling-organelles/{smp}.prefilter-organelles.tmp.ids",
        unclassified_ids=config["rdir"]
        + "/prefilter-profiling-organelles/{smp}_unclassified.txt",
        k_organelles_out_report=(
            config["rdir"]
            + "/prefilter-profiling-organelles/{smp}.prefilter-organelles.report.txt"
        ),
        k_confidence=config["k_confidence_organelles"],
        k_rtl=config["k_rtl_organelles"],
        k_min_hit_group=config["k_min_hit_group"],
        k_extra_parms_organelles=config["k_extra_parms_organelles"],
        label="{smp}",
        rdir=config["rdir"] + "/prefilter-profiling-organelles",
        wdir=config["wdir"],
        k_out_fq=config["rdir"]
        + "/prefilter-profiling-organelles/{smp}.prefilter-coarse.fq",
        header="label,name,taxonomy_id,taxonomy_lvl,kraken_assigned_reads,added_reads,new_est_reads,fraction_total_reads,tax_string",
    log:
        config["rdir"]
        + "/logs/prefilter-profiling-organelles/{smp}.prefilter-organelles.log",
    benchmark:
        (
            config["rdir"]
            + "/benchmarks/prefilter-profiling-organelles/{smp}.prefilter-organelles.bmk"
        )
    conda:
        "../envs/prefilter.yaml"
    message:
        """--- High-res taxonomic profiling with prefilter."""
    shell:
        """
        set -x
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        rm -rf {params.k_out_fq} {output.k_out} {params.rdir}/{params.label}*gz {params.read_ids_tmp}

        # Check if we have read euk_ids and root_ids
        if [ ! -s {input.euk_ids} ] && [ ! -s {input.root_ids} ] && [ ! -s {input.unclassified_ids} ]; then
            echo "Cannot find euk_ids or root_ids"
            touch {output.k_out}
            touch {output.k_organelles_out_tax}
            touch {output.euk_ids}
            touch {output.root_ids}
            touch {output.unclassified_ids}
            exit 0
        fi

        # Extract reads
        if [[ -s {input.euk_ids} ]]; then
            zcat {input.euk_ids} >  {params.read_ids_tmp}
        fi
        if [[ -s {input.root_ids} ]]; then
            zcat {input.root_ids} >> {params.read_ids_tmp}
        fi
        if [[ -s {input.unclassified_ids} ]]; then
            zcat {input.unclassified_ids} >> {params.read_ids_tmp}
        fi

        {params.seqkit_bin} grep \
            -f <(sort -u --parallel {threads} -S20% {params.read_ids_tmp}) \
            -j {threads} \
            -o {params.k_out_fq} {input.reads}

        # run prefilter
        {params.prefilter_k2_bin} \
            --memory-mapping \
            --db {params.db} \
            --threads {threads} \
            --output {params.k_out} \
            --confidence {params.k_confidence} \
            --minimum-hit-groups {params.k_min_hit_group} \
            --report {params.k_organelles_out_report} \
            {params.k_extra_parms_organelles} \
            {params.k_out_fq} >> {log} 2>&1

        # filter prefilter output
        {params.conifer_bin} -i {params.k_out} \
            -d {params.prefilter_k2_db_organelles_taxo} \
            -b -a > {params.k_out_filt}

        # separate bewteen classified (k_out_filt_tmp) and unclassified (unclassified_ids)
        {params.awk_bin} -vFS="\t" -vR={params.k_rtl} -vA={params.k_out_filt_tmp} -vB={params.unclassified_ids} '{{if ($1 == "C" && $7 > R){{print $2"\t"$3 > A}}else{{print $2 > B}}}}' {params.k_out_filt} 

        # get tax domain
        {params.taxonkit_bin} lineage -j {threads} -i 2 --data-dir {params.prefilter_k2_db_organelles_taxonomy} {params.k_out_filt_tmp} \
            | {params.awk_bin} -vFS="\t" -vL={params.label} '{{print L"\t"$0}}' \
            | gzip > {output.k_organelles_out_tax} 

        zcat {output.k_organelles_out_tax} \
            | {params.awk_bin} -vFS=$'\t' -vL={params.label} '{{split($4, a, ";"); gsub("d__", "", a[1]); if(a[1]==""){{a[1]="Unclassified"}}; print $2 > L"_"a[1]".txt"}}'

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

        # Create files to be combined with arc/bac/viruses. Will contain all organelle classified reads, root

        pigz -p {threads} {params.k_out} {params.label}_*.txt {params.k_organelles_out_report}
        rm -rf {params.k_out_filt_tmp} {params.k_out_fq} {params.k_out_filt} {params.read_ids_tmp}
        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
