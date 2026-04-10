rule agp_damage_annotate_kegg:
    """Annotate damage on KEGG protein alignments"""
    input:
        hits=config["rdir"] + "/agp-kegg/{smp}.hits.tsv.gz",
    output:
        read_damage=config["rdir"] + "/agp-kegg-damage/{smp}.read_damage.tsv",
        protein_damage=config["rdir"] + "/agp-kegg-damage/{smp}.protein_damage.tsv",
        sites=config["rdir"] + "/agp-kegg-damage/{smp}.sites.tsv",
        anvio_ko=config["rdir"] + "/agp-kegg-damage/{smp}.anvio_ko.tsv.gz",
    threads: config.get("damage_threads", 8)
    params:
        agp_bin=config["agp_bin"],
        hits_noheader=config["rdir"] + "/agp-kegg-damage/{smp}.hits_noheader.tsv",
        anvio_ko_tmp=config["rdir"] + "/agp-kegg-damage/{smp}.anvio_ko.tsv",
        prior_ancient=config.get("prior_ancient", 0.10),
        min_reads=config.get("damage_min_reads", 3),
        min_breadth=config.get("damage_min_breadth", 0.10),
        min_depth=config.get("damage_min_depth", 0),
        preset=config.get("damage_preset", ""),
        auto_calibrate="--auto-calibrate-spurious" if config.get("damage_auto_calibrate", True) else "",
        kegg_gene_list=config["kegg_gene_list"],
        annotation_source=config.get("annotation_source", "KEGGgenes"),
    log:
        config["rdir"] + "/logs/agp-kegg-damage/{smp}.damage-annotate.log",
    benchmark:
        config["rdir"] + "/benchmarks/agp-kegg-damage/{smp}.damage-annotate.bmk"
    message:
        """--- AGP damage-annotate KEGG: {wildcards.smp}"""
    shell:
        """
        zcat {input.hits} > {params.hits_noheader}

        {params.agp_bin} damage-annotate \
            -i {params.hits_noheader} \
            -o {output.read_damage} \
            --protein-summary {output.protein_damage} \
            --sites {output.sites} \
            --anvio-ko {params.anvio_ko_tmp} \
            --map {params.kegg_gene_list} \
            --annotation-source {params.annotation_source} \
            --prior-ancient {params.prior_ancient} \
            --min-reads {params.min_reads} \
            --min-breadth {params.min_breadth} \
            --min-depth {params.min_depth} \
            $([ -n "{params.preset}" ] && echo "--preset {params.preset}") \
            {params.auto_calibrate} \
            -t {threads} \
            -v 2>&1 | tee {log}

        pigz -p {threads} -c {params.anvio_ko_tmp} > {output.anvio_ko}
        rm -f {params.hits_noheader} {params.anvio_ko_tmp}
        """


rule agp_damage_annotate_dbcan:
    """Annotate damage on dbCAN protein alignments"""
    input:
        hits=config["rdir"] + "/agp-dbcan/{smp}.hits.tsv.gz",
    output:
        read_damage=config["rdir"] + "/agp-dbcan-damage/{smp}.read_damage.tsv",
        protein_damage=config["rdir"] + "/agp-dbcan-damage/{smp}.protein_damage.tsv",
        sites=config["rdir"] + "/agp-dbcan-damage/{smp}.sites.tsv",
    threads: config.get("damage_threads", 8)
    params:
        agp_bin=config["agp_bin"],
        hits_noheader=config["rdir"] + "/agp-dbcan-damage/{smp}.hits_noheader.tsv",
        prior_ancient=config.get("prior_ancient", 0.10),
        min_reads=config.get("damage_min_reads", 3),
        min_breadth=config.get("damage_min_breadth", 0.10),
        min_depth=config.get("damage_min_depth", 0),
        preset=config.get("damage_preset", ""),
        auto_calibrate="--auto-calibrate-spurious" if config.get("damage_auto_calibrate", True) else "",
    log:
        config["rdir"] + "/logs/agp-dbcan-damage/{smp}.damage-annotate.log",
    benchmark:
        config["rdir"] + "/benchmarks/agp-dbcan-damage/{smp}.damage-annotate.bmk"
    message:
        """--- AGP damage-annotate dbCAN: {wildcards.smp}"""
    shell:
        """
        zcat {input.hits} > {params.hits_noheader}

        {params.agp_bin} damage-annotate \
            -i {params.hits_noheader} \
            -o {output.read_damage} \
            --protein-summary {output.protein_damage} \
            --sites {output.sites} \
            --prior-ancient {params.prior_ancient} \
            --min-reads {params.min_reads} \
            --min-breadth {params.min_breadth} \
            --min-depth {params.min_depth} \
            $([ -n "{params.preset}" ] && echo "--preset {params.preset}") \
            {params.auto_calibrate} \
            -t {threads} \
            -v 2>&1 | tee {log}

        rm -f {params.hits_noheader}
        """


rule agp_damage_annotate_viral:
    """Annotate damage on viral protein alignments"""
    input:
        hits=config["rdir"] + "/agp-viral/{smp}.hits.tsv.gz",
    output:
        read_damage=config["rdir"] + "/agp-viral-damage/{smp}.read_damage.tsv",
        protein_damage=config["rdir"] + "/agp-viral-damage/{smp}.protein_damage.tsv",
        sites=config["rdir"] + "/agp-viral-damage/{smp}.sites.tsv",
    threads: config.get("damage_threads", 8)
    params:
        agp_bin=config["agp_bin"],
        hits_noheader=config["rdir"] + "/agp-viral-damage/{smp}.hits_noheader.tsv",
        prior_ancient=config.get("prior_ancient", 0.10),
        min_reads=config.get("damage_min_reads", 3),
        min_breadth=config.get("damage_min_breadth", 0.10),
        min_depth=config.get("damage_min_depth", 0),
        preset=config.get("damage_preset", ""),
        auto_calibrate="--auto-calibrate-spurious" if config.get("damage_auto_calibrate", True) else "",
    log:
        config["rdir"] + "/logs/agp-viral-damage/{smp}.damage-annotate.log",
    benchmark:
        config["rdir"] + "/benchmarks/agp-viral-damage/{smp}.damage-annotate.bmk"
    message:
        """--- AGP damage-annotate viral: {wildcards.smp}"""
    shell:
        """
        zcat {input.hits} > {params.hits_noheader}

        {params.agp_bin} damage-annotate \
            -i {params.hits_noheader} \
            -o {output.read_damage} \
            --protein-summary {output.protein_damage} \
            --sites {output.sites} \
            --prior-ancient {params.prior_ancient} \
            --min-reads {params.min_reads} \
            --min-breadth {params.min_breadth} \
            --min-depth {params.min_depth} \
            $([ -n "{params.preset}" ] && echo "--preset {params.preset}") \
            {params.auto_calibrate} \
            -t {threads} \
            -v 2>&1 | tee {log}

        rm -f {params.hits_noheader}
        """


rule agp_filter_damaged_kegg:
    """Filter KEGG results to only damaged proteins for anvi'o"""
    input:
        protein_damage=config["rdir"] + "/agp-kegg-damage/{smp}.protein_damage.tsv",
        anvio=config["rdir"] + "/agp-kegg-damage/{smp}.anvio_ko.tsv.gz",
    output:
        damaged_anvio=config["rdir"] + "/agp-kegg-damaged/{smp}.damaged_group-abundances-anvio.tsv.gz",
        damaged_list=config["rdir"] + "/agp-kegg-damaged/{smp}.damaged_proteins.txt",
    threads: 1
    params:
        damage_threshold=config.get("damage_threshold", 0.5),
    log:
        config["rdir"] + "/logs/agp-kegg-damaged/{smp}.filter-damaged.log",
    message:
        """--- Filter damaged KEGG: {wildcards.smp}"""
    run:
        import pandas as pd
        import gzip

        # Read damage scores
        damage = pd.read_csv(input.protein_damage, sep='\t')

        # Filter to damaged proteins (p_protein_damaged >= threshold)
        damaged_ids = damage[damage['p_protein_damaged'] >= params.damage_threshold]['protein_id'].tolist()

        # Save damaged protein list
        with open(output.damaged_list, 'w') as f:
            for pid in damaged_ids:
                f.write(f"{pid}\n")

        # Filter anvi'o file (gene_callers_id column matches protein_id)
        anvio = pd.read_csv(input.anvio, sep='\t')
        damaged_anvio = anvio[anvio['gene_callers_id'].isin(damaged_ids)]
        damaged_anvio.to_csv(output.damaged_anvio, sep='\t', index=False, compression='gzip')

        print(f"Damaged proteins: {len(damaged_ids)} / {len(damage)} ({100*len(damaged_ids)/len(damage):.1f}%)")
