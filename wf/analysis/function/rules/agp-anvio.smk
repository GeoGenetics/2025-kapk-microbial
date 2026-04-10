rule agp_anvio_kegg_damaged:
    """Run anvi-estimate-metabolism on damage-filtered KEGG hits"""
    input:
        kegg_abun=config["rdir"] + "/agp-kegg-damaged/{smp}.damaged_group-abundances-anvio.tsv.gz",
    output:
        hits_out=config["rdir"] + "/agp-kegg-anvio-damaged/{smp}_hits.txt",
        hits_modules_out=config["rdir"] + "/agp-kegg-anvio-damaged/{smp}_module_paths.txt",
        steps_modules_out=config["rdir"] + "/agp-kegg-anvio-damaged/{smp}_module_steps.txt",
        modules_out=config["rdir"] + "/agp-kegg-anvio-damaged/{smp}_modules.txt",
    threads: config.get("anvio_threads", 4)
    params:
        rdir=config["rdir"] + "/agp-kegg-anvio-damaged",
        wdir=config["wdir"],
        label="{smp}",
        anvi_estimate_extra_parms=config.get("anvi_estimate_extra_parms", ""),
    conda:
        "../envs/kegg-xfilter-anvio.yaml"
    log:
        config["rdir"] + "/logs/agp-kegg-anvio-damaged/{smp}.anvio.log",
    benchmark:
        config["rdir"] + "/benchmarks/agp-kegg-anvio-damaged/{smp}.anvio.bmk"
    message:
        """--- Anvi'o metabolism (damaged): {wildcards.smp}"""
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}
        rm -rf {output.hits_out} {output.hits_modules_out} {output.modules_out} {output.steps_modules_out}

        if [ ! -s {input.kegg_abun} ]; then
            echo "KEGG abundance file is empty"
            touch {output.hits_out}
            touch {output.hits_modules_out}
            touch {output.modules_out}
            touch {output.steps_modules_out}
            cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
            exit 0
        fi

        anvi-estimate-metabolism \
            --add-coverage \
            --enzymes-txt {input.kegg_abun} \
            --output-modes modules,module_paths,module_steps,hits \
            -O {params.label} \
            --include-kos-not-in-kofam \
            {params.anvi_estimate_extra_parms} \
            2>&1 | tee {log}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """


rule agp_anvio_kegg_all:
    """Run anvi-estimate-metabolism on all KEGG hits (not filtered by damage)"""
    input:
        kegg_abun=config["rdir"] + "/agp-kegg-damage/{smp}.anvio_ko.tsv.gz",
    output:
        hits_out=config["rdir"] + "/agp-kegg-anvio/{smp}_hits.txt",
        hits_modules_out=config["rdir"] + "/agp-kegg-anvio/{smp}_module_paths.txt",
        steps_modules_out=config["rdir"] + "/agp-kegg-anvio/{smp}_module_steps.txt",
        modules_out=config["rdir"] + "/agp-kegg-anvio/{smp}_modules.txt",
    threads: config.get("anvio_threads", 4)
    params:
        rdir=config["rdir"] + "/agp-kegg-anvio",
        wdir=config["wdir"],
        label="{smp}",
        anvi_estimate_extra_parms=config.get("anvi_estimate_extra_parms", ""),
    conda:
        "../envs/kegg-xfilter-anvio.yaml"
    log:
        config["rdir"] + "/logs/agp-kegg-anvio/{smp}.anvio.log",
    benchmark:
        config["rdir"] + "/benchmarks/agp-kegg-anvio/{smp}.anvio.bmk"
    message:
        """--- Anvi'o metabolism (all): {wildcards.smp}"""
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}
        rm -rf {output.hits_out} {output.hits_modules_out} {output.modules_out} {output.steps_modules_out}

        if [ ! -s {input.kegg_abun} ]; then
            echo "KEGG abundance file is empty"
            touch {output.hits_out}
            touch {output.hits_modules_out}
            touch {output.modules_out}
            touch {output.steps_modules_out}
            cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
            exit 0
        fi

        anvi-estimate-metabolism \
            --add-coverage \
            --enzymes-txt {input.kegg_abun} \
            --output-modes modules,module_paths,module_steps,hits \
            -O {params.label} \
            --include-kos-not-in-kofam \
            {params.anvi_estimate_extra_parms} \
            2>&1 | tee {log}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
