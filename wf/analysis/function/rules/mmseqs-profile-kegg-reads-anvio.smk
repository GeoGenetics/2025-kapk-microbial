rule mmseqs_profile_kegg_reads_anvio:
    input:
        kegg_abun=(
            config["rdir"]
            + "/kegg-reads/{smp}.profile.kegg_group-abundances-anvio.tsv.gz"
        ),
    output:
        hits_out=config["rdir"] + "/kegg-reads-anvio/{smp}_hits.txt",
        hits_modules_out=(config["rdir"] + "/kegg-reads-anvio/{smp}_module_paths.txt"),
        steps_modules_out=(config["rdir"] + "/kegg-reads-anvio/{smp}_module_steps.txt"),
        modules_out=config["rdir"] + "/kegg-reads-anvio/{smp}_modules.txt",
    threads: config["anvio_threads"]
    params:
        rdir=config["rdir"] + "/kegg-reads-anvio",
        wdir=config["wdir"],
        label="{smp}",
        anvi_estimate_extra_parms=config["anvi_estimate_extra_parms"],
    conda:
        "../envs/kegg-xfilter-anvio.yaml"
    log:
        config["rdir"] + "/logs/kegg-reads-anvio/{smp}.kegg-reads.log",
    benchmark:
        config["rdir"] + "/benchmarks/kegg-reads-anvio/{smp}.kegg-reads.bmk"
    message:
        """--- Annotate KEGGs on reads."""
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        rm -rf {output.hits_out} {output.hits_modules_out} {output.modules_out}

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
            {params.anvi_estimate_extra_parms}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
