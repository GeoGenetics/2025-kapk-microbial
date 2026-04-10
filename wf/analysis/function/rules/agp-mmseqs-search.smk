rule agp_mmseqs_kegg:
    """Search AGP proteins against KEGG using MMseqs2 with VTML20"""
    input:
        faa=config["rdir"] + "/agp-predict/{smp}.faa",
    output:
        out=config["rdir"] + "/agp-kegg/{smp}.hits.tsv.gz",
    threads: config["mmseqs2_threads"]
    params:
        mmseqs_bin=config["mmseqs_bin"],
        kegg_db=config["kegg_db"],
        vtml20=config["vtml20_matrix"],
        out=config["rdir"] + "/agp-kegg/{smp}.hits.tsv",
        rdir=config["rdir"] + "/agp-kegg",
        wdir=config["wdir"],
        tmp="{smp}-tmp",
        min_seq_id=config.get("agp_min_seq_id", "0.5"),
    conda:
        "../envs/kegg-xfilter.yaml"
    log:
        config["rdir"] + "/logs/agp-kegg/{smp}.mmseqs.log",
    benchmark:
        config["rdir"] + "/benchmarks/agp-kegg/{smp}.mmseqs.bmk"
    message:
        """--- MMseqs2 KEGG search: {wildcards.smp}"""
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}
        rm -rf {params.tmp}

        {params.mmseqs_bin} easy-search {input.faa} \
            {params.kegg_db} \
            {params.out} \
            {params.tmp} \
            --min-length 12 \
            -e 10.0 \
            --min-seq-id {params.min_seq_id} \
            -c 0.65 \
            --cov-mode 2 \
            --format-mode 0 \
            --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qaln,taln \
            --comp-bias-corr 0 \
            --mask 0 \
            --exact-kmer-matching 1 \
            --sub-mat "aa:{params.vtml20}" \
            --seed-sub-mat "aa:{params.vtml20}" \
            -s 2 \
            -k 6 \
            --spaced-kmer-pattern 11011101 \
            --max-seqs 10000 \
            --max-rejected 10 \
            --threads {threads} \
            --remove-tmp-files 0 \
            --use-all-table-starts 1 \
            2>&1 | tee {log}

        pigz -p {threads} {params.out}
        rm -rf {params.tmp}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """


rule agp_mmseqs_dbcan:
    """Search AGP proteins against dbCAN (CAZy) using MMseqs2 with VTML20"""
    input:
        faa=config["rdir"] + "/agp-predict/{smp}.faa",
    output:
        out=config["rdir"] + "/agp-dbcan/{smp}.hits.tsv.gz",
    threads: config["mmseqs2_threads"]
    params:
        mmseqs_bin=config["mmseqs_bin"],
        dbcan_db=config["dbcan_db"],
        vtml20=config["vtml20_matrix"],
        out=config["rdir"] + "/agp-dbcan/{smp}.hits.tsv",
        rdir=config["rdir"] + "/agp-dbcan",
        wdir=config["wdir"],
        tmp="{smp}-tmp",
        min_seq_id=config.get("agp_min_seq_id", "0.5"),
    conda:
        "../envs/kegg-xfilter.yaml"
    log:
        config["rdir"] + "/logs/agp-dbcan/{smp}.mmseqs.log",
    benchmark:
        config["rdir"] + "/benchmarks/agp-dbcan/{smp}.mmseqs.bmk"
    message:
        """--- MMseqs2 dbCAN search: {wildcards.smp}"""
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}
        rm -rf {params.tmp}

        {params.mmseqs_bin} easy-search {input.faa} \
            {params.dbcan_db} \
            {params.out} \
            {params.tmp} \
            --min-length 12 \
            -e 10.0 \
            --min-seq-id {params.min_seq_id} \
            -c 0.65 \
            --cov-mode 2 \
            --format-mode 0 \
            --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qaln,taln \
            --comp-bias-corr 0 \
            --mask 0 \
            --exact-kmer-matching 1 \
            --sub-mat "aa:{params.vtml20}" \
            --seed-sub-mat "aa:{params.vtml20}" \
            -s 2 \
            -k 6 \
            --spaced-kmer-pattern 11011101 \
            --max-seqs 10000 \
            --max-rejected 10 \
            --threads {threads} \
            --remove-tmp-files 0 \
            --use-all-table-starts 1 \
            2>&1 | tee {log}

        pigz -p {threads} {params.out}
        rm -rf {params.tmp}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """


rule agp_mmseqs_viral:
    """Search AGP proteins against viral protein database using MMseqs2 with VTML20"""
    input:
        faa=config["rdir"] + "/agp-predict/{smp}.faa",
    output:
        out=config["rdir"] + "/agp-viral/{smp}.hits.tsv.gz",
    threads: config["mmseqs2_threads"]
    params:
        mmseqs_bin=config["mmseqs_bin"],
        viral_db=config["viral_db"],
        vtml20=config["vtml20_matrix"],
        out=config["rdir"] + "/agp-viral/{smp}.hits.tsv",
        rdir=config["rdir"] + "/agp-viral",
        wdir=config["wdir"],
        tmp="{smp}-tmp",
        min_seq_id=config.get("agp_min_seq_id", "0.5"),
    conda:
        "../envs/kegg-xfilter.yaml"
    log:
        config["rdir"] + "/logs/agp-viral/{smp}.mmseqs.log",
    benchmark:
        config["rdir"] + "/benchmarks/agp-viral/{smp}.mmseqs.bmk"
    message:
        """--- MMseqs2 viral search: {wildcards.smp}"""
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}
        rm -rf {params.tmp}

        {params.mmseqs_bin} easy-search {input.faa} \
            {params.viral_db} \
            {params.out} \
            {params.tmp} \
            --min-length 12 \
            -e 10.0 \
            --min-seq-id {params.min_seq_id} \
            -c 0.65 \
            --cov-mode 2 \
            --format-mode 0 \
            --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qaln,taln \
            --comp-bias-corr 0 \
            --mask 0 \
            --exact-kmer-matching 1 \
            --sub-mat "aa:{params.vtml20}" \
            --seed-sub-mat "aa:{params.vtml20}" \
            -s 2 \
            -k 6 \
            --spaced-kmer-pattern 11011101 \
            --max-seqs 10000 \
            --max-rejected 10 \
            --threads {threads} \
            --remove-tmp-files 0 \
            --use-all-table-starts 1 \
            2>&1 | tee {log}

        pigz -p {threads} {params.out}
        rm -rf {params.tmp}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
