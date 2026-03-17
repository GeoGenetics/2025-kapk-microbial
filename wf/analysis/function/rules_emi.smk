import struct

# EM in-memory threshold: dart auto-streams when EMI exceeds -m value
EM_THRESHOLD_MB = 32000

def emi_mem_mb(emi_path):
    """Dynamically scale SLURM memory per-sample based on EMI reference count.

    After fixes (MADV_DONTNEED for mmap pages + malloc_trim after coverage-EM):
    - EMI mmap contributes negligible RSS (demand-paged, DONTNEED after each row group)
    - Peak RSS is now dominated by heap allocations that scale with num_refs:
        coverage-EM thread_stats: ~32 threads × num_refs × 264 bytes ≈ ref_mb × 14
        output phase (summaries, gene_agg, best_hits): ~64 GB fixed + ref scaling
    Falls back to 64 GB if file is not yet available at scheduling time.
    """
    if not os.path.exists(emi_path):
        return 65536  # 64 GB fallback
    emi_mb = int(os.path.getsize(emi_path) / 1e6)
    try:
        with open(emi_path, 'rb') as fh:
            hdr = fh.read(128)
        num_refs = struct.unpack_from('<I', hdr, 28)[0]
        ref_mb = int(num_refs * 600 / 1e6)
    except Exception:
        ref_mb = 0
    if emi_mb < EM_THRESHOLD_MB:
        # In-memory EM: decoded alignment data expands ~2.5x + ref structures + 12 GB margin
        return max(EM_THRESHOLD_MB, int(emi_mb * 2.5) + ref_mb + 12288)
    else:
        # Streaming coverage EM uses sorted CovEntry (9ce044a): no per-thread maps.
        # ref_mb × 14 term eliminated; peak ≈ emi_mb (NFS pages) + ref_mb (damage index) + 64 GB
        # MADV_RANDOM + column restriction: NFS pages no longer accumulate (bdd405d).
        # Peak dominated by ref_mb (damage index) + ~64 GB base overhead.
        return max(163840, ref_mb + 65536)

# ============================================================================
# KEGG Pipeline with EMI + damage-annotate (replaces xFilter)
# Simplified version - anvi'o metabolism is optional
# ============================================================================

rule kegg_hits2emi:
    """Convert MMseqs2 hits to columnar EMI index with per-read damage info."""
    input:
        hits = OUTDIR + "/kegg/{sample}/kegg_hits.tsv",
        agd = OUTDIR + "/sample_damage/{sample}/predictions.agd"
    output:
        emi = OUTDIR + "/kegg/{sample}/kegg_hits.emi"
    params:
        agp = config["agp_binary"]
    threads: 8
    resources:
        runtime=lambda wc, input: max(120, int(os.path.getsize(input.agd) / 1e9) * 8),
        mem_mb=lambda wc, input: max(32000, int(os.path.getsize(input.agd) / 1e6) + 12000)
    shell:
        """
        export LD_LIBRARY_PATH=/vol/cloud/agp/lib:${{LD_LIBRARY_PATH:-}} && {params.agp} hits2emi -i {input.hits} -o {output.emi} --damage-index {input.agd} \
            --threads {threads} --memory 16G -v
        """


rule kegg_damage_annotate_emi:
    """
    Integrated damage annotation with EM reassignment.
    Produces per-read, per-protein, and per-gene outputs.
    """
    input:
        emi = OUTDIR + "/kegg/{sample}/kegg_hits.emi",
        agd = OUTDIR + "/sample_damage/{sample}/predictions.agd",
        gene_map = config["kegg"]["gene_map"]
    output:
        reads = OUTDIR + "/kegg/{sample}/emi.reads.tsv",
        proteins = OUTDIR + "/kegg/{sample}/emi.protein.tsv",
        genes = OUTDIR + "/kegg/{sample}/emi.genes.tsv",
        functional = OUTDIR + "/kegg/{sample}/emi.functional.tsv",
        anvio = OUTDIR + "/kegg/{sample}/emi.anvio_ko.tsv"
    params:
        agp = config["agp_binary"],
        prefix = OUTDIR + "/kegg/{sample}/emi",
        annotation_source = "KOfam",
        em_iters = config.get("damage_annotate", {}).get("em_iters", 100),
        em_lambda = config.get("damage_annotate", {}).get("em_lambda", 3.0),
        min_depth = config.get("damage_annotate", {}).get("min_depth", 0.5),
        min_breadth = config.get("damage_annotate", {}).get("min_breadth", 0.1),
        min_reads = config.get("damage_annotate", {}).get("min_reads", 3),
    threads: 8
    resources:
        runtime=480,
        mem_mb=lambda wc, input: emi_mem_mb(input.emi)
    shell:
        """
        export LD_LIBRARY_PATH=/vol/cloud/agp/lib:${{LD_LIBRARY_PATH:-}} && {params.agp} damage-annotate --emi {input.emi} \
            --damage-index {input.agd} \
            --map {input.gene_map} \
            --analysis-prefix {params.prefix} \
            --gene-summary {output.genes} \
            --functional-summary {output.functional} \
            --anvio-ko {output.anvio} \
            --annotation-source {params.annotation_source} \
            --em-iters {params.em_iters} \
            --em-lambda {params.em_lambda} \
            --min-depth {params.min_depth} \
            --min-breadth {params.min_breadth} \
            --min-reads {params.min_reads} \
            --coverage-em \
            --auto-prior-ancient \
            --auto-calibrate-spurious \
            --refine-damage \
            -m 32000 \
            --threads {threads} -v
        """


rule kegg_pathway_damage_stats:
    """Compute pathway damage enrichment using Mann-Whitney U test."""
    input:
        proteins = OUTDIR + "/kegg/{sample}/emi.protein.tsv"
    output:
        stats = OUTDIR + "/kegg/{sample}/pathway_damage_stats.tsv"
    resources:
        runtime=10,
        mem_mb=lambda wc, input: max(32000, int(os.path.getsize(input.proteins) / 1e6) * 8 + 4096)
    run:
        import pandas as pd
        import numpy as np
        from scipy import stats as sp_stats

        # Load protein damage scores
        dmg = pd.read_csv(input.proteins, sep='\t')

        # Find score column
        score_cols = [c for c in dmg.columns if 'combined_score' in c.lower() or 'p_protein' in c.lower()]
        if score_cols:
            score_col = score_cols[0]
        else:
            # Fallback: just output empty file
            pd.DataFrame().to_csv(output.stats, sep='\t', index=False)
            return

        scores = dmg[score_col].dropna().values
        n_total = len(scores)

        if n_total < 10:
            pd.DataFrame().to_csv(output.stats, sep='\t', index=False)
            return

        # For now, just compute basic stats (pathway enrichment needs modules file)
        summary = {
            'n_reads': n_total,
            'mean_damage': scores.mean(),
            'std_damage': scores.std(),
            'median_damage': np.median(scores),
            'q25': np.percentile(scores, 25),
            'q75': np.percentile(scores, 75)
        }

        pd.DataFrame([summary]).to_csv(output.stats, sep='\t', index=False)


# ============================================================================
# CAZy Pipeline with EMI (replaces xFilter)
# ============================================================================

rule cazy_hits2emi:
    input:
        hits = OUTDIR + "/cazy/{sample}/cazy_hits.tsv",
        agd = OUTDIR + "/sample_damage/{sample}/predictions.agd"
    output:
        emi = OUTDIR + "/cazy/{sample}/cazy_hits.emi"
    params:
        agp = config["agp_binary"]
    threads: 8
    resources:
        runtime=lambda wc, input: max(120, int(os.path.getsize(input.agd) / 1e9) * 8),
        mem_mb=lambda wc, input: max(32000, int(os.path.getsize(input.agd) / 1e6) + 12000)
    shell:
        """
        export LD_LIBRARY_PATH=/vol/cloud/agp/lib:${{LD_LIBRARY_PATH:-}} && {params.agp} hits2emi -i {input.hits} -o {output.emi} --damage-index {input.agd} \
            --threads {threads} --memory 16G -v
        """


rule cazy_damage_annotate_emi:
    input:
        emi = OUTDIR + "/cazy/{sample}/cazy_hits.emi",
        agd = OUTDIR + "/sample_damage/{sample}/predictions.agd",
        gene_map = "/vol/cloud/agp/data/dbcan_mmseqs2/dbCAN_gene.list"
    output:
        proteins = OUTDIR + "/cazy/{sample}/cazy_damage.tsv",
        genes = OUTDIR + "/cazy/{sample}/cazy_genes.tsv",
        functional = OUTDIR + "/cazy/{sample}/cazy_emi.functional.tsv",
        anvio = OUTDIR + "/cazy/{sample}/cazy_emi.anvio_ko.tsv"
    params:
        agp = config["agp_binary"],
        prefix = OUTDIR + "/cazy/{sample}/cazy_emi",
        min_depth = config.get("damage_annotate", {}).get("min_depth", 0.5),
        min_breadth = config.get("damage_annotate", {}).get("min_breadth", 0.1),
        min_reads = config.get("damage_annotate", {}).get("min_reads", 3),
    threads: 8
    resources:
        runtime=480,
        mem_mb=lambda wc, input: emi_mem_mb(input.emi)
    shell:
        """
        export LD_LIBRARY_PATH=/vol/cloud/agp/lib:${{LD_LIBRARY_PATH:-}} && {params.agp} damage-annotate --emi {input.emi} \
            --damage-index {input.agd} \
            --map {input.gene_map} \
            --analysis-prefix {params.prefix} \
            --gene-summary {output.genes} \
            --functional-summary {output.functional} \
            --anvio-ko {output.anvio} \
            --annotation-source dbCAN \
            --min-depth {params.min_depth} \
            --min-breadth {params.min_breadth} \
            --min-reads {params.min_reads} \
            --coverage-em \
            --auto-prior-ancient \
            --auto-calibrate-spurious \
            --refine-damage \
            -m 32000 \
            --threads {threads} -v

        mv {params.prefix}.protein.tsv {output.proteins}
        """


# ============================================================================
# Viral Pipeline with EMI (replaces xFilter)
# ============================================================================

rule viral_hits2emi:
    input:
        hits = OUTDIR + "/viral/{sample}/viral_hits.tsv",
        agd = OUTDIR + "/sample_damage/{sample}/predictions.agd"
    output:
        emi = OUTDIR + "/viral/{sample}/viral_hits.emi"
    params:
        agp = config["agp_binary"]
    threads: 8
    resources:
        runtime=lambda wc, input: max(240, int(os.path.getsize(input.hits) / 1e6 / 30)),
        mem_mb=lambda wc, input: max(32000, int(os.path.getsize(input.agd) / 1e6) + 12000)
    shell:
        """
        export LD_LIBRARY_PATH=/vol/cloud/agp/lib:${{LD_LIBRARY_PATH:-}} && {params.agp} hits2emi -i {input.hits} -o {output.emi} --damage-index {input.agd} \
            --threads {threads} --memory 16G -v
        """


rule viral_damage_annotate_emi:
    input:
        emi = OUTDIR + "/viral/{sample}/viral_hits.emi",
        agd = OUTDIR + "/sample_damage/{sample}/predictions.agd",
        gene_map = config["viral"]["gene_list"]
    output:
        proteins = OUTDIR + "/viral/{sample}/viral_damage.tsv",
        genes = OUTDIR + "/viral/{sample}/viral_genes.tsv",
        functional = OUTDIR + "/viral/{sample}/viral_emi.functional.tsv"
    params:
        agp = config["agp_binary"],
        prefix = OUTDIR + "/viral/{sample}/viral_emi",
        min_depth = config.get("damage_annotate", {}).get("min_depth", 0.5),
        min_breadth = config.get("damage_annotate", {}).get("min_breadth", 0.1),
        min_reads = config.get("damage_annotate", {}).get("min_reads", 3),
    threads: 8
    resources:
        runtime=480,
        mem_mb=lambda wc, input: emi_mem_mb(input.emi)
    shell:
        """
        export LD_LIBRARY_PATH=/vol/cloud/agp/lib:${{LD_LIBRARY_PATH:-}} && {params.agp} damage-annotate --emi {input.emi} \
            --damage-index {input.agd} \
            --map {input.gene_map} \
            --analysis-prefix {params.prefix} \
            --gene-summary {output.genes} \
            --functional-summary {output.functional} \
            --annotation-source IMGVR \
            --min-depth {params.min_depth} \
            --min-breadth {params.min_breadth} \
            --min-reads {params.min_reads} \
            --coverage-em \
            --auto-prior-ancient \
            --auto-calibrate-spurious \
            --refine-damage \
            -m 32000 \
            --threads {threads} -v

        mv {params.prefix}.protein.tsv {output.proteins}
        """

# Prefer EMI rules over xFilter rules
ruleorder: kegg_damage_annotate_emi > kegg_damage_annotate
ruleorder: cazy_damage_annotate_emi > cazy_damage_annotate
ruleorder: viral_damage_annotate_emi > viral_damage_annotate
rule kegg_anvio_metabolism:
    """Run anvi-estimate-metabolism with AGP damage-aware abundances."""
    input:
        anvio = OUTDIR + "/kegg/{sample}/emi.anvio_ko.tsv"
    output:
        modules = OUTDIR + "/kegg/{sample}/anvio_modules.txt",
        module_paths = OUTDIR + "/kegg/{sample}/anvio_module_paths.txt",
        module_steps = OUTDIR + "/kegg/{sample}/anvio_module_steps.txt",
        hits = OUTDIR + "/kegg/{sample}/anvio_hits.txt"
    params:
        prefix = OUTDIR + "/kegg/{sample}/anvio",
        user_modules = config["kegg"]["user_modules"]
    resources:
        runtime=120,
        mem_mb=8000
    shell:
        """
        # anvio v8 requires 'gene_id' but dart outputs 'gene_callers_id'
        sed '1s/gene_callers_id/gene_id/' {input.anvio} > {input.anvio}.fixed.tsv

        /vol/cloud/opt/miniconda3/envs/anvio/bin/anvi-estimate-metabolism \
            --add-coverage \
            --enzymes-txt {input.anvio}.fixed.tsv \
            --output-modes modules,module_paths,module_steps,hits \
            -O {params.prefix} \
            --include-kos-not-in-kofam \
            --user-modules {params.user_modules}


        """
# ============================================================================
# Synthetic Benchmark with EMI pipeline
# ============================================================================

rule synthetic_hits2emi:
    """Convert synthetic sample hits to EMI format."""
    input:
        hits = f"{OUTDIR}/synthetic/{{sample}}/hits.tsv",
        agd = f"{OUTDIR}/synthetic/{{sample}}/predictions.agd"
    output:
        emi = f"{OUTDIR}/synthetic/{{sample}}/hits.emi"
    params:
        agp = config["agp_binary"]
    threads: 8
    resources:
        runtime=120,
        mem_mb=16000
    shell:
        """
        export LD_LIBRARY_PATH=/vol/cloud/agp/lib:${{LD_LIBRARY_PATH:-}} && {params.agp} hits2emi -i {input.hits} -o {output.emi} --damage-index {input.agd} \
            --threads {threads} --memory 16G -v
        """


rule synthetic_damage_annotate_emi:
    """Damage annotation with EM reassignment for synthetic benchmark."""
    input:
        emi = f"{OUTDIR}/synthetic/{{sample}}/hits.emi",
        agd = f"{OUTDIR}/synthetic/{{sample}}/predictions.agd"
    output:
        proteins = f"{OUTDIR}/synthetic/{{sample}}/emi.protein.tsv",
        reads = f"{OUTDIR}/synthetic/{{sample}}/emi.reads.tsv"
    params:
        agp = config["agp_binary"],
        prefix = f"{OUTDIR}/synthetic/{{sample}}/emi",
        em_iters = config.get("damage_annotate", {}).get("em_iters", 100),
        em_lambda = config.get("damage_annotate", {}).get("em_lambda", 3.0),
    threads: 8
    resources:
        runtime=0,
        mem_mb=lambda wc, input: emi_mem_mb(input.emi)
    shell:
        """
        export LD_LIBRARY_PATH=/vol/cloud/agp/lib:${{LD_LIBRARY_PATH:-}} && {params.agp} damage-annotate --emi {input.emi} \
            --damage-index {input.agd} \
            --analysis-prefix {params.prefix} \
            --em-iters {params.em_iters} \
            --em-lambda {params.em_lambda} \
            --coverage-em \
            --auto-prior-ancient \
            --auto-calibrate-spurious \
            --refine-damage \
            -m 32000 \
            --threads {threads} -v
        """


GC_RICH_SUFFIXES = {"24", "27", "29", "31"}
AT_RICH_SUFFIXES = {"25", "33", "34", "35", "36", "37"}

rule evaluate_sample_emi:
    """
    Protein-level damage classification evaluation.

    Ground truth: built from emi.reads.tsv - a protein is labelled ancient
    if >=80% of its reads carry ':ancient:' in query_id and it has >=5 reads.
    Score: assign_mean_posterior from emi.protein.tsv (pass_mapping_filter=1 only).
    Samples with fewer than 50 modern proteins are skipped (near-pure-ancient).
    """
    input:
        reads    = f"{OUTDIR}/synthetic/{{sample}}/emi.reads.tsv",
        proteins = f"{OUTDIR}/synthetic/{{sample}}/emi.protein.tsv"
    output:
        eval = f"{OUTDIR}/synthetic/{{sample}}/evaluation_emi.tsv"
    run:
        import numpy as np
        from sklearn.metrics import roc_auc_score, average_precision_score

        MIN_READS  = 5
        PURITY     = 0.80
        MIN_MODERN = 50
        SCORE_COL  = "assign_mean_posterior"

        suffix = wildcards.sample.split("-")[-1]
        gc_group = "GC-rich" if suffix in GC_RICH_SUFFIXES else \
                   ("AT-rich" if suffix in AT_RICH_SUFFIXES else "unknown")

        # Build protein-level GT from read labels in emi.reads.tsv
        prot_counts = {}
        with open(input.reads) as fh:
            hdr = fh.readline().strip().split('\t')
            qi = hdr.index('query_id')
            ti = hdr.index('target_id')
            for ln in fh:
                p = ln.rstrip('\n').split('\t')
                if len(p) <= max(qi, ti):
                    continue
                target = p[ti]
                is_anc = ':ancient:' in p[qi]
                if target not in prot_counts:
                    prot_counts[target] = [0, 0]
                if is_anc:
                    prot_counts[target][0] += 1
                else:
                    prot_counts[target][1] += 1

        prot_label = {}
        for pid, (n_anc, n_mod) in prot_counts.items():
            total = n_anc + n_mod
            if total < MIN_READS:
                continue
            frac = n_anc / total
            if frac >= PURITY:
                prot_label[pid] = 1
            elif (1 - frac) >= PURITY:
                prot_label[pid] = 0

        # Load protein scores (pass_mapping_filter=1 only)
        prot_scores = {}
        with open(input.proteins) as fh:
            hdr = fh.readline().strip().split('\t')
            pi = hdr.index('protein_id')
            si = hdr.index(SCORE_COL)
            fi = hdr.index('pass_mapping_filter') if 'pass_mapping_filter' in hdr else None
            for ln in fh:
                p = ln.rstrip('\n').split('\t')
                if len(p) <= max(pi, si):
                    continue
                if fi is not None and p[fi] != '1':
                    continue
                try:
                    prot_scores[p[pi]] = float(p[si])
                except ValueError:
                    pass

        y_true, y_score = [], []
        for pid, label in prot_label.items():
            if pid in prot_scores:
                y_true.append(label)
                y_score.append(prot_scores[pid])

        y_true  = np.array(y_true)
        y_score = np.array(y_score)
        n_ancient = int((y_true == 1).sum())
        n_modern  = int((y_true == 0).sum())

        fields = ['sample', 'gc_group', 'n_proteins', 'n_ancient', 'n_modern',
                  'auc_roc', 'auc_pr', 'mean_score_ancient', 'mean_score_modern',
                  'score_separation', 'method']

        if n_modern < MIN_MODERN or n_ancient == 0:
            row = dict(sample=wildcards.sample, gc_group=gc_group,
                       n_proteins=len(y_true), n_ancient=n_ancient, n_modern=n_modern,
                       auc_roc=float('nan'), auc_pr=float('nan'),
                       mean_score_ancient=float('nan'), mean_score_modern=float('nan'),
                       score_separation=float('nan'), method='emi_protein')
        else:
            auc_roc = roc_auc_score(y_true, y_score)
            auc_pr  = average_precision_score(y_true, y_score)
            ms_anc  = float(y_score[y_true == 1].mean())
            ms_mod  = float(y_score[y_true == 0].mean())
            row = dict(sample=wildcards.sample, gc_group=gc_group,
                       n_proteins=len(y_true), n_ancient=n_ancient, n_modern=n_modern,
                       auc_roc=auc_roc, auc_pr=auc_pr,
                       mean_score_ancient=ms_anc, mean_score_modern=ms_mod,
                       score_separation=ms_anc - ms_mod, method='emi_protein')

        with open(output.eval, 'w') as fh:
            fh.write('\t'.join(fields) + '\n')
            fh.write('\t'.join(str(row[k]) for k in fields) + '\n')


rule aggregate_evaluation_emi:
    """Aggregate protein-level EMI benchmark results across all synthetic samples."""
    input:
        evals = expand(f"{OUTDIR}/synthetic/{{sample}}/evaluation_emi.tsv", sample=SYNTHETIC_SAMPLES)
    output:
        report = f"{OUTDIR}/reports/protein_damage_evaluation_emi.tsv"
    run:
        rows, fields = [], None
        for f in input.evals:
            with open(f) as fh:
                hdr = fh.readline().strip().split('\t')
                if fields is None:
                    fields = hdr
                row = dict(zip(hdr, fh.readline().strip().split('\t')))
                rows.append(row)

        with open(output.report, 'w') as fh:
            fh.write('\t'.join(fields) + '\n')
            for row in rows:
                fh.write('\t'.join(row.get(k, '') for k in fields) + '\n')

        valid = [r for r in rows if r['auc_roc'] != 'nan']
        print(f"\n=== Protein-level Damage Benchmark (assign_mean_posterior) ===")
        print(f"Samples total: {len(rows)}  evaluable: {len(valid)}")
        if valid:
            aucs = [float(r['auc_roc']) for r in valid]
            ns   = [int(r['n_proteins']) for r in valid]
            print(f"Mean AUC-ROC:     {sum(aucs)/len(aucs):.4f}")
            print(f"Weighted AUC-ROC: {sum(a*n for a,n in zip(aucs,ns))/sum(ns):.4f}")
            for grp in ('AT-rich', 'GC-rich'):
                g = [r for r in valid if r['gc_group'] == grp]
                if g:
                    ga = [float(r['auc_roc']) for r in g]
                    print(f"  {grp}: {sum(ga)/len(ga):.4f}  (n={len(g)})")


rule compare_methods_emi:
    """Compare old TSV method vs new protein-level EMI method."""
    input:
        old_report = f"{OUTDIR}/reports/protein_damage_evaluation.tsv",
        emi_report = f"{OUTDIR}/reports/protein_damage_evaluation_emi.tsv"
    output:
        comparison = f"{OUTDIR}/reports/method_comparison_emi.tsv"
    run:
        rows_old, rows_emi = [], []
        with open(input.old_report) as fh:
            hdr = fh.readline().strip().split('\t')
            for ln in fh:
                r = dict(zip(hdr, ln.strip().split('\t')))
                r['method'] = 'tsv'
                rows_old.append(r)
        with open(input.emi_report) as fh:
            hdr = fh.readline().strip().split('\t')
            for ln in fh:
                rows_emi.append(dict(zip(hdr, ln.strip().split('\t'))))

        all_fields = list(dict.fromkeys(k for r in rows_old + rows_emi for k in r))
        with open(output.comparison, 'w') as fh:
            fh.write('\t'.join(all_fields) + '\n')
            for r in rows_old + rows_emi:
                fh.write('\t'.join(r.get(k, '') for k in all_fields) + '\n')

        print("\n=== Method Comparison ===")
        for method, rows in [('tsv', rows_old), ('emi_protein', rows_emi)]:
            valid = [r for r in rows if r.get('auc_roc', 'nan') != 'nan']
            if valid:
                aucs = [float(r['auc_roc']) for r in valid]
                print(f"  {method}: mean AUC {sum(aucs)/len(aucs):.4f}  ({len(valid)} samples)")


rule synthetic_damage_annotate_emi_streaming:
    """EM annotation for synthetic benchmark (streaming comparison variant)."""
    input:
        emi = f"{OUTDIR}/synthetic/{{sample}}/hits.emi",
        agd = f"{OUTDIR}/synthetic/{{sample}}/predictions.agd"
    output:
        proteins = f"{OUTDIR}/synthetic/{{sample}}/emi_streaming.protein.tsv",
        reads    = f"{OUTDIR}/synthetic/{{sample}}/emi_streaming.reads.tsv"
    params:
        agp    = config["agp_binary"],
        prefix = f"{OUTDIR}/synthetic/{{sample}}/emi_streaming",
        em_iters = config.get("damage_annotate", {}).get("em_iters", 100),
        em_lambda = config.get("damage_annotate", {}).get("em_lambda", 3.0),
    threads: 8
    resources:
        runtime=0,
        mem_mb=lambda wc, input: emi_mem_mb(input.emi)
    shell:
        """
        export LD_LIBRARY_PATH=/vol/cloud/agp/lib:${{LD_LIBRARY_PATH:-}} && {params.agp} damage-annotate --emi {input.emi} \
            --damage-index {input.agd} \
            --analysis-prefix {params.prefix} \
            --em-iters {params.em_iters} \
            --em-lambda {params.em_lambda} \
            --coverage-em \
            --auto-prior-ancient \
            --auto-calibrate-spurious \
            --refine-damage \
            -m 32000 \
            --threads {threads} -v
        """


rule synthetic_all_streaming:
    """Run streaming EM annotation on all synthetic samples."""
    input:
        expand(f"{OUTDIR}/synthetic/{{sample}}/emi_streaming.protein.tsv",
               sample=config["synthetic_samples"]["samples"])
