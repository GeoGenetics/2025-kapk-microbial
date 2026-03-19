suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(openxlsx)
    library(janitor)
})

setwd("/maps/projects/fernandezguerra/apps/repos/2025-kapk-microbial/analysis")

# Sample metadata
kapk_cdata     <- read_tsv("./data/cdata/KapK-cdata-manuscript-20221211.tsv",
                           show_col_types = FALSE)
samples_to_keep <- read_tsv("./data/cdata/KapK-cdata-manuscript-20221211-labels-nobloom.tsv",
                             show_col_types = FALSE)

sample_meta_dart <- kapk_cdata |>
    filter(figure_names %in% samples_to_keep$label) |>
    mutate(
        sample      = substr(label, 1, 10),
        short_label = paste(site, sub(".*-", "", figure_names), sep = "_")
    ) |>
    select(sample, short_label, label)

# ── Build KEGG per-sample module table from raw per-sample files ──────────────
#
# For each sample:
#   anvio_modules.txt  → module metadata + completeness + coverage
#   emi.functional.tsv → KO-level damage stats (level == "group")
#   anvio_hits.txt     → KO → module(s) mapping (modules_with_enzyme)
#
# KO stats are aggregated to module level via the KO→module mapping.

build_kegg_sample <- function(sample_dir) {
    sample <- basename(sample_dir)

    # Module completeness + metadata + coverage
    anvio <- fread(file.path(sample_dir, "anvio_modules.txt"),
                   showProgress = FALSE)
    anvio_mod <- anvio[, .(
        module,
        module_name    = module_name,
        module_class   = module_class,
        module_category    = module_category,
        module_subcategory = module_subcategory,
        stepwise_completeness = stepwise_module_completeness,
        pathwise_completeness = pathwise_module_completeness,
        prop_unique_enzymes   = proportion_unique_enzymes_present,
        avg_coverage   = emi_anvio_ko_tsv_fixed_tsv_avg_coverage,
        avg_detection  = emi_anvio_ko_tsv_fixed_tsv_avg_detection
    )]

    # KO → module mapping (one row per KO-module pair)
    hits <- fread(file.path(sample_dir, "anvio_hits.txt"),
                  showProgress = FALSE)
    ko_mod <- hits[modules_with_enzyme != ""][
        , .(module = unlist(strsplit(modules_with_enzyme, ","), use.names = FALSE)),
        by = .(ko = enzyme)
    ] |> unique()

    # KO-level damage stats
    emi <- fread(file.path(sample_dir, "emi.functional.tsv"),
                 showProgress = FALSE)
    ko_dmg <- emi[, .(
        ko              = function_id,
        n_genes, n_reads, n_ancient, ancient_frac,
        ci_low, ci_high,
        mean_posterior,
        n_damaged_genes, damaged_gene_frac, damage_enrichment,
        coverage_mean, avg_identity, avg_read_length
    )]

    # Join KO damage with KO→module mapping
    merged <- merge(ko_mod, ko_dmg, by = "ko", all.x = FALSE)

    # Aggregate KO stats to module level (weighted by n_genes)
    mod_dmg <- merged[, .(
        n_ko              = .N,
        n_genes           = sum(n_genes),
        n_reads           = sum(n_reads),
        n_ancient         = sum(n_ancient),
        ancient_frac      = weighted.mean(ancient_frac, n_genes, na.rm = TRUE),
        ci_low            = weighted.mean(ci_low, n_genes, na.rm = TRUE),
        ci_high           = weighted.mean(ci_high, n_genes, na.rm = TRUE),
        mean_posterior    = weighted.mean(mean_posterior, n_genes, na.rm = TRUE),
        n_damaged_genes   = sum(n_damaged_genes),
        damaged_gene_frac = weighted.mean(damaged_gene_frac, n_genes, na.rm = TRUE),
        damage_enrichment = weighted.mean(damage_enrichment, n_genes, na.rm = TRUE),
        coverage_mean     = weighted.mean(coverage_mean, n_genes, na.rm = TRUE),
        avg_identity      = weighted.mean(avg_identity, n_genes, na.rm = TRUE),
        avg_read_length   = weighted.mean(avg_read_length, n_genes, na.rm = TRUE)
    ), by = module]

    # Join module completeness with damage stats
    result <- merge(anvio_mod, mod_dmg, by = "module", all.x = TRUE)
    result[, sample := sample]
    result
}

kegg_dirs <- list.dirs("./results/functional_agp/kegg",
                        recursive = FALSE, full.names = TRUE)
kegg_dirs <- kegg_dirs[file.exists(file.path(kegg_dirs, "anvio_modules.txt"))]

cat("Building KEGG table from", length(kegg_dirs), "samples...\n")
kegg_all <- rbindlist(lapply(kegg_dirs, build_kegg_sample), fill = TRUE) |>
    as_tibble() |>
    inner_join(sample_meta_dart, by = "sample") |>
    select(
        short_label, label_orig = label,
        module, module_name, module_class, module_category, module_subcategory,
        stepwise_completeness, pathwise_completeness, prop_unique_enzymes,
        avg_coverage, avg_detection,
        n_ko, n_genes, n_reads, n_ancient, ancient_frac,
        ci_low, ci_high,
        mean_posterior, n_damaged_genes, damaged_gene_frac, damage_enrichment,
        coverage_mean, avg_identity, avg_read_length
    )

sup_table_5_s1 <- kegg_all
sup_table_5_s2 <- kegg_all |> filter(pathwise_completeness >= 0.75)

# ── Build CAZy per-sample table from raw per-sample files ─────────────────────
#
# cazy_emi.functional.tsv is already at family level (function_id = family)

build_cazy_sample <- function(sample_dir) {
    sample <- basename(sample_dir)
    f <- file.path(sample_dir, "cazy_emi.functional.tsv")
    if (!file.exists(f)) return(NULL)
    d <- fread(f, showProgress = FALSE)
    d[, sample := sample]
    d
}

cazy_dirs <- list.dirs("./results/functional_agp/cazy",
                        recursive = FALSE, full.names = TRUE)

cat("Building CAZy table from", length(cazy_dirs), "samples...\n")
cazy_all <- rbindlist(lapply(cazy_dirs, build_cazy_sample), fill = TRUE) |>
    as_tibble() |>
    rename(family = function_id) |>
    inner_join(sample_meta_dart, by = "sample") |>
    select(
        short_label, label_orig = label,
        family, db,
        n_genes, n_reads, n_ancient, n_modern, ancient_frac,
        ci_low, ci_high,
        mean_posterior, n_damaged_genes, damaged_gene_frac, damage_enrichment,
        coverage_mean, avg_identity, avg_read_length
    )

sup_table_5_s3 <- cazy_all
sup_table_5_s4 <- cazy_all |> filter(mean_posterior >= 0.7)

# ── Write sup_table_5 ─────────────────────────────────────────────────────────

wb5 <- createWorkbook()
addWorksheet(wb5, "S13 - KEGG detected")
addWorksheet(wb5, "S14 - KEGG complete")
addWorksheet(wb5, "S15 - CAZy all")
addWorksheet(wb5, "S16 - CAZy authenticated")

writeData(wb5, "S13 - KEGG detected",  sup_table_5_s1 |> clean_names(case = "sentence"))
writeData(wb5, "S14 - KEGG complete",  sup_table_5_s2 |> clean_names(case = "sentence"))
writeData(wb5, "S15 - CAZy all",       sup_table_5_s3 |> clean_names(case = "sentence"))
writeData(wb5, "S16 - CAZy authenticated", sup_table_5_s4 |> clean_names(case = "sentence"))

saveWorkbook(wb5, "../supp-tab-v2/sup_table_5.xlsx", overwrite = TRUE)
cat("Saved sup_table_5.xlsx\n")
cat("  S13 KEGG detected:", nrow(sup_table_5_s1), "rows,", ncol(sup_table_5_s1), "cols\n")
cat("  S14 KEGG complete:", nrow(sup_table_5_s2), "rows\n")
cat("  S15 CAZy all:     ", nrow(sup_table_5_s3), "rows,", ncol(sup_table_5_s3), "cols\n")
cat("  S16 CAZy auth:    ", nrow(sup_table_5_s4), "rows\n")

# ── sup_table_6: Viral ────────────────────────────────────────────────────────

viral_files   <- list.files("./results/functional_agp/viral",
                            pattern = "^viral_emi\\.functional\\.tsv$",
                            recursive = TRUE, full.names = TRUE)
viral_samples <- basename(dirname(viral_files))

viral_emi_raw <- rbindlist(mapply(function(f, s) {
    d <- fread(f, showProgress = FALSE)[level == "group"]
    d[, sample := s]
    d
}, viral_files, viral_samples, SIMPLIFY = FALSE), fill = TRUE) |>
    as_tibble() |>
    rename(reference = function_id) |>
    inner_join(sample_meta_dart, by = "sample")

sup_table_6_s1 <- viral_emi_raw |>
    filter(mean_posterior >= 0.7) |>
    select(short_label, label_orig = label, reference, n_genes, n_reads,
           n_ancient, ancient_frac, mean_posterior, coverage_mean, avg_identity)

imgvr <- fread("./data/cdata/IMGVR_all_Sequence_information-high_confidence.tsv",
               showProgress = FALSE,
               select = c("UVIG", "Ecosystem classification",
                          "Gene content (total genes;cds;tRNA;geNomad marker)")) |>
    clean_names() |>
    mutate(n_cds = as.integer(sub("^[^;]+;([^;]+);.*", "\\1",
        gene_content_total_genes_cds_t_rna_ge_nomad_marker))) |>
    select(reference = uvig, n_cds, ecosystem = ecosystem_classification)

sup_table_6_s2 <- imgvr |> filter(reference %in% sup_table_6_s1$reference)

wb6 <- createWorkbook()
addWorksheet(wb6, "S17 - Viral references")
addWorksheet(wb6, "S18 - IMGVR annotation")

writeData(wb6, "S17 - Viral references", sup_table_6_s1 |> clean_names(case = "sentence"))
writeData(wb6, "S18 - IMGVR annotation", sup_table_6_s2 |> clean_names(case = "sentence"))

saveWorkbook(wb6, "../supp-tab-v2/sup_table_6.xlsx", overwrite = TRUE)
cat("Saved sup_table_6.xlsx\n")
cat("  S17 Viral references:", nrow(sup_table_6_s1), "rows\n")
cat("  S18 IMGVR annotation:", nrow(sup_table_6_s2), "rows\n")
