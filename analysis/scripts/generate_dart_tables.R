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

# ── sup_table_5: KEGG + CAZy ─────────────────────────────────────────────────

kegg_dart <- read_tsv("./results/functional_agp/kegg_module_damage.tsv",
                      show_col_types = FALSE) |>
    inner_join(sample_meta_dart, by = "sample")

sup_table_5_s1 <- kegg_dart |>
    select(short_label, label_orig = label, module, module_name, module_class,
           module_category, module_subcategory, avg_coverage, n_proteins,
           n_damaged, damage_rate, mean_p_damaged, mean_d_aa)

sup_table_5_s2 <- sup_table_5_s1 |> filter(mean_p_damaged >= 0.7)

cazy_dart <- read_tsv("./results/functional_agp/cazy_family_damage.tsv",
                      show_col_types = FALSE) |>
    inner_join(sample_meta_dart, by = "sample")

sup_table_5_s3 <- cazy_dart |>
    select(short_label, label_orig = label, family, cazy_class, deg_group,
           n_proteins, n_damaged, damage_rate, mean_p_damaged, mean_d_aa,
           coverage_mean)

sup_table_5_s4 <- sup_table_5_s3 |> filter(mean_p_damaged >= 0.7)

wb5 <- createWorkbook()
addWorksheet(wb5, "S13 - KEGG all")
addWorksheet(wb5, "S14 - KEGG authenticated")
addWorksheet(wb5, "S15 - CAZy all")
addWorksheet(wb5, "S16 - CAZy authenticated")

writeData(wb5, "S13 - KEGG all",           sup_table_5_s1 |> clean_names(case = "sentence"))
writeData(wb5, "S14 - KEGG authenticated", sup_table_5_s2 |> clean_names(case = "sentence"))
writeData(wb5, "S15 - CAZy all",           sup_table_5_s3 |> clean_names(case = "sentence"))
writeData(wb5, "S16 - CAZy authenticated", sup_table_5_s4 |> clean_names(case = "sentence"))

saveWorkbook(wb5, "../supp-tab-v2/sup_table_5.xlsx", overwrite = TRUE)
cat("Saved sup_table_5.xlsx\n")
cat("  S13 KEGG all:          ", nrow(sup_table_5_s1), "rows\n")
cat("  S14 KEGG authenticated:", nrow(sup_table_5_s2), "rows\n")
cat("  S15 CAZy all:          ", nrow(sup_table_5_s3), "rows\n")
cat("  S16 CAZy authenticated:", nrow(sup_table_5_s4), "rows\n")

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
               showProgress = FALSE) |>
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
