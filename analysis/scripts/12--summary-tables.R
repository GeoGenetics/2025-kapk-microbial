library(tidyverse)
library(readxl)
library(janitor)
library(openxlsx)
library(data.table)

# ── Metadata ──────────────────────────────────────────────────────────────────

kapk_cdata       <- read_tsv("./data/cdata/KapK-cdata-manuscript-20221211.tsv")
kapk_cdata_agg   <- read_tsv("./data/cdata/KapK-cdata-agg-manuscript-20221211.tsv")
samples_to_keep  <- read_tsv("./data/cdata/KapK-cdata-manuscript-20221211-labels-nobloom.tsv")
label_to_file    <- read_tsv("./data/cdata/KapK-label-to-file-20221211.tsv")
control_samples  <- read_tsv("./data/cdata/KapK_samples-controls-20221211.tsv")
dmg_thresholds   <- read_tsv("./results/damage/dmg_thresholds.tsv")

tax_info <- fread("./data/taxonomy/hires-organelles-viruses-arctic.tax.tsv",
                  col.names = c("reference", "tax_string"), showProgress = FALSE) |>
    as_tibble() |>
    mutate(
        domain  = str_extract(tax_string, "d__[^;]+"),
        phylum  = str_extract(tax_string, "p__[^;]+"),
        class   = str_extract(tax_string, "c__[^;]+"),
        order   = str_extract(tax_string, "o__[^;]+"),
        family  = str_extract(tax_string, "f__[^;]+"),
        genus   = str_extract(tax_string, "g__[^;]+"),
        species = str_extract(tax_string, "s__[^;]+")
    ) |>
    select(-tax_string)

sample_meta_dart <- kapk_cdata |>
    filter(figure_names %in% samples_to_keep$label) |>
    mutate(
        sample      = substr(label, 1, 10),
        short_label = paste(site, sub(".*-", "", figure_names), sep = "_")
    ) |>
    select(sample, short_label, label)

# ── Stats ─────────────────────────────────────────────────────────────────────

initial_stats <- read_tsv("./data/stats/all.stats-initial-summary.tsv.gz") |>
    select(label, initial_num_reads = num_seqs)

derep_stats <- read_tsv("./data/stats/all.stats-derep-summary.tsv.gz") |>
    select(label, derep_num_reads = num_seqs)

reads_mapped <- read_tsv("./data/cdata/all.tp-classified-reads.tsv",
                         col_names = c("file", "reads_mapped")) |>
    inner_join(label_to_file)

# ── Supplementary Table 1 ─────────────────────────────────────────────────────
# S1: Sample information
# S2: Control samples
# S3: Uplift model data
# S5: Re-extraction metadata
# S6: Re-extraction taxonomy
# S7: Lipid biomarkers

sup_table_1_s1 <- kapk_cdata |>
    inner_join(initial_stats) |>
    inner_join(derep_stats) |>
    inner_join(reads_mapped) |>
    left_join(samples_to_keep |> mutate(used_in_analysis = "Yes") |> rename(figure_names = label)) |>
    mutate(used_in_analysis = ifelse(is.na(used_in_analysis), "No", "Yes")) |>
    inner_join(kapk_cdata) |>
    select(-label) |>
    mutate(label = paste0(str_split_i(figure_names, "_", i = 1), "_",
                          str_split_i(figure_names, "-", i = -1))) |>
    select(label, label_orig, member_unit, site, figure_names,
           initial_num_reads, derep_num_reads, reads_mapped, used_in_analysis)

control_stats <- initial_stats |>
    inner_join(derep_stats) |>
    inner_join(reads_mapped) |>
    inner_join(control_samples) |>
    select(label_orig, initial_num_reads, derep_num_reads, reads_mapped)

sup_table_1_s1 <- sup_table_1_s1 |>
    bind_rows(control_stats) |>
    mutate(used_in_analysis = case_when(
        grepl("EC", label_orig) ~ "Extraction control",
        grepl("LB", label_orig) ~ "Library blank",
        TRUE ~ used_in_analysis
    ))

# Control taxonomy (S2)
tax_data <- read_tsv("./data/taxonomy/tp-mapping-filtered.summary.tsv.gz") |>
    filter(breadth >= 0.01) |>
    mutate(abundance = ifelse(tax_abund_tad == 0, tax_abund_read, tax_abund_tad)) |>
    inner_join(tax_info)

dmg_local <- read_csv("./data/taxonomy/tp-mdmg.weight-1.csv.gz") |>
    rename(reference = tax_id) |>
    select(label, reference, damage, significance)

sup_table_1_s2 <- control_samples |>
    select(label_orig, label) |>
    inner_join(tax_data) |>
    inner_join(dmg_local) |>
    select(-label)

# Uplift model (S3) — pre-computed TSV
sup_table_1_s3 <- read_tsv("./manuscript/tables/sup_table_1_s3.tsv")

# Re-extraction (S5, S6)
sup_table_1_s5 <- read_tsv("./manuscript/tables/kapk_reextractions-cdata.tsv")
sup_table_1_s6 <- read_tsv("./manuscript/tables/kapk_reextractions-taxonomy.tsv")

# Biomarkers (S7)
sup_table_1_s7 <- read_tsv("./data/biomarkers/kapk-20230925-biomarkers-cleaned.tsv")

sup_table_1 <- createWorkbook()
addWorksheet(sup_table_1, "S1 - Sample information")
addWorksheet(sup_table_1, "S2 - Control samples")
addWorksheet(sup_table_1, "S3 - Uplift model data")
addWorksheet(sup_table_1, "S5 - Re-extraction metadata")
addWorksheet(sup_table_1, "S6 - Re-extraction taxonomy")
addWorksheet(sup_table_1, "S7 - Lipid biomarkers")

writeData(sup_table_1, "S1 - Sample information",   sup_table_1_s1 |> clean_names(case = "sentence"))
writeData(sup_table_1, "S2 - Control samples",       sup_table_1_s2 |> clean_names(case = "sentence"))
writeData(sup_table_1, "S3 - Uplift model data",     sup_table_1_s3 |> clean_names(case = "sentence"))
writeData(sup_table_1, "S5 - Re-extraction metadata",sup_table_1_s5 |> clean_names(case = "sentence"))
writeData(sup_table_1, "S6 - Re-extraction taxonomy",sup_table_1_s6 |> clean_names(case = "sentence"))
writeData(sup_table_1, "S7 - Lipid biomarkers",      sup_table_1_s7 |> clean_names(case = "sentence"))

saveWorkbook(sup_table_1, file = "../supp-tab-v2/sup_table_1.xlsx", overwrite = TRUE)

# ── Supplementary Table 2 ─────────────────────────────────────────────────────
# S8: Source samples
# S9: Source contributions

source_cdata <- read_tsv("./data/sourcetracker/cdata/kapk-biomes__combined.tsv")

st_labels <- read_tsv("./results/sourcetracker/st-biome-gm.map") |>
    clean_names() |>
    filter(source_sink == "source") |>
    pull(number_sample_id)

sup_table_2_s1 <- source_cdata |> filter(run_accession %in% st_labels)

st_results <- read_tsv("./results/sourcetracker/mixing_proportions.txt") |>
    rename(figure_names = "#SampleID")

st_results_long <- st_results |>
    pivot_longer(cols = -figure_names, names_to = "biome", values_to = "proportion")

sup_table_2_s3 <- kapk_cdata |>
    select(label = figure_names) |>
    inner_join(kapk_cdata_agg) |>
    select(short_label, figure_names = label) |>
    inner_join(st_results_long) |>
    rename(label_orig = figure_names)

sup_table_2 <- createWorkbook()
addWorksheet(sup_table_2, "S8 - Source samples")
addWorksheet(sup_table_2, "S9 - Source contributions")

writeData(sup_table_2, "S8 - Source samples",      sup_table_2_s1 |> clean_names(case = "sentence"))
writeData(sup_table_2, "S9 - Source contributions", sup_table_2_s3 |> clean_names(case = "sentence"))

saveWorkbook(sup_table_2, file = "../supp-tab-v2/sup_table_2.xlsx", overwrite = TRUE)

# ── Supplementary Table 3 ─────────────────────────────────────────────────────
# S10: GEM + tg2g MAG mapping
# S11: Woodcroft mapping

mg_mapping <- read_tsv("./data/mag-distribution/GEM-20220926__tg2g-20220926-tp-mapping-filtered.summary.tsv.gz") |>
    filter(breadth >= 0.01) |>
    inner_join(
        read_csv("./data/mag-distribution/GEM-20220926__tg2g-20220926-tp-mdmg.weight-1.csv.gz") |>
            filter(label %in% kapk_cdata$label) |>
            rename(reference = tax_id) |>
            filter(significance > dmg_thresholds$signf, damage >= dmg_thresholds$damage)
    )

mg_mapping <- kapk_cdata |>
    select(label, figure_names) |>
    inner_join(mg_mapping) |>
    select(-label) |>
    rename(label_orig = figure_names)

sup_table_3_s1 <- kapk_cdata_agg |>
    select(short_label, label) |>
    distinct() |>
    inner_join(mg_mapping) |>
    rename(label_orig = label)

woodcroft_mapping <- read_tsv("./data/mag-distribution/woodcroft2018-tp-mapping-filtered.summary.tsv.gz") |>
    filter(breadth >= 0.01) |>
    inner_join(
        read_csv("./data/mag-distribution/woodcroft2018-tp-mdmg.weight-1.csv.gz") |>
            filter(label %in% kapk_cdata$label) |>
            rename(reference = tax_id) |>
            filter(significance > dmg_thresholds$signf, damage >= dmg_thresholds$damage)
    )

woodcroft_mapping <- kapk_cdata |>
    select(label, figure_names) |>
    inner_join(woodcroft_mapping) |>
    select(-label) |>
    rename(label_orig = figure_names)

sup_table_3_s2 <- kapk_cdata_agg |>
    select(short_label, label) |>
    distinct() |>
    inner_join(woodcroft_mapping) |>
    rename(label_orig = label)

sup_table_3 <- createWorkbook()
addWorksheet(sup_table_3, "S10 - GEM mapping")
addWorksheet(sup_table_3, "S11 - Woodcroft mapping")

writeData(sup_table_3, "S10 - GEM mapping",      sup_table_3_s1 |> clean_names(case = "sentence"))
writeData(sup_table_3, "S11 - Woodcroft mapping", sup_table_3_s2 |>
              select(names(sup_table_3_s1)) |> clean_names(case = "sentence"))

saveWorkbook(sup_table_3, file = "../supp-tab-v2/sup_table_3.xlsx", overwrite = TRUE)

# ── Supplementary Table 7 (Briggs) ────────────────────────────────────────────

sup_table_7_s1 <- read_tsv("./manuscript/tables/briggs_rl-percid.tsv.gz")
sup_table_7_s2 <- read_tsv("./manuscript/tables/briggs_rl-ancient-isrecal1.tsv.gz")

sup_table_7 <- createWorkbook()
addWorksheet(sup_table_7, "S1 - Briggs RL Percid")
addWorksheet(sup_table_7, "S2 - Briggs RL Ancient Isrecal1")

writeData(sup_table_7, "S1 - Briggs RL Percid",          sup_table_7_s1)
writeData(sup_table_7, "S2 - Briggs RL Ancient Isrecal1", sup_table_7_s2)

saveWorkbook(sup_table_7, file = "../supp-tab-v2/sup_table_7.xlsx", overwrite = TRUE)

# ── Supplementary Table 8 (Simulations) ──────────────────────────────────────

sim_cdata <- kapk_cdata_agg |>
    select(short_label, label) |>
    inner_join(
        kapk_cdata |>
            inner_join(label_to_file) |>
            select(label, figure_names, file) |>
            select(-label, label = figure_names, Sample = file)
    ) |>
    select(Label = short_label, label_orig = label, Sample)

sup_table_8_s1 <- sim_cdata |>
    inner_join(read_tsv("./manuscript/tables/bowtie2-benchmarks.tsv")) |>
    select(-Sample)

sup_table_8_s2 <- sim_cdata |>
    inner_join(read_tsv("./manuscript/tables/saturation.tsv") |> rename(Sample = Label)) |>
    select(-Sample)

sup_table_8_s3 <- read_tsv("./manuscript/tables/simulations-detection-stats-composition.tsv")
sup_table_8_s4 <- read_tsv("./manuscript/tables/simulations-abundance-mdae.tsv")
sup_table_8_s5 <- read_tsv("./manuscript/tables/simulations-damage-mdae.tsv")

sup_table_8 <- createWorkbook()
addWorksheet(sup_table_8, "S18 - Bowtie2 benchmarks")
addWorksheet(sup_table_8, "S19 - Saturation")
addWorksheet(sup_table_8, "S20 - Simulations composition")
addWorksheet(sup_table_8, "S21 - Simulations abundance")
addWorksheet(sup_table_8, "S22 - Simulations damage")

writeData(sup_table_8, "S18 - Bowtie2 benchmarks",     sup_table_8_s1 |> clean_names(case = "sentence"))
writeData(sup_table_8, "S19 - Saturation",              sup_table_8_s2 |> clean_names(case = "sentence"))
writeData(sup_table_8, "S20 - Simulations composition", sup_table_8_s3 |> clean_names(case = "sentence"))
writeData(sup_table_8, "S21 - Simulations abundance",   sup_table_8_s4 |> clean_names(case = "sentence"))
writeData(sup_table_8, "S22 - Simulations damage",      sup_table_8_s5 |> clean_names(case = "sentence"))

saveWorkbook(sup_table_8, file = "../supp-tab-v2/sup_table_8.xlsx", overwrite = TRUE)
