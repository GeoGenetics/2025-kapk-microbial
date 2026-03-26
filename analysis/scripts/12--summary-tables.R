library(tidyverse)
library(readxl)

# Supp table 1
# - Sheet 1: Samples used and the information associated with them which can be
#   found linked to Kjær et al.
#   - This table includes the sample name, number of original reads, dereplicated reads, mapped
# - Sheet 2: Control data and the taxonomic composition of the control samples
# - Sheet 3: The uplift model data
# - Sheet 4: Taxonomic composition of the samples and damage statistics
# - Sheet 5: Re-extraction data
# - Sheet 6: Methane budget estimates

# Control data


# Sample data
kapk_cdata <- read_tsv("./data/cdata/KapK-cdata-manuscript-20221211.tsv")
kapk_cdata_agg <- read_tsv("./data/cdata/KapK-cdata-agg-manuscript-20221211.tsv")
samples_to_keep <- read_tsv("./data/cdata/KapK-cdata-manuscript-20221211-labels-nobloom.tsv")
label_to_file <- read_tsv("./data/cdata/KapK-label-to-file-20221211.tsv")

initial_stats <- read_tsv("./data/stats/all.stats-initial-summary.tsv.gz") |>
    select(label, initial_num_reads = num_seqs)


derep_stats <- read_tsv("./data/stats/all.stats-derep-summary.tsv.gz") |>
    select(label, derep_num_reads = num_seqs)


# Get number of reads mapped
reads_mapped <- read_tsv("./data/cdata/all.tp-classified-reads.tsv", col_names = c("file", "reads_mapped")) |>
    inner_join(label_to_file)

sup_table_1_s1 <- kapk_cdata |>
    inner_join(initial_stats) |>
    inner_join(derep_stats) |>
    inner_join(reads_mapped) |>
    left_join(samples_to_keep |> mutate(used_in_analysis = "Yes") |> rename(figure_names = label)) |>
    mutate(used_in_analysis = ifelse(is.na(used_in_analysis), "No", "Yes")) |>
    inner_join(kapk_cdata) |>
    select(-label) |>
    mutate(label = paste0(str_split_i(figure_names, "_", i = 1), "_", str_split_i(figure_names, "-", i = -1))) |>
    select(label, label_orig, member_unit, site, figure_names, initial_num_reads, derep_num_reads, reads_mapped, used_in_analysis)

control_samples <- read_tsv("./data/cdata/KapK_samples-controls-20221211.tsv")
# Get some stats

control_stats <- initial_stats |>
    inner_join(derep_stats) |>
    inner_join(reads_mapped) |>
    inner_join(control_samples) |>
    select(label_orig, initial_num_reads, derep_num_reads, reads_mapped)

sup_table_1_s1 <- sup_table_1_s1 |>
    bind_rows(control_stats) |>
    mutate(
        used_in_analysis = case_when(
            grepl("EC", label_orig) ~ "Extraction control",
            grepl("LB", label_orig) ~ "Library blank",
            TRUE ~ used_in_analysis
        )
    )


tax_data <- read_tsv("./data/taxonomy/profiling/standard/k1000/94/tp-mapping-filtered.summary.tsv.gz") |>
    filter(breadth >= 0.01) |>
    mutate(abundance = ifelse(tax_abund_tad == 0, tax_abund_read, tax_abund_tad)) |>
    inner_join(tax_info)

dmg_local <- read_csv("./data/taxonomy/dmg/standard/k1000/94/local/tp-mdmg.weight-1.csv.gz") |>
    rename(reference = tax_id) |>
    select(label, reference, damage, significance)


sup_table_1_s2 <- control_samples |>
    select(label_orig, label) |>
    inner_join(tax_data) |>
    inner_join(dmg_local) |>
    select(-label)



# uplift model data
sup_table_1_s3 <- readxl::read_xlsx("./data/cdata/KapK_uplift_model.xlsx") |>
    clean_names() |>
    select(-x5) |>
    setNames(c("site", "member_unit", "top_masl", "max_sample_depth", "sample_depth_2M", "sample_depth_today", "linear_model", "time_0_masl_ka")) |>
    clean_names() |>
    select(-linear_model)




# taxonomic information and damage statistics

tax_data <- read_tsv("./results/taxonomy/tp-mapping-filtered.nocontam.tax.tsv.gz") |>
    inner_join(kapk_cdata |> select(label, figure_names)) |>
    inner_join(tax_info) |>
    filter(figure_names %in% kapk_cdata_agg$label)

dmg_local <- read_csv("./data/taxonomy/dmg/standard/k1000/94/local/tp-mdmg.weight-1.csv.gz") |>
    rename(reference = tax_id) |>
    inner_join(kapk_cdata |> select(label, figure_names)) |>
    select(-sample) |>
    filter(figure_names %in% kapk_cdata_agg$label) |>
    select(label, figure_names, reference, damage, significance)

sup_table_1_s4 <- kapk_cdata_agg |>
    select(short_label, figure_names = label) |>
    inner_join(tax_data) |>
    inner_join(dmg_local) |>
    select(-label) |>
    rename(label_orig = figure_names, label = short_label)

# Re-extraction data
sup_table_1_s5 <- read_tsv("./manuscript/tables/kapk_reextractions-cdata.tsv")

sup_table_1_s6 <- read_tsv("./manuscript/tables/kapk_reextractions-taxonomy.tsv")

sup_table_1_s7 <- read_tsv("./data/biomarkers/kapk-20230925-biomarkers-cleaned.tsv")

# sourcetracker data
st_labels <- read_tsv("./results/sourcetracker/st-biome-gm.map") |>
    clean_names() |>
    filter(source_sink == "source") |>
    pull(number_sample_id)

sup_table_2_s1 <- source_cdata |> filter(run_accession %in% st_labels)
sup_table_2_s2 <- read_tsv("./data/sourcetracker/standard/k1000/95/tp-mapping-filtered.summary.tsv.gz") |>
    inner_join(tax_info) |>
    filter(label %in% st_labels)

st_results <- read_tsv("./results/sourcetracker/mixing_proportions.txt") |>
    rename(figure_names = "#SampleID")

st_results_sd <- read_tsv("./results/sourcetracker/mixing_proportions_stds.txt") |>
    rename(figure_names = "#SampleID")

st_results_long <- st_results |>
    pivot_longer(
        cols = -figure_names,
        names_to = "biome",
        values_to = "proportion"
    )

sup_table_2_s3 <- kapk_cdata |>
    select(label = figure_names) |>
    inner_join(kapk_cdata_agg) |>
    select(short_label, figure_names = label) |>
    inner_join(st_results_long) |>
    rename(label_orig = figure_names)





# Other data
kapk_cdata |> select(label, figure_names)
mg_mapping <- read_tsv(file = "./data/mag-distribution/reads-dmg-distribution-mags-summary/GEM-20220926__tg2g-20220926/damaged/tp-mapping-filtered.summary.tsv.gz") |>
    filter(breadth >= 0.01) |>
    inner_join(read_csv(file = "./data/mag-distribution/reads-dmg-distribution-mags-dmg/GEM-20220926__tg2g-20220926/damaged/tp-mdmg.weight-1.csv.gz") |>
        filter(label %in% kapk_cdata$label) |>
        rename(reference = tax_id) |>
        filter(significance > dmg_thresholds$signf, damage >= dmg_thresholds$damage))

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


woodcroft_mapping <- read_tsv("./data/mag-distribution/reads-dmg-distribution-mags-summary/woodcroft2018/damaged/tp-mapping-filtered.summary.tsv.gz") |>
    filter(breadth >= 0.01) |>
    inner_join(
        read_csv("./data/mag-distribution/reads-dmg-distribution-mags-dmg/woodcroft2018/damaged/tp-mdmg.weight-1.csv.gz") |>
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

# Sample lookup: 10-char sample folder → short_label + full label
sample_meta_dart <- kapk_cdata |>
    filter(figure_names %in% samples_to_keep$label) |>
    mutate(
        sample      = substr(label, 1, 10),
        short_label = paste(site, sub(".*-", "", figure_names), sep = "_")
    ) |>
    select(sample, short_label, label)

# KEGG — DART/AGP outputs (all modules)
kegg_dart <- read_tsv("./results/functional_agp/kegg_module_damage.tsv") |>
    inner_join(sample_meta_dart, by = "sample")

sup_table_5_s1 <- kegg_dart |>
    select(short_label, label_orig = label, module, module_name, module_class,
           module_category, module_subcategory, avg_coverage, n_proteins,
           n_damaged, damage_rate, mean_p_damaged, mean_d_aa)

# KEGG — DART-authenticated (mean posterior >= 0.7)
sup_table_5_s2 <- sup_table_5_s1 |>
    filter(mean_p_damaged >= 0.7)

# CAZy — DART/AGP outputs (all families)
cazy_dart <- read_tsv("./results/functional_agp/cazy_family_damage.tsv") |>
    inner_join(sample_meta_dart, by = "sample")

sup_table_5_s3 <- cazy_dart |>
    select(short_label, label_orig = label, family, cazy_class, deg_group,
           n_proteins, n_damaged, damage_rate, mean_p_damaged, mean_d_aa,
           coverage_mean)

# CAZy — DART-authenticated (mean posterior >= 0.7)
sup_table_5_s4 <- sup_table_5_s3 |>
    filter(mean_p_damaged >= 0.7)

# Viruses — DART/AGP outputs aggregated across samples
viral_files <- list.files("./results/functional_agp/viral",
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
    janitor::clean_names() |>
    mutate(n_cds = as.integer(sub("^[^;]+;([^;]+);.*", "\\1",
        gene_content_total_genes_cds_t_rna_ge_nomad_marker))) |>
    select(reference = uvig, n_cds, ecosystem = ecosystem_classification)

sup_table_6_s2 <- imgvr |>
    filter(reference %in% sup_table_6_s1$reference)

# Briggs
sup_table_7_s1 <- read_tsv("./manuscript/tables/briggs_rl-percid.tsv.gz")
sup_table_7_s2 <- read_tsv("./manuscript/tables/briggs_rl-ancient-isrecal1.tsv.gz")


# Simulations
sim_cdata <- kapk_cdata_agg |>
    select(short_label, label) |>
    inner_join(
        kapk_cdata |>
            inner_join(label_to_file) |>
            select(label, figure_names, file) |>
            select(-label, label = figure_names, Sample = file)
    ) |>
    select(Label = short_label, label_orig = label, Sample)
sup_table_8_s1 <- read_tsv("./manuscript/tables/bowtie2-benchmarks.tsv")
sup_table_8_s1 <- sim_cdata |>
    inner_join(sup_table_8_s1) |>
    select(-Sample)
sup_table_8_s2 <- read_tsv("./manuscript/tables/saturation.tsv") |>
    rename(Sample = Label) |>
    inner_join(sim_cdata)
sup_table_8_s2 <- sim_cdata |>
    inner_join(sup_table_8_s2) |>
    select(-Sample)
sup_table_8_s3 <- read_tsv("./manuscript/tables/simulations-detection-stats-composition.tsv")
sup_table_8_s4 <- read_tsv("./manuscript/tables/simulations-abundance-mdae.tsv")
sup_table_8_s5 <- read_tsv("./manuscript/tables/simulations-damage-mdae.tsv")


library(openxlsx)

# Create the xlsx file
sup_table_1 <- createWorkbook()
addWorksheet(sup_table_1, "S1 - Sample information")
addWorksheet(sup_table_1, "S2 - Control samples")
addWorksheet(sup_table_1, "S3 - Uplift model data")
# addWorksheet(sup_table_1, "S4 - Taxonomy and damage")
addWorksheet(sup_table_1, "S5 - Re-extraction metadata")
addWorksheet(sup_table_1, "S6 - Re-extraction taxonomy")
addWorksheet(sup_table_1, "S7 - Lipid biomarkers")

writeData(sup_table_1, sheet = "S1 - Sample information", x = sup_table_1_s1 |> janitor::clean_names(case = "sentence"))
writeData(sup_table_1, sheet = "S2 - Control samples", x = sup_table_1_s2 |> janitor::clean_names(case = "sentence"))
writeData(sup_table_1, sheet = "S3 - Uplift model data", x = sup_table_1_s3 |> janitor::clean_names(case = "sentence"))
# writeData(sup_table_1, sheet = "S4 - Taxonomy and damage", x = sup_table_1_s4 |> head(50000))
writeData(sup_table_1, sheet = "S5 - Re-extraction metadata", x = sup_table_1_s5 |> janitor::clean_names(case = "sentence"))
writeData(sup_table_1, sheet = "S6 - Re-extraction taxonomy", x = sup_table_1_s6 |> janitor::clean_names(case = "sentence"))
writeData(sup_table_1, sheet = "S7 - Lipid biomarkers", x = sup_table_1_s7 |> janitor::clean_names(case = "sentence"))

saveWorkbook(sup_table_1, file = "./manuscript/supplementary/sup_table_1.xlsx", overwrite = TRUE)


# write_tsv(sup_table_1_s1, "./manuscript/tables/sup_table_1_s1.tsv")
# write_tsv(sup_table_1_s2, "./manuscript/tables/sup_table_1_s2.tsv")
# write_tsv(sup_table_1_s3, "./manuscript/tables/sup_table_1_s3.tsv")
write_tsv(sup_table_1_s4 |> janitor::clean_names(case = "sentence"), "./manuscript/tables/sup_table_1_s4.tsv")
# write_tsv(sup_table_1_s5, "./manuscript/tables/sup_table_1_s5.tsv")
# write_tsv(sup_table_1_s6, "./manuscript/tables/sup_table_1_s6.tsv")




sup_table_2 <- createWorkbook()
addWorksheet(sup_table_2, "S8 - Source samples")
addWorksheet(sup_table_2, "S9 - Source contributions")

writeData(sup_table_2, sheet = "S8 - Source samples", x = sup_table_2_s1 |> janitor::clean_names(case = "sentence"))
writeData(sup_table_2, sheet = "S9 - Source contributions", x = sup_table_2_s3 |> janitor::clean_names(case = "sentence"))

saveWorkbook(sup_table_2, file = "./manuscript/supplementary/sup_table_2.xlsx", overwrite = TRUE)

sup_table_3 <- createWorkbook()
addWorksheet(sup_table_3, "S10 - GEM mapping")
addWorksheet(sup_table_3, "S11 - Woodcroft mapping")


writeData(sup_table_3, sheet = "S10 - GEM mapping", x = sup_table_3_s1 |> janitor::clean_names(case = "sentence"))
writeData(sup_table_3, sheet = "S11 - Woodcroft mapping", x = sup_table_3_s2 |> select(names(sup_table_3_s1)) |> janitor::clean_names(case = "sentence"))

saveWorkbook(sup_table_3, file = "./manuscript/supplementary/sup_table_3.xlsx", overwrite = TRUE)

sup_table_5 <- createWorkbook()
addWorksheet(sup_table_5, "S12 - KEGG all")
addWorksheet(sup_table_5, "S13 - KEGG authenticated")
addWorksheet(sup_table_5, "S14 - CAZy all")
addWorksheet(sup_table_5, "S15 - CAZy authenticated")

writeData(sup_table_5, sheet = "S12 - KEGG all",           x = sup_table_5_s1 |> janitor::clean_names(case = "sentence"))
writeData(sup_table_5, sheet = "S13 - KEGG authenticated", x = sup_table_5_s2 |> janitor::clean_names(case = "sentence"))
writeData(sup_table_5, sheet = "S14 - CAZy all",           x = sup_table_5_s3 |> janitor::clean_names(case = "sentence"))
writeData(sup_table_5, sheet = "S15 - CAZy authenticated", x = sup_table_5_s4 |> janitor::clean_names(case = "sentence"))

saveWorkbook(sup_table_5, file = "./manuscript/supplementary/sup_table_5.xlsx", overwrite = TRUE)

sup_table_6 <- createWorkbook()
addWorksheet(sup_table_6, "S16 - Viral references")
addWorksheet(sup_table_6, "S17 - IMGVR annotation")

writeData(sup_table_6, sheet = "S16 - Viral references",  x = sup_table_6_s1 |> janitor::clean_names(case = "sentence"))
writeData(sup_table_6, sheet = "S17 - IMGVR annotation",  x = sup_table_6_s2 |> janitor::clean_names(case = "sentence"))

saveWorkbook(sup_table_6, file = "./manuscript/supplementary/sup_table_6.xlsx", overwrite = TRUE)


sup_table_7 <- createWorkbook()
addWorksheet(sup_table_7, "S1 - Briggs RL Percid")
addWorksheet(sup_table_7, "S2 - Briggs RL Ancient Isrecal1")

writeData(sup_table_7, sheet = "S1 - Briggs RL Percid", x = sup_table_7_s1)
writeData(sup_table_7, sheet = "S2 - Briggs RL Ancient Isrecal1", x = sup_table_7_s2)

saveWorkbook(sup_table_7, file = "./manuscript/supplementary/sup_table_7.xlsx", overwrite = TRUE)


sup_table_8 <- createWorkbook()
addWorksheet(sup_table_8, "S18 - Bowtie2 benchmarks")
addWorksheet(sup_table_8, "S19 - Saturation")
addWorksheet(sup_table_8, "S20 - Simulations composition")
addWorksheet(sup_table_8, "S21 - Simulations abundance")
addWorksheet(sup_table_8, "S22 - Simulations damage")

writeData(sup_table_8, sheet = "S18 - Bowtie2 benchmarks", x = sup_table_8_s1 |> janitor::clean_names(case = "sentence"))
writeData(sup_table_8, sheet = "S19 - Saturation", x = sup_table_8_s2 |> janitor::clean_names(case = "sentence"))
writeData(sup_table_8, sheet = "S20 - Simulations composition", x = sup_table_8_s3 |> janitor::clean_names(case = "sentence"))
writeData(sup_table_8, sheet = "S21 - Simulations abundance", x = sup_table_8_s4 |> janitor::clean_names(case = "sentence"))
writeData(sup_table_8, sheet = "S22 - Simulations damage", x = sup_table_8_s5 |> janitor::clean_names(case = "sentence"))

saveWorkbook(sup_table_8, file = "./manuscript/supplementary/sup_table_8.xlsx", overwrite = TRUE)
