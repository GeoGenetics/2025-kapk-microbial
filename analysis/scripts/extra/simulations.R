library(tidyverse)
library(purrr)
library(ggpubr)
library(rstatix)
library(unixtools)
library(lvplot)
library(janitor)
library(Metrics)
library(broom)
library(showtext)
showtext_auto()
set.tempdir("/maps/projects/fernandezguerra/scratch/tmp")

source("./scripts/lib-sims.R")
# Objectives:
# 1. Use the simulated data to evaluate the mapping strategy and the filtering strategy at species level
# 2. Use the simulated data to evaluate the damage model at species level

k_colors <- c(
    k50 = "#1f77b4",
    k100 = "#ff7f0e",
    k250 = "#2ca02c",
    k500 = "#d62728",
    k1000 = "#9467bd"
)

# read taxonomic annotations
tax_info <- read_tsv("./data/sims/cdata/hires-organelles-viruses.tax.tsv",
    col_names = c("reference", "tax_string"), show_col_types = FALSE
) |>
    separate(
        col = tax_string,
        sep = ";",
        into = c(
            "domain",
            "lineage",
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species"
        )
    )

# Read baseline data
baseline_files <- list.files(path = "./data/sims/baseline-v4", pattern = "tsv", full.names = TRUE)
baseline_data <- map_dfr_progress(baseline_files, read_tsv, show_col_types = FALSE) |>
    separate(Taxon, into = c("Taxon", "reference"), sep = "----") |>
    select(
        reference,
        label = Community,
        abundance = Perc_rel_abund
    ) |>
    mutate(abundance = abundance / 100) |>
    inner_join(tax_info)

baseline_data |>
    ggplot(aes(x = label, y = abundance)) +
    geom_boxplot() +
    scale_y_log10(labels = scales::percent) +
    coord_flip()

# Read simulated data
sim_files <- list.files(
    path = "./data/sims/taxonomic-profiling-summary",
    pattern = "tp-mapping-filtered.summary.tsv.gz", full.names = TRUE, recursive = TRUE
)
sim_data <- map_dfr_progress(sim_files, get_sim_data) |>
    inner_join(tax_info)

# Calculate the number of false positives and false negatives species
# and  the number of true positives and true negatives species
# for each simulated dataset

# False positive species that is not present in the baseline data.
# False negative as the failure to detect a species present in the baseline data
# True positive species that is present in the baseline data.


# Precision = true positives /(true positives + false positives)
# Recall = true positives /(true positives + false negatives)
# F1 = (2 * precision * recall)/(precision + recall)
# F0.5 = ((1+0.5^2) * precision * recall)/((0.5^2 * precision) + recall)

# Let's see how the different breadth thresholds affect the results
bf <- c(0, 0.01, 0.02, 0.05, 0.1)

baseline_data_species <- baseline_data |>
    select(label, species) |>
    distinct()

# True positives
true_positives <- map_dfr_progress(bf, function(X) {
    sim_data |>
        filter(breadth >= X) |>
        select(k, wf, id, label, species) |>
        distinct() |>
        inner_join(baseline_data_species, by = c("label", "species")) |>
        mutate(breadth_filter = X)
})
# False positives
false_positives <- map_dfr_progress(bf, function(X) {
    sim_data |>
        filter(breadth >= X) |>
        select(k, wf, id, label, species) |>
        distinct() |>
        anti_join(baseline_data_species, by = c("label", "species")) |>
        mutate(breadth_filter = X)
})

# Get the different combinations of k, id, wf that we used for the testing
k <- sim_data |>
    select(k) |>
    distinct() |>
    pull(k)

id <- sim_data |>
    select(id) |>
    distinct() |>
    pull(id)

wf <- sim_data |>
    select(wf) |>
    distinct() |>
    pull(wf)
# create all combinations of k, id, wf
sim_data_combinations <- expand.grid(k, id, wf, bf) |>
    as_tibble() |>
    rename(k = Var1, id = Var2, wf = Var3, breadth_filter = Var4)

# False negatives
false_negatives <- pbmcapply::pbmclapply(seq_len(nrow(sim_data_combinations)),
    function(x) {
        data <- sim_data_combinations[x, ]
        k <- data$k
        id <- data$id
        wf <- data$wf
        breadth_filter <- data$breadth_filter
        s <- sim_data |>
            filter(k == k, id == id, wf == wf, breadth >= breadth_filter) |>
            select(label, species) |>
            distinct()
        baseline_data_species |>
            anti_join(s, by = c("label", "species")) |>
            mutate(k = k, id = id, wf = wf, label = label, species = species, breadth_filter = breadth_filter)
    },
    mc.cores = 16
)

false_negatives <- bind_rows(false_negatives) |>
    group_by(k, wf, id, breadth_filter, label) |>
    count(name = "false_negatives") |>
    ungroup() |>
    complete(
        label = true_positives$label,
        k = true_positives$k,
        wf = true_positives$wf,
        id = true_positives$id,
        breadth_filter = true_positives$breadth_filter,
        fill = list(false_negatives = 0)
    )



detection_stats <- true_positives |>
    group_by(k, wf, id, breadth_filter, label) |>
    count(name = "true_positives") |>
    inner_join(false_positives |>
        group_by(k, wf, id, breadth_filter, label) |>
        count(name = "false_positives")) |>
    inner_join(false_negatives) |>
    mutate(
        precision = true_positives / (true_positives + false_positives),
        recall = true_positives / (true_positives + false_negatives),
        F1 = (2 * precision * recall) / (precision + recall),
        F05 = ((1 + 0.5^2) * precision * recall) / ((0.5^2 * precision) + recall)
    )

plot_detection_stats <- function(wf_mode = "raw") {
    detection_stats |>
        select(k, wf, id, breadth_filter, label, precision, recall, F1, F05) |>
        group_by(k, wf, id, breadth_filter) |>
        ungroup() |>
        pivot_longer(cols = c(precision, recall, F1, F05), names_to = "metric", values_to = "value") |>
        mutate(
            k = fct_relevel(k, c("k50", "k100", "k250", "k500", "k1000")),
            metric = fct_relevel(metric, "precision", "recall", "F1", "F05")
        ) |>
        filter(wf == wf_mode) |>
        ggplot(aes(x = id, y = value, fill = k)) +
        geom_hline(yintercept = 0.8, linetype = "dashed") +
        geom_violin(position = position_dodge(width = 0.8), alpha = 0.5, linewidth = 0.1) +
        geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.5, width = .2, linewidth = 0.1) +
        theme_bw() +
        theme(
            text = element_text(size = 14),
            legend.position = "top",
            legend.title = element_blank(),
        ) +
        xlab("%ANI filtering") +
        ylab("Value") +
        facet_grid(breadth_filter ~ metric) +
        scale_fill_manual(values = k_colors)
    fname <- paste0("./manuscript/figures/simulations/detection-stats-", wf_mode, ".pdf")
    ggsave(fname, width = 12, height = 6)
}



detection_stats |>
    select(k, wf, id, breadth_filter, label, precision, recall, F1, F05) |>
    group_by(k, wf, id, breadth_filter) |>
    ungroup() |>
    pivot_longer(cols = c(precision, recall, F1, F05), names_to = "metric", values_to = "value") |>
    mutate(
        k = fct_relevel(k, c("k50", "k100", "k250", "k500", "k1000")),
        metric = fct_relevel(metric, "precision", "recall", "F1", "F05")
    ) |>
    filter(wf == "standard") |>
    ggplot(aes(x = id, y = value, fill = k)) +
    geom_hline(yintercept = 0.8, linetype = "dashed") +
    geom_violin(position = position_dodge(width = 0.8), alpha = 0.5, linewidth = 0.1) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.5, width = .2, linewidth = 0.1) +
    theme_bw() +
    theme(
        text = element_text(size = 12),
        legend.position = "top",
        legend.title = element_blank(),
    ) +
    xlab("%ANI filtering") +
    ylab("Value") +
    facet_grid(breadth_filter ~ metric) +
    scale_fill_manual(values = k_colors)


plot_detection_stats("raw")
plot_detection_stats("standard")

# We will keep the results for k1000 and breadth_filter = 0.01

# Which is the best id for each metric?
# Do some stats to compare the differences using the wilcoxon rank sum test
my_comparisons <- list(c(1, 6))
detection_stats |>
    select(k, wf, id, breadth_filter, label, precision, recall, F1, F05) |>
    group_by(k, wf, id, breadth_filter) |>
    ungroup() |>
    pivot_longer(cols = c(precision, recall, F1, F05), names_to = "metric", values_to = "value") |>
    filter(k == "k1000", breadth_filter == 0.01) |>
    mutate(
        id = as.factor(id),
        metric = fct_relevel(metric, "precision", "recall", "F1", "F05")
    ) |>
    ggplot(aes(x = id, y = value)) +
    geom_violin() +
    geom_boxplot(width = 0.1) +
    geom_pwc(
        aes(group = id),
        tip.length = 0.01,
        method = "wilcox.test",
        label = "{p.adj.format} {p.adj.signif}",
        method.args = list(
            alternative = "two.sided",
            paired = FALSE
        ),
        p.adjust.method = "BH",
        p.adjust.by = "group",
        hide.ns = TRUE
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 14),
    ) +
    xlab("%ANI filtering") +
    ylab("Value") +
    facet_grid(wf ~ metric) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, NA))


stats_composition <- detection_stats |>
    select(k, wf, id, breadth_filter, label, precision, recall, F1, F05) |>
    pivot_longer(cols = c(precision, recall, F1, F05), names_to = "metric", values_to = "value") |>
    # filter(k == "k1000", breadth_filter == 0.01) |>
    group_by(k, wf, id, breadth_filter, metric) |>
    summarise(
        mean = mean(value),
        sd = sd(value),
        n = n(),
        median = median(value),
    ) |>
    ungroup() |>
    filter(wf == "standard") |>
    select(-wf)

write_tsv(stats_composition, "./manuscript/tables/simulations-detection-stats-composition.tsv")


detection_stats |>
    select(k, wf, id, breadth_filter, label, precision, recall, F1, F05) |>
    group_by(k, wf, id, breadth_filter) |>
    ungroup() |>
    pivot_longer(cols = c(precision, recall, F1, F05), names_to = "metric", values_to = "value") |>
    filter(k == "k1000", wf == "standard", breadth_filter == 0.01) |>
    mutate(
        id = as.factor(id),
    ) |>
    group_by(metric) |>
    pairwise_wilcox_test(value ~ id, alternative = "two.sided", p.adjust.method = "BH", paired = FALSE) |>
    knitr::kable()

# It seems that the best id for the profiling is 95%, although the difference is not significant with 94%


# What about the abundances?
# We will do it easy and calculate the MAE and correlations between the baseline and the simulated data

# First, we need to get the abundance of the baseline
# We will aggregate at the species level
baseline_data_abundance <- baseline_data |>
    select(species, label, abundance) |>
    group_by(species, label) |>
    summarise(prop_baseline = sum(abundance)) |>
    ungroup() |>
    select(species, label, prop_baseline)

# Now, we need to get the abundance of the simulated data
# We will aggregate at the species level
sim_data_abundance <- pbmcapply::pbmclapply(bf, function(X) {
    sim_data |>
        filter(breadth >= X) |>
        mutate(abundance = ifelse(tax_abund_tad == 0, tax_abund_read, tax_abund_tad)) |>
        select(k, wf, id, label, species, abundance) |>
        group_by(k, wf, id, label, species) |>
        summarise(abundance = round_half_up(mean(abundance))) |>
        ungroup() |>
        group_by(k, wf, id, label) |>
        mutate(total_abundance = sum(abundance)) |>
        ungroup() |>
        mutate(prop_simdata = abundance / total_abundance) |>
        select(k, wf, id, label, species, prop_simdata) |>
        mutate(breadth_filter = X)
}, mc.cores = 16) |>
    bind_rows()



# Let's checl if there are differences between the raw and the standard workflow
sim_data_abundance_raw <- sim_data_abundance |>
    filter(wf == "raw") |>
    select(k, id, label, species, breadth_filter, prop_simdata_raw = prop_simdata, -wf)

sim_data_abundance_std <- sim_data_abundance |>
    filter(wf == "standard") |>
    select(k, id, label, species, breadth_filter, prop_simdata_std = prop_simdata, -wf)

# Now, we can join the data and calculate the MDAE
sim_data_abundance_raw_vs_std_mdae <- sim_data_abundance_raw |>
    inner_join(sim_data_abundance_std) |>
    group_by(k, id, label, breadth_filter) |>
    summarize(
        mdae = Metrics::mdae(prop_simdata_raw * 100, prop_simdata_std * 100)
    ) |>
    ungroup() |>
    mutate(k = fct_relevel(k, c("k50", "k100", "k250", "k500", "k1000")))

ggplot(sim_data_abundance_raw_vs_std_mdae, aes(x = id, y = mdae, fill = k)) +
    geom_violin(position = position_dodge(width = 0.8), alpha = 0.5) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.5, width = .2) +
    theme_bw() +
    theme(
        text = element_text(size = 14),
        legend.position = "top",
        legend.title = element_blank(),
    ) +
    xlab("%ANI filtering") +
    ylab("Median absolute error") +
    facet_grid(breadth_filter ~ .) +
    scale_fill_manual(values = k_colors)
ggsave("./manuscript/figures/simulations-abundance-mdae.pdf", width = 12, height = 6)


sim_data_abundance_raw_vs_std_mdae |>
    group_by(k, id, breadth_filter) |>
    summarize(
        mean = mean(mdae),
        sd = sd(mdae),
        n = n(),
        median = median(mdae),
    ) |>
    write_tsv("./manuscript/tables/simulations-abundance-mdae.tsv")



# There are almost no differences between the raw and the standard workflow

# Let's check how correlated are the results of the raw and the standard workflow
sim_data_abundance_raw_vs_std_cor <- sim_data_abundance_raw |>
    inner_join(sim_data_abundance_std) |>
    group_by(k, id, label, breadth_filter) |>
    do(tidy(Hmisc::rcorr(.$prop_simdata_raw, .$prop_simdata_std, type = "spearman"))) |>
    ungroup() |>
    mutate(
        p.value.adj.bh = p.adjust(p.value, method = "BH"),
        p.value.adj.holm = p.adjust(p.value, method = "holm"),
        p.value.adj.fdr = p.adjust(p.value, method = "fdr")
    )

ggplot(sim_data_abundance_raw_vs_std_cor, aes(x = id, y = estimate, fill = k)) +
    geom_violin(position = position_dodge(width = 0.6), alpha = 0.5) +
    geom_boxplot(position = position_dodge(width = 0.6), outlier.size = 0.5, width = .1) +
    theme_bw() +
    theme(
        text = element_text(size = 14),
        legend.position = "top",
        legend.title = element_blank(),
    ) +
    xlab("%ANI filtering") +
    ylab("Spearman correlation (ρ)") +
    facet_grid(breadth_filter ~ .)

# Results are almost identical between the raw and the standard workflow
# Let's visuallize the correlation between the raw and the standard workflow at 0.01 breadth filter
sim_data_abundance_raw |>
    inner_join(sim_data_abundance_std) |>
    filter(breadth_filter == 0.01) |>
    mutate(k = fct_relevel(k, c("k50", "k100", "k250", "k500", "k1000"))) |>
    ggplot(aes(x = prop_simdata_raw, y = prop_simdata_std)) +
    geom_hdr() +
    # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
    scale_fill_viridis_c(direction = -1, option = "inferno") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey", alpha = 0.5) +
    facet_grid(id ~ k) +
    scale_x_log10(expand = c(0, 0)) +
    scale_y_log10(expand = c(0, 0)) +
    theme_bw() +
    theme(
        legend.position = "top",
        text = element_text(size = 14),
        legend.key.width = unit(1, "cm"),
    ) +
    xlab("Abundance raw") +
    ylab("Abundance standard")


# We will use the standard workflow for the rest of the analysis
# Let's check the correlation between the baseline and the simulated data

# We will start by looking at the median absolute error between the baseline and the simulated data
sim_data_abundance_std |>
    inner_join(baseline_data_abundance) |>
    group_by(k, id, label, breadth_filter) |>
    summarize(
        mdae = Metrics::mdae(prop_baseline * 100, prop_simdata_std * 100)
    ) |>
    ungroup() |>
    mutate(k = fct_relevel(k, c("k50", "k100", "k250", "k500", "k1000"))) |>
    ggplot(aes(x = id, y = mdae, fill = k)) +
    geom_violin(position = position_dodge(width = 0.8), alpha = 0.5) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.5, width = .2) +
    theme_bw() +
    theme(
        text = element_text(size = 14),
        legend.position = "top",
        legend.title = element_blank(),
    ) +
    xlab("%ANI filtering") +
    ylab("Median absolute error") +
    facet_grid(breadth_filter ~ .)

# The median absolute error is very low, so the difference between the baseline and the simulated data is very small

# Let's check the correlation between the baseline and the simulated data
sim_data_abundance_cor <- sim_data_abundance_std |>
    inner_join(baseline_data_abundance) |>
    group_by(k, id, label, breadth_filter) |>
    do(tidy(Hmisc::rcorr(.$prop_baseline, .$prop_simdata_std, type = "spearman"))) |>
    ungroup() |>
    mutate(
        p.value.adj.bh = p.adjust(p.value, method = "BH"),
        p.value.adj.holm = p.adjust(p.value, method = "holm"),
        p.value.adj.fdr = p.adjust(p.value, method = "fdr")
    )

sim_data_abundance_cor |>
    mutate(k = fct_relevel(k, c("k50", "k100", "k250", "k500", "k1000"))) |>
    ggplot(aes(x = id, y = estimate, fill = k)) +
    geom_violin(position = position_dodge(width = 0.8), alpha = 0.5) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.5, width = .2) +
    # geom_jitter(position = position_dodge(), shape = 21, color = "black") +
    theme_bw() +
    theme(
        text = element_text(size = 14),
        legend.position = "top",
        legend.title = element_blank(),
    ) +
    xlab("%ANI filtering") +
    ylab("Spearman Correlation (ρ)") +
    facet_grid(breadth_filter ~ .)

# At 0.01 breadth filter, the correlation is very high

# Let's explore the damage estimates
# Prepare the baseline data
kapk_cdata <- read_tsv("./data/cdata/KapK-cdata-manuscript-20221211.tsv") |>
    mutate(label = str_sub(label, 1, 10)) |>
    select(sample = label, label = figure_names)

baseline_dmg_files <- list.files(path = "./data/sims/baseline", pattern = ".tp-mdmg.weight-1.csv.gz", full.names = TRUE)

baseline_dmg_data <- map_dfr_progress(baseline_dmg_files, function(X) {
    read_csv(X, show_col_types = FALSE) |>
        select(sample, reference = tax_id, damage, N_reads, significance)
}) |>
    inner_join(kapk_cdata) |>
    select(-sample) |>
    inner_join(baseline_data)

# Read the simulated data
sim_dmg_files <- list.files(path = "./data/sims/taxonomic-profiling-dmg", pattern = ".tp-mdmg.weight-1.csv.gz", full.names = TRUE, recursive = TRUE)
sim_dmg_files <- sim_dmg_files[grepl("local", sim_dmg_files)]

sim_data_dmg <- pbmcapply::pbmclapply(sim_dmg_files, get_dmg_sim_data, mc.cores = 16) |>
    bind_rows()

sim_data_dmg_bf <- map_dfr_progress(bf, function(X) {
    sim_data_dmg |>
        filter(significance >= 2) |>
        inner_join(sim_data, by = c("label", "reference", "k", "id", "wf")) |>
        filter(breadth >= X) |>
        rename(damage_sim = damage) |>
        ungroup() |>
        mutate(breadth_filter = X)
})

# We will do two different evaluations, at the reference level and at the species level

# Let's start comparing the damage estimates at the reference level,
# Let's check the differences between both workflows
sim_data_dmg_bf_raw <- sim_data_dmg_bf |>
    filter(wf == "raw") |>
    select(reference, k, id, label, breadth_filter, damage_sim_raw = damage_sim, n_alns, n_reads)

sim_data_dmg_bf_std <- sim_data_dmg_bf |>
    filter(wf == "standard") |>
    select(reference, k, id, label, breadth_filter, damage_sim_std = damage_sim, n_alns, n_reads)

sim_data_dmg_bf_raw |>
    inner_join(sim_data_dmg_bf_std) |>
    mutate(
        dmg_diff = ((damage_sim_std) - (damage_sim_raw)),
        k = fct_relevel(k, c("k50", "k100", "k250", "k500", "k1000"))
    ) |>
    ggplot(aes(x = id, y = dmg_diff, fill = k)) +
    geom_violin(position = position_dodge(width = 0.6), alpha = 0.5, linewidth = 0.5) +
    geom_boxplot(position = position_dodge(width = 0.6), outlier.shape = NA, width = .2, linewidth = 0.1) +
    theme_bw() +
    theme(
        text = element_text(size = 14),
        legend.position = "top",
        legend.title = element_blank(),
    ) +
    xlab("%ANI filtering") +
    ylab("Standard-Raw damage difference") +
    facet_grid(breadth_filter ~ .)

# The standard workflow produces slightly higher damage estimates than the raw workflow
sim_data_dmg_bf_raw |>
    inner_join(sim_data_dmg_bf_std) |>
    filter(k == "k1000") |>
    ggplot(aes(x = damage_sim_raw, y = damage_sim_std)) +
    geom_hdr() +
    # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
    scale_fill_viridis_c(direction = -1, option = "inferno") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey", alpha = 0.5) +
    theme_bw() +
    theme(
        legend.position = "top",
        text = element_text(size = 14),
        legend.key.width = unit(1, "cm"),
    ) +
    xlab("Damage raw") +
    ylab("Damage standard") +
    facet_grid(breadth_filter ~ id)

# We will go with the standard workflow and k1000
# Let's check the differences between the baseline and the simulated data
sim_data_dmg_bf_std |>
    inner_join(baseline_dmg_data |> select(label, reference, damage, N_reads) |> rename(damage_base = damage)) |>
    mutate(
        dmg_diff = ((damage_sim_std) - (damage_base)),
        k = fct_relevel(k, c("k50", "k100", "k250", "k500", "k1000"))
    ) |>
    ggplot(aes(x = id, y = dmg_diff, fill = k)) +
    geom_violin(position = position_dodge(width = 0.6), alpha = 0.5, linewidth = 0.5) +
    geom_boxplot(position = position_dodge(width = 0.6), outlier.shape = NA, width = .2, linewidth = 0.1) +
    theme_bw() +
    theme(
        text = element_text(size = 14),
        legend.position = "top",
        legend.title = element_blank(),
    ) +
    xlab("%ANI filtering") +
    ylab("Standard-Baseline damage difference") +
    facet_grid(breadth_filter ~ .)

# The simulated damaged estimates are slightly lower than the baseline damage estimates

sim_data_dmg_bf_std |>
    inner_join(baseline_dmg_data |> select(label, reference, damage, N_reads) |> rename(damage_base = damage)) |>
    filter(k == "k1000") |>
    ggplot(aes(x = damage_base, y = damage_sim_std)) +
    # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
    geom_hdr() +
    scale_fill_viridis_c(direction = -1, option = "inferno") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey", alpha = 0.5) +
    theme_bw() +
    theme(
        legend.position = "top",
        text = element_text(size = 14),
        legend.key.width = unit(1, "cm"),
    ) +
    xlab("Damage baseline") +
    ylab("Damage standard") +
    facet_grid(breadth_filter ~ id)


# Let's calculate the MDAE
sim_data_dmg_bf_std |>
    inner_join(baseline_dmg_data |> select(label, reference, damage, N_reads) |> rename(damage_base = damage)) |>
    group_by(k, id, label, breadth_filter) |>
    summarise(
        mdae = Metrics::mae(damage_base, damage_sim_std)
    ) |>
    ungroup() |>
    mutate(k = fct_relevel(k, c("k50", "k100", "k250", "k500", "k1000"))) |>
    ggplot(aes(x = id, y = mdae, fill = k)) +
    geom_violin(position = position_dodge(width = 0.8), alpha = 0.5) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.5, width = .2) +
    # geom_jitter(position = position_dodge(), shape = 21, color = "black") +
    theme_bw() +
    theme(
        text = element_text(size = 14),
        legend.position = "top",
        legend.title = element_blank(),
    ) +
    xlab("%ANI filtering") +
    ylab("Median Absolute Error") +
    facet_grid(breadth_filter ~ .) +
    scale_fill_manual(values = k_colors)
ggsave("./manuscript/figures/simulations-damage-mdae.pdf", width = 12, height = 6)

sim_data_dmg_bf_std |>
    inner_join(baseline_dmg_data |> select(label, reference, damage, N_reads) |> rename(damage_base = damage)) |>
    group_by(k, id, label, breadth_filter) |>
    summarise(
        mdae = Metrics::mae(damage_base, damage_sim_std)
    ) |>
    ungroup() |>
    group_by(k, id, breadth_filter) |>
    summarize(
        mean = mean(mdae),
        sd = sd(mdae),
        n = n(),
        median = median(mdae),
    ) |>
    write_tsv("./manuscript/tables/simulations-damage-mdae.tsv")

# As before the estimates are slightly different, but not too much



# Let's check the correlation between the baseline and the simulated damage estimates

sim_data_dmg_bf_std_cor <- sim_data_dmg_bf_std |>
    inner_join(baseline_dmg_data |> select(label, reference, damage, N_reads) |> rename(damage_base = damage)) |>
    group_by(k, id, label, breadth_filter) |>
    do(tidy(Hmisc::rcorr(.$damage_base, .$damage_sim_std, type = "pearson"))) |>
    ungroup() |>
    mutate(
        p.value.adj.bh = p.adjust(p.value, method = "BH"),
        p.value.adj.holm = p.adjust(p.value, method = "holm"),
        p.value.adj.fdr = p.adjust(p.value, method = "fdr")
    )

sim_data_dmg_bf_std_cor |>
    mutate(k = fct_relevel(k, c("k50", "k100", "k250", "k500", "k1000"))) |>
    ggplot(aes(x = id, y = estimate, fill = k)) +
    geom_violin(position = position_dodge(width = 0.8), alpha = 0.5) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.5, width = .2) +
    # geom_jitter(position = position_dodge(), shape = 21, color = "black") +
    theme_bw() +
    theme(
        text = element_text(size = 14),
        legend.position = "top",
        legend.title = element_blank(),
    ) +
    xlab("%ANI filtering") +
    ylab("Spearman Correlation (ρ)") +
    facet_grid(breadth_filter ~ .) +
    scale_y_continuous(limits = c(0, 1))

# The correlations are pretty high

# Let's repeat the same analysis but at the species level
baseline_dmg_data_sp <- baseline_dmg_data |>
    group_by(label, species) |>
    summarise(damage_base = mean(damage)) |>
    ungroup()

sim_data_dmg_bf_raw_sp <- sim_data_dmg_bf |>
    filter(wf == "raw") |>
    select(reference, k, id, label, breadth_filter, damage_sim_raw = damage_sim, species, n_alns, n_reads) |>
    group_by_all() |>
    summarise(
        damage_sim_raw = mean(damage_sim_raw),
        n_alns = janitor::round_half_up(mean(n_alns)),
        n_reads = janitor::round_half_up(mean(n_reads))
    ) |>
    ungroup()

sim_data_dmg_bf_std_sp <- sim_data_dmg_bf |>
    filter(wf == "standard") |>
    select(reference, k, id, label, breadth_filter, damage_sim_std = damage_sim, species, n_alns, n_reads) |>
    group_by_all() |>
    summarise(
        damage_sim_std = mean(damage_sim_std),
        n_alns = janitor::round_half_up(mean(n_alns)),
        n_reads = janitor::round_half_up(mean(n_reads))
    ) |>
    ungroup()

sim_data_dmg_bf_raw_sp |>
    inner_join(sim_data_dmg_bf_std_sp) |>
    mutate(
        dmg_diff = ((damage_sim_std) - (damage_sim_raw)),
        k = fct_relevel(k, c("k50", "k100", "k250", "k500", "k1000"))
    ) |>
    ggplot(aes(x = id, y = dmg_diff, fill = k)) +
    geom_violin(position = position_dodge(width = 0.6), alpha = 0.5, linewidth = 0.5) +
    geom_boxplot(position = position_dodge(width = 0.6), outlier.shape = NA, width = .2, linewidth = 0.1) +
    theme_bw() +
    theme(
        text = element_text(size = 14),
        legend.position = "top",
        legend.title = element_blank(),
    ) +
    xlab("%ANI filtering") +
    ylab("Standard-Raw damage difference") +
    facet_grid(breadth_filter ~ .)

# The standard workflow produces slightly higher damage estimates than the raw workflow

sim_data_dmg_bf_raw_sp |>
    inner_join(sim_data_dmg_bf_std_sp) |>
    filter(k == "k1000") |>
    ggplot(aes(x = damage_sim_raw, y = damage_sim_std)) +
    geom_hdr() +
    scale_fill_viridis_c(direction = -1, option = "inferno") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey", alpha = 0.5) +
    theme_bw() +
    theme(
        legend.position = "top",
        text = element_text(size = 14),
        legend.key.width = unit(1, "cm"),
    ) +
    xlab("Damage raw") +
    ylab("Damage standard") +
    facet_grid(breadth_filter ~ id)

# We will go with the standard workflow and k1000
# Let's check the differences between the baseline and the simulated data
sim_data_dmg_bf_std_sp |>
    inner_join(baseline_dmg_data_sp |> select(label, species, damage_base)) |>
    mutate(dmg_diff = ((damage_sim_std) - (damage_base))) |>
    mutate(
        k = fct_relevel(k, c("k50", "k100", "k250", "k500", "k1000"))
    ) |>
    ggplot(aes(x = id, y = dmg_diff, fill = k)) +
    geom_violin(position = position_dodge(width = 0.9), alpha = 0.5, linewidth = 0.5) +
    geom_boxplot(position = position_dodge(width = 0.9), outlier.shape = NA, width = .2, linewidth = 0.1) +
    theme_bw() +
    theme(
        text = element_text(size = 14),
        legend.position = "top",
        legend.title = element_blank(),
    ) +
    xlab("%ANI filtering") +
    ylab("Standard-Baseline damage difference") +
    facet_grid(breadth_filter ~ .)

sim_data_dmg_bf_std_sp |>
    inner_join(baseline_dmg_data |> select(label, species, damage_base = damage)) |>
    filter(k == "k1000") |>
    mutate(dmg_diff = ((damage_sim_std) - (damage_base))) |>
    ggplot(aes(y = dmg_diff, x = damage_base)) +
    geom_hdr() +
    scale_fill_viridis_c(direction = -1, option = "inferno") +
    # geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey", alpha = 0.5) +
    # geom_hline(yintercept = 0.15, linetype = "dashed", linewidth = 0.5, color = "grey", alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey", alpha = 0.5) +
    theme_bw() +
    theme(
        legend.position = "top",
        text = element_text(size = 14),
        legend.key.width = unit(1, "cm"),
    ) +
    xlab("Baseline damage") +
    ylab("Standard-baseline damage difference") +
    facet_grid(breadth_filter ~ id)



# The simulated damaged estimates are slightly lower than the baseline damage estimates

# Let's calculate the MDAE
sim_data_dmg_bf_std_sp |>
    inner_join(baseline_dmg_data_sp |> select(label, species, damage_base)) |>
    group_by(k, id, label, breadth_filter) |>
    summarise(
        mdae = Metrics::mae(damage_base, damage_sim_std)
    ) |>
    ungroup() |>
    mutate(k = fct_relevel(k, c("k50", "k100", "k250", "k500", "k1000"))) |>
    ggplot(aes(x = id, y = mdae, fill = k)) +
    geom_violin(position = position_dodge(width = 0.9), alpha = 0.5) +
    geom_boxplot(position = position_dodge(width = 0.9), outlier.size = 0.5, width = .2) +
    # geom_jitter(position = position_dodge(), shape = 21, color = "black") +
    theme_bw() +
    theme(
        text = element_text(size = 14),
        legend.position = "top",
        legend.title = element_blank(),
    ) +
    xlab("%ANI filtering") +
    ylab("Median Absolute Error") +
    facet_grid(breadth_filter ~ .)

# As before the estimates are slightly different, but not too much

# Let's check the correlation between the baseline and the simulated damage estimates

sim_data_dmg_bf_std_cor_sp <- sim_data_dmg_bf_std_sp |>
    inner_join(baseline_dmg_data_sp |> select(label, species, damage_base)) |>
    group_by(k, id, label, breadth_filter) |>
    do(tidy(Hmisc::rcorr(.$damage_base, .$damage_sim_std, type = "pearson"))) |>
    ungroup() |>
    mutate(
        p.value.adj.bh = p.adjust(p.value, method = "BH"),
        p.value.adj.holm = p.adjust(p.value, method = "holm"),
        p.value.adj.fdr = p.adjust(p.value, method = "fdr")
    )

sim_data_dmg_bf_std_cor_sp |>
    mutate(k = fct_relevel(k, c("k50", "k100", "k250", "k500", "k1000"))) |>
    ggplot(aes(x = id, y = estimate, fill = k)) +
    geom_violin(position = position_dodge(width = 0.8), alpha = 0.5) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.5, width = .2) +
    # geom_jitter(position = position_dodge(), shape = 21, color = "black") +
    theme_bw() +
    theme(
        text = element_text(size = 14),
        legend.position = "top",
        legend.title = element_blank(),
    ) +
    xlab("%ANI filtering") +
    ylab("Spearman Correlation (ρ)") +
    facet_grid(breadth_filter ~ .) +
    scale_y_continuous(limits = c(0, 1))

sim_data_dmg_bf_std_sp |>
    inner_join(baseline_dmg_data_sp |> select(label, species, damage_base)) |>
    filter(k == "k1000") |>
    ggplot(aes(x = damage_base, y = damage_sim_std)) +
    geom_hdr() +
    scale_fill_viridis_c(direction = -1, option = "inferno") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey", alpha = 0.5) +
    geom_hline(yintercept = 0.15, linetype = "dashed", linewidth = 0.5, color = "grey", alpha = 0.5) +
    geom_vline(xintercept = 0.15, linetype = "dashed", linewidth = 0.5, color = "grey", alpha = 0.5) +
    theme_bw() +
    theme(
        legend.position = "top",
        text = element_text(size = 14),
        legend.key.width = unit(1, "cm"),
    ) +
    xlab("Damage standard") +
    ylab("Damage baseline") +
    facet_grid(breadth_filter ~ id)



# The correlation is slightly lower than the baseline

# Let's check the difference between the baseline and the simulated damage estimates and the number of alignments
sim_data_dmg_bf_std |>
    inner_join(baseline_dmg_data |> select(label, reference, damage_base = damage)) |>
    filter(k == "k1000") |>
    mutate(dmg_diff = ((damage_sim_std) - (damage_base))) |>
    ggplot(aes(y = n_alns, x = dmg_diff)) +
    geom_hdr() +
    scale_fill_viridis_c(direction = -1, option = "inferno") +
    # geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey", alpha = 0.5) +
    # geom_hline(yintercept = 0.15, linetype = "dashed", linewidth = 0.5, color = "grey", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey", alpha = 0.5) +
    theme_bw() +
    theme(
        legend.position = "top",
        text = element_text(size = 14),
        legend.key.width = unit(1, "cm"),
    ) +
    ylab("Number of alignments") +
    xlab("Standard-baseline damage difference") +
    facet_grid(breadth_filter ~ id) +
    scale_y_log10()

sim_data_dmg_bf_std_sp |>
    inner_join(baseline_dmg_data |> select(label, species, damage_base = damage)) |>
    filter(k == "k1000") |>
    mutate(dmg_diff = ((damage_sim_std) - (damage_base))) |>
    ggplot(aes(y = n_alns, x = dmg_diff)) +
    stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
    scale_fill_viridis_c(direction = -1, option = "inferno") +
    # geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey", alpha = 0.5) +
    # geom_hline(yintercept = 0.15, linetype = "dashed", linewidth = 0.5, color = "grey", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey", alpha = 0.5) +
    theme_bw() +
    theme(
        legend.position = "top",
        text = element_text(size = 14),
        legend.key.width = unit(1, "cm"),
    ) +
    ylab("Number of alignments") +
    xlab("Standard-baseline damage difference") +
    facet_grid(breadth_filter ~ id) +
    scale_y_log10()
