library(tidyverse)
library(janitor)
library(ggthemr)
library(ggpubr)
library(phytools)
library(lvplot)
library(showtext)

# Constants
DATA_DIR <- "./data"
RESULTS_DIR <- "./results"
MIN_SEQUENCES <- 10e6
CONTROL_REFS_PATH <- file.path(RESULTS_DIR, "controls/controls.references2remove.tsv.gz")

# Initialize text rendering
showtext_auto()

#' Load custom library functions
load_custom_lib <- function() {
    source("libs/lib.R")
}

#' Load all metadata files
load_metadata <- function() {
    list(
        cdata = read_tsv(file.path(DATA_DIR, "cdata/KapK-cdata-manuscript-20221211.tsv")),
        cdata_agg = read_tsv(file.path(DATA_DIR, "cdata/KapK-cdata-agg-manuscript-20221211.tsv"))
    )
}

#' Load and process taxonomic annotations
load_taxonomy <- function() {
    read_tsv(
        file.path(DATA_DIR, "taxonomy/hires-organelles-viruses-arctic.tax.tsv"),
        col_names = c("reference", "tax_string")
    ) |>
        separate(
            col = tax_string,
            sep = ";",
            into = c(
                "domain", "lineage", "kingdom", "phylum",
                "class", "order", "family", "genus",
                "species", "strain"
            )
        )
}

#' Load and process statistics for different stages
load_stats <- function() {
    stats_dir <- file.path(DATA_DIR, "stats")
    list(
        initial = read_tsv(file.path(stats_dir, "all.stats-initial-summary.tsv.gz")),
        extension = read_tsv(file.path(stats_dir, "all.stats-extension-summary.tsv.gz")),
        derep = read_tsv(file.path(stats_dir, "all.stats-derep-summary.tsv.gz"))
    )
}

#' Calculate average read length differences
calculate_read_differences <- function(initial_stats, extension_stats) {
    initial_stats |>
        select(initial_avg_len = avg_len, label) |>
        inner_join(extension_stats |>
            select(extension_avg_len = avg_len, label)) |>
        mutate(avg_len_diff = extension_avg_len - initial_avg_len)
}

#' Create sequence distribution plot
plot_sequence_distribution <- function(initial_stats, derep_stats, kapk_cdata) {
    initial_stats |>
        select(num_seqs, label) |>
        mutate(step = "initial") |>
        bind_rows(derep_stats |>
            select(num_seqs, label) |>
            mutate(step = "derep")) |>
        group_by(label) |>
        arrange(desc(num_seqs), .by_group = TRUE) |>
        mutate(
            ratio = 1 - (num_seqs / lag(num_seqs)),
            diff = abs(num_seqs - lag(num_seqs))
        ) |>
        ungroup() |>
        group_by(label) |>
        arrange(desc(num_seqs), .by_group = TRUE) |>
        mutate(diff1 = lag(num_seqs) - diff) |>
        ungroup() |>
        filter(!is.na(diff1)) |>
        select(label, duplicated = diff, unique = diff1) |>
        pivot_longer(
            cols = c(duplicated, unique),
            names_to = "type", values_to = "num_seqs"
        ) |>
        inner_join(kapk_cdata, by = c("label" = "label")) |>
        group_by(figure_names, type, site, member_unit, site_rnk) |>
        summarise(num_seqs = sum(num_seqs), .groups = "drop") |>
        mutate(
            figure_names = fct_reorder(figure_names, site_rnk),
            member_unit = fct_relevel(member_unit, c("B3", "B2", "B1"))
        ) |>
        ggplot(aes(x = figure_names, y = num_seqs, fill = type)) +
        geom_col(position = "stack", color = "black") +
        theme_bw() +
        labs(
            x = "Sample",
            y = "Number of sequences"
        ) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16),
            legend.position = "top",
            legend.title = element_blank(),
            text = element_text(size = 24)
        ) +
        scale_fill_manual(values = c("duplicated" = "#e6001b", "unique" = "#464144")) +
        facet_grid(~member_unit, scales = "free_x", space = "free") +
        scale_y_continuous(labels = scales::label_number(suffix = "M", scale = 1e-6))
}

#' Process and filter taxonomic data
process_taxonomic_data <- function(kapk_cdata, control_references) {
    tax_data <- read_tsv(file.path(DATA_DIR, "taxonomy/tp-mapping-filtered.summary.tsv.gz")) |>
        filter(label %in% kapk_cdata$label) |>
        filter(!reference %in% control_references) |>
        filter(breadth >= 0.01) |>
        mutate(abundance = ifelse(tax_abund_tad == 0, tax_abund_read, tax_abund_tad))

    write_tsv(tax_data, file.path(RESULTS_DIR, "taxonomy/tp-mapping-filtered.nocontam.tax.tsv.gz"))
    return(tax_data)
}

#' Filter samples based on minimum sequence threshold
filter_samples <- function(tax_data, derep_stats, kapk_cdata) {
    derep_stats_filtered <- derep_stats |>
        inner_join(kapk_cdata) |>
        mutate(
            step = "derep",
            included = ifelse(num_seqs >= MIN_SEQUENCES, "Included", "Excluded")
        )

    tax_data_filtered <- tax_data |>
        filter(label %in% (derep_stats_filtered |>
            filter(included == "Included") |>
            pull(label)))

    write_tsv(
        tax_data_filtered,
        file.path(RESULTS_DIR, "taxonomy/tp-mapping-filtered.nocontam.10M.tax.tsv.gz")
    )
    return(tax_data_filtered)
}

#' Aggregate filtered taxonomic data
aggregate_taxonomic_data <- function(tax_data_filtered, kapk_cdata) {
    tax_data_agg <- tax_data_filtered |>
        select(
            label, reference, abundance, breadth,
            read_ani_mean, coverage_mean_trunc, coverage_mean
        ) |>
        inner_join(kapk_cdata |> select(label, figure_names)) |>
        group_by(figure_names, reference) |>
        summarise(
            abundance = janitor::round_half_up(mean(abundance)),
            read_ani_mean = mean(read_ani_mean),
            coverage_mean_trunc = mean(coverage_mean_trunc),
            breadth = mean(breadth),
            coverage_mean = mean(coverage_mean),
            .groups = "drop"
        ) |>
        rename(label = figure_names)

    write_tsv(
        tax_data_agg,
        file.path(RESULTS_DIR, "taxonomy/tp-mapping-filtered.nocontam.10M.agg.tax.tsv.gz")
    )
    return(tax_data_agg)
}

#' Main execution function
main <- function() {
    # Load custom library
    load_custom_lib()

    # Load metadata
    metadata <- load_metadata()
    kapk_cdata <- metadata$cdata
    kapk_cdata_agg <- metadata$cdata_agg

    # Load taxonomy
    tax_info <- load_taxonomy()

    # Count samples
    sample_counts <- kapk_cdata |> count(label)

    # Load control references
    control_references <- read_tsv(CONTROL_REFS_PATH)

    # Load statistics
    stats <- load_stats()

    # Calculate read differences
    avg_read_diff <- calculate_read_differences(stats$initial, stats$extension)

    # Create sequence distribution plot
    seq_dist_plot <- plot_sequence_distribution(stats$initial, stats$derep, kapk_cdata)

    # Process taxonomic data
    tax_data <- process_taxonomic_data(kapk_cdata, control_references)

    # Filter samples
    tax_data_filtered <- filter_samples(tax_data, stats$derep, kapk_cdata)

    # Aggregate taxonomic data
    tax_data_agg <- aggregate_taxonomic_data(tax_data_filtered, kapk_cdata)
}

# Run the main function
main()
