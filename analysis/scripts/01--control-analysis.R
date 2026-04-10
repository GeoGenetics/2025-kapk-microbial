library(tidyverse)

DATA_DIR <- "./data"
RESULTS_DIR <- "./results/controls"
TAXONOMY_DIR <- file.path(DATA_DIR, "taxonomy")
STATS_DIR <- file.path(DATA_DIR, "stats")

#' Load and process taxonomic annotations
load_taxonomy <- function() {
    read_tsv(file.path(TAXONOMY_DIR, "hires-organelles-viruses-arctic.tax.tsv"),
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

#' Load control samples data
load_control_samples <- function() {
    read_tsv(file.path(DATA_DIR, "cdata/KapK_samples-controls-20221211.tsv"))
}

#' Load statistics for a given stage
load_stage_stats <- function(stage, control_labels) {
    read_tsv(file.path(STATS_DIR, sprintf("all.stats-%s-summary.tsv.gz", stage))) |>
        filter(label %in% control_labels)
}

#' Calculate and display sequence statistics
calculate_sequence_stats <- function(initial_stats, derep_stats, control_samples) {
    initial_stats |>
        rename(num_seqs_initial = num_seqs, sum_len_initial = sum_len) |>
        inner_join(derep_stats |>
            select(label, num_seqs_derep = num_seqs, sum_len_derep = sum_len)) |>
        inner_join(control_samples) |>
        mutate(redundancy = 100 * (1 - (num_seqs_derep / num_seqs_initial))) |>
        select(
            label_orig, num_seqs_initial, num_seqs_derep,
            redundancy, sum_len_initial, sum_len_derep
        ) |>
        knitr::kable()
}

#' Create sequence distribution plot
plot_sequence_distribution <- function(initial_stats, derep_stats, control_samples) {
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
        inner_join(control_samples) |>
        select(label_orig, duplicated = diff, unique = diff1) |>
        pivot_longer(
            cols = c(duplicated, unique),
            names_to = "type", values_to = "num_seqs"
        ) |>
        ungroup() |>
        ggplot(aes(x = label_orig, y = num_seqs, fill = type)) +
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
        scale_y_continuous(labels = scales::label_number(suffix = "M", scale = 1e-6))
}

#' Process taxonomic data
process_taxonomic_data <- function(control_samples, tax_info) {
    read_tsv(file.path(TAXONOMY_DIR, "tp-mapping-filtered.summary.tsv.gz")) |>
        filter(breadth >= 0.01) |>
        inner_join(control_samples |> select(label, label_orig)) |>
        mutate(abundance = ifelse(tax_abund_tad == 0, tax_abund_read, tax_abund_tad)) |>
        inner_join(tax_info)
}

#' Process and save control taxa data
process_control_taxa <- function(tax_data) {
    control_taxa_nreads <- tax_data |>
        select(domain:species, reference, n_reads) |>
        group_by(across(domain:reference)) |>
        summarise(n_reads = sum(n_reads), n = n(), .groups = "drop") |>
        arrange(desc(n_reads))

    # Save results
    control_taxa_nreads |>
        write_tsv(file.path(RESULTS_DIR, "controls.taxa_nreads.tsv.gz"))

    control_taxa_nreads |>
        pull(reference) |>
        unique() |>
        enframe(name = "reference") |>
        write_tsv(file.path(RESULTS_DIR, "controls.references2remove.tsv.gz"))

    control_taxa_nreads |>
        top_n(10, n_reads) |>
        knitr::kable()

    return(control_taxa_nreads)
}

#' Create damage plot
plot_damage <- function(control_samples, tax_data) {
    read_csv(file.path(TAXONOMY_DIR, "tp-mdmg.weight-1.csv.gz")) |>
        inner_join(control_samples |> select(label, label_orig)) |>
        rename(reference = tax_id) |>
        select(label, label_orig, reference, damage, significance) |>
        inner_join(control_samples) |>
        inner_join(tax_data) |>
        ggplot(aes(x = significance, y = damage, fill = label_orig, size = n_reads)) +
        geom_point(shape = 21) +
        theme_bw() +
        theme(text = element_text(size = 14)) +
        labs(
            x = "Significance",
            y = "Damage"
        ) +
        scale_size_continuous(
            range = c(1, 5),
            trans = "sqrt",
            name = "# reads",
            label = scales::comma
        ) +
        guides(fill = guide_legend(
            title = "Control",
            override.aes = list(size = 5)
        ))
}

# Main execution
main <- function() {
    # Load initial data
    control_samples <- load_control_samples()
    tax_info <- load_taxonomy()

    # Load statistics
    initial_stats <- load_stage_stats("initial", control_samples$label)
    extension_stats <- load_stage_stats("extension", control_samples$label)
    derep_stats <- load_stage_stats("derep", control_samples$label)

    # Calculate and display sequence statistics
    calculate_sequence_stats(initial_stats, derep_stats, control_samples)

    # Create sequence distribution plot
    plot_sequence_distribution(initial_stats, derep_stats, control_samples)

    # Process taxonomic data
    tax_data <- process_taxonomic_data(control_samples, tax_info)

    # Process and save control taxa data
    control_taxa_nreads <- process_control_taxa(tax_data)

    # Create damage plot
    plot_damage(control_samples, tax_data)
}

# Run the main function
main()
