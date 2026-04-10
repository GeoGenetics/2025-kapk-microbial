library(tidyverse)
library(lvplot)
library(showtext)
library(unixtools)

# Load custom library
source("libs/lib.R")

# Constants
DATA_DIR <- "./data"
RESULTS_DIR <- "./results"
SIGNIFICANCE_THRESHOLD <- 2
BREADTH_FILTER <- 0.01
DOMAINS <- c("d__Eukaryota", "d__Bacteria", "d__Archaea", "d__Viruses")

# Initialize settings
showtext_auto()
set.tempdir("/maps/projects/fernandezguerra/scratch/tmp")

#' Load and process metadata and taxonomy files
load_initial_data <- function() {
    # Load metadata
    kapk_cdata <- read_tsv(
        file.path(DATA_DIR, "cdata/KapK-cdata-manuscript-20221211.tsv"),
        show_col_types = FALSE
    )
    kapk_cdata_agg <- read_tsv(
        file.path(DATA_DIR, "cdata/KapK-cdata-agg-manuscript-20221211.tsv"),
        show_col_types = FALSE
    )

    # Load filtered taxonomy data
    tax_data <- read_tsv(
        file.path(RESULTS_DIR, "taxonomy/tp-mapping-filtered.nocontam.10M.tax.tsv.gz"),
        show_col_types = FALSE
    ) |>
        inner_join(kapk_cdata)

    # Load taxonomic annotations
    tax_info <- read_tsv(
        file.path(DATA_DIR, "taxonomy/hires-organelles-viruses-arctic.tax.tsv"),
        col_names = c("reference", "tax_string"),
        show_col_types = FALSE
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

    # Load damage data
    dmg_local <- read_csv(file.path(DATA_DIR, "taxonomy/tp-mdmg.weight-1.csv.gz")) |>
        filter(label %in% kapk_cdata$label)

    list(
        kapk_cdata = kapk_cdata,
        kapk_cdata_agg = kapk_cdata_agg,
        tax_data = tax_data,
        tax_info = tax_info,
        dmg_local = dmg_local
    )
}

#' Create and save damage thresholds
save_damage_thresholds <- function(thresholds) {
    tibble(
        domain = "d__Eukaryota",
        rank = NA,
        mode = "local",
        damage = thresholds$damage,
        signf = thresholds$signf,
        n_reads = thresholds$n_reads
    ) |>
        write_tsv(file.path(RESULTS_DIR, "damage/dmg_thresholds.tsv"))
}

#' Analyze Eukaryotic damage distribution
analyze_euk_damage <- function(dmg_local, tax_data, tax_info, signf) {
    dmg_local_euk <- dmg_local |>
        filter(significance >= signf) |>
        mutate(reference = tax_id) |>
        inner_join(tax_data) |>
        inner_join(tax_info) |>
        filter(domain == "d__Eukaryota")

    # Calculate damage statistics
    dmg_local_euk_stats <- dmg_local_euk |>
        pull(damage) |>
        lsum(l = 8)

    dmg_threshold <- dmg_local_euk_stats |>
        filter(letter == "C") |>
        pull(lower)

    # Calculate read number statistics
    dmg_local_euk_nreads <- dmg_local_euk |>
        filter(damage >= dmg_threshold) |>
        pull(n_reads) |>
        lsum(l = 8)

    nreads_threshold <- dmg_local_euk_nreads |>
        filter(letter == "Z") |>
        pull(lower)

    # Create thresholds object
    thresholds <- list(
        damage = dmg_threshold,
        n_reads = nreads_threshold,
        signf = signf
    )

    # Save thresholds file
    save_damage_thresholds(thresholds)

    list(
        damage_threshold = dmg_threshold,
        nreads_threshold = nreads_threshold,
        euk_data = dmg_local_euk,
        thresholds = thresholds
    )
}

#' Plot read number distribution
plot_reads_distribution <- function(dmg_local, tax_info, tax_data, thresholds) {
    dmg_local |>
        filter(significance >= thresholds$signf) |>
        select(label, tax_id) |>
        rename(reference = tax_id) |>
        inner_join(tax_info) |>
        inner_join(tax_data) |>
        mutate(domain = fct_relevel(domain, rev(DOMAINS))) |>
        ggplot(aes(x = domain, y = n_reads)) +
        geom_lv(aes(fill = after_stat(LV)),
            size = 0.5,
            width.method = "height",
            color = "#404040",
            width = 0.5,
            alpha = 1
        ) +
        geom_hline(
            yintercept = thresholds$n_reads,
            linetype = "dashed",
            color = "red"
        ) +
        coord_flip() +
        theme_bw() +
        labs(y = "Number of reads", x = "") +
        theme(
            legend.position = "top",
            legend.title = element_blank(),
            text = element_text(size = 16)
        ) +
        guides(fill = guide_legend(nrow = 2)) +
        scale_y_log10(label = scales::comma)
}

#' Plot damage distribution
plot_damage_distribution <- function(dmg_local, tax_info, tax_data, thresholds) {
    dmg_local |>
        select(label, tax_id, damage, significance) |>
        rename(reference = tax_id) |>
        inner_join(tax_info) |>
        inner_join(tax_data) |>
        filter(
            significance >= thresholds$signf,
            n_reads >= thresholds$n_reads
        ) |>
        mutate(domain = fct_relevel(domain, rev(DOMAINS))) |>
        ggplot(aes(x = domain, y = damage)) +
        geom_lv(aes(fill = after_stat(LV)),
            size = 0.5,
            width.method = "height",
            color = "#404040",
            width = 0.5,
            alpha = 1
        ) +
        geom_hline(
            yintercept = thresholds$damage,
            linetype = "dashed",
            color = "red"
        ) +
        coord_flip() +
        theme_bw() +
        labs(y = "Damage", x = "") +
        theme(
            legend.position = "top",
            legend.title = element_blank(),
            text = element_text(size = 12)
        ) +
        guides(fill = guide_legend(nrow = 2))
}

#' Process and filter damage data
process_damage_data <- function(dmg_local, tax_data, thresholds) {
    dmg_local |>
        mutate(reference = tax_id) |>
        inner_join(tax_data) |>
        mutate(
            is_dmg = ifelse(
                (significance >= thresholds$signf &
                    damage >= thresholds$damage &
                    n_reads >= thresholds$n_reads),
                "Damaged",
                "Non-damaged"
            )
        )
}

#' Plot reference counts by damage status
plot_reference_counts <- function(dmg_data, tax_info) {
    dmg_data |>
        inner_join(tax_info) |>
        filter(domain != "d__Eukaryota") |>
        select(label = figure_names, reference, site_rnk, member_unit, is_dmg) |>
        distinct() |>
        group_by(label, site_rnk, member_unit, is_dmg) |>
        count(sort = TRUE) |>
        ungroup() |>
        arrange(desc(n)) |>
        mutate(
            label = fct_reorder(label, site_rnk),
            member_unit = fct_relevel(member_unit, c("B3", "B2", "B1"))
        ) |>
        ggplot(aes(label, n, fill = is_dmg)) +
        geom_col(position = position_dodge()) +
        facet_grid(~member_unit, space = "free", scales = "free_x") +
        theme_bw() +
        labs(y = "Number of references", x = "") +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "top",
            legend.title = element_blank(),
            text = element_text(size = 16)
        ) +
        guides(fill = guide_legend(nrow = 1)) +
        scale_fill_manual(values = c("Non-damaged" = "#1f78b4", "Damaged" = "#e31a1c"))
}

#' Calculate taxonomic proportions
calculate_taxonomic_proportions <- function(dmg_data, tax_info) {
    dmg_data |>
        inner_join(tax_info) |>
        filter(domain != "d__Eukaryota") |>
        select(
            label = figure_names,
            member_unit,
            is_dmg,
            abundance,
            site_rnk,
            domain,
            species
        ) |>
        group_by(label, member_unit, is_dmg, site_rnk, domain, species) |>
        summarise(
            abundance = janitor::round_half_up(gm_mean(abundance)),
            .groups = "drop"
        ) |>
        group_by(label, member_unit) |>
        mutate(abundance = abundance / sum(abundance)) |>
        ungroup() |>
        mutate(
            label = fct_reorder(label, site_rnk),
            member_unit = fct_relevel(member_unit, c("B3", "B2", "B1")),
            # Ensure domain is a factor in the correct order
            domain = factor(domain, levels = DOMAINS[-1]) # Remove Eukaryota
        )
}

#' Plot taxonomic proportions
plot_taxonomic_proportions <- function(tax_props) {
    tax_props |>
        # Ensure domain is present and grouped correctly
        group_by(label, member_unit, is_dmg, site_rnk, domain) |>
        summarise(abundance = sum(abundance), .groups = "drop") |>
        ggplot(aes(x = label, y = abundance, fill = domain)) +
        geom_col(position = "stack", color = "black", linewidth = 0.3, width = 1) +
        geom_hline(yintercept = 0.5, linetype = 2) +
        facet_grid(is_dmg ~ member_unit, space = "free", scales = "free_x") +
        theme_bw() +
        labs(y = "Proportion", x = "") +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "top",
            legend.title = element_blank(),
            text = element_text(size = 16)
        ) +
        guides(fill = guide_legend(nrow = 1)) +
        scale_fill_manual(
            values = c(
                "d__Archaea" = "#9D443C",
                "d__Bacteria" = "#74AFB8",
                "d__Viruses" = "#EECC66"
            )
        ) +
        scale_y_continuous(labels = scales::percent_format())
}

#' Process and save species-level aggregations
process_species_level <- function(tax_data_filt, dmg_local_filt, tax_info,
                                  kapk_cdata, thresholds) {
    # Aggregate taxonomy data at species level
    tax_data_filt_sp <- tax_data_filt |>
        inner_join(tax_info) |>
        group_by(label, species, domain) |>
        summarise(
            abundance = janitor::round_half_up(gm_mean(abundance)),
            read_ani_mean = mean(read_ani_mean),
            coverage_mean_trunc = mean(coverage_mean_trunc),
            breadth = mean(breadth),
            coverage_mean = mean(coverage_mean),
            .groups = "drop"
        ) |>
        inner_join(kapk_cdata |> select(label, figure_names)) |>
        select(-label) |>
        rename(label = figure_names)

    write_tsv(
        tax_data_filt_sp,
        file.path(RESULTS_DIR, "taxonomy/tp-mapping-filtered.nocontam.10M.agg.sp.tax.tsv.gz")
    )

    # Aggregate damage data at species level
    dmg_local_filt_sp <- dmg_local_filt |>
        rename(reference = tax_id) |>
        inner_join(tax_data_filt |> select(reference, label, n_reads)) |>
        inner_join(tax_info) |>
        group_by(label, species, domain) |>
        summarise(
            damage = mean(damage),
            n_reads_agg = janitor::round_half_up(mean(n_reads)),
            significance = mean(significance),
            .groups = "drop"
        ) |>
        mutate(
            is_dmg = if_else(
                damage >= thresholds$damage &
                    significance >= thresholds$signf &
                    n_reads_agg >= thresholds$n_reads,
                "Damaged",
                "Non-damaged"
            )
        ) |>
        inner_join(kapk_cdata |> select(label, figure_names)) |>
        select(-label) |>
        rename(label = figure_names)

    write_csv(
        dmg_local_filt_sp,
        file.path(RESULTS_DIR, "damage/tp-mdmg.weight-1.local.10M.agg.sp.csv.gz")
    )
}

#' Filter and save results based on damage proportions
save_filtered_results <- function(tax_data, dmg_local, tax_proportions, kapk_cdata) {
    # Identify samples to keep
    samples_to_keep <- tax_proportions |>
        group_by(label, is_dmg) |>
        summarise(abundance = sum(abundance), .groups = "drop") |>
        filter(abundance >= 0.5, is_dmg == "Damaged")

    write_tsv(
        samples_to_keep,
        file.path(DATA_DIR, "cdata/KapK-cdata-manuscript-20221211-labels-nobloom.tsv")
    )

    # Get labels to keep
    labels_to_keep <- kapk_cdata |>
        filter(figure_names %in% samples_to_keep$label)

    # Filter and save taxonomy data
    tax_data_filt <- tax_data |>
        filter(label %in% labels_to_keep$label)
    write_tsv(
        tax_data_filt,
        file.path(RESULTS_DIR, "taxonomy/tp-mapping-filtered.nocontam.10M.nobloom.tax.tsv.gz")
    )

    # Filter and save damage data
    dmg_local_filt <- dmg_local |>
        filter(label %in% labels_to_keep$label)
    write_csv(
        dmg_local_filt,
        file.path(RESULTS_DIR, "damage/tp-mdmg.weight-1.local.10M.nobloom.csv.gz")
    )

    list(
        tax_data_filt = tax_data_filt,
        dmg_local_filt = dmg_local_filt,
        samples_to_keep = samples_to_keep
    )
}

#' Main execution function
main <- function() {
    # Load initial data
    data <- load_initial_data()

    # Analyze Eukaryotic damage
    euk_analysis <- analyze_euk_damage(
        data$dmg_local,
        data$tax_data,
        data$tax_info,
        SIGNIFICANCE_THRESHOLD
    )

    # Create visualization plots
    reads_plot <- plot_reads_distribution(
        data$dmg_local,
        data$tax_info,
        data$tax_data,
        euk_analysis$thresholds
    )

    damage_plot <- plot_damage_distribution(
        data$dmg_local,
        data$tax_info,
        data$tax_data,
        euk_analysis$thresholds
    )

    # Process damage data
    dmg_processed <- process_damage_data(
        data$dmg_local,
        data$tax_data,
        euk_analysis$thresholds
    )

    # Create reference counts plot
    ref_counts_plot <- plot_reference_counts(
        dmg_processed,
        data$tax_info
    )

    # Calculate taxonomic proportions and create plot
    tax_props <- calculate_taxonomic_proportions(dmg_processed, data$tax_info)
    tax_props_plot <- plot_taxonomic_proportions(tax_props)

    # Save filtered results
    filtered_results <- save_filtered_results(
        data$tax_data,
        data$dmg_local,
        tax_props,
        data$kapk_cdata
    )

    # Process species-level aggregations
    process_species_level(
        filtered_results$tax_data_filt,
        filtered_results$dmg_local_filt,
        data$tax_info,
        data$kapk_cdata,
        euk_analysis$thresholds
    )

    # Return all results
    list(
        plots = list(
            reads = reads_plot,
            damage = damage_plot,
            ref_counts = ref_counts_plot,
            tax_proportions = tax_props_plot
        ),
        data = list(
            euk_analysis = euk_analysis,
            processed_damage = dmg_processed,
            tax_proportions = tax_props,
            filtered_results = filtered_results
        )
    )
}

# Run the main function and store results
results <- main()
