# Required libraries
library(tidyverse)
library(janitor)
library(ggthemr)
library(ggpubr)
library(phytools)
library(lvplot)
library(phyloseq)
library(ggdensity)
library(showtext)

# Load custom library
source("libs/lib.R")

# Initialize settings
showtext_auto()

# Constants
DATA_DIR <- "./data"
RESULTS_DIR <- "./results"
ABUNDANCE_THRESHOLD <- 0.01
COLOR_MAPPING <- c(
    "d__Archaea" = "#9D443C",
    "d__Bacteria" = "#74AFB8",
    "d__Viruses" = "#EECC66",
    "d__Eukaryota" = "#009E73"
)

#' Load and process initial data
#' @return list containing required data frames
load_initial_data <- function() {
    # Load metadata
    kapk_cdata <- read_tsv("./data/cdata/KapK-cdata-manuscript-20221211.tsv")
    kapk_cdata_agg <- read_tsv("./data/cdata/KapK-cdata-agg-manuscript-20221211.tsv")

    # Load taxonomic annotations
    tax_info <- read_tsv("./data/taxonomy/hires-organelles-viruses-arctic.tax.tsv",
        col_names = c("reference", "tax_string")
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
                "species",
                "strain"
            )
        )

    # Load aggregated taxonomy data
    tax_data_agg <- read_tsv("./results/taxonomy/tp-mapping-filtered.nocontam.10M.agg.nobloom.tax.tsv.gz")
    tax_data <- read_tsv("./results/taxonomy/tp-mapping-filtered.nocontam.10M.nobloom.tax.tsv.gz")

    # Load damage data
    dmg_local_agg <- read_csv("./results/damage/tp-mdmg.weight-1.local.10M.agg.nobloom.csv.gz")
    dmg_local <- read_csv("./results/damage/tp-mdmg.weight-1.local.10M.nobloom.csv.gz")

    # Load and process derep stats
    derep_stats <- read_tsv("./data/stats/all.stats-derep-summary.tsv.gz") |>
        inner_join(kapk_cdata) |>
        select(label = figure_names, member_unit, num_seqs, site_rnk, avg_len) |>
        group_by(label, member_unit, site_rnk) |>
        summarise(
            num_seqs = janitor::round_half_up(mean(num_seqs)),
            avg_len = mean(avg_len)
        ) |>
        ungroup() |>
        mutate(step = "derep") |>
        filter(label %in% (tax_data_agg |> distinct() |> pull(label)))

    # Load and process mapped reads data
    label_to_file <- read_tsv("./data/cdata/KapK-label-to-file-20221211.tsv")
    reads_mapped <- read_tsv("./data/cdata/all.tp-classified-reads.tsv",
        col_names = c("file", "reads_mapped")
    ) |>
        inner_join(label_to_file) |>
        inner_join(kapk_cdata |> select(label, figure_names)) |>
        select(-label, -file) |>
        rename(label = figure_names)

    # Process derep stats with mapped reads
    derep_stats <- derep_stats |>
        inner_join(reads_mapped) |>
        mutate(
            prop_mapped = (reads_mapped / num_seqs),
            not_mapped = num_seqs - reads_mapped
        )

    list(
        kapk_cdata = kapk_cdata,
        kapk_cdata_agg = kapk_cdata_agg,
        tax_info = tax_info,
        tax_data_agg = tax_data_agg,
        tax_data = tax_data,
        dmg_local_agg = dmg_local_agg,
        dmg_local = dmg_local,
        derep_stats = derep_stats
    )
}
#' Process taxonomic data
#' @param tax_data_agg Aggregated taxonomy data
#' @param tax_info Taxonomy information
#' @param dmg_local_agg Aggregated damage data
#' @param kapk_cdata_agg Aggregated sample metadata
#' @return list containing processed taxonomic data
process_taxonomic_data <- function(tax_data_agg, tax_info, dmg_local_agg, kapk_cdata_agg) {
    # Combine taxonomy and damage data
    tax_data_agg <- tax_data_agg |>
        inner_join(dmg_local_agg)

    # Process viruses data
    tax_data_agg_viruses <- tax_data_agg |>
        inner_join(tax_info) |>
        filter(domain == "d__Viruses")

    # Process species-level data
    tax_data_agg_dmg_sp <- tax_data_agg |>
        inner_join(tax_info) |>
        filter(domain != "d__Eukaryota") |>
        select(
            label = figure_names,
            member_unit, is_dmg, abundance,
            site_rnk, domain, species,
            damage, read_ani_mean, breadth
        ) |>
        group_by(label, member_unit, is_dmg, site_rnk, domain, species) |>
        summarise(
            abundance = janitor::round_half_up(gm_mean(abundance)),
            breadth = mean(breadth),
            damage = mean(damage),
            read_ani_mean = mean(read_ani_mean)
        ) |>
        ungroup()

    # Get short label order
    short_label_order <- kapk_cdata_agg |>
        filter(label %in% (tax_data_agg_dmg_sp |> distinct() |> pull(label))) |>
        mutate(label = fct_reorder(label, -site_rnk)) |>
        pull(label) |>
        levels() |>
        enframe(value = "label") |>
        inner_join(kapk_cdata_agg |> select(label, short_label)) |>
        pull(short_label)

    list(
        tax_data_agg = tax_data_agg,
        tax_data_agg_viruses = tax_data_agg_viruses,
        tax_data_agg_dmg_sp = tax_data_agg_dmg_sp,
        short_label_order = short_label_order
    )
} #' Process taxonomic data
#' @param tax_data_agg Aggregated taxonomy data
#' @param tax_info Taxonomy information
#' @param dmg_local_agg Aggregated damage data
#' @param kapk_cdata_agg Aggregated sample metadata
#' @return list containing processed taxonomic data
process_taxonomic_data <- function(tax_data_agg, tax_info, dmg_local_agg, kapk_cdata_agg) {
    # Combine taxonomy and damage data
    tax_data_agg <- tax_data_agg |>
        inner_join(dmg_local_agg)

    # Process viruses data
    tax_data_agg_viruses <- tax_data_agg |>
        inner_join(tax_info) |>
        filter(domain == "d__Viruses")

    # Process species-level data
    tax_data_agg_dmg_sp <- tax_data_agg |>
        inner_join(tax_info) |>
        filter(domain != "d__Eukaryota") |>
        inner_join(kapk_cdata_agg) |>
        select(
            label, short_label, domain, lineage, kingdom, phylum, class, order,
            family, genus, species, site, member_unit, site_rnk, abundance, is_dmg,
            damage, read_ani_mean, breadth
        ) |>
        group_by(
            label, short_label, domain, lineage, kingdom, phylum, class, order,
            family, genus, species, site, member_unit, site_rnk, is_dmg
        ) |>
        summarise(
            abundance = janitor::round_half_up(gm_mean(abundance)),
            damage = mean(damage),
            read_ani_mean = mean(read_ani_mean),
            breadth = mean(breadth)
        ) |>
        ungroup()

    # Get short label order
    short_label_order <- kapk_cdata_agg |>
        filter(label %in% (tax_data_agg_dmg_sp |> distinct() |> pull(label))) |>
        mutate(label = fct_reorder(label, -site_rnk)) |>
        pull(label) |>
        levels() |>
        enframe(value = "label") |>
        inner_join(kapk_cdata_agg |> select(label, short_label)) |>
        pull(short_label)

    list(
        tax_data_agg = tax_data_agg,
        tax_data_agg_viruses = tax_data_agg_viruses,
        tax_data_agg_dmg_sp = tax_data_agg_dmg_sp,
        short_label_order = short_label_order
    )
}

#' Process genus-level data for damaged taxa
#' @param tax_data_agg_dmg_sp Species-level damage data
#' @param abun_thresh Abundance threshold
#' @return list containing genus-level data
process_genus_level_data <- function(tax_data_agg_dmg_sp, abun_thresh = ABUNDANCE_THRESHOLD) {
    # Process family level data
    tax_data_agg_family <- tax_data_agg_dmg_sp |>
        filter(domain %in% c("d__Archaea", "d__Bacteria")) |>
        select(label, domain, phylum, class, order, family, abundance, is_dmg) |>
        group_by(label, domain, phylum, class, order, family, is_dmg) |>
        summarise(abundance = sum(abundance)) |>
        ungroup() |>
        group_by(label) |>
        mutate(total_abundance = sum(abundance)) |>
        ungroup() |>
        mutate(abundance = abundance / total_abundance)

    # Get damaged families
    dmg_family <- tax_data_agg_family |>
        filter(abundance > abun_thresh, is_dmg == "Damaged") |>
        select(domain, phylum, class, order, family) |>
        distinct()

    # Process genus level for damaged families
    tax_data_agg_dmg_genus <- tax_data_agg_dmg_sp |>
        filter(domain %in% c("d__Archaea", "d__Bacteria")) |>
        select(label, domain, phylum, class, order, family, genus, abundance, is_dmg) |>
        group_by(label, domain, phylum, class, order, family, genus, is_dmg) |>
        summarise(abundance = sum(abundance)) |>
        ungroup() |>
        group_by(label) |>
        mutate(total_abundance = sum(abundance)) |>
        ungroup() |>
        mutate(abundance = abundance / total_abundance) |>
        filter(family %in% dmg_family$family, is_dmg == "Damaged")

    list(
        family_data = tax_data_agg_family,
        dmg_family = dmg_family,
        dmg_genus = tax_data_agg_dmg_genus
    )
}

#' Process tree data and save results
#' @param dmg_genus Damaged genus data
#' @param tree_path Path to GTDB tree
#' @param output_path Output path for filtered tree
#' @return Filtered tree object
process_tree_data <- function(dmg_genus, tree_path, output_path) {
    # Read and process tree
    gtdb_tree <- read_gtdb_tree(tree_path)
    gtdb_tree$node.label <- NULL
    gtdb_tree$tip.label <- gsub('"', "", gtdb_tree$tip.label)
    gtdb_tree$tip.label <- gsub(" ", "_", gtdb_tree$tip.label)

    # Filter tree
    to_drop <- setdiff(gtdb_tree$tip.label, dmg_genus$genus)
    gtdb_tree_filt <- ape::drop.tip(gtdb_tree, to_drop)
    gtdb_tree_filt$tip.label <- gsub("g__", "", gtdb_tree_filt$tip.label)

    # Save filtered tree
    write.tree(gtdb_tree_filt, file = output_path)

    gtdb_tree_filt
}

#' Create phyloseq object
#' @param tax_data_agg Aggregated taxonomy data
#' @param tax_info Taxonomy information
#' @param kapk_cdata_agg Aggregated sample metadata
#' @return phyloseq object
create_phyloseq_object <- function(tax_data_agg, tax_info, kapk_cdata_agg) {
    # Filter damaged taxa
    tax_filt <- tax_data_agg |>
        filter(is_dmg == "Damaged") |>
        inner_join(tax_info)

    # Process bacteria and archaea data
    tax_filt_ba <- tax_filt |>
        filter(domain %in% c("d__Bacteria", "d__Archaea")) |>
        group_by(label, species) |>
        summarise(
            abundance = janitor::round_half_up(gm_mean(abundance))
        ) |>
        ungroup()

    # Create OTU table
    tax_df <- tax_filt_ba |>
        select(species, label, abundance) |>
        pivot_wider(
            names_from = "label",
            values_from = abundance,
            values_fill = 0
        ) |>
        as.data.frame() |>
        column_to_rownames("species")

    # Create sample data
    kapk_cdata_df <- kapk_cdata_agg |>
        mutate(label_cl = label) |>
        as.data.frame() |>
        column_to_rownames("label_cl")

    # Create taxonomy table
    taxonomy_df <- tax_info |>
        filter(domain %in% c("d__Bacteria", "d__Archaea")) |>
        select(species, domain, phylum, class, order, family, genus) |>
        distinct() |>
        mutate(reference = species) |>
        distinct() |>
        as.data.frame() |>
        column_to_rownames("reference")

    # Create phyloseq object
    phyloseq(
        otu_table(tax_df, taxa_are_rows = TRUE),
        tax_table(as.matrix(taxonomy_df)),
        sample_data(kapk_cdata_df)
    )
}

#' Create reference counts plot
#' @param tax_data_agg Aggregated taxonomy data
#' @param kapk_cdata_agg Aggregated sample metadata
#' @param short_label_order Order of short labels
#' @return ggplot object
plot_reference_counts <- function(tax_data_agg, kapk_cdata_agg, short_label_order) {
    tax_data_agg |>
        select(label, reference, is_dmg) |>
        group_by(label, is_dmg) |>
        count() |>
        ungroup() |>
        inner_join(kapk_cdata_agg) |>
        mutate(
            short_label = fct_relevel(short_label, short_label_order),
            member_unit = fct_relevel(member_unit, c("B3", "B2", "B1"))
        ) |>
        ggplot(aes(x = short_label, y = n)) +
        geom_col(color = "black", linewidth = 0.3, width = 1) +
        theme_bw() +
        xlab("Sample") +
        ylab("Number of references") +
        theme(
            legend.position = "top",
            legend.title = element_blank(),
            text = element_text(size = 10),
            strip.background = element_blank(),
            strip.text.y = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        facet_grid(member_unit ~ is_dmg, scales = "free", space = "free_y") +
        ggpubr::rotate()
}

#' Create dereplicated sequences plot
#' @param derep_stats Dereplicated statistics
#' @param kapk_cdata_agg Aggregated sample metadata
#' @param short_label_order Order of short labels
#' @return ggplot object
plot_derep_seqs <- function(derep_stats, kapk_cdata_agg, short_label_order) {
    derep_stats |>
        inner_join(kapk_cdata_agg |> select(label, short_label)) |>
        select(label, short_label, member_unit, site_rnk, num_seqs, step) |>
        mutate(
            label = fct_reorder(label, -site_rnk),
            member_unit = fct_relevel(member_unit, c("B3", "B2", "B1")),
            short_label = fct_relevel(short_label, short_label_order)
        ) |>
        ggplot(aes(label, num_seqs)) +
        geom_col(
            color = "black", linewidth = 0.3, width = 1,
            alpha = 1, fill = "#66767A"
        ) +
        theme_bw() +
        xlab("Sample") +
        ylab("# reads") +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            strip.background = element_blank(),
            strip.text.y = element_blank(),
            legend.position = "none",
            legend.title = element_blank(),
            text = element_text(size = 10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        facet_grid(member_unit ~ step, scales = "free", space = "free") +
        scale_y_continuous(labels = scales::label_number(suffix = "M", scale = 1e-6)) +
        ggpubr::rotate()
}

#' Create mapped sequences proportion plot
#' @param derep_stats Dereplicated statistics
#' @param kapk_cdata_agg Aggregated sample metadata
#' @param short_label_order Order of short labels
#' @return ggplot object
plot_mapped_proportions <- function(derep_stats, kapk_cdata_agg, short_label_order) {
    derep_stats |>
        inner_join(kapk_cdata_agg |> select(label, short_label)) |>
        select(label, short_label, member_unit, site_rnk, prop_mapped, step) |>
        mutate(
            label = fct_reorder(label, -site_rnk),
            member_unit = fct_relevel(member_unit, c("B3", "B2", "B1")),
            short_label = fct_relevel(short_label, short_label_order)
        ) |>
        ggplot(aes(label, prop_mapped)) +
        geom_col(
            color = "black", linewidth = 0.3, width = 1,
            alpha = 1, fill = "#C75F44"
        ) +
        theme_bw() +
        xlab("Sample") +
        ylab("Reads mapped") +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            strip.background = element_blank(),
            legend.position = "none",
            legend.title = element_blank(),
            text = element_text(size = 10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        facet_grid(member_unit ~ step, scales = "free", space = "free") +
        scale_y_continuous(labels = scales::percent) +
        ggpubr::rotate()
}

#' Create domain proportions plot
#' @param tax_data_agg_dmg_sp Species-level damage data
#' @param short_label_order Order of short labels
#' @return ggplot object
plot_domain_proportions <- function(tax_data_agg_dmg_sp, short_label_order) {
    tax_props <- tax_data_agg_dmg_sp |>
        group_by(short_label, domain, site, member_unit, site_rnk, is_dmg) |>
        summarise(abundance = sum(abundance)) |>
        ungroup() |>
        group_by(short_label) |>
        mutate(
            total_abundance = sum(abundance),
            abundance = abundance / total_abundance
        ) |>
        ungroup() |>
        mutate(
            short_label = fct_relevel(short_label, short_label_order),
            member_unit = fct_relevel(member_unit, c("B3", "B2", "B1")),
            domain = fct_relevel(domain, c("d__Viruses", "d__Bacteria", "d__Archaea"))
        )

    ggplot(tax_props, aes(x = short_label, y = abundance, fill = domain)) +
        geom_col(position = "stack", color = "black", linewidth = 0.3, width = 1) +
        theme_bw() +
        xlab("Sample") +
        ylab("Proportion") +
        theme(
            legend.position = "top",
            legend.title = element_blank(),
            text = element_text(size = 10),
            strip.background = element_blank(),
            strip.text.y = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        facet_grid(member_unit ~ is_dmg, scales = "free_y", space = "free") +
        scale_y_continuous(labels = scales::percent) +
        scale_fill_manual(values = COLOR_MAPPING) +
        ggpubr::rotate()
}

#' Create HDR plots
#' @param tax_data_agg Aggregated taxonomy data
#' @param tax_info Taxonomy information
#' @param kapk_cdata_agg Aggregated sample metadata
#' @param dmg_thresholds Damage thresholds
#' @return list of ggplot objects
create_hdr_plots <- function(tax_data_agg, tax_info, kapk_cdata_agg, dmg_thresholds) {
    # Prepare data for HDR plots
    plot_data <- tax_data_agg |>
        inner_join(tax_info) |>
        inner_join(kapk_cdata_agg) |>
        filter(domain != "d__Eukaryota") |>
        select(label, domain, abundance, species, damage, read_ani_mean, breadth) |>
        group_by(label, domain, species) |>
        summarise(
            abundance = janitor::round_half_up(gm_mean(abundance)),
            breadth = mean(breadth),
            damage = mean(damage),
            read_ani_mean = mean(read_ani_mean),
            .groups = "drop"
        ) |>
        group_by(label) |>
        mutate(total_abundance = sum(abundance)) |>
        ungroup() |>
        mutate(
            abundance = abundance / total_abundance,
            domain = fct_relevel(domain, c("d__Viruses", "d__Bacteria", "d__Archaea"))
        )

    # Abundance HDR plot
    hdr_abun <- plot_data |>
        ggplot(aes(y = abundance, x = damage)) +
        annotate("rect",
            xmin = dmg_thresholds$damage, xmax = Inf,
            ymin = 0, ymax = Inf,
            fill = "#F2B705", alpha = 0.2
        ) +
        annotate("rect",
            xmin = 0, xmax = dmg_thresholds$damage,
            ymin = 0, ymax = Inf,
            fill = "#1C7085", alpha = 0.2
        ) +
        geom_hdr() +
        geom_vline(
            xintercept = dmg_thresholds$damage,
            linetype = "dashed", color = "#B04035", alpha = 0.5
        ) +
        theme_bw() +
        labs(y = "Proportion", x = "Damage") +
        theme(
            legend.position = "top",
            legend.title = element_blank(),
            text = element_text(size = 12),
            strip.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.text.x = element_blank()
        ) +
        facet_grid(~domain) +
        scale_y_log10(
            labels = scales::percent_format(accuracy = 1),
            breaks = c(0.01, 0.1, 1),
            limits = c(3.7e-05, 1)
        )

    # Breadth HDR plot
    hdr_breadth <- plot_data |>
        ggplot(aes(y = breadth, x = damage)) +
        annotate("rect",
            xmin = dmg_thresholds$damage, xmax = Inf,
            ymin = 0, ymax = Inf,
            fill = "#F2B705", alpha = 0.2
        ) +
        annotate("rect",
            xmin = 0, xmax = dmg_thresholds$damage,
            ymin = 0, ymax = Inf,
            fill = "#1C7085", alpha = 0.2
        ) +
        geom_hdr() +
        geom_vline(
            xintercept = dmg_thresholds$damage,
            linetype = "dashed", color = "#B04035", alpha = 0.5
        ) +
        theme_bw() +
        labs(y = "Detection", x = "Damage") +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(size = 12),
            legend.position = "top",
            legend.title = element_blank(),
            text = element_text(size = 12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.text.x = element_blank()
        ) +
        facet_grid(~domain) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1))

    # ANI HDR plot
    hdr_ani <- plot_data |>
        ggplot(aes(y = read_ani_mean / 100, x = damage)) +
        annotate("rect",
            xmin = dmg_thresholds$damage, xmax = Inf,
            ymin = -Inf, ymax = Inf,
            fill = "#F2B705", alpha = 0.2
        ) +
        annotate("rect",
            xmin = 0, xmax = dmg_thresholds$damage,
            ymin = -Inf, ymax = Inf,
            fill = "#1C7085", alpha = 0.2
        ) +
        geom_hdr() +
        geom_vline(
            xintercept = dmg_thresholds$damage,
            linetype = "dashed", color = "#B04035", alpha = 0.5
        ) +
        theme_bw() +
        labs(y = "ANI", x = "Damage") +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(size = 12),
            legend.position = "top",
            legend.title = element_blank(),
            text = element_text(size = 12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        facet_grid(~domain) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1))

    list(
        abundance = hdr_abun,
        breadth = hdr_breadth,
        ani = hdr_ani
    )
}

#' Main execution function
#' @return list containing all results and plots
main <- function() {
    # Load initial data
    data <- load_initial_data()

    # Process taxonomic data
    tax_processed <- process_taxonomic_data(
        data$tax_data_agg,
        data$tax_info,
        data$dmg_local_agg,
        data$kapk_cdata_agg
    )

    # Create reference counts plot
    ref_counts_plot <- plot_reference_counts(
        tax_processed$tax_data_agg,
        data$kapk_cdata_agg,
        tax_processed$short_label_order
    )

    # Create sequence plots
    derep_plot <- plot_derep_seqs(
        data$derep_stats,
        data$kapk_cdata_agg,
        tax_processed$short_label_order
    )

    mapped_plot <- plot_mapped_proportions(
        data$derep_stats,
        data$kapk_cdata_agg,
        tax_processed$short_label_order
    )

    # Create domain proportions plot
    domain_plot <- plot_domain_proportions(
        tax_processed$tax_data_agg_dmg_sp,
        tax_processed$short_label_order
    )

    # Process genus-level data
    genus_data <- process_genus_level_data(
        tax_processed$tax_data_agg_dmg_sp
    )

    # Process tree data
    gtdb_tree_filt <- process_tree_data(
        genus_data$dmg_genus,
        file.path(DATA_DIR, "taxonomy/gtdb-r207_genus.tree"),
        file.path(RESULTS_DIR, "taxonomy/anvio/gtdb-r207_dmg-genus_filt.tree")
    )

    # Load damage thresholds
    dmg_thresholds <- read_tsv(
        file.path(RESULTS_DIR, "damage/dmg_thresholds.tsv")
    )

    # Create HDR plots
    hdr_plots <- create_hdr_plots(
        tax_processed$tax_data_agg,
        data$tax_info,
        data$kapk_cdata_agg,
        dmg_thresholds
    )

    # Create phyloseq object
    ps_obj <- create_phyloseq_object(
        tax_processed$tax_data_agg,
        data$tax_info,
        data$kapk_cdata_agg
    )

    # Save phyloseq object
    saveRDS(ps_obj, file.path(RESULTS_DIR, "taxonomy/kapk_ps_ba_gm.rds"))

    # Process and save anvio data
    tax_data_agg_dmg_genus <- genus_data$dmg_genus |>
        mutate(genus = gsub("g__", "", genus)) |>
        filter(genus %in% gtdb_tree_filt$tip.label) |>
        mutate(item_name = genus) |>
        select(item_name, label, domain, phylum, class, order, family, genus, abundance) |>
        mutate(abundance = 100 * abundance) |>
        pivot_wider(
            names_from = label,
            values_from = abundance,
            values_fill = list(abundance = 0)
        ) |>
        mutate(B1 = "B1", B2 = "B2", B3 = "B3")

    write_tsv(
        tax_data_agg_dmg_genus,
        file.path(RESULTS_DIR, "taxonomy/anvio/gtdb_genus_filt_dmg-abundance.tsv")
    )

    # Combine plots for output
    combined_plot <- ggpubr::ggarrange(
        domain_plot, derep_plot, mapped_plot,
        ncol = 3, nrow = 1,
        common.legend = TRUE,
        legend = "none",
        align = "h",
        heights = c(1, 1, 1),
        widths = c(0.67, 0.16, 0.19)
    )

    hdr_combined_plot <- ggpubr::ggarrange(
        hdr_plots$ani,
        hdr_plots$breadth,
        hdr_plots$abundance,
        ncol = 1, nrow = 3,
        align = "hv",
        legend = "none"
    )

    # Return results
    list(
        data = list(
            processed_tax = tax_processed,
            genus_data = genus_data,
            tree = gtdb_tree_filt,
            phyloseq = ps_obj,
            anvio_data = tax_data_agg_dmg_genus
        ),
        plots = list(
            reference_counts = ref_counts_plot,
            derep_sequences = derep_plot,
            mapped_proportions = mapped_plot,
            domain_proportions = domain_plot,
            hdr_plots = hdr_plots,
            combined = combined_plot,
            hdr_combined = hdr_combined_plot
        )
    )
}

# Run main function and store results
showtext::showtext_opts(dpi = 300)
results <- main()
