# Required libraries
library(tidyverse)
library(showtext)
library(ggdensity)
library(ggpol)

# Initialize settings
source("./libs/lib.R")
showtext_auto()

#' Clean vector names
#' @param x Vector to clean
#' @param refactor Whether to return as factor
#' @return Cleaned vector
clean_vec <- function(x, refactor = FALSE) {
    require(magrittr, quietly = TRUE)

    if (!(is.character(x) || is.factor(x))) {
        return(x)
    }

    x_is_factor <- is.factor(x)
    old_names <- as.character(x)

    new_names <- old_names |>
        gsub("'", "", .) |>
        gsub("\"", "", .) |>
        gsub("%", "percent", .) |>
        gsub("^[ ]+", "", .) |>
        make.names(.) |>
        gsub("[.]+", "_", .) |>
        gsub("[_]+", "_", .) |>
        tolower(.) |>
        gsub("_$", "", .)

    if (x_is_factor && refactor) factor(new_names) else new_names
}

#' Load required data for analysis
#' @return list of loaded data
load_data <- function() {
    # Load metadata
    kapk_cdata <- read_tsv("./data/cdata/KapK-cdata-manuscript-20221211.tsv",
        show_col_types = FALSE
    )

    kapk_cdata_agg <- read_tsv("./data/cdata/KapK-cdata-agg-manuscript-20221211.tsv",
        show_col_types = FALSE
    )

    label_nobloom <- read_tsv("./data/cdata/KapK-cdata-manuscript-20221211-labels-nobloom.tsv",
        show_col_types = FALSE
    )

    # Load damage thresholds
    dmg_thresholds <- read_tsv("./results/damage/dmg_thresholds.tsv",
        show_col_types = FALSE
    )

    # Load taxonomic data
    gem_tax_info <- read_tsv("./data/cdata/GEM-20220926-tax-gtdbtk-complete.tsv",
        show_col_types = FALSE
    ) |>
        select(reference = accession, tax_string = taxon) |>
        separate(
            col = tax_string,
            sep = ";",
            into = c(
                "domain", "lineage", "kingdom", "phylum",
                "class", "order", "family", "genus", "species"
            ),
            fill = "right"
        )

    tg2g_tax_info <- read_tsv("./data/cdata/tg2g-20220926-tax-gtdbtk-complete.tsv",
        show_col_types = FALSE
    ) |>
        select(reference = accession, tax_string = taxon) |>
        separate(
            col = tax_string,
            sep = ";",
            into = c(
                "domain", "lineage", "kingdom", "phylum",
                "class", "order", "family", "genus", "species"
            ),
            fill = "right"
        )

    mg_tax_info <- bind_rows(gem_tax_info, tg2g_tax_info)

    # Load GEM metadata
    gem_metadata <- read_tsv("./data/cdata/GEM-20220926_genome-metadata.tsv",
        show_col_types = FALSE
    ) |>
        select(reference = genome_id, habitat)

    # Load Woodcroft data
    woodcroft_tax_info <- read_tsv("./data/cdata/woodcroft2018-tax-gtdbtk-complete.tsv",
        show_col_types = FALSE
    ) |>
        select(reference = accession, tax_string = taxon) |>
        separate(
            col = tax_string,
            sep = ";",
            into = c(
                "domain", "lineage", "kingdom", "phylum",
                "class", "order", "family", "genus", "species"
            ),
            fill = "right"
        )

    # Load Woodcroft metadata
    woodcroft_metadata <- read_tsv("./data/cdata/woodcroft2018-metadata.txt",
        show_col_types = FALSE
    )

    list(
        kapk_cdata = kapk_cdata,
        kapk_cdata_agg = kapk_cdata_agg,
        label_nobloom = label_nobloom,
        dmg_thresholds = dmg_thresholds,
        gem_tax_info = gem_tax_info,
        mg_tax_info = mg_tax_info,
        gem_metadata = gem_metadata,
        woodcroft_tax_info = woodcroft_tax_info,
        woodcroft_metadata = woodcroft_metadata
    )
}

#' Process MAG data
#' @param data Loaded data list
#' @return processed MAG data
process_mag_data <- function(data) {
    # Load and process mapping data
    mg_mapping <- read_tsv("./data/mag-distribution/GEM-20220926__tg2g-20220926-tp-mapping-filtered.summary.tsv.gz",
        show_col_types = FALSE
    ) |>
        filter(breadth >= 0.01) |>
        mutate(abundance = ifelse(tax_abund_tad == 0, tax_abund_read, tax_abund_tad))

    # Process metadata
    mg_metadata <- mg_mapping |>
        select(reference) |>
        distinct() |>
        filter(!reference %in% data$gem_tax_info$reference) |>
        mutate(habitat = "Glacier") |>
        bind_rows(data$gem_metadata)

    # Load damage data
    mg_dmg <- read_csv("./data/mag-distribution/GEM-20220926__tg2g-20220926-tp-mdmg.weight-1.csv.gz",
        show_col_types = FALSE
    ) |>
        filter(label %in% data$kapk_cdata$label) |>
        select(label, reference = tax_id, damage, significance) |>
        filter(
            significance > data$dmg_thresholds$signf,
            damage >= data$dmg_thresholds$damage
        )

    # Process species-level data
    mg_mapping_filt_sp <- mg_mapping |>
        filter(breadth >= 0.01, n_reads >= data$dmg_thresholds$n_reads) |>
        inner_join(mg_dmg) |>
        inner_join(data$kapk_cdata |> select(label, figure_names)) |>
        inner_join(data$mg_tax_info) |>
        inner_join(mg_metadata) |>
        group_by(figure_names, domain, species, habitat) |>
        summarise(
            abundance = janitor::round_half_up(gm_mean(abundance)),
            damage = mean(damage),
            read_ani_mean = mean(read_ani_mean),
            breadth = mean(breadth)
        ) |>
        ungroup() |>
        rename(label = figure_names) |>
        filter(label %in% data$label_nobloom$label)

    mg_mapping_filt_sp
}

#' Process Woodcroft data
#' @param data Loaded data list
#' @return processed Woodcroft data
process_woodcroft_data <- function(data) {
    # Load and process mapping data
    woodcroft_mapping <- read_tsv("./data/mag-distribution/woodcroft2018-tp-mapping-filtered.summary.tsv.gz",
        show_col_types = FALSE
    ) |>
        filter(breadth >= 0.01) |>
        mutate(abundance = ifelse(tax_abund_tad == 0, tax_abund_read, tax_abund_tad)) |>
        inner_join(data$woodcroft_tax_info) |>
        separate(reference, c("acc", "sample"), sep = ".1_", remove = FALSE) |>
        select(-acc) |>
        inner_join(data$woodcroft_metadata |> select(sample, habitat))

    # Load damage data
    woodcroft_dmg <- read_csv("./data/mag-distribution/woodcroft2018-tp-mdmg.weight-1.csv.gz",
        show_col_types = FALSE
    ) |>
        filter(label %in% data$kapk_cdata$label) |>
        select(label, reference = tax_id, damage, significance) |>
        filter(
            significance > data$dmg_thresholds$signf,
            damage >= data$dmg_thresholds$damage
        )

    # Process final data
    woodcroft_mapping_filt_sp <- woodcroft_mapping |>
        filter(breadth >= 0.01, n_reads >= data$dmg_thresholds$n_reads) |>
        inner_join(woodcroft_dmg) |>
        inner_join(data$kapk_cdata) |>
        select(-reference) |>
        group_by(
            figure_names, domain, lineage, kingdom, phylum, class,
            order, family, genus, species, habitat
        ) |>
        summarise(
            abundance = janitor::round_half_up(gm_mean(abundance)),
            damage = mean(damage),
            read_ani_mean = mean(read_ani_mean),
            breadth = mean(breadth)
        ) |>
        ungroup() |>
        rename(label = figure_names) |>
        filter(label %in% data$label_nobloom$label)

    woodcroft_mapping_filt_sp
}

#' Process data for habitat metrics plot
#' @param mag_data Processed MAG data
#' @param kapk_cdata_agg Aggregated metadata
#' @return list containing plot data and top habitats
process_habitat_metrics_data <- function(mag_data, kapk_cdata_agg) {
    # Get top habitats
    top_habitats <- mag_data |>
        mutate(habitat = clean_vec(habitat)) |>
        filter(breadth > 0.5) |>
        group_by(habitat) |>
        summarise(
            max_breadth = max(breadth),
            n = n()
        ) |>
        filter(n > 1) |>
        arrange(-n) |>
        head(5) |>
        mutate(rnk = row_number())

    # Process plot data
    plot_data <- mag_data |>
        inner_join(kapk_cdata_agg) |>
        select(
            label, domain, species, damage, breadth, abundance,
            member_unit, read_ani_mean, habitat
        ) |>
        group_by(label) |>
        mutate(prop = abundance / sum(abundance)) |>
        ungroup() |>
        mutate(
            habitat = clean_vec(habitat),
            habitat = ifelse(habitat %in% top_habitats$habitat, habitat, "Other")
        )

    list(
        data = plot_data,
        top_habitats = top_habitats
    )
}

#' Create habitat metrics plot
#' @param data Processed habitat data
#' @param top_habitats Top habitats data
#' @return ggplot object
plot_habitat_metrics <- function(data, top_habitats) {
    data |>
        select(species, domain, damage,
            Detection = breadth,
            prop, habitat, ANI = read_ani_mean
        ) |>
        mutate(ANI = ANI / 100) |>
        pivot_longer(
            cols = c(ANI, Detection),
            names_to = "metric",
            values_to = "value"
        ) |>
        filter(habitat != "Other") |>
        mutate(habitat = fct_relevel(habitat, top_habitats$habitat)) |>
        ggplot(aes(x = damage, y = value, fill = domain)) +
        geom_hdr() +
        geom_point(shape = 21, color = "black", alpha = 0.8) +
        facet_grid(habitat ~ metric, scales = "free") +
        theme_bw() +
        theme(
            legend.position = "top",
            text = element_text(size = 10),
            strip.background = element_blank(),
            axis.title.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        ylab("Metric") +
        xlab("Damage") +
        guides(fill = guide_legend(override.aes = list(size = 3)), name = "") +
        scale_fill_manual(values = c("d__Archaea" = "#9D443C", "d__Bacteria" = "#74AFB8")) +
        scale_y_continuous(labels = scales::percent) +
        scale_x_continuous(limits = c(0, 0.65)) +
        ggpubr::rotate()
}

#' Create Woodcroft plot
#' @param data Processed Woodcroft data
#' @param kapk_cdata_agg Aggregated metadata
#' @return ggplot object
plot_woodcroft <- function(data, kapk_cdata_agg) {
    data |>
        inner_join(kapk_cdata_agg) |>
        group_by(label) |>
        mutate(
            total = sum(abundance),
            habitat = case_when(
                habitat == "fen" ~ "Fen",
                habitat == "palsa" ~ "Palsa",
                habitat == "bog" ~ "Bog"
            )
        ) |>
        ungroup() |>
        mutate(
            perc = abundance / total,
            member_unit = fct_relevel(member_unit, c("B3", "B2", "B1")),
            habitat = fct_relevel(habitat, rev(c("Palsa", "Bog", "Fen"))),
            label = fct_reorder(label, -site_rnk)
        ) |>
        ggplot(aes(y = damage, x = perc, fill = domain)) +
        geom_hdr() +
        geom_point(shape = 21, color = "black", alpha = 0.8) +
        facet_grid(habitat ~ .) +
        scale_fill_manual(values = c("d__Archaea" = "#9D443C", "d__Bacteria" = "#74AFB8")) +
        theme_bw() +
        theme(
            legend.position = "top",
            text = element_text(size = 10),
            strip.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        xlab("Proportion") +
        ylab("Damage") +
        scale_x_continuous(labels = scales::percent) +
        guides(fill = guide_legend(override.aes = list(size = 3)), name = "") +
        ggpubr::rotate()
}

#' Main execution function
#' @return list containing results and plots
main <- function() {
    # Load data
    data <- load_data()

    # Process MAG data
    mag_data <- process_mag_data(data)

    # Process habitat metrics data
    habitat_results <- process_habitat_metrics_data(
        mag_data,
        data$kapk_cdata_agg
    )

    # Create habitat metrics plot
    habitat_plot <- plot_habitat_metrics(
        habitat_results$data,
        habitat_results$top_habitats
    )

    # Process Woodcroft data
    woodcroft_data <- process_woodcroft_data(data)

    # Create Woodcroft plot
    woodcroft_plot <- plot_woodcroft(
        woodcroft_data,
        data$kapk_cdata_agg
    )

    # Return results
    list(
        plots = list(
            habitat_metrics = habitat_plot,
            woodcroft = woodcroft_plot
        ),
        data = list(
            habitat_metrics = habitat_results$data,
            habitat_top = habitat_results$top_habitats,
            woodcroft = woodcroft_data,
            mag = mag_data
        )
    )
}

# Run main function and store results
results <- main()
