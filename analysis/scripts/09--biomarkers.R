# Required libraries
library(tidyverse)
library(readxl)
library(showtext)
library(microViz)
library(ggtext)

# Initialize settings
showtext_auto()

# Data Loading Functions ----

#' Load sample metadata
#' @param path Path to metadata file
#' @return Tibble with cleaned metadata
load_metadata <- function(path) {
    read_tsv(path, show_col_types = FALSE)
}

#' Load biomarker data
#' @param path Path to biomarker Excel file
#' @return Tibble with raw biomarker data
load_biomarker_data <- function(path) {
    read_xlsx(path) %>%
        janitor::clean_names() %>%
        rename(short_label = sample_name) %>%
        select(-unit, -sample_code, -x4)
}

# Data Processing Functions ----

#' Process biomarker data
#' @param raw_data Raw biomarker data
#' @param combined_biomarker Biomarker metadata with short names and environments
#' @return Processed biomarker data tibble
process_biomarker_data <- function(raw_data, combined_biomarker) {
    raw_data %>%
        pivot_longer(
            cols = -short_label,
            names_to = "short_name",
            values_to = "abundance"
        ) %>%
        mutate(
            short_name = ifelse(short_name %in% c("unsaturated_bht_25", "unsaturated_bht_26"),
                "unsaturated_bht", short_name
            ),
            short_name = ifelse(short_name %in% c("me_adenosylhopane_hg_di_me_8", "me_adenosylhopane_hg_di_me_9"),
                "me_adenosylhopane_hg_di_me", short_name
            )
        ) %>%
        group_by(short_label, short_name) %>%
        summarise(abundance = sum(abundance), .groups = "drop") %>%
        group_by(short_label) %>%
        mutate(total = sum(abundance)) %>%
        ungroup() %>%
        mutate(abundance = abundance / total) %>%
        inner_join(combined_biomarker, by = "short_name")
}

#' Create biomarker metadata
#' @return Tibble with biomarker names, short names, environments, and colors
create_biomarker_metadata <- function() {
    biomarkers <- c(
        "adenosylhopane", "adenosylhopane<sub><sub>HG-Me<sub><sub>", "adenosylhopane<sub><sub>HG-diMe<sub><sub>",
        "Me-adenosylhopane<sub><sub>HG-diMe<sub><sub>", "diMe-adenosylhopane<sub><sub>HG-diMe<sub><sub>",
        "Me-adenosylhopane", "Me-adenosylhopane<sub><sub>HG-Me<sub><sub>", "adenosylhopane type-2",
        "adenosylhopane type-2 deMe", "aminopentol", "unsaturated aminopentol", "aminotetrol",
        "aminotriol", "BHhexol", "BHpentol", "unsaturated BHpentol", "BHT",
        "BHT-CE", "anhydro BHT", "unsaturated BHT", "Me-BHT"
    )
    short_names <- c(
        "adenosylhopane", "adenosylhopane_hg_me", "adenosylhopane_hg_di_me",
        "me_adenosylhopane_hg_di_me", "di_me_adenosylhopane_hg_di_me",
        "me_adenosylhopane", "me_adenosylhopane_hg_me", "adenosylhopane_type_2",
        "adenosylhopane_type_2_de_me", "aminopentol", "unsaturated_aminopentol",
        "aminotetrol", "aminotriol", "b_hhexol", "b_hpentol", "unsaturated_bhpentol",
        "bht", "bht_ce", "anhydro_bht", "unsaturated_bht", "me_bht"
    )
    envs <- c(rep("Soil", 9), rep("Methanotrophs", 3), rep("Various sources", 9))

    # Define colors
    soil_colors <- c(
        "#882d17", "#654522", "#848482", "#8db600", "#c2b280",
        "#2b3d26", "#875692", "#b3446c", "#008856"
    )
    methanotrophs_colors <- c("#f3c300", "#f6a600", "#e25822")
    various_sources_colors <- c(
        "#e68fac", "#a1caf1", "#0067a5", "lightgrey",
        "#f38400", "#dcd300", "#be0032", "#604e97", "#f99379"
    )
    colors <- c(soil_colors, methanotrophs_colors, various_sources_colors)

    tibble(biomarker = biomarkers, short_name = short_names, env = envs, color = colors)
}

#' Get short label order
#' @param kapk_cdata_agg Aggregated sample metadata
#' @param biomarker_data Processed biomarker data
#' @return Character vector of ordered short labels
get_short_label_order <- function(kapk_cdata_agg, biomarker_data) {
    kapk_cdata_agg %>%
        filter(short_label %in% unique(biomarker_data$short_label)) %>%
        mutate(label = fct_reorder(label, -site_rnk)) %>%
        pull(label) %>%
        levels() %>%
        enframe(value = "label") %>%
        inner_join(kapk_cdata_agg %>% select(label, short_label), by = "label") %>%
        pull(short_label)
}

# Visualization Functions ----

#' Plot biomarker proportions
#' @param biomarker_data Processed biomarker data
#' @param kapk_cdata_agg Aggregated sample metadata
#' @param short_label_order Ordered short labels
#' @param colors Named vector of colors
#' @return ggplot object
plot_biomarker_proportions <- function(biomarker_data, kapk_cdata_agg, short_label_order, colors) {
    biomarker_data %>%
        inner_join(kapk_cdata_agg, by = "short_label") %>%
        mutate(
            short_label = fct_relevel(short_label, short_label_order),
            member_unit = fct_relevel(member_unit, c("B3", "B2", "B1")),
            biomarker = fct_relevel(biomarker, rev(names(colors))),
            env = fct_relevel(env, c("Soil", "Methanotrophs", "Various sources"))
        ) %>%
        ggplot(aes(x = short_label, y = abundance, fill = biomarker)) +
        geom_col(color = "black", linewidth = 0.1, width = 1, alpha = 0.8) +
        theme_bw() +
        labs(x = "Sample", y = "Proportion") +
        theme(
            # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5),
            legend.position = "top",
            legend.text = element_markdown(),
            legend.title = element_blank(),
            text = element_text(size = 10),
            strip.background = element_blank(),
            strip.text.y = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        scale_fill_manual(values = colors) +
        facet_grid(member_unit ~ env, space = "free_y", scale = "free") +
        scale_y_continuous(labels = scales::percent) +
        guides(fill = guide_legend(keywidth = 1, keyheight = 1, reverse = TRUE)) +
        ggpubr::rotate()
}

# Main Execution Function ----

#' Main execution function
#' @return List containing results and plots
main <- function() {
    # Load metadata
    kapk_cdata <- load_metadata("./data/cdata/KapK-cdata-manuscript-20221211.tsv")
    kapk_cdata_agg <- load_metadata("./data/cdata/KapK-cdata-agg-manuscript-20221211.tsv")

    # Load and process biomarker data
    raw_biomarker_data <- load_biomarker_data("./data/biomarkers/kapk-20230925-biomarkers.xlsx")
    combined_biomarker <- create_biomarker_metadata()
    biomarker_data <- process_biomarker_data(raw_biomarker_data, combined_biomarker)

    # Get short label order
    short_label_order <- get_short_label_order(kapk_cdata_agg, biomarker_data)

    # Define colors using the full set of biomarkers
    colors <- setNames(combined_biomarker$color, combined_biomarker$biomarker)

    # Generate figure
    figure <- plot_biomarker_proportions(biomarker_data, kapk_cdata_agg, short_label_order, colors)

    # Prepare results
    tables <- list(
        kapk_cdata = kapk_cdata,
        kapk_cdata_agg = kapk_cdata_agg,
        raw_biomarker_data = raw_biomarker_data,
        combined_biomarker = combined_biomarker,
        biomarker_data = biomarker_data
    )

    figures <- list(
        biomarker_proportions = figure
    )

    # Return results
    list(
        tables = tables,
        figures = figures
    )
}

# Run main function and store results
results <- main()
