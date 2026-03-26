library(tidyverse)
library(readxl)
library(janitor)
library(showtext)
library(ggpubr)

# Initialize settings
source("libs/lib.R")
showtext_auto()


#' Load sample metadata from TSV
load_kapk_cdata <- function(path) {
    read_tsv(path, show_col_types = FALSE)
}

#' Load sample metadata from Excel (paper)
load_kapk_cdata_paper <- function(path) {
    readxl::read_xlsx(path) %>%
        clean_names()
}

#' Load extraction metadata from Excel
load_kapk_cdata_extractions <- function(path) {
    readxl::read_xlsx(path) %>%
        clean_names()
}

#' Load initial stats from TSV
load_initial_stats <- function(path) {
    read_tsv(path, show_col_types = FALSE)
}

#' Load derep stats from TSV
load_derep_stats <- function(path) {
    read_tsv(path, show_col_types = FALSE)
}

#' Load new initial stats from multiple files
load_initial_stats_new <- function(dir) {
    files <- list.files(dir, pattern = "*initial*", full.names = TRUE)
    map_dfr(files, function(X) {
        fname <- basename(X)
        label <- str_replace(fname, ".stats-initial.txt", "")
        read_tsv(X, show_col_types = FALSE) %>%
            mutate(label = label)
    })
}

#' Load new derep stats from multiple files
load_derep_stats_new <- function(dir) {
    files <- list.files(dir, pattern = "*derep*", full.names = TRUE)
    map_dfr(files, function(X) {
        fname <- basename(X)
        label <- str_replace(fname, ".stats-derep.txt", "")
        read_tsv(X, show_col_types = FALSE) %>%
            mutate(label = label)
    })
}

#' Load taxonomic annotations
load_tax_info <- function(path) {
    read_tsv(path, col_names = c("reference", "tax_string"), show_col_types = FALSE) %>%
        separate(
            col = tax_string,
            sep = ";",
            into = c("domain", "lineage", "kingdom", "phylum", "class", "order", "family", "genus", "species", "strain")
        )
}

#' Load new taxonomic data from multiple files
load_tax_data_new <- function(dir) {
    files <- list.files(dir, pattern = "*stats-filtered.tsv.gz", full.names = TRUE)
    map_dfr(files, function(X) {
        fname <- basename(X)
        label <- str_replace(fname, ".dedup_stats-filtered.tsv.gz", "")
        read_tsv(X, show_col_types = FALSE) %>%
            mutate(label = label)
    })
}

#' Load new damage data from multiple files
load_dmg_data_new <- function(dir) {
    files <- list.files(dir, pattern = "*.tp-mdmg.weight-1.csv.gz", full.names = TRUE)
    map_dfr(files, function(X) {
        fname <- basename(X)
        label <- str_replace(fname, ".tp-mdmg.weight-1.csv.gz", "")
        read_csv(X, show_col_types = FALSE) %>%
            mutate(label = label)
    })
}

#' Load damage thresholds
load_dmg_thresholds <- function(path) {
    read_tsv(path, show_col_types = FALSE)
}


#' Process extraction metadata
process_extractions <- function(kapk_cdata_extractions, kapk_cdata_paper, kapk_cdata) {
    kapk_cdata_extractions %>%
        inner_join(
            kapk_cdata_paper %>%
                filter(method == "Shotgun") %>%
                mutate(
                    label = basename(collapsed_read_files_path_29bp_and_q_29_only_collapsed_reads),
                    label = gsub("\\..*", "", label)
                ) %>%
                select(label_orig = label, sample_id) %>%
                distinct()
        ) %>%
        inner_join(kapk_cdata) %>% # No explicit 'by' to match original implicit join
        mutate(short_label = paste0(str_split_i(figure_names, "_", i = 1), "_", str_split_i(figure_names, "-", i = -1))) %>%
        select(-label_orig, -label) %>%
        distinct()
}

#' Process new KapK metadata
process_kapk_cdata_new <- function(initial_stats_new) {
    initial_stats_new %>%
        select(label) %>%
        separate(label, into = c("sample"), sep = "_", extra = "drop", remove = FALSE) %>%
        mutate(protocol = gsub("-", "", str_extract(sample, "-(F|M)$"))) %>%
        mutate(sample = gsub("-F$", "", sample)) %>%
        mutate(sample = gsub("-M$", "", sample)) %>%
        mutate(sample = gsub("KapK-", "", sample)) %>%
        mutate(file_name = paste("KapK", sample, sep = "-")) %>%
        filter(protocol == "F")
}

#' Process extraction stats
process_extraction_stats <- function(kapk_cdata_extractions, initial_stats_new, derep_stats_new, kapk_cdata_new) {
    kapk_cdata_extractions %>%
        inner_join(
            initial_stats_new %>%
                select(label, initial_num_reads = num_seqs) %>%
                inner_join(kapk_cdata_new, by = "label") %>%
                select(file_name, initial_num_reads) %>%
                inner_join(
                    derep_stats_new %>%
                        select(label, derep_num_reads = num_seqs) %>%
                        inner_join(kapk_cdata_new, by = "label") %>%
                        select(file_name, derep_num_reads),
                    by = "file_name"
                ),
            by = "file_name"
        ) %>%
        select(short_label, member_unit, site, names(kapk_cdata_extractions), initial_num_reads, derep_num_reads) %>%
        select(-file_name, -figure_names, -site_rnk)
}

#' Process taxonomic data
process_tax_data <- function(tax_data_new, dmg_data_new, tax_info, dmg_thresholds) {
    dmg_data <- dmg_data_new %>%
        select(label, reference = tax_id, damage, significance)

    tax_data_new %>%
        inner_join(dmg_data, by = c("label", "reference")) %>%
        filter(
            breadth >= 0.01,
            significance >= dmg_thresholds$signf,
            n_reads >= 100
        ) %>%
        mutate(abundance = ifelse(tax_abund_tad == 0, tax_abund_read, tax_abund_tad)) %>%
        inner_join(tax_info, by = "reference") %>%
        filter(domain != "d__Eukaryota") %>%
        select(label, abundance, damage, domain, lineage, kingdom, phylum, class, order, family, genus, species) %>%
        group_by(label, domain, lineage, kingdom, phylum, class, order, family, genus, species) %>%
        summarise(abundance = janitor::round_half_up(gm_mean(abundance)), damage = gm_mean(damage), .groups = "drop") %>%
        mutate(
            is_dmg = ifelse(damage >= dmg_thresholds$damage, "Damaged", "Non damaged")
        )
}

#' Process re-extraction taxonomy
process_reextraction_taxonomy <- function(kapk_cdata_new, kapk_cdata_extractions, tax_data_new, tax_info, dmg_data_new) {
    kapk_cdata_new %>%
        inner_join(kapk_cdata_extractions %>% select(short_label, file_name), by = "file_name") %>%
        select(short_label, file_name, label) %>%
        inner_join(tax_data_new, by = "label") %>%
        mutate(abundance = ifelse(tax_abund_tad == 0, tax_abund_read, tax_abund_tad)) %>%
        inner_join(tax_info, by = "reference") %>%
        inner_join(dmg_data_new %>% rename(reference = tax_id) %>% select(label, reference, damage, significance), by = c("label", "reference")) %>%
        select(-label, -file_name)
}


#' Define common plotting theme
get_common_theme <- function(size = 10) {
    theme_bw() +
        theme(
            legend.position = "top",
            legend.title = element_blank(),
            text = element_text(size = size),
            strip.background = element_blank()
        )
}

#' Plot re-extraction proportions
plot_reextraction_proportions <- function(tax_data, kapk_cdata_new, kapk_cdata_extractions) {
    tax_data %>%
        inner_join(kapk_cdata_new, by = "label") %>%
        inner_join(kapk_cdata_extractions %>% select(file_name, short_label, reason), by = "file_name") %>%
        group_by(short_label, domain, is_dmg, reason) %>%
        summarise(abundance = sum(abundance), .groups = "drop") %>%
        group_by(short_label) %>%
        mutate(abundance = abundance / sum(abundance)) %>%
        ungroup() %>%
        mutate(
            domain = fct_relevel(domain, c("d__Viruses", "d__Bacteria", "d__Archaea"))
        ) %>%
        ggplot(aes(x = short_label, y = abundance, fill = domain)) +
        geom_col(position = "stack", color = "black", linewidth = 0.3, width = 1) +
        labs(x = "Sample", y = "Proportion") +
        facet_grid(reason ~ is_dmg, scales = "free_y", space = "free") +
        scale_y_continuous(labels = scales::percent) +
        scale_fill_manual(values = c("d__Archaea" = "#9D443C", "d__Bacteria" = "#74AFB8", "d__Viruses" = "#EECC66", "d__Eukaryota" = "#009E73")) +
        ggpubr::rotate() +
        get_common_theme()
}


#' Main execution function
main <- function() {
    # Load metadata
    kapk_cdata <- load_kapk_cdata("./data/cdata/KapK-cdata-manuscript-20221211.tsv")
    kapk_cdata_paper <- load_kapk_cdata_paper("./data/cdata/41586_2022_5453_MOESM3_ESM.xlsx")
    kapk_cdata_extractions <- load_kapk_cdata_extractions("./data/cdata/kapk-extraction-eDNA_dataforpaper.xlsx")

    # Process extractions
    kapk_cdata_extractions <- process_extractions(kapk_cdata_extractions, kapk_cdata_paper, kapk_cdata)

    # Load stats
    initial_stats <- load_initial_stats("./data/stats/all.stats-initial-summary.tsv.gz")
    derep_stats <- load_derep_stats("./data/stats/all.stats-derep-summary.tsv.gz")
    initial_stats_new <- load_initial_stats_new("./data/re-extractions/stats/")
    derep_stats_new <- load_derep_stats_new("./data/re-extractions/stats/")

    # Process new KapK metadata
    kapk_cdata_new <- process_kapk_cdata_new(initial_stats_new)

    # Process extraction stats
    extraction_stats <- process_extraction_stats(kapk_cdata_extractions, initial_stats_new, derep_stats_new, kapk_cdata_new)

    # Load taxonomic data
    tax_info <- load_tax_info("./data/taxonomy/hires-organelles-viruses-arctic.tax.tsv")
    tax_data_new <- load_tax_data_new("./data/re-extractions/taxonomy")
    dmg_data_new <- load_dmg_data_new("./data/re-extractions/taxonomy")
    dmg_thresholds <- load_dmg_thresholds("./results/damage/dmg_thresholds.tsv")

    # Process taxonomic data
    tax_data <- process_tax_data(tax_data_new, dmg_data_new, tax_info, dmg_thresholds)

    # Process re-extraction taxonomy
    reextraction_taxonomy <- process_reextraction_taxonomy(kapk_cdata_new, kapk_cdata_extractions, tax_data_new, tax_info, dmg_data_new)

    # Generate plot
    reextraction_plot <- plot_reextraction_proportions(tax_data, kapk_cdata_new, kapk_cdata_extractions)

    # Prepare results
    tables <- list(
        kapk_cdata = kapk_cdata,
        kapk_cdata_paper = kapk_cdata_paper,
        kapk_cdata_extractions = kapk_cdata_extractions,
        initial_stats = initial_stats,
        derep_stats = derep_stats,
        initial_stats_new = initial_stats_new,
        derep_stats_new = derep_stats_new,
        kapk_cdata_new = kapk_cdata_new,
        extraction_stats = extraction_stats,
        tax_info = tax_info,
        tax_data_new = tax_data_new,
        dmg_data_new = dmg_data_new,
        dmg_thresholds = dmg_thresholds,
        tax_data = tax_data,
        reextraction_taxonomy = reextraction_taxonomy
    )

    figures <- list(
        kapk_reextraction = reextraction_plot
    )

    # Return results
    list(
        tables = tables,
        figures = figures
    )
}

# Run main function and store results
results <- main()