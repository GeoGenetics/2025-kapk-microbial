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
library(unikn)
library(readxl)

# Initialize settings
source("./libs/lib.R")
source("./libs/lib-vir.R")
showtext_auto()

# Data Loading Functions ----

#' Load sample metadata
#' @param path Path to metadata file
#' @return Tibble with cleaned metadata
load_metadata <- function(path) {
    read_tsv(path, show_col_types = FALSE)
}

#' Load taxonomic information
#' @param path Path to taxonomy file
#' @return Tibble with parsed taxonomic data
load_taxonomy <- function(path) {
    read_tsv(path, col_names = c("reference", "tax_string"), show_col_types = FALSE) |>
        separate(tax_string,
            sep = ";",
            into = c(
                "domain", "lineage", "kingdom", "phylum", "class",
                "order", "family", "genus", "species", "strain"
            )
        )
}

#' Load abundance data
#' @param path Path to abundance file
#' @param col_names Optional column names
#' @return Tibble with abundance data
load_abundance <- function(path, col_names = NULL) {
    read_tsv(path, show_col_types = FALSE)
}

#' Load damage data
#' @param path Path to damage file
#' @return Tibble with damage data
load_damage <- function(path) {
    read_csv(path, show_col_types = FALSE)
}

#' Load IMGVR metadata
#' @param path Path to IMGVR file
#' @return Tibble with cleaned IMGVR metadata
load_imgvr <- function(path) {
    imgvr_cdata <- read_tsv(path) |>
        rename(reference = UVIG) |>
        janitor::clean_names() |>
        mutate(taxon_oid = str_split(reference, "_") |> map_chr(3))
    jgi_url <- "https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid="

    # Let's add a column with the JGI url
    imgvr_cdata_1 <- imgvr_cdata |>
        filter(!grepl("GVMAG", taxon_oid))
    imgvr_cdata_2 <- imgvr_cdata |>
        filter(grepl("GVMAG", taxon_oid)) |>
        mutate(taxon_oid = str_split(taxon_oid, "-") |> map_chr(3))
    bind_rows(imgvr_cdata_1, imgvr_cdata_2) |>
        mutate(jgi_url = paste0(jgi_url, taxon_oid))
}

#' Load amino acid data
#' @param agg_path Path to aggregated AA results
#' @param hits_path Path to AA hits
#' @param ngenes_path Path to gene counts
#' @return List of tibbles with AA data
load_aa_data <- function(agg_path, hits_path, ngenes_path) {
    viral_ngenes <- read_tsv(ngenes_path,
        col_names = c("reference", "n_genes_total"),
        show_col_types = FALSE
    )
    aa_results <- read_tsv(agg_path, show_col_types = FALSE) |>
        rename(reference = group) |>
        inner_join(viral_ngenes) |>
        mutate(prot_coverage = n_genes / n_genes_total)
    aa_hits <- read_tsv(hits_path, show_col_types = FALSE) |>
        rename(gene = reference, reference = group)
    list(results = aa_results, hits = aa_hits, ngenes = viral_ngenes)
}

#' Load protein annotations
#' @param pfam_path Path to PFAM results
#' @param phrog_path Path to PHROG results
#' @param phrog_cdata_path Path to PHROG metadata
#' @return List of filtered annotation tibbles
load_protein_annotations <- function(pfam_path, phrog_path, phrog_cdata_path) {
    pfam_cols <- c(
        "gene", "pfam_acc", "probability", "e-value", "Score", "Cols",
        "q_start", "q_stop", "t_start", "t_stop", "q_len", "t_len",
        "q_cov", "t_cov", "pfam_name", "pfam_description"
    )
    phrog_cols <- c(
        "gene", "phrog_acc", "probability", "e-value", "Score", "Cols",
        "q_start", "q_stop", "t_start", "t_stop", "q_len", "t_len",
        "q_cov", "t_cov", "name", "description"
    )

    pfam <- read_tsv(pfam_path, col_names = pfam_cols, show_col_types = FALSE) |>
        filter(probability >= 90, q_cov >= 0.4 | t_cov >= 0.4) |>
        group_by(gene) |>
        arrange(desc(probability), desc(q_cov), desc(t_cov)) |>
        slice(1) |>
        ungroup()

    phrog <- read_tsv(phrog_path, col_names = phrog_cols, show_col_types = FALSE) |>
        filter(probability >= 90, q_cov >= 0.4 | t_cov >= 0.4) |>
        group_by(gene) |>
        arrange(desc(probability), desc(q_cov), desc(t_cov)) |>
        slice(1) |>
        ungroup() |>
        inner_join(read_csv(phrog_cdata_path, show_col_types = FALSE) |>
            clean_names() |>
            rename(phrog_acc = number_phrog))

    list(pfam = pfam, phrog = phrog)
}

# Data Processing Functions ----

#' Process dereplication stats
#' @param derep_path Path to derep stats
#' @param kapk_cdata Sample metadata
#' @param reads_mapped_path Path to reads mapped data
#' @param label_to_file_path Path to label-to-file mapping
#' @param tax_data_agg Aggregated taxonomic data
#' @return Tibble with processed derep stats
process_derep_stats <- function(derep_path, kapk_cdata, reads_mapped_path, label_to_file_path, tax_data_agg) {
    label_to_file <- read_tsv(label_to_file_path, show_col_types = FALSE)
    reads_mapped <- read_tsv(reads_mapped_path,
        col_names = c("file", "reads_mapped"),
        show_col_types = FALSE
    ) |>
        inner_join(label_to_file) |>
        inner_join(kapk_cdata |> select(label, figure_names)) |>
        select(-label, -file) |>
        rename(label = figure_names)

    read_tsv(derep_path, show_col_types = FALSE) |>
        inner_join(kapk_cdata) |>
        select(label = figure_names, member_unit, num_seqs, site_rnk, avg_len) |>
        group_by(label, member_unit, site_rnk) |>
        summarise(
            num_seqs = janitor::round_half_up(mean(num_seqs)),
            avg_len = mean(avg_len),
            .groups = "drop"
        ) |>
        mutate(step = "derep") |>
        filter(label %in% unique(tax_data_agg$label)) |>
        inner_join(reads_mapped) |>
        mutate(
            prop_mapped = reads_mapped / num_seqs,
            not_mapped = num_seqs - reads_mapped
        )
}

#' Process viral data
#' @param tax_data_agg Aggregated taxonomic data
#' @param tax_info Taxonomic information
#' @param dmg_local_agg Aggregated damage data
#' @param kapk_cdata_agg Aggregated sample metadata
#' @param imgvr_cdata IMGVR metadata
#' @return List of processed viral data
process_viral_data <- function(tax_data_agg, tax_info, dmg_local_agg, kapk_cdata_agg, imgvr_cdata) {
    tax_data_agg <- tax_data_agg |> inner_join(dmg_local_agg)
    tax_data_agg_viruses <- tax_data_agg |>
        inner_join(tax_info) |>
        filter(domain == "d__Viruses")

    non_damaged <- tax_data_agg_viruses |>
        filter(is_dmg == "Non-damaged", breadth >= 0.1) |>
        arrange(desc(coverage_mean)) |>
        inner_join(imgvr_cdata |> select(reference, ecosystem_classification, host_taxonomy_prediction, jgi_url)) |>
        inner_join(kapk_cdata_agg)

    damaged <- tax_data_agg_viruses |>
        filter(is_dmg == "Damaged", breadth >= 0.1) |>
        arrange(desc(coverage_mean)) |>
        inner_join(imgvr_cdata |> select(reference, ecosystem_classification, host_taxonomy_prediction, jgi_url)) |>
        inner_join(kapk_cdata_agg)

    list(non_damaged = non_damaged, damaged = damaged, all_viruses = tax_data_agg_viruses)
}

#' Process amino acid analysis
#' @param aa_data List of amino acid data
#' @param protein_annotations List of protein annotations
#' @param kapk_cdata Sample metadata
#' @param tax_info Taxonomic information
#' @param imgvr_cdata IMGVR metadata
#' @return Processed amino acid results
process_aa_analysis <- function(aa_data, protein_annotations, kapk_cdata, tax_info, imgvr_cdata) {
    pfam_agg <- aa_data$hits |>
        distinct(reference, gene, label) |>
        inner_join(protein_annotations$pfam) |>
        group_by(reference, label) |>
        summarise(
            pfam_n_annotations = sum(!is.na(pfam_acc)),
            pfam_acc = paste0(pfam_acc, collapse = "|"),
            pfam_name = paste0(pfam_name, collapse = "|"),
            pfam_description = paste0(pfam_description, collapse = "|"),
            .groups = "drop"
        ) |>
        inner_join(kapk_cdata_filt) |>
        select(-label, -label_orig) |>
        rename(label = figure_names)

    phrog_agg <- aa_data$hits |>
        distinct(reference, gene, label) |>
        inner_join(protein_annotations$phrog) |>
        group_by(reference, label) |>
        summarise(
            phrog_n_annotations = sum(!is.na(phrog_acc)),
            phrog_acc = paste0(phrog_acc, collapse = "|"),
            phrog_category = paste(category, collapse = "|"),
            phrog_annotation = paste(annotation, collapse = "|"),
            phrog_host_domain = paste(host_domain, collapse = "|"),
            .groups = "drop"
        ) |>
        inner_join(kapk_cdata_filt) |>
        select(-label, -label_orig) |>
        rename(label = figure_names)

    aa_results_imgvr <- aa_data$results |>
        filter(prot_coverage >= 0.1, grepl("IMGVR", reference)) |>
        inner_join(imgvr_cdata |> select(reference, ecosystem_classification, host_taxonomy_prediction, jgi_url)) |>
        inner_join(tax_info) |>
        inner_join(kapk_cdata_filt) |>
        select(-label, -label_orig) |>
        rename(label = figure_names)

    aa_results_non_imgvr <- aa_data$results |>
        filter(prot_coverage >= 0.1, !grepl("IMGVR", reference)) |>
        inner_join(kapk_cdata_filt) |>
        select(-label, -label_orig) |>
        rename(label = figure_names)

    bind_rows(aa_results_imgvr, aa_results_non_imgvr) |>
        mutate(rank = case_when(
            prot_coverage >= 0.25 & n_genes > 1 ~ "green",
            prot_coverage >= 0.2 & prot_coverage < 0.25 & n_genes > 1 ~ "yellow",
            prot_coverage < 0.2 ~ "red",
            n_genes == 1 & phylum == "p__Cressdnaviricota" & prot_coverage >= 0.2 ~ "blue",
            n_genes == 1 & prot_coverage >= 0.2 ~ "grey",
            TRUE ~ "NA"
        )) |>
        left_join(pfam_agg) |>
        left_join(phrog_agg) |>
        select(reference, label, rank, prot_coverage,
            n_proteins_detected = n_genes,
            n_proteins_genome = n_genes_total, pfam_n_annotations, phrog_n_annotations,
            coverage_mean, coverage_stdev, coverage_median, coverage_sum,
            avg_read_length, stdev_read_length, avg_identity, stdev_identity,
            pfam_acc, pfam_name, pfam_description, phrog_acc, phrog_category,
            phrog_annotation, phrog_host_domain, ecosystem_classification,
            host_taxonomy_prediction, jgi_url, domain, lineage, kingdom, phylum,
            class, order, family, genus, species, strain, member_unit, site, site_rnk
        )
}

# Visualization Functions ----

#' Define common plotting theme
#' @return ggplot theme object
get_common_theme <- function(size = 12) {
    theme_bw() +
        theme(
            legend.position = "none",
            legend.title = element_blank(),
            text = element_text(size = size),
            strip.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        )
}

#' Plot nucleotide vs amino acid damage and breadth
#' @param nt_aa_data Combined NT/AA data
#' @param dmg_thresholds Damage thresholds
#' @return ggplot object
plot_nt_aa_damage_breadth <- function(nt_aa_data, dmg_thresholds) {
    colors <- c("#AFCCB8", "#C94A6B", "#4097AA")
    names(colors) <- c("B3", "B2", "B1")

    ggplot(nt_aa_data, aes(x = damage, y = breadth, fill = member_unit)) +
        geom_vline(xintercept = dmg_thresholds$damage, linetype = "dashed", color = "#B04035", alpha = 0.5) +
        geom_point(alpha = 0.8, size = 1.5, shape = 21) +
        facet_grid(~is_aa) +
        scale_fill_manual(values = colors) +
        labs(x = "Damage", y = "Detection") +
        scale_y_continuous(labels = scales::percent) +
        get_common_theme()
}

#' Plot ecosystem ANIr
#' @param nt_aa_data Combined NT/AA data
#' @param ecosystems List of ecosystems to include
#' @return ggplot object
plot_ecosystem_anir <- function(nt_aa_data, ecosystems) {
    com_colors <- c("#F2B705", "#1C7085")
    names(com_colors) <- c("Pioneer", "Permafrost")

    ecosystem_nref <- nt_aa_data |>
        filter(breadth >= 0.1) |>
        distinct(reference, ecosystem_classification) |>
        count(ecosystem_classification, name = "n_ref") |>
        arrange(desc(n_ref))

    nt_aa_data |>
        filter(breadth >= 0.1, ecosystem_classification %in% ecosystems) |>
        group_by(ecosystem_classification) |>
        mutate(n = n()) |>
        ungroup() |>
        inner_join(ecosystem_nref) |>
        mutate(read_ani_mean = read_ani_mean / 100) |>
        mutate(ecosystem_classification = fct_relevel(ecosystem_classification, rev(ecosystems))) |>
        mutate(is_dmg = if_else(is_dmg == "Damaged", "Pioneer", "Permafrost")) |>
        ggplot(aes(ecosystem_classification, read_ani_mean, fill = is_dmg)) +
        geom_point(alpha = 0.2, position = position_jitter(width = 0.1, seed = 3922), shape = 21, size = 1) +
        geom_pointrange(
            stat = "summary", fun = median,
            fun.min = ~ quantile(.x, 0.25), fun.max = ~ quantile(.x, 0.75),
            shape = 21, size = 0.3, position = position_dodge(width = 0.6), linewidth = 0.5
        ) +
        scale_fill_manual(values = com_colors) +
        ggpubr::rotate() +
        labs(x = "", y = "ANIr") +
        scale_y_continuous(labels = scales::percent) +
        get_common_theme(size = 11)
}

#' Plot methanogenic host counts
#' @param aa_results_nored_methano Methanogenic amino acid results
#' @param methano_cdata Methanogenic metadata
#' @return ggplot object
plot_methano_hosts <- function(aa_results_nored_methano, methano_cdata) {
    aa_results_nored_methano |>
        inner_join(methano_cdata) |>
        group_by(host) |>
        summarise(n = n(), .groups = "drop") |>
        mutate(host = fct_reorder(host, n, .desc = TRUE)) |>
        ggplot(aes(host, n)) +
        geom_col() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        labs(x = "", y = "References") +
        get_common_theme()
}

#' Plot nucleotide vs amino acid phyla distribution
#' @param nt_aa_data Combined NT/AA data
#' @param aa_results_nored_filt Filtered amino acid results
#' @param vir_phyla_order Ordered viral phyla
#' @return ggplot object
plot_nt_aa_phyla <- function(nt_aa_data, aa_results_nored_filt, vir_phyla_order) {
    vir_colors <- c("#7faec2", "#f0cfb9", "#f7e7cf", "#f8b3b5", "#87809c", "#94b18a")
    names(vir_colors) <- vir_phyla_order

    aa_results_nored_filt_nrefs <- aa_results_nored_filt |>
        distinct(reference, phylum) |>
        count(phylum, name = "n_ref", sort = TRUE)

    nt_aa_data |>
        distinct(reference, phylum) |>
        count(phylum, name = "n_ref") |>
        mutate(class = "NT") |>
        bind_rows(aa_results_nored_filt_nrefs |> mutate(class = "AA")) |>
        group_by(class) |>
        mutate(perc = n_ref / sum(n_ref)) |>
        ungroup() |>
        mutate(phylum = fct_relevel(phylum, rev(vir_phyla_order))) |>
        ggplot(aes(class, perc, fill = phylum)) +
        geom_col(width = 0.5, color = "black", linewidth = 0.2, alpha = 0.8) +
        scale_fill_manual(values = vir_colors) +
        ggpubr::rotate() +
        scale_y_continuous(labels = scales::percent) +
        labs(x = "", y = "") +
        get_common_theme(size = 11) +
        theme(legend.position = "top")
}

#' Plot amino acid hits (depth, identity, and functional categories)
#' @param aa_hits_filt Filtered amino acid hits
#' @param aa_results_nored_filt_annotated Annotated amino acid results
#' @param vir_phyla_order Ordered viral phyla
#' @return ggplot object (combined)
plot_aa_hits <- function(aa_hits_filt, aa_results_nored_filt_annotated, vir_phyla_order) {
    colors <- c("#AFCCB8", "#C94A6B", "#4097AA")
    names(colors) <- c("B3", "B2", "B1")
    func_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#999999", "#CC79A7")
    names(func_colors) <- c(
        "Other", "Virus structural protein", "DNA/RNA transactions/processing",
        "Lysis", "TA/RM systems", "AMG and other", "Unknown"
    )


    aa_hits_filt_nprots <- aa_hits_filt |>
        inner_join(aa_results_nored_filt_annotated) |>
        select(label, reference, gene, member_unit) |>
        group_by(label, reference, member_unit) |>
        count(name = "n_prots") |>
        ungroup() |>
        inner_join(viral_ngenes) |>
        mutate(prot_cov = n_prots / n_genes_total)

    # Calculate average identity
    aa_hits_filt_avg_identity <- aa_hits_filt |>
        inner_join(aa_results_nored_filt_annotated) |>
        select(label, reference, avg_identity, member_unit) |>
        group_by(label, reference, member_unit) |>
        summarise(avg_identity = mean(avg_identity)) |>
        ungroup()


    aa_hits_filt_avg_depth <- aa_hits_filt |>
        inner_join(aa_results_nored_filt_annotated) |>
        inner_join(aa_hits_filt |> select(label, reference, phylum) |> distinct()) |>
        select(label, reference, depth_mean, member_unit) |>
        group_by(label, reference, member_unit) |>
        summarise(depth = janitor::round_half_up(gm_mean(depth_mean))) |>
        ungroup()


    # Depth plot
    p00 <- aa_hits_filt_avg_depth |>
        group_by(label, member_unit, reference) |>
        summarise(depth = sum(depth)) |>
        ungroup() |>
        group_by(label, member_unit) |>
        mutate(total_depth = sum(depth)) |>
        ungroup() |>
        mutate(prop_depth = depth / total_depth) |>
        inner_join(aa_hits_filt |> select(label, reference, phylum) |> distinct()) |>
        mutate(phylum = gsub("p__", "", phylum)) |>
        mutate(phylum = fct_relevel(phylum, rev(gsub("p__", "", vir_phyla_order)))) |>
        mutate(org_cat = ifelse(phylum == "Methanogenic viruses", "B", "A")) |>
        ggplot(aes(x = phylum, y = prop_depth, fill = member_unit)) +
        geom_point(alpha = 0.2, position = position_jitter(width = 0.1, seed = 3922), shape = 21, size = 1) +
        geom_pointrange(
            stat = "summary", fun = median,
            fun.min = ~ quantile(.x, 0.25), fun.max = ~ quantile(.x, 0.75),
            shape = 21, size = 0.3, position = position_dodge(width = 0.6), linewidth = 0.5
        ) +
        facet_grid(org_cat ~ ., scales = "free_y", space = "free_y") +
        scale_fill_manual(values = colors) +
        ggpubr::rotate() +
        labs(x = "", y = "Relative abundance") +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1), trans = "sqrt") +
        get_common_theme() +
        theme(strip.text.y = element_blank(), axis.text.y = element_blank())

    # Identity plot
    p01 <- aa_hits_filt_avg_identity |>
        inner_join(aa_hits_filt |> select(label, reference, phylum) |> distinct()) |>
        filter(phylum %in% vir_phyla) |>
        mutate(phylum = gsub("p__", "", phylum)) |>
        mutate(phylum = fct_relevel(phylum, rev(gsub("p__", "", vir_phyla_order)))) |>
        mutate(org_cat = ifelse(phylum == "Methanogenic viruses", "B", "A")) |>
        ggplot(aes(x = phylum, y = avg_identity, fill = member_unit)) +
        geom_point(alpha = 0.2, position = position_jitter(width = 0.1, seed = 3922), shape = 21, size = 1) +
        geom_pointrange(
            stat = "summary", fun = median,
            fun.min = ~ quantile(.x, 0.25), fun.max = ~ quantile(.x, 0.75),
            shape = 21, size = 0.3, position = position_dodge(width = 0.6), linewidth = 0.5
        ) +
        facet_grid(org_cat ~ ., scales = "free_y", space = "free_y") +
        scale_fill_manual(values = colors) +
        ggpubr::rotate() +
        labs(x = "", y = "Similarity (AA)") +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        get_common_theme() +
        theme(strip.text.y = element_blank(), axis.text.y = element_blank())

    # Functional categories plot
    p02 <- aa_hits_filt |>
        inner_join(aa_results_nored_filt_annotated) |>
        group_by(phylum, functional_category) |>
        summarise(n = n(), .groups = "drop") |>
        group_by(phylum) |>
        mutate(prop = n / sum(n)) |>
        ungroup() |>
        mutate(phylum = gsub("p__", "", phylum), phylum = fct_relevel(phylum, rev(gsub("p__", "", vir_phyla_order)))) |>
        mutate(functional_category = fct_reorder(functional_category, prop, mean)) |>
        mutate(org_cat = if_else(phylum == "Methanogenic viruses", "B", "A")) |>
        ggplot(aes(phylum, prop, fill = functional_category)) +
        geom_col(width = 1, color = "black", linewidth = 0.2, alpha = 0.8) +
        facet_grid(org_cat ~ ., scales = "free_y", space = "free_y") +
        scale_fill_manual(values = func_colors) +
        ggpubr::rotate() +
        labs(x = "", y = "Proportion of proteins") +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        get_common_theme() +
        theme(legend.position = "top", strip.text.y = element_blank(), axis.text.y = element_blank())

    ggarrange(p00, p01, p02, align = "hv", ncol = 3, nrow = 1)
}

#' Plot amino acid ecosystem similarity
#' @param aa_results_nored_filt Filtered amino acid results
#' @param ecosystems List of ecosystems to include
#' @param vir_phyla_order Ordered viral phyla
#' @return ggplot object
plot_aa_ecosystem_similarity <- function(aa_results_nored_filt, ecosystems, vir_phyla_order) {
    vir_colors <- c("#7faec2", "#f0cfb9", "#f7e7cf", "#f8b3b5", "#87809c", "#94b18a")
    names(vir_colors) <- vir_phyla_order

    aa_results_nored_filt |>
        filter(ecosystem_classification %in% ecosystems) |>
        mutate(
            ecosystem_classification = fct_relevel(ecosystem_classification, rev(ecosystems)),
            phylum = fct_relevel(phylum, vir_phyla_order)
        ) |>
        ggplot(aes(x = ecosystem_classification, y = avg_identity, fill = phylum)) +
        geom_point(alpha = 0.2, position = position_jitter(width = 0.1, seed = 3922), shape = 21, size = 1) +
        geom_pointrange(
            stat = "summary", fun = median,
            fun.min = ~ quantile(.x, 0.25), fun.max = ~ quantile(.x, 0.75),
            shape = 21, size = 0.3, position = position_dodge(width = 0.6), linewidth = 0.5
        ) +
        scale_fill_manual(values = vir_colors) +
        ggpubr::rotate() +
        labs(x = "", y = "Similarity (AA)") +
        scale_y_continuous(labels = scales::percent) +
        get_common_theme(size = 11)
}

# Main Execution Function ----

#' Main execution function
#' @return List containing results and plots
main <- function() {
    # Load metadata
    kapk_cdata <- load_metadata("./data/cdata/KapK-cdata-manuscript-20221211.tsv")

    kapk_cdata_agg <- load_metadata("./data/cdata/KapK-cdata-agg-manuscript-20221211.tsv")
    tax_info <- load_taxonomy("./data/taxonomy/hires-organelles-viruses-arctic.tax.tsv")

    # Load abundance and damage data
    tax_data_agg <- load_abundance("./results/taxonomy/tp-mapping-filtered.nocontam.10M.agg.nobloom.tax.tsv.gz")
    tax_data <- load_abundance("./results/taxonomy/tp-mapping-filtered.nocontam.10M.nobloom.tax.tsv.gz")
    dmg_local_agg <- load_damage("./results/damage/tp-mdmg.weight-1.local.10M.agg.nobloom.csv.gz")
    dmg_thresholds <- load_abundance("./results/damage/dmg_thresholds.tsv")

    kapk_cdata_filt <- kapk_cdata |>
        filter(label %in% (tax_data$label |> unique()))

    # Load derep stats
    derep_stats <- process_derep_stats(
        derep_path = "./data/stats/all.stats-derep-summary.tsv.gz",
        kapk_cdata = kapk_cdata,
        reads_mapped_path = "./data/cdata/all.tp-classified-reads.tsv",
        label_to_file_path = "./data/cdata/KapK-label-to-file-20221211.tsv",
        tax_data_agg = tax_data_agg
    )

    # Load IMGVR and AA data
    imgvr_cdata <- load_imgvr("./data/cdata/IMGVR_all_Sequence_information-high_confidence.tsv")
    aa_data <- load_aa_data(
        agg_path = "./data/taxonomy/viruses-aa/IMGVR4_derepG-archaea-profiling.group-abundances-agg.tsv.gz",
        hits_path = "./data/taxonomy/viruses-aa/IMGVR4_derepG-archaea-profiling.group-abundances.tsv.gz",
        ngenes_path = "./data/taxonomy/IMGVR_all_proteins-high_confidence_derepG-archaeal.ngenes.tsv"
    )
    protein_annotations <- load_protein_annotations(
        pfam_path = "./data/taxonomy/viruses-aa/pfamA_res.tsv.gz",
        phrog_path = "./data/taxonomy/viruses-aa/phrog_res.tsv.gz",
        phrog_cdata_path = "./data/taxonomy/viruses-aa/PHROG_index.csv"
    )

    # Process data
    viral_data <- process_viral_data(tax_data_agg, tax_info, dmg_local_agg, kapk_cdata_agg, imgvr_cdata)
    aa_results <- process_aa_analysis(aa_data, protein_annotations, kapk_cdata_filt, tax_info, imgvr_cdata)

    # Process methanogenic data
    methano_cdata <- read_tsv("./data/taxonomy/viruses-aa/methano_cdata.tsv",
        col_names = c("tax_class", "reference", "host")
    )

    aa_results_nored <- aa_results |> filter(rank != "red")
    aa_results_nored_methano <- aa_results_nored |>
        filter(rank %in% c("green", "yellow", "blue"), is.na(phylum)) |>
        mutate(
            domain = if_else(reference %in% c("vir291", "vir303"), "d__Viruses", domain),
            lineage = if_else(reference %in% c("vir291", "vir303"), "l__Duplodnaviria", lineage),
            kingdom = if_else(reference %in% c("vir291", "vir303"), "k__Heunggongvirae", kingdom),
            phylum = if_else(reference %in% c("vir291", "vir303"), "p__Uroviricota", phylum),
            class = if_else(reference %in% c("vir291", "vir303"), "c__Caudoviricetes", class),
            phylum = if_else(!(reference %in% c("vir291", "vir303")), "Methanogenic viruses", phylum)
        )

    # Process AA hits and annotations
    aa_results_nored_filt <- aa_results_nored |>
        bind_rows(aa_results_nored_methano) |>
        filter(rank %in% c("green", "yellow", "blue"), phylum %in% vir_phyla)
    aa_hits_filt <- aa_data$hits |>
        inner_join(kapk_cdata |> select(label, figure_names)) |>
        select(-label) |>
        rename(label = figure_names) |>
        inner_join(aa_results_nored_filt |> distinct(reference, label, phylum, member_unit))

    pfam_hits <- aa_hits_filt |>
        distinct(reference, gene, label) |>
        left_join(protein_annotations$pfam) |>
        select(label, gene, reference, pfam_acc, pfam_name, pfam_description)
    phrog_hits <- aa_hits_filt |>
        distinct(reference, gene, label) |>
        left_join(protein_annotations$phrog) |>
        select(label, gene, reference, phrog_acc,
            phrog_host_domain = host_domain,
            phrog_category = category, phrog_annotation = annotation
        )
    aa_results_nored_filt_annotated <- pfam_hits |>
        inner_join(phrog_hits) |>
        select(-label) |>
        distinct()

    unannotated <- aa_results_nored_filt_annotated |>
        filter(is.na(pfam_acc) & is.na(phrog_acc)) |>
        select(gene) |>
        mutate(annotation = NA, functional_category = "Unknown") |>
        distinct()
    vir_annotations <- aa_results_nored_filt_annotated |>
        filter(phrog_annotation %in% virus_structural_protein) |>
        mutate(functional_category = "Virus structural protein") |>
        rename(annotation = phrog_annotation) |>
        select(gene, annotation, functional_category) |>
        bind_rows(
            aa_results_nored_filt_annotated |>
                filter(phrog_annotation %in% dna_rna_transactions_processing) |>
                mutate(functional_category = "DNA/RNA transactions/processing") |>
                rename(annotation = phrog_annotation) |>
                select(gene, annotation, functional_category),
            aa_results_nored_filt_annotated |>
                filter(phrog_annotation %in% lysis) |>
                mutate(functional_category = "Lysis") |>
                rename(annotation = phrog_annotation) |>
                select(gene, annotation, functional_category),
            aa_results_nored_filt_annotated |>
                filter(phrog_annotation %in% ta_rm_systems) |>
                mutate(functional_category = "TA/RM systems") |>
                rename(annotation = phrog_annotation) |>
                select(gene, annotation, functional_category),
            aa_results_nored_filt_annotated |>
                filter(phrog_annotation %in% amg_and_other) |>
                mutate(functional_category = "AMG and other") |>
                rename(annotation = phrog_annotation) |>
                select(gene, annotation, functional_category),
            unannotated
        )
    aa_results_nored_filt_annotated <- aa_results_nored_filt_annotated |>
        left_join(vir_annotations) |>
        mutate(functional_category = if_else(is.na(functional_category), "Unknown", functional_category)) |>
        left_join(read_xlsx("./data/taxonomy/viruses-aa/vir-pfam-annot-refined.xlsx") |>
            select(-n, category_new = Refined)) |>
        mutate(functional_category = if_else(!is.na(category_new) & functional_category == "Unknown",
            category_new, functional_category
        ))

    # Prepare NT/AA data
    nt_aa_data <- viral_data$all_viruses |>
        inner_join(imgvr_cdata |> select(reference, ecosystem_classification, jgi_url)) |>
        filter(breadth >= 0.1) |>
        inner_join(kapk_cdata_agg) |>
        distinct(
            reference, label, is_dmg, breadth, damage, phylum, member_unit,
            ecosystem_classification, jgi_url, read_ani_mean
        ) |>
        filter(phylum %in% vir_phyla) |>
        left_join(aa_results_nored_filt |> distinct(reference, label, avg_identity) |> mutate(is_aa = "NT/AA")) |>
        mutate(is_aa = if_else(is.na(is_aa), "NT", is_aa))

    # Define ecosystems and phyla order
    ecosystems <- c(
        "Environmental;Terrestrial;Soil;Wetlands",
        "Environmental;Aquatic;Marine;Coastal",
        "Environmental;Aquatic;Freshwater;River",
        "Environmental;Aquatic;Freshwater;Lake",
        "Environmental;Aquatic;Marine;Wetlands",
        "Environmental;Aquatic;Freshwater;Ice",
        "Environmental;Aquatic;Non-marine Saline and Alkaline;Saline",
        "Environmental;Aquatic;Deep subsurface;Groundwater"
    )
    vir_phyla_order <- aa_results_nored_filt |>
        distinct(reference, phylum) |>
        count(phylum, name = "n_ref", sort = TRUE) |>
        pull(phylum)

    # Generate figures
    figures <- list(
        nt_aa_damage_breadth = plot_nt_aa_damage_breadth(nt_aa_data, dmg_thresholds),
        ecosystem_anir = plot_ecosystem_anir(nt_aa_data, ecosystems),
        methano_hosts = plot_methano_hosts(aa_results_nored_methano, methano_cdata),
        nt_aa_phyla = plot_nt_aa_phyla(nt_aa_data, aa_results_nored_filt, vir_phyla_order),
        aa_hits = plot_aa_hits(aa_hits_filt, aa_results_nored_filt_annotated, vir_phyla_order),
        aa_ecosystem_similarity = plot_aa_ecosystem_similarity(aa_results_nored_filt, ecosystems, vir_phyla_order)
    )

    # Prepare tables
    tables <- list(
        derep_stats = derep_stats,
        viral_non_damaged = viral_data$non_damaged,
        viral_damaged = viral_data$damaged,
        viral_all = viral_data$all_viruses,
        aa_results = aa_results,
        aa_results_nored = aa_results_nored,
        aa_results_nored_filt = aa_results_nored_filt,
        aa_hits_filt = aa_hits_filt,
        aa_results_nored_filt_annotated = aa_results_nored_filt_annotated,
        nt_aa_data = nt_aa_data
    )

    # Return results
    list(
        tables = tables,
        figures = figures
    )
}

# Run main function and store results
results <- main()
