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

# Initialize settings
source("libs/lib.R")
source("./libs/lib-vir.R")
showtext_auto()

#' Process viral taxonomy and damage data
#' @return Processed viral taxonomy data
process_viral_taxonomy <- function() {
    # Read metadata
    kapk_cdata <- read_tsv("./data/cdata/KapK-cdata-manuscript-20221211.tsv")
    kapk_cdata_agg <- read_tsv("./data/cdata/KapK-cdata-agg-manuscript-20221211.tsv")

    # Read taxonomy
    tax_info <- read_tsv("./data/taxonomy/db/hires-organelles-viruses-arctic.tax.tsv",
        col_names = c("reference", "tax_string")
    ) |>
        separate(
            col = tax_string,
            sep = ";",
            into = c(
                "domain", "lineage", "kingdom", "phylum", "class",
                "order", "family", "genus", "species", "strain"
            )
        )

    # Read taxonomic profiles
    tax_data_agg <- read_tsv("./results/taxonomy/tp-mapping-filtered.nocontam.10M.agg.nobloom.tax.tsv.gz")
    tax_data <- read_tsv("./results/taxonomy/tp-mapping-filtered.nocontam.10M.nobloom.tax.tsv.gz")

    # Read damage data
    dmg_local_agg <- read_csv("./results/damage/tp-mdmg.weight-1.local.10M.agg.nobloom.csv.gz")
    dmg_local <- read_csv("./results/damage/tp-mdmg.weight-1.local.10M.nobloom.csv.gz")

    # Process derep stats
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

    # Get mapped reads
    label_to_file <- read_tsv("./data/cdata/KapK-label-to-file-20221211.tsv")
    reads_mapped <- read_tsv("./data/cdata/all.tp-classified-reads.tsv",
        col_names = c("file", "reads_mapped")
    ) |>
        inner_join(label_to_file) |>
        inner_join(kapk_cdata |> select(label, figure_names)) |>
        select(-label, -file) |>
        rename(label = figure_names)

    # Update derep stats
    derep_stats <- derep_stats |>
        inner_join(reads_mapped) |>
        mutate(
            prop_mapped = (reads_mapped / num_seqs),
            not_mapped = num_seqs - reads_mapped
        )

    # Join damage data
    tax_data_agg <- tax_data_agg |>
        inner_join(dmg_local_agg)

    # Get viral data
    tax_data_agg_viruses <- tax_data_agg |>
        inner_join(tax_info) |>
        filter(domain == "d__Viruses")

    list(
        kapk_cdata = kapk_cdata,
        kapk_cdata_agg = kapk_cdata_agg,
        tax_data = tax_data,
        tax_data_agg = tax_data_agg,
        tax_data_agg_viruses = tax_data_agg_viruses,
        tax_info = tax_info,
        derep_stats = derep_stats
    )
}

#' Process IMGVR data
#' @return Processed IMGVR data
process_imgvr_data <- function() {
    # Read IMGVR data
    imgvr_cdata <- read_tsv("./data/cdata/IMGVR_all_Sequence_information-high_confidence.tsv") |>
        rename(reference = UVIG) |>
        janitor::clean_names() |>
        mutate(taxon_oid = str_split(reference, "_") |> map_chr(3))

    # Add JGI URL
    jgi_url <- "https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid="
    imgvr_cdata_1 <- imgvr_cdata |> filter(!grepl("GVMAG", taxon_oid))
    imgvr_cdata_2 <- imgvr_cdata |>
        filter(grepl("GVMAG", taxon_oid)) |>
        mutate(taxon_oid = str_split(taxon_oid, "-") |> map_chr(3))

    bind_rows(imgvr_cdata_1, imgvr_cdata_2) |>
        mutate(jgi_url = paste0(jgi_url, taxon_oid))
}

#' Process amino acid data
#' @return List containing processed AA data
process_aa_data <- function() {
    # Read AA data
    viral_ngenes <- read_tsv("./data/taxonomy/db/IMGVR_all_proteins-high_confidence_derepG-archaeal.ngenes.tsv",
        col_names = c("reference", "n_genes_total")
    )

    aa_results <- read_tsv("/projects/fernandezguerra/apps/repos/kapk-microbial-analysis/data/taxonomy/viruses-aa/IMGVR4_derepG-archaea-profiling.group-abundances-agg.tsv.gz") |>
        rename(reference = group) |>
        inner_join(viral_ngenes) |>
        mutate(prot_coverage = n_genes / n_genes_total)

    aa_hits <- read_tsv("/projects/fernandezguerra/apps/repos/kapk-microbial-analysis/data/taxonomy/viruses-aa/IMGVR4_derepG-archaea-profiling.group-abundances.tsv.gz") |>
        rename(gene = reference, reference = group)

    list(
        viral_ngenes = viral_ngenes,
        aa_results = aa_results,
        aa_hits = aa_hits
    )
}

#' Process annotation data
#' @return List containing processed annotation data
process_annotation_data <- function() {
    # Read PFAM data
    pfamA_res <- read_tsv("/projects/fernandezguerra/apps/repos/kapk-microbial-analysis/data/taxonomy/viruses-aa/pfamA_res.tsv.gz",
        col_names = c(
            "gene", "pfam_acc", "probability", "e-value", "Score",
            "Cols", "q_start", "q_stop", "t_start", "t_stop", "q_len",
            "t_len", "q_cov", "t_cov", "pfam_name", "pfam_description"
        )
    )

    # Read PHROG data
    phrog_res <- read_tsv("/projects/fernandezguerra/apps/repos/kapk-microbial-analysis/data/taxonomy/viruses-aa/phrog_res.tsv.gz",
        col_names = c(
            "gene", "phrog_acc", "probability", "e-value", "Score",
            "Cols", "q_start", "q_stop", "t_start", "t_stop", "q_len",
            "t_len", "q_cov", "t_cov", "name", "description"
        )
    )

    phrog_cdata <- read_csv("/projects/fernandezguerra/apps/repos/kapk-microbial-analysis/data/taxonomy/viruses-aa/PHROG_index.csv") |>
        janitor::clean_names()

    # Filter results
    pfamA_res_filt <- pfamA_res |>
        filter(probability >= 90, (q_cov >= 0.4 | t_cov >= 0.4)) |>
        group_by(gene) |>
        arrange(desc(probability), desc(q_cov), desc(t_cov), .by_group = T) |>
        filter(row_number() == 1) |>
        ungroup()

    phrog_res_filt <- phrog_res |>
        filter(probability >= 90, (q_cov >= 0.4 | t_cov >= 0.4)) |>
        group_by(gene) |>
        arrange(desc(probability), desc(q_cov), desc(t_cov), .by_group = T) |>
        filter(row_number() == 1) |>
        ungroup() |>
        inner_join(phrog_cdata |> rename(phrog_acc = number_phrog))

    list(
        pfamA_res = pfamA_res,
        phrog_res = phrog_res,
        phrog_cdata = phrog_cdata,
        pfamA_res_filt = pfamA_res_filt,
        phrog_res_filt = phrog_res_filt
    )
}

#' Process and analyze viral data
#' @return List containing results and plots
analyze_viral_data <- function() {
    # Process data
    tax_data <- process_viral_taxonomy()
    imgvr_data <- process_imgvr_data()
    aa_data <- process_aa_data()
    annot_data <- process_annotation_data()

    # Filter KAPK data
    kapk_cdata_filt <- tax_data$kapk_cdata |>
        filter(label %in% (tax_data$tax_data$label |> unique()))

    # Process damaged and non-damaged viruses
    non_damaged_viruses <- tax_data$tax_data_agg_viruses |>
        filter(is_dmg == "Non-damaged", breadth >= 0.1) |>
        arrange(desc(coverage_mean)) |>
        inner_join(imgvr_data |> select(reference, ecosystem_classification, host_taxonomy_prediction, jgi_url)) |>
        inner_join(tax_data$kapk_cdata_agg)

    damaged_viruses <- tax_data$tax_data_agg_viruses |>
        filter(is_dmg == "Damaged", breadth >= 0.1) |>
        arrange(desc(coverage_mean)) |>
        inner_join(imgvr_data |> select(reference, ecosystem_classification, host_taxonomy_prediction, jgi_url)) |>
        inner_join(tax_data$kapk_cdata_agg)

    # Create results directory
    dir.create("./results/taxonomy", recursive = TRUE, showWarnings = FALSE)

    # Save results
    write_tsv(non_damaged_viruses, "./results/taxonomy/non-damaged-viruses-nt.tsv")
    write_tsv(damaged_viruses, "./results/taxonomy/damaged-viruses-nt.tsv")

    list(
        tax_data = tax_data,
        imgvr_data = imgvr_data,
        aa_data = aa_data,
        annot_data = annot_data,
        non_damaged_viruses = non_damaged_viruses,
        damaged_viruses = damaged_viruses
    )
}

# Run analysis
results <- analyze_viral_data()
