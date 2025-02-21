# Required libraries
library(tidyverse)
library(biomformat)
library(phyloseq)
library(showtext)
library(ggpol)

# Initialize settings
showtext_auto()
source("./libs/lib.R")

#' Load and process taxonomic data
#' @return list containing taxonomy data
load_taxonomy_data <- function() {
    # Get taxonomic data
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

    tax_info
}

#' Load and process source data
#' @param tax_info Taxonomy information
#' @return list containing source data
load_source_data <- function(tax_info) {
    # Load source metadata
    source_cdata <- read_tsv("./data/sourcetracker/cdata/kapk-biomes-download.txt",
        show_col_types = FALSE
    )

    source_cdata_short <- source_cdata |>
        select(label = run_accession, biome, biome_class, biome_subclass) |>
        mutate(SourceSink = "source") |>
        select(label, SourceSink, biome, biome_class, biome_subclass)

    # Load and process source abundance data
    source_data <- read_tsv("./data/sourcetracker/tp-mapping-filtered.summary.tsv.gz",
        show_col_types = FALSE
    ) |>
        filter(breadth >= 0.01) |>
        mutate(abundance = ifelse(tax_abund_tad == 0, tax_abund_read, tax_abund_tad)) |>
        select(label, reference, abundance)

    source_data_sp <- source_data |>
        inner_join(tax_info, by = "reference") |>
        filter(domain %in% c("d__Bacteria", "d__Archaea")) |>
        group_by(label, domain, species) |>
        summarize(
            abundance = janitor::round_half_up(gm_mean(abundance)),
            .groups = "drop"
        ) |>
        select(-domain)

    source_taxonomy <- tax_info |>
        filter(domain %in% c("d__Bacteria", "d__Archaea")) |>
        select(-reference) |>
        distinct() |>
        mutate(reference = species) |>
        filter(species %in% (source_data_sp |> select(species) |> distinct() |> pull(species)))

    list(
        cdata = source_cdata,
        cdata_short = source_cdata_short,
        data = source_data,
        data_sp = source_data_sp,
        taxonomy = source_taxonomy
    )
}

#' Load and process sink data
#' @return list containing sink data
load_sink_data <- function() {
    # Read sink data
    sink_data <- readRDS("./results/taxonomy/kapk_ps_ba_gm.rds")

    sink_data_sp <- speedyseq::psmelt(sink_data) |>
        as_tibble() |>
        select(species = OTU, abundance = Abundance, label = Sample)

    sink_taxonomy <- speedyseq::psmelt(sink_data) |>
        as_tibble() |>
        filter(domain %in% c("d__Bacteria", "d__Archaea")) |>
        select(reference = OTU, domain, phylum, class, order, family, genus, species) |>
        distinct()

    sink_cdata <- speedyseq::psmelt(sink_data) |>
        as_tibble() |>
        select(label = Sample, biome = member_unit, biome_class = member_unit, biome_subclass = member_unit) |>
        distinct() |>
        mutate(SourceSink = "sink")

    list(
        data = sink_data,
        data_sp = sink_data_sp,
        taxonomy = sink_taxonomy,
        cdata = sink_cdata
    )
}

#' Create phyloseq object from source and sink data
#' @param source_data_sp Source species data
#' @param sink_data_sp Sink species data
#' @param tax_info Taxonomy information
#' @param source_cdata_short Source metadata
#' @param sink_cdata Sink metadata
#' @return phyloseq object
create_phyloseq_object <- function(source_data_sp, sink_data_sp, tax_info, source_cdata_short, sink_cdata) {
    # Create OTU table
    st_df <- bind_rows(source_data_sp, sink_data_sp) |>
        pivot_wider(names_from = species, values_from = abundance, values_fill = 0) |>
        as.data.frame() |>
        column_to_rownames("label")

    # Create taxonomy table
    st_taxonomy <- tax_info |>
        filter(species %in% colnames(st_df)) |>
        filter(domain %in% c("d__Bacteria", "d__Archaea")) |>
        select(species, domain, phylum, class, order, family, genus) |>
        distinct() |>
        mutate(reference = species) |>
        distinct() |>
        as.data.frame() |>
        column_to_rownames("reference")

    # Create sample data
    st_cdata <- bind_rows(source_cdata_short, sink_cdata) |>
        as.data.frame() |>
        mutate(clabel = label) |>
        column_to_rownames("clabel")

    # Create phyloseq object
    phyloseq(
        otu_table(st_df, taxa_are_rows = FALSE),
        tax_table(as.matrix(st_taxonomy)),
        sample_data(st_cdata)
    )
}

#' Filter and export phyloseq object for SourceTracker
#' @param st_ps Phyloseq object
#' @return filtered phyloseq object and exported files
process_sourcetracker_data <- function(st_ps) {
    # Filter phyloseq object
    st_ps_filt <- get_st(st_ps, nsites = 0.01, vcoeff = 3)

    # Create and export BIOM file
    st_biom <- biomformat::make_biom(
        data = t(as((otu_table(st_ps_filt, taxa_are_rows = FALSE)), "matrix")),
        matrix_element_type = "int"
    )
    biomformat::write_biom(st_biom, biom_file = "./results/sourcetracker/st-biome-class-gm.biom")

    # Export mapping files for different biome levels
    # Biome class
    st_ps_filt |>
        speedyseq::psmelt() |>
        as_tibble() |>
        select(Sample, SourceSink, Env = biome_class) |>
        distinct() |>
        rename("#SampleID" = Sample) |>
        write_tsv("./results/sourcetracker/st-biome-class-gm.map")

    # Biome subclass
    st_ps_filt |>
        speedyseq::psmelt() |>
        as_tibble() |>
        select(Sample, SourceSink, Env = biome_subclass) |>
        distinct() |>
        rename("#SampleID" = Sample) |>
        write_tsv("./results/sourcetracker/st-biome-subclass-gm.map")

    # Biome
    st_ps_filt |>
        speedyseq::psmelt() |>
        as_tibble() |>
        select(Sample, SourceSink, Env = biome) |>
        distinct() |>
        rename("#SampleID" = Sample) |>
        write_tsv("./results/sourcetracker/st-biome-gm.map")

    st_ps_filt
}

#' Load and process SourceTracker results
#' @return list containing processed results
load_sourcetracker_results <- function() {
    # Read results
    st_results <- read_tsv("./results/sourcetracker/mixing_proportions.txt",
        show_col_types = FALSE
    ) |>
        rename(figure_names = "#SampleID")

    st_results_sd <- read_tsv("./results/sourcetracker/mixing_proportions_stds.txt",
        show_col_types = FALSE
    ) |>
        rename(figure_names = "#SampleID")

    # Convert to long format
    st_results_long <- st_results |>
        pivot_longer(
            cols = -figure_names,
            names_to = "biome",
            values_to = "proportion"
        )

    st_result_sd_long <- st_results_sd |>
        pivot_longer(
            cols = -figure_names,
            names_to = "biome",
            values_to = "std"
        )

    list(
        results = st_results,
        results_sd = st_results_sd,
        results_long = st_results_long,
        results_sd_long = st_result_sd_long
    )
}

#' Process feature contributions
#' @return processed feature contributions
process_feature_contributions <- function() {
    cont_files <- list.files(
        "./results/sourcetracker",
        pattern = "feature_table.txt",
        full.names = TRUE
    )

    map_dfr(cont_files, function(X) {
        label <- basename(X)
        label <- gsub(".feature_table.txt", "", label)
        read_tsv(X, show_col_types = FALSE) |>
            rename(biome = "...1") |>
            pivot_longer(-biome, names_to = "species", values_to = "contribution") |>
            filter(contribution > 0) |>
            mutate(label = label)
    })
}

#' Process contribution data
#' @param contribution_data Feature contribution data
#' @param tax_info Taxonomy information
#' @return top contributions by genus
get_top_contributions <- function(contribution_data, tax_info) {
    contribution_data |>
        select(-label) |>
        inner_join(tax_info |> select(species, genus)) |>
        group_by(biome, genus) |>
        summarise(contribution = mean(contribution)) |>
        arrange(desc(contribution), .by_group = TRUE) |>
        mutate(rank = row_number()) |>
        filter(rank <= 10, contribution > 10)
}

#' Plot source biome data
#' @param source_biome_data Source biome data
#' @return ggplot object
plot_source_biomes <- function(source_biome_data) {
    source_biome_data |>
        select(label, biome_subclass) |>
        distinct() |>
        count(biome_subclass) |>
        mutate(biome_subclass = reorder(biome_subclass, n, decreasing = TRUE)) |>
        ggplot(aes(x = biome_subclass, y = n)) +
        geom_col(color = "black", width = 1, linewidth = 0.5, fill = "#D6604D") +
        geom_text(aes(label = n), vjust = -0.3, hjust = 0.5) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            text = element_text(size = 12)
        ) +
        labs(x = "Biome subclass", y = "Number of samples") +
        xlab("")
}

#' Plot biome classification results
#' @param st_results_long Long format SourceTracker results
#' @param st_result_sd_long Long format SourceTracker standard deviations
#' @param sink_data Sink data
#' @return ggplot object
plot_biome_classification <- function(st_results_long, st_result_sd_long, sink_data) {
    # Get top biomes
    top_biomes <- st_results_long |>
        group_by(biome) |>
        summarise(mean = max(proportion)) |>
        arrange(desc(mean)) |>
        ungroup() |>
        top_n(4)

    # Prepare data for plotting
    data <- st_results_long |>
        filter(biome %in% top_biomes$biome) |>
        inner_join(st_result_sd_long) |>
        inner_join(sink_data |>
            speedyseq::psmelt() |>
            select(Sample, member_unit, site_rnk) |>
            distinct() |>
            rename(figure_names = Sample)) |>
        mutate(biome = gsub("root:Environmental:", "", biome)) |>
        mutate(
            figure_names = fct_reorder(figure_names, site_rnk),
            member_unit = fct_relevel(member_unit, rev(c("B3", "B2", "B1"))),
            biome = fct_reorder(biome, proportion, max),
            type = ifelse(biome == "Unknown", "Unknown", "Known")
        )

    # Create plot
    ggplot(data, aes(x = biome, y = proportion, fill = member_unit)) +
        geom_point(alpha = 0.2, position = position_jitter(width = 0.1, seed = 3922)) +
        geom_pointrange(
            mapping = aes(x = biome, y = proportion),
            stat = "summary",
            fun.min = function(z) {
                quantile(z, 0.25)
            },
            fun.max = function(z) {
                quantile(z, 0.75)
            },
            fun = median,
            shape = 21,
            size = 0.5
        ) +
        ggpubr::rotate() +
        facet_grid(type ~ member_unit, space = "free", scale = "free_y") +
        theme_bw() +
        theme(
            axis.title.y = element_blank(),
            text = element_text(size = 10),
            strip.background = element_blank(),
            legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        ylab("Proportion") +
        scale_y_sqrt(labels = scales::percent) +
        scale_fill_manual(values = c("#C2E4EF", "#EAECCC", "#FEDA8B"))
}

#' Plot decOM results
#' @param decom_results DecOM results
#' @param label_to_file Label to file mapping
#' @param kapk_cdata Sample metadata
#' @param st_results SourceTracker results
#' @return processed decOM results
process_decom_results <- function(decom_results, label_to_file, kapk_cdata, st_results) {
    decom_results |>
        select(file = Sink, starts_with("p_")) |>
        pivot_longer(-file, names_to = "biome", values_to = "proportion") |>
        inner_join(label_to_file) |>
        inner_join(kapk_cdata) |>
        filter(figure_names %in% st_results$figure_names) |>
        mutate(source = gsub("^p_", "", biome)) |>
        group_by(biome) |>
        summarise(mean = mean(proportion)) |>
        arrange(desc(mean)) |>
        ungroup()
}

#' Plot biome distributions by class and subclass
#' @param source_biome_data Source biome data
#' @return list of plots
plot_biome_distributions <- function(source_biome_data) {
    # Biome class plot
    biome_class_plot <- source_biome_data |>
        select(label, biome_class) |>
        distinct() |>
        count(biome_class) |>
        mutate(biome_class = reorder(biome_class, n, decreasing = TRUE)) |>
        ggplot(aes(x = biome_class, y = n)) +
        geom_col() +
        geom_text(aes(label = n), vjust = -0.25, hjust = 0.5) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(x = "Biome class", y = "Number of samples")

    # Biome subclass plot (original plot)
    biome_subclass_plot <- source_biome_data |>
        select(label, biome_subclass) |>
        distinct() |>
        count(biome_subclass) |>
        mutate(biome_subclass = reorder(biome_subclass, n, decreasing = TRUE)) |>
        ggplot(aes(x = biome_subclass, y = n)) +
        geom_col(color = "black", width = 1, linewidth = 0.5, fill = "#D6604D") +
        geom_text(aes(label = n), vjust = -0.3, hjust = 0.5) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            text = element_text(size = 12)
        ) +
        labs(x = "Biome subclass", y = "Number of samples") +
        xlab("")

    # Biome plot
    biome_plot <- source_biome_data |>
        select(label, biome) |>
        distinct() |>
        count(biome) |>
        mutate(biome = reorder(biome, n, decreasing = TRUE)) |>
        ggplot(aes(x = biome, y = n)) +
        geom_col() +
        geom_text(aes(label = n), vjust = -0.25, hjust = 0.5) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(x = "Biome", y = "Number of samples")

    list(
        class = biome_class_plot,
        subclass = biome_subclass_plot,
        biome = biome_plot
    )
}

#' Main execution function
#' @return list containing results and all plots
main <- function() {
    # Load initial data
    tax_info <- load_taxonomy_data()
    source_data <- load_source_data(tax_info)
    sink_data <- load_sink_data()

    # Create and process phyloseq object
    st_ps <- create_phyloseq_object(
        source_data$data_sp,
        sink_data$data_sp,
        tax_info,
        source_data$cdata_short,
        sink_data$cdata
    )

    # Filter and export for SourceTracker
    st_ps_filt <- process_sourcetracker_data(st_ps)

    # Load SourceTracker results
    st_results <- load_sourcetracker_results()

    # Process source biome data
    source_biome_data <- st_ps_filt |>
        speedyseq::psmelt() |>
        as_tibble() |>
        distinct() |>
        filter(SourceSink == "source") |>
        inner_join(source_data$cdata |>
            rename(label = run_accession) |>
            select(label, biome, biome_class, biome_subclass))

    # Create all biome distribution plots
    biome_plots <- plot_biome_distributions(source_biome_data)

    # Create biome classification plot
    biome_class_plot <- plot_biome_classification(
        st_results$results_long,
        st_results$results_sd_long,
        sink_data$data
    )

    # Process feature contributions
    contribution_data <- process_feature_contributions()
    top_contributions <- get_top_contributions(contribution_data, tax_info)

    # Load additional metadata
    decom_results <- read_csv("./results/sourcetracker/decOM_output-dmg-k29.csv",
        show_col_types = FALSE
    )
    label_to_file <- read_tsv("./data/cdata/KapK-label-to-file-20221211.tsv",
        show_col_types = FALSE
    )
    kapk_cdata <- read_tsv("./data/cdata/KapK-cdata-manuscript-20221211.tsv",
        show_col_types = FALSE
    )

    # Process decOM results
    decom_processed <- process_decom_results(
        decom_results,
        label_to_file,
        kapk_cdata,
        st_results$results
    )

    # Return results
    list(
        data = list(
            source_data = source_data,
            sink_data = sink_data,
            st_results = st_results,
            contribution_data = contribution_data,
            decom_results = decom_processed,
            top_contributions = top_contributions,
            source_biome_data = source_biome_data
        ),
        plots = list(
            biome_distributions = biome_plots,
            biome_classification = biome_class_plot
        ),
        phyloseq = list(
            raw = st_ps,
            filtered = st_ps_filt
        )
    )
}

# Run main function and store results
results <- main()

# If plots need to be saved, they can be saved afterward:
# ggsave("biome-class.pdf", results$plots$biome_distributions$class)
# ggsave("biome-classification.pdf", results$plots$biome_classification)
