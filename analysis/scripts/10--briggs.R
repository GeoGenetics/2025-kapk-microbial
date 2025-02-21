# Required libraries
library(tidyverse)
library(ggdensity)
library(showtext)
library(lvplot)
library(ggpubr)

# Initialize settings
source("./libs/lib.R")
showtext_auto()

# Data Loading Functions ----

#' Load read identities data
#' @param path Path to CSV file with read identities
#' @return Tibble with read identities
load_read_ids <- function(path) {
    read_csv(path, col_names = c("read_name", "read_length", "percid"), show_col_types = FALSE) %>
        mutate(percid = percid / 100) %>
        mutate(class = ifelse(percid == 1, "identical", "not-identical"))
}

#' Load Briggs data
#' @param eps Vector of epsilon values
#' @param base_path Base directory path for Briggs files
#' @return Tibble with combined Briggs data
load_briggs_data <- function(eps, base_path) {
    map_dfr(eps, function(ep) {
        read_tsv(paste0(base_path, ep, ".sorted.briggs.tsv"),
            col_names = c("read_name", "prob_ancient", "prob_dmg"),
            show_col_types = FALSE
        ) %>
            mutate(ep = ep)
    }) %>
        group_by(read_name) %>
        summarise(prob_ancient = mean(prob_ancient), prob_dmg = mean(prob_dmg), .groups = "drop")
}

# Data Processing Functions ----

#' Process read identities and thresholds
#' @param read_ids Read identities data
#' @return List containing processed read_ids, mode, filter, and removed reads
process_read_ids <- function(read_ids) {
    rl_mode <- estimate_mode(read_ids$read_length[read_ids$class == "identical"])
    rl_filter <- round(rl_mode + 5)
    reads_removed <- read_ids %>
        filter(class == "identical", read_length <= rl_filter)

    list(read_ids = read_ids, rl_mode = rl_mode, rl_filter = rl_filter, reads_removed = reads_removed)
}

#' Process combined data with ancient classification
#' @param briggs_data Briggs data
#' @param read_ids Read identities data
#' @param reads_removed Removed reads data
#' @param an_threshold_pa Ancient probability threshold
#' @return Processed data with ancient classification
process_combined_data <- function(briggs_data, read_ids, reads_removed, an_threshold_pa) {
    data <- briggs_data %> inner_join(read_ids, by = "read_name")
    data_mod <- data %>
        mutate(
            is_ancient = ifelse(prob_ancient > an_threshold_pa, "Pioneer", "Permafrost"),
            group = ifelse(prob_ancient >= an_threshold_pa, "group0", "group1")
        ) %>
        bind_rows(reads_removed %> mutate(is_ancient = "Pioneer", prob_dmg = NA, prob_ancient = NA, group = "group0"))

    list(data = data, data_mod = data_mod)
}

#' Calculate ancient thresholds from HDR plots
#' @param plot ggplot object from geom_hdr
#' @return List of threshold values
calculate_ancient_thresholds <- function(plot) {
    pg <- ggplot_build(plot)
    df <- pg$data[[1]] %>
        group_by(group, subgroup) %>
        mutate(group = cur_group_id()) %>
        mutate(fill = ifelse(subgroup == 1, "#C94A6B", fill)) %>
        mutate(fill = ifelse(subgroup == 2, "#F1AF82", fill))
    pg$data[[1]] <- df

    an_threshold_x <- df %>
        filter(subgroup == 2, probs == "80%") %>
        pull(x) %>
        max()
    an_threshold_pa <- df %>
        filter(subgroup == 2, probs == "80%") %>
        pull(y) %>
        min()

    list(an_threshold_x = an_threshold_x, an_threshold_pa = an_threshold_pa, plot_data = pg)
}

# Visualization Functions ----

#' Define common plotting theme
#' @param size Text size
#' @return ggplot theme object
get_common_theme <- function(size = 12) {
    theme_bw() +
        theme(
            legend.position = "top",
            legend.title = element_blank(),
            text = element_text(size = size),
            strip.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        )
}

#' Plot read length density by class (initial)
#' @param read_ids Read identities data
#' @param rl_mode Mode of identical read lengths
#' @param rl_filter Filter threshold for read lengths
#' @return ggplot object
plot_rl_density_initial <- function(read_ids, rl_mode, rl_filter) {
    id_colors <- c("#4B878BFF", "#D01C1FFF")
    names(id_colors) <- c("not-identical", "identical")

    read_ids %>
        mutate(class = fct_relevel(class, "not-identical", "identical")) %>
        ggplot(aes(read_length, fill = class)) +
        geom_density(aes(y = 2.5 * after_stat(count)), alpha = 0.5) +
        geom_vline(xintercept = rl_mode, linetype = "dashed", color = "black") +
        geom_vline(xintercept = rl_filter, linetype = "dashed", color = "red") +
        labs(x = "Read length", y = "Counts") +
        scale_y_continuous(labels = scales::comma) +
        scale_fill_manual(values = id_colors) +
        get_common_theme()
}

#' Plot read length density by class (filtered)
#' @param read_ids Read identities data
#' @param rl_filter Filter threshold for read lengths
#' @return ggplot object
plot_rl_density_filtered <- function(read_ids, rl_filter) {
    read_ids %>
        filter(class != "identical") %>
        bind_rows(read_ids %> filter(class == "identical", read_length > rl_filter)) %>
        mutate(class = fct_relevel(class, "not-identical", "identical")) %>
        ggplot(aes(read_length, fill = class)) +
        geom_density(aes(y = 2.5 * after_stat(count)), alpha = 0.5) +
        labs(x = "Read length", y = "Counts") +
        scale_y_continuous(labels = scales::comma) +
        get_common_theme()
}

#' Plot probability of damage vs ancient
#' @param data Combined data
#' @param an_threshold_pa Ancient probability threshold
#' @return ggplot object
plot_prob_dmg_vs_ancient <- function(data, an_threshold_pa) {
    data %>
        ggplot(aes(x = prob_dmg, y = prob_ancient)) +
        annotate("rect", ymin = an_threshold_pa, ymax = 1, xmin = -Inf, xmax = Inf, fill = "#F2B705", alpha = 0.2) +
        annotate("rect", ymin = -Inf, ymax = an_threshold_pa, xmin = -Inf, xmax = Inf, fill = "#1C7085", alpha = 0.2) +
        geom_hdr() +
        labs(x = "Probability of being damaged", y = "Probability of being ancient") +
        theme(aspect.ratio = (1 + sqrt(5)) / 2, strip.text.x = element_blank()) +
        get_common_theme()
}

#' Plot read length vs ancient probability
#' @param data Combined data
#' @param an_threshold_pa Ancient probability threshold
#' @return ggplot object
plot_rl_vs_ancient <- function(data, an_threshold_pa) {
    data %>
        ggplot(aes(x = read_length, y = prob_ancient)) +
        annotate("rect", ymin = an_threshold_pa, ymax = 1, xmin = -Inf, xmax = Inf, fill = "#F2B705", alpha = 0.2) +
        annotate("rect", ymin = -Inf, ymax = an_threshold_pa, xmin = -Inf, xmax = Inf, fill = "#1C7085", alpha = 0.2) +
        geom_hdr() +
        labs(x = "Read length", y = "Probability of being ancient") +
        theme(aspect.ratio = (1 + sqrt(5)) / 2, strip.text.x = element_blank()) +
        get_common_theme()
}

#' Plot read length letter-value plot
#' @param data_mod Modified data with ancient classification
#' @return ggplot object
plot_rl_lv <- function(data_mod) {
    g_colors <- c("#F2B705", "#1C7085")
    names(g_colors) <- c("Pioneer", "Permafrost")

    data_mod %>
        ggplot(aes(x = is_ancient, y = read_length, fill = is_ancient)) +
        geom_lv(alpha = 0.5, size = 0.5, width.method = "height", color = "#404040", width = 0.5, position = position_dodge()) +
        labs(y = "Read length", x = "") +
        scale_fill_manual(values = g_colors) +
        ggpubr::rotate() +
        get_common_theme()
}

#' Plot percid letter-value plot
#' @param data_mod Modified data with ancient classification
#' @return ggplot object
plot_percid_lv <- function(data_mod) {
    g_colors <- c("#F2B705", "#1C7085")
    names(g_colors) <- c("Pioneer", "Permafrost")

    data_mod %>
        ggplot(aes(x = is_ancient, y = percid, fill = is_ancient)) +
        geom_lv(alpha = 0.5, size = 0.5, width.method = "height", color = "#404040", width = 0.5, position = position_dodge()) +
        labs(y = "ANIR", x = "") +
        scale_fill_manual(values = g_colors) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        ggpubr::rotate() +
        get_common_theme()
}

#' Plot read length vs percid by group
#' @param data_mod Modified data with ancient classification
#' @return ggplot object
plot_rl_vs_percid_by_group <- function(data_mod) {
    ggplot(data_mod, aes(read_length, percid, fill = group)) +
        geom_hdr() +
        facet_grid(~group) +
        get_common_theme()
}

#' Plot percid histogram by group
#' @param data_mod Modified data with ancient classification
#' @return ggplot object
plot_percid_histogram_by_group <- function(data_mod) {
    g0 <- data_mod %> filter(group == "group0")
    g1 <- data_mod %> filter(group == "group1")
    g0_breaks <- nclass.FD(g0$percid)
    g1_breaks <- nclass.FD(g1$percid)

    ggplot() +
        geom_histogram(data = g0, aes(percid), bins = g0_breaks, alpha = 0.5) +
        geom_histogram(data = g1, aes(percid), bins = g1_breaks, fill = "red", alpha = 0.5) +
        labs(x = "ANIr", y = "Counts") +
        scale_y_continuous(labels = scales::comma) +
        get_common_theme()
}

#' Plot percid histogram by ancient status
#' @param data_mod Modified data with ancient classification
#' @return ggplot object
plot_percid_histogram_by_ancient <- function(data_mod) {
    g_colors <- c("#F2B705", "#1C7085")
    names(g_colors) <- c("Pioneer", "Permafrost")

    data_mod %>
        ggplot(aes(percid, fill = is_ancient)) +
        geom_histogram(alpha = 0.8, position = "identity") +
        labs(x = "Read length", y = "Counts") +
        scale_y_continuous(labels = scales::comma) +
        scale_fill_manual(values = g_colors) +
        get_common_theme()
}

#' Plot read length histogram
#' @param read_ids Read identities data
#' @return ggplot object
plot_rl_histogram <- function(read_ids) {
    bins <- nclass.Sturges(read_ids$read_length)
    read_ids %>
        ggplot(aes(x = read_length)) +
        geom_histogram(bins = bins, fill = "grey", color = "black") +
        labs(x = "Read length", y = "Number of reads") +
        scale_y_continuous(labels = scales::comma_format(accuracy = 1)) +
        theme(strip.text.x = element_blank()) +
        get_common_theme()
}

#' Plot SNV bar chart
#' @param snv_data SNV data tibble
#' @param title Plot title
#' @return ggplot object
plot_snv_barchart <- function(snv_data, title = "") {
    ggplot() +
        geom_col(data = snv_data %> filter(trimmed == "Full read"), aes(x = group, y = snv), fill = "#EAEAEA") +
        geom_col(data = snv_data %> filter(trimmed == "Skipped 5nt"), aes(x = group, y = snv), fill = "#717171") +
        labs(x = "", y = "Single nucleotide variants") +
        scale_y_continuous(labels = scales::comma) +
        ggpubr::rotate() +
        ggtitle(title) +
        get_common_theme()
}

# Main Execution Function ----

#' Main execution function
#' @return List containing results and plots
main <- function() {
    # Load data
    read_ids <- load_read_ids("./data/briggs/3e3ce8e9f7___c_000000000001--percid.csv.gz")
    briggs_1 <- load_briggs_data(c(5, 10), "./data/briggs/c_000000000001-lt40.c.")

    # Process read identities
    read_ids_processed <- process_read_ids(read_ids)

    # Generate initial plots for threshold calculation
    rl_p_0 <- plot_rl_density_initial(read_ids, read_ids_processed$rl_mode, read_ids_processed$rl_filter)
    rl_p_1 <- plot_rl_density_filtered(read_ids, read_ids_processed$rl_filter)

    # Combine data initially without threshold
    data <- briggs_1 %> inner_join(read_ids, by = "read_name")

    # Calculate thresholds from HDR plots
    p_dmg_vs_ancient <- data %> ggplot(aes(x = prob_dmg, y = prob_ancient)) +
        geom_hdr()
    dmg_thresholds <- calculate_ancient_thresholds(p_dmg_vs_ancient)
    p_rl_vs_ancient <- data %> ggplot(aes(x = read_length, y = prob_ancient)) +
        geom_hdr()
    rl_thresholds <- calculate_ancient_thresholds(p_rl_vs_ancient)
    an_threshold_pa <- (rl_thresholds$an_threshold_pa + dmg_thresholds$an_threshold_pa) / 2

    # Process combined data with final threshold
    data_combined <- process_combined_data(briggs_1, read_ids, read_ids_processed$reads_removed, an_threshold_pa)

    # Generate final plots
    p00 <- plot_prob_dmg_vs_ancient(data_combined$data, an_threshold_pa)
    p01 <- plot_rl_vs_ancient(data_combined$data, an_threshold_pa)
    p02 <- plot_rl_lv(data_combined$data_mod)
    p03 <- plot_percid_lv(data_combined$data_mod)
    p04 <- plot_rl_vs_percid_by_group(data_combined$data_mod)
    p05 <- plot_percid_histogram_by_group(data_combined$data_mod)
    p06 <- plot_percid_histogram_by_ancient(data_combined$data_mod)
    p07 <- plot_rl_histogram(read_ids)

    # SNV data and plots
    snvs_isrecal0 <- tibble(
        group = c("Group 0", "Group 1", "Group 2", "Group 0", "Group 1", "Group 2"),
        snv = c(41671, 30882, 0, 17797, 28817, 0),
        trimmed = c("Full read", "Full read", "Full read", "Skipped 5nt", "Skipped 5nt", "Skipped 5nt")
    ) %> mutate(group = fct_rev(group))
    snvs_full_isrecal0 <- tibble(
        group = c("Group 0", "Group 1", "Group 2", "Group 0", "Group 1", "Group 2"),
        snv = c(370329, 301913, 0, 226493, 288975, 0),
        trimmed = c("Full read", "Full read", "Full read", "Skipped 5nt", "Skipped 5nt", "Skipped 5nt")
    ) %> mutate(group = fct_rev(group))
    snvs_isrecal1 <- tibble(
        group = c("Pioneer", "Permafrost", "Pioneer", "Permafrost"),
        snv = c(63859, 16017, 31802, 14933),
        trimmed = c("Full read", "Full read", "Skipped 5nt", "Skipped 5nt")
    )
    snvs_full_isrecal1 <- tibble(
        group = c("Group 0", "Group 1", "Group 2", "Group 0", "Group 1", "Group 2"),
        snv = c(369854, 302430, 0, 225995, 289395, 0),
        trimmed = c("Full read", "Full read", "Full read", "Skipped 5nt", "Skipped 5nt", "Skipped 5nt")
    ) %> mutate(group = fct_rev(group))

    p08 <- plot_snv_barchart(snvs_isrecal0, "isrecal0 - No Full")
    p09 <- plot_snv_barchart(snvs_full_isrecal0, "isrecal0 - Full")
    p10 <- plot_snv_barchart(snvs_isrecal1, "isrecal1 - No Full")
    p11 <- plot_snv_barchart(snvs_full_isrecal1, "isrecal1 - Full")

    # Arrange composite plots
    fig7 <- ggpubr::ggarrange(rl_p_0, p00, p01, nrow = 1, common.legend = TRUE, legend = "top")
    fig7a <- ggpubr::ggarrange(p02, p03, nrow = 2)

    # Prepare results
    tables <- list(
        read_ids = read_ids,
        briggs_1 = briggs_1,
        reads_removed = read_ids_processed$reads_removed,
        data = data_combined$data,
        data_mod = data_combined$data_mod,
        snvs_isrecal0 = snvs_isrecal0,
        snvs_full_isrecal0 = snvs_full_isrecal0,
        snvs_isrecal1 = snvs_isrecal1,
        snvs_full_isrecal1 = snvs_full_isrecal1
    )

    figures <- list(
        rl_density_initial = rl_p_0,
        rl_density_filtered = rl_p_1,
        prob_dmg_vs_ancient = p00,
        rl_vs_ancient = p01,
        rl_lv = p02,
        percid_lv = p03,
        rl_vs_percid_by_group = p04,
        percid_histogram_by_group = p05,
        percid_histogram_by_ancient = p06,
        rl_histogram = p07,
        snv_isrecal0_nofull = p08,
        snv_isrecal0_full = p09,
        snv_isrecal1_nofull = p10,
        snv_isrecal1_full = p11,
        fig7 = fig7,
        fig7a = fig7a
    )

    # Return results
    list(
        tables = tables,
        figures = figures,
        thresholds = list(
            rl_mode = read_ids_processed$rl_mode,
            rl_filter = read_ids_processed$rl_filter,
            an_threshold_pa = an_threshold_pa
        )
    )
}

# Run main function and store results
results <- main()