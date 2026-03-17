# 07a--agp-functional.R
# Functional profiling using DART/AGP damage data
#
# Direct port of 07--functional-profiling.R, replacing:
#   data/function/kegg/{all,damaged}/kegg-modules-summary.tsv.gz
#   data/function/dbcan/{all,damaged}/dbcan.group-abundances-agg.tsv.gz
# With per-sample AGP outputs:
#   results/functional_agp/kegg/{sample}/anvio_modules.txt
#   results/functional_agp/kegg/{sample}/emi.protein.slim.tsv.gz
#   results/functional_agp/cazy/{sample}/cazy_emi.functional.tsv
#   results/functional_agp/cazy/{sample}/cazy_damage.slim.tsv.gz
#
# The "Damaged" condition uses DART mean_p_damaged >= DAMAGE_THRESHOLD (default 0.5)
# instead of a separate damaged-reads coverage file.

library(tidyverse)
library(data.table)
library(janitor)
library(showtext)
library(ggh4x)

source("./libs/lib.R")
showtext_auto()


# =============================================================================
# Constants (identical to 07--functional-profiling.R)
# =============================================================================

UNIT_COLORS <- c("B1" = "#4097AA", "B2" = "#C94A6B", "B3" = "#AFCCB8")

CELLULOSE_FAMILIES <- c("GH5","GH9","3.2.1.4","GH51","GH6","GH7","GH48","3.2.1.91")

HEMICELLULOSE_FAMILIES <- c(
    "GH5","GH8","GH10","GH11","GH43","3.2.1.8","GH3","GH30",
    "GH39","GH52","GH54","GH116","GH120","3.2.1.37","GH67","GH115",
    "3.2.1.139","CE1","CE2","CE3","CE4","CE5","CE6","CE7","CE12","3.1.1.72"
)

FERMENTATION_MODULES <- c(
    "ethanol fermentation", "acetatogenesis",
    "lactate fermentation", "propionate fermentation"
)

METHANE_MODULES_PATHWAY <- c(
    "Methanogenesis, CO2 => methane",
    "Methanogenesis, acetate => methane",
    "Methanogenesis, methanol => methane",
    "Methanogenesis, methylamine/dimethylamine/trimethylamine => methane"
)

METHANE_MODULES_WOODCROFT <- c("acetoclastic methanogenesis", "hydrogenotrophic methanogenesis")

DATA_ALL_ORDER <- c(
    "Hemicellulose xylan degradation",
    "Cellulose degradation",
    "Methanotrophy",
    "Methylotrophic methanogenesis (methanol)",
    "Methylotrophic methanogenesis (methylamine)",
    "Acetoclastic methanogenesis",
    "Hydrogenotrophic methanogenesis",
    "Ethanol fermentation",
    "Propionate fermentation",
    "Lactate fermentation",
    "Acetogen"
)

# Dotplot pathway facet order (top → bottom), mirroring heatmap:
#   Central carbon block (top): Acetatogenesis → Ethanol
#   Methane block (bottom):     Hydrogenotrophic → Methanotrophy
DOTPLOT_PATHWAY_ORDER <- c(
    "Acetatogenesis",
    "Acetogenesis",
    "Lactate fermentation",
    "Propionate fermentation",
    "Ethanol fermentation",
    "Hydrogenotrophic methanogenesis",
    "Acetoclastic methanogenesis",
    "Methylotrophic methanogenesis",
    "Methanotrophy"
)

# DART 3-class thresholds (from src/cli/cmd_damage_annotate.cpp defaults)
DAMAGED_THRESHOLD <- 0.60  # p_protein_damaged >= this → damaged (ancient)
# NOTE: In the slim TSVs, DART returns p=0.5 (its Bayesian prior) when no damage
# signal is detected. The slim filter is p_protein_damaged >= 0.5, so the data is
# bimodal: 0.5 (no signal) and ~1.0 (high confidence damage). There are NO proteins
# with 0.25 < p < 0.60 in the slim data. We set MODERN_THRESHOLD = 0.50 so that
# p=0.5 proteins (no-signal/undetermined) serve as the reference class for Fisher's test.
MODERN_THRESHOLD  <- 0.50  # p_protein_damaged <= this → reference (no damage signal)
# (nothing falls in the 0.50–0.60 ambiguous zone in slim data; zone is structurally empty)
DAMAGE_THRESHOLD  <- 0.70  # module-level mean posterior >= this → bordered tile in heatmap
                           # matches the reference-level filter used in the virome analysis

# Figure modules shown in the heatmap (highlighted in enrichment plots)
FIGURE_MODULES_KEGG <- c(
    FERMENTATION_MODULES,
    METHANE_MODULES_PATHWAY,
    METHANE_MODULES_WOODCROFT,
    "Methane oxidation, methanotroph, methane => formaldehyde",
    "Acetogen"
)

# Key enzyme KOs for the carbon cycling pathways shown in Fig. 5.
# pathway / class / pathway_order / class_order mirror DATA_ALL_ORDER and the
# class levels used in plot_functional_heatmap so both figures share the same
# top-to-bottom row ordering.
KEY_KOS_TABLE <- data.table(
    ko = c(
        "K10944", "K16157", "K16158",          # pmoA/mmoXY – methanotrophy
        "K14081", "K14080",                    # mtaBA   – methylotrophic MG
        "K00193", "K00194", "K00197",          # cdhACE  – acetoclastic MG
        "K00399", "K00401", "K00402",          # mcrAGB  – hydrogenotrophic MG
        "K13953", "K00121", "K00128",          # adhE/adhP/aldB – ethanol ferm.
        "K00625", "K01026", "K18118",          # MCEE/pct/pccB – propionate ferm.
        "K00016",                              # ldh     – lactate ferm.
        "K01938", "K01491",                    # fhs/folD – Wood-Ljungdahl
        "K00925"                               # ack     – acetatogenesis
    ),
    gene = c(
        "pmoA/mmoX", "mmoX", "mmoY",
        "mtaB", "mtaA",
        "cdhA", "cdhC", "cdhE",
        "mcrA", "mcrG", "mcrB",
        "adhE", "adhP", "aldB",
        "MCEE", "pct", "pccB",
        "ldh",
        "fhs", "folD",
        "ack"
    ),
    pathway = c(
        rep("Methanotrophy", 3),
        rep("Methylotrophic methanogenesis", 2),
        rep("Acetoclastic methanogenesis", 3),
        rep("Hydrogenotrophic methanogenesis", 3),
        rep("Ethanol fermentation", 3),
        rep("Propionate fermentation", 3),
        "Lactate fermentation",
        rep("Acetogenesis", 2),
        "Acetatogenesis"
    ),
    class = c(
        rep("Methane metabolism", 11),
        rep("Central carbon metabolism", 10)
    ),
    pathway_order = c(
        rep(1L, 3), rep(2L, 2), rep(3L, 3), rep(4L, 3),  # methane metabolism
        rep(5L, 3), rep(6L, 3), 7L, rep(8L, 2), 9L        # central carbon
    ),
    class_order = c(rep(1L, 11), rep(2L, 10)),
    ko_idx      = 1L:21L  # row position in table (used for within-pathway gene ordering)
)


# =============================================================================
# DART Data Loading
# =============================================================================

#' Discover sample folders
discover_samples <- function(db_dir) {
    list.dirs(db_dir, full.names = FALSE, recursive = FALSE)
}

#' Load KEGG module completeness for all samples
#'
#' Replaces kegg-modules-summary.tsv.gz. Returns one row per sample per module
#' with columns: label (= full MD5), module, module_name, module_class,
#' module_category, module_subcategory, stepwise_module_completeness, avg_coverage.
load_kegg_modules <- function(kegg_dir, sample_meta) {
    samples <- intersect(discover_samples(kegg_dir), sample_meta$sample)

    dt <- rbindlist(lapply(samples, function(s) {
        path <- file.path(kegg_dir, s, "anvio_modules.txt")
        if (!file.exists(path)) { warning("Missing: ", path); return(NULL) }
        d <- fread(path, showProgress = FALSE)
        cov_col <- grep("avg_coverage$", names(d), value = TRUE)[1]
        if (cov_col != "avg_coverage") setnames(d, cov_col, "avg_coverage")

        # --- DART EM ancient coverage per module ---
        # Load per-KO ancient read counts from DART EM output
        func_path <- file.path(kegg_dir, s, "emi.functional.tsv")
        ko_ancient <- if (file.exists(func_path)) {
            f <- fread(func_path, showProgress = FALSE)[level == "group"]
            f[, .(ko = function_id, ancient_cov = coverage_mean * n_ancient / n_reads)]
        } else {
            data.table(ko = character(0), ancient_cov = numeric(0))
        }

        # Explode enzyme_hits_in_module → unique KOs per module
        mods_ne <- d[!is.na(enzyme_hits_in_module) & enzyme_hits_in_module != ""]
        if (nrow(mods_ne) > 0) {
            ko_map <- mods_ne[, .(ko = unique(trimws(unlist(strsplit(enzyme_hits_in_module, ","))))),
                               by = module]
            ko_map <- ko_map[ko != ""]
            ko_map <- merge(ko_map, ko_ancient, by = "ko", all.x = TRUE)
            mod_anc <- ko_map[, .(ancient_coverage = mean(ancient_cov, na.rm = TRUE)), by = module]
        } else {
            mod_anc <- data.table(module = character(0), ancient_coverage = numeric(0))
        }

        d <- merge(
            d[, .(module, module_name, module_class, module_category, module_subcategory,
                  stepwise_module_completeness, stepwise_module_is_complete, avg_coverage)],
            mod_anc, by = "module", all.x = TRUE
        )
        d[, sample := s]
        d
    }), fill = TRUE)

    merge(dt, sample_meta[, .(sample, label)], by = "sample")
}


#' Compute per-module DART damage stats
#'
#' For each (sample, module), computes mean_p_damaged by joining
#' anvio_modules.txt gene lists with emi.protein.slim.tsv.gz damage scores.
#'
#' @param threshold p_protein_damaged cutoff to call a protein damaged
#' @return data.table: sample, module, module_name, module_class, mean_p_damaged, label
compute_kegg_damage <- function(kegg_dir, sample_meta, threshold = 0.7) {
    samples <- intersect(discover_samples(kegg_dir), sample_meta$sample)

    dt <- rbindlist(lapply(samples, function(s) {
        mod_path <- file.path(kegg_dir, s, "anvio_modules.txt")
        dmg_path <- file.path(kegg_dir, s, "emi.protein.slim.tsv.gz")
        if (!file.exists(mod_path) || !file.exists(dmg_path)) return(NULL)

        mods <- fread(mod_path, showProgress = FALSE)
        cov_col <- grep("avg_coverage$", names(mods), value = TRUE)[1]
        if (cov_col != "avg_coverage") setnames(mods, cov_col, "avg_coverage")

        # Explode gene_caller_ids_in_module → one row per protein
        mods_ne <- mods[!is.na(gene_caller_ids_in_module) & gene_caller_ids_in_module != ""]
        if (nrow(mods_ne) == 0) return(NULL)

        gene_map <- mods_ne[, .(
            gene_caller_id = trimws(unlist(strsplit(gene_caller_ids_in_module, ",")))
        ), by = .(module, module_name, module_class, avg_coverage)]
        gene_map <- gene_map[gene_caller_id != ""]

        prot <- fread(dmg_path, showProgress = FALSE)
        prot <- prot[pass_mapping_filter == 1]

        joined <- merge(gene_map,
                        prot[, .(protein_id, p_protein_damaged)],
                        by.x = "gene_caller_id", by.y = "protein_id")

        joined[, .(
            mean_p_damaged = mean(p_protein_damaged, na.rm = TRUE),
            avg_coverage   = unique(avg_coverage)[1]
        ), by = .(module, module_name, module_class)][, sample := s]
    }), fill = TRUE)

    merge(dt, sample_meta[, .(sample, label)], by = "sample")
}


#' Load CAZy family abundance for all samples
#'
#' Replaces dbcan.group-abundances-agg.tsv.gz.
#' Returns: label (MD5), group (family), mean (coverage_mean).
load_cazy_abundances <- function(cazy_dir, sample_meta) {
    samples <- intersect(discover_samples(cazy_dir), sample_meta$sample)

    dt <- rbindlist(lapply(samples, function(s) {
        path <- file.path(cazy_dir, s, "cazy_emi.functional.tsv")
        if (!file.exists(path)) { warning("Missing: ", path); return(NULL) }
        d <- fread(path, showProgress = FALSE)[level == "group"]
        d[, .(sample = s, group = function_id, mean = coverage_mean)]
    }), fill = TRUE)

    merge(dt, sample_meta[, .(sample, label)], by = "sample")
}


#' Compute per-family DART damage stats
#'
#' @return data.table: sample, group (family), mean, mean_p_damaged, label
compute_cazy_damage <- function(cazy_dir, sample_meta, threshold = 0.7) {
    samples <- intersect(discover_samples(cazy_dir), sample_meta$sample)

    dt <- rbindlist(lapply(samples, function(s) {
        abund_path <- file.path(cazy_dir, s, "cazy_emi.functional.tsv")
        dmg_path   <- file.path(cazy_dir, s, "cazy_damage.slim.tsv.gz")
        if (!file.exists(abund_path) || !file.exists(dmg_path)) return(NULL)

        abund <- fread(abund_path, showProgress = FALSE)[level == "group",
                       .(group = function_id, mean = coverage_mean)]

        prot <- fread(dmg_path, showProgress = FALSE)[pass_mapping_filter == 1]
        prot[, group := sub("^.*\\|", "", protein_id)]

        prot_agg <- prot[, .(mean_p_damaged = mean(p_protein_damaged, na.rm = TRUE)),
                          by = group]

        merged <- merge(abund, prot_agg, by = "group", all.x = TRUE)
        merged[, sample := s]
        merged
    }), fill = TRUE)

    merge(dt, sample_meta[, .(sample, label)], by = "sample")
}


#' Load per-protein module assignments with DART damage scores
#'
#' Returns protein-level data (not aggregated): sample, label, module,
#' module_name, module_class, gene_caller_id, p_protein_damaged.
#' Used as input for enrichment tests.
load_protein_module_data <- function(kegg_dir, sample_meta) {
    samples <- intersect(discover_samples(kegg_dir), sample_meta$sample)

    dt <- rbindlist(lapply(samples, function(s) {
        mod_path <- file.path(kegg_dir, s, "anvio_modules.txt")
        dmg_path <- file.path(kegg_dir, s, "emi.protein.slim.tsv.gz")
        if (!file.exists(mod_path) || !file.exists(dmg_path)) return(NULL)

        mods <- fread(mod_path, showProgress = FALSE)
        mods_ne <- mods[!is.na(gene_caller_ids_in_module) & gene_caller_ids_in_module != ""]
        if (nrow(mods_ne) == 0) return(NULL)

        gene_map <- mods_ne[, .(
            gene_caller_id = trimws(unlist(strsplit(gene_caller_ids_in_module, ",")))
        ), by = .(module, module_name, module_class)]
        gene_map <- gene_map[gene_caller_id != ""]

        prot <- fread(dmg_path, showProgress = FALSE)[pass_mapping_filter == 1]

        merged <- merge(gene_map,
                        prot[, .(protein_id, p_protein_damaged)],
                        by.x = "gene_caller_id", by.y = "protein_id")
        merged[, sample := s]
        merged
    }), fill = TRUE)

    merge(dt, sample_meta[, .(sample, label)], by = "sample")
}


# =============================================================================
# build_functional_info  (mirrors get_functional_info from original)
# =============================================================================

#' Build functional info list from DART data tables
#'
#' Produces the same list structure as get_functional_info():
#'   list(polysac_deg, central_carbon_meta, methane_metabolism, signature_modules)
#'
#' @param kegg_dt  data.table from load_kegg_modules() (or damage-filtered subset)
#' @param cazy_dt  data.table from load_cazy_abundances() (or damage-filtered subset)
build_functional_info <- function(kegg_dt, cazy_dt) {
    # Polysaccharide degradation (CAZy)
    polysac_deg <- rbind(
        cazy_dt[group %in% CELLULOSE_FAMILIES][, process := "Cellulose degradation"],
        cazy_dt[group %in% HEMICELLULOSE_FAMILIES][, process := "Hemicellulose xylan degradation"]
    )

    # Fermentation / central carbon (Woodcroft2018)
    woodcroft <- kegg_dt[module_class == "Woodcroft2018"]
    fermentation <- woodcroft[module_name %in% FERMENTATION_MODULES]
    central_carbon_meta <- fermentation[,
        .(annotation = module, abundance = avg_coverage, label, process = module_name)]

    # Methane metabolism
    pathway <- kegg_dt[module_class == "Pathway modules"]
    borrel  <- copy(kegg_dt[module_class == "Borrell2022"])
    borrel[, module_class := "Borrel2023"]

    methane_metabolism <- rbind(
        pathway[module_name %in% METHANE_MODULES_PATHWAY][,
            .(annotation = module, label, process = module_name,
              abundance = avg_coverage, module_class)],
        woodcroft[module_name %in% METHANE_MODULES_WOODCROFT][,
            .(annotation = module, label, process = module_name,
              abundance = avg_coverage, module_class)],
        pathway[module_name == "Methane oxidation, methanotroph, methane => formaldehyde"][,
            .(annotation = module, label, process = module_name,
              abundance = avg_coverage, module_class)],
        borrel[, .(annotation = module, label, process = module_name,
                   abundance = avg_coverage, module_class)],
        fill = TRUE
    )

    # Signature modules (contains Acetogen)
    signature_modules <- kegg_dt[module_class == "Signature modules"]

    list(
        polysac_deg         = polysac_deg,
        central_carbon_meta = central_carbon_meta,
        methane_metabolism  = methane_metabolism,
        signature_modules   = signature_modules
    )
}


# =============================================================================
# process_functional_data  (identical to 07--functional-profiling.R)
# =============================================================================

process_functional_data <- function(func_all, func_dmg, kapk_cdata, sample_order = NULL) {
    # Signature paths (Acetogen)
    signature_paths <- as_tibble(func_all$signature_modules) %>%
        mutate(type = "All") %>%
        bind_rows(as_tibble(func_dmg$signature_modules) %>% mutate(type = "Damaged")) %>%
        filter(module_category == "Module set",
               module_subcategory == "Metabolic capacity") %>%
        select(label, abundance = avg_coverage, type, process = module_name) %>%
        filter(process == "Acetogen") %>%
        mutate(prank = 2, class = "Central carbon metabolism")

    # Methane metabolism
    methane_metabolism <- as_tibble(func_all$methane_metabolism) %>%
        mutate(type = "All") %>%
        bind_rows(as_tibble(func_dmg$methane_metabolism) %>% mutate(type = "Damaged")) %>%
        filter(
            module_class == "Borrel2023" |
            process == "Methane oxidation, methanotroph, methane => formaldehyde" |
            process == "hydrogenotrophic methanogenesis" |
            process == "Methanogenesis, methanol => methane" |
            process == "Methanogenesis, methylamine/dimethylamine/trimethylamine => methane"
        ) %>%
        select(label, abundance, type, process) %>%
        mutate(
            process = case_when(
                process == "Methane oxidation, methanotroph, methane => formaldehyde"
                    ~ "Methanotrophy",
                process == "Methanogenesis, methanol => methane"
                    ~ "Methylotrophic methanogenesis (methanol)",
                process == "Methanogenesis, methylamine/dimethylamine/trimethylamine => methane"
                    ~ "Methylotrophic methanogenesis (methylamine)",
                TRUE ~ paste0(toupper(substr(process, 1, 1)), substr(process, 2, nchar(process)))
            ),
            prank = 3,
            class = "Methane metabolism"
        )

    # Central carbon metabolism
    central_carbon_meta <- as_tibble(func_all$central_carbon_meta) %>%
        mutate(type = "All") %>%
        bind_rows(as_tibble(func_dmg$central_carbon_meta) %>% mutate(type = "Damaged")) %>%
        select(label, abundance, type, process) %>%
        mutate(
            process = paste0(toupper(substr(process, 1, 1)), substr(process, 2, nchar(process))),
            prank   = 2,
            class   = "Central carbon metabolism"
        ) %>%
        filter(process != "Acetatogenesis")

    # Polysaccharide degradation
    polysac_deg <- as_tibble(func_all$polysac_deg) %>%
        mutate(type = "All") %>%
        bind_rows(as_tibble(func_dmg$polysac_deg) %>% mutate(type = "Damaged")) %>%
        select(label, abundance = mean, type, process) %>%
        filter(!grepl("Xylose", process)) %>%
        mutate(prank = 1, class = "Polysaccharide degradation")

    # Combine
    combined <- bind_rows(polysac_deg, central_carbon_meta, methane_metabolism, signature_paths)

    data_all <- combined %>%
        filter(type == "All") %>%
        droplevels() %>%
        complete(label, nesting(process, class), type, fill = list(abundance = NA)) %>%
        inner_join(kapk_cdata, by = "label") %>%
        select(label = short_label, member_unit, process, class, abundance, type, site_rnk) %>%
        ungroup() %>%
        mutate(
            label       = if (!is.null(sample_order))
                              factor(label, levels = sample_order)
                          else
                              fct_reorder(label, site_rnk),
            member_unit = fct_relevel(member_unit, c("B1", "B2", "B3")),
            class       = fct_relevel(class, c("Polysaccharide degradation",
                                               "Central carbon metabolism",
                                               "Methane metabolism")),
            process     = fct_relevel(process, DATA_ALL_ORDER)
        )

    data_dmg <- combined %>%
        filter(type == "Damaged") %>%
        droplevels() %>%
        inner_join(kapk_cdata, by = "label") %>%
        mutate(
            short_label  = if (!is.null(sample_order))
                               factor(short_label, levels = sample_order)
                           else
                               fct_reorder(short_label, site_rnk),
            member_unit  = fct_relevel(member_unit, c("B1", "B2", "B3")),
            class        = fct_relevel(class, c("Polysaccharide degradation",
                                                "Central carbon metabolism",
                                                "Methane metabolism")),
            process      = fct_relevel(process, DATA_ALL_ORDER)
        ) %>%
        select(label = short_label, member_unit, process, class, abundance, type)

    list(all = data_all, dmg = data_dmg)
}


# =============================================================================
# plot_functional_heatmap  (identical to 07--functional-profiling.R)
# =============================================================================

plot_functional_heatmap <- function(data_all, data_dmg) {
    # Capture factor levels before any reassignment.
    label_lvls   <- levels(data_all$label)
    process_lvls <- levels(data_all$process)
    class_lvls   <- levels(data_all$class)
    unit_lvls    <- levels(data_all$member_unit)

    restore_levels <- function(df) {
        df %>% mutate(
            label       = factor(as.character(label),       levels = label_lvls),
            process     = factor(as.character(process),     levels = process_lvls),
            class       = factor(as.character(class),       levels = class_lvls),
            member_unit = factor(as.character(member_unit), levels = unit_lvls)
        )
    }

    # Aggregate to one row per (label, process): multiple CAZy families share a process label.
    data_all <- data_all %>%
        group_by(label, member_unit, process, class, site_rnk) %>%
        summarize(
            abundance = if (all(is.na(abundance))) NA_real_ else mean(abundance, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        restore_levels()

    data_dmg <- data_dmg %>%
        group_by(label, member_unit, process, class) %>%
        summarize(
            abundance = if (all(is.na(abundance))) NA_real_ else mean(abundance, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        restore_levels()

    # Only annotate damage where coverage is visible (>= scale minimum of 1x).
    # Use character comparison to avoid any factor-level mismatch in the join.
    visible_keys <- data_all %>%
        filter(!is.na(abundance) & abundance >= 1) %>%
        transmute(lbl = as.character(label), proc = as.character(process))
    data_dmg <- data_dmg %>%
        mutate(lbl = as.character(label), proc = as.character(process)) %>%
        semi_join(visible_keys, by = c("lbl", "proc")) %>%
        select(-lbl, -proc)

    p <- ggplot(data_all, aes(x = label, y = process, fill = abundance)) +
        geom_tile(color = "grey60", linewidth = 0.15) +
        geom_tile(
            data      = data_dmg,
            aes(x = label, y = process),
            fill      = NA,
            color     = "#323232",
            linewidth = 0.9
        ) +
        facet_grid(class ~ member_unit, scales = "free", space = "free", switch = "x",
                   labeller = labeller(class = label_wrap_gen(12))) +
        theme_bw(base_size = 9) +
        theme(
            axis.text.x            = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
            axis.text.y            = element_text(size = 8),
            axis.ticks             = element_blank(),
            panel.border           = element_blank(),
            panel.spacing          = unit(1.5, "mm"),
            legend.position        = "top",
            legend.justification   = "left",
            legend.direction       = "horizontal",
            legend.title           = element_text(size = 8),
            legend.text            = element_text(size = 8),
            legend.key.height      = unit(3, "mm"),
            legend.key.width       = unit(8, "mm"),
            legend.background      = element_blank(),
            legend.margin          = margin(t = 0, b = 1, l = 0, r = 0),
            legend.box.spacing     = unit(0, "mm"),
            strip.text.x           = element_text(size = 9, face = "bold"),
            strip.text.y.right     = element_text(size = 8, angle = 0, hjust = 0),
            strip.background.x     = element_blank(),
            strip.background.y     = element_blank(),
            strip.placement        = "outside",
            plot.margin            = margin(t = 2, r = 8, b = 2, l = 2, unit = "mm")
        ) +
        scale_fill_gradientn(
            colors   = c("#D7E7F0", "#F7E2B5", "#F1AF82", "#E67B70", "#C94A6B", "#9B2471"),
            breaks   = c(1, 2, 4, 8, 16, 32),
            labels   = c("1", "2", "4", "8", "16", "32"),
            limits   = c(1, 32),
            guide    = guide_colorbar(title.position = "top", title.hjust = 0,
                                      barwidth = unit(32, "mm"), barheight = unit(3, "mm")),
            na.value = "#ffffff",
            trans    = "log2",
            name     = "Average coverage (log2)"
        ) +
        coord_fixed(ratio = 1) +
        xlab("") + ylab("")

    # Build gtable and color the bottom (switched) column strips by depth unit.
    # In ggplot2 4.x the strip grob structure:
    #   gtable$grobs[[strip_idx]] (gTree) -> grobs[[1]] -> children[[1]] (zeroGrob background)
    g <- ggplotGrob(p)
    strip_idx <- which(grepl("^strip-b-", g$layout$name))
    strip_idx <- strip_idx[order(g$layout$l[strip_idx])]
    unit_order <- levels(data_all$member_unit)
    for (i in seq_along(strip_idx)) {
        fill_col <- alpha(UNIT_COLORS[unit_order[i]], 0.35)
        rect_grob <- grid::rectGrob(gp = grid::gpar(fill = fill_col, col = NA))
        g$grobs[[strip_idx[i]]]$grobs[[1]]$children[[1]] <- rect_grob
    }

    # Position the top guide-box just above the panels, shifted right of B1.
    guide_idx_top <- which(g$layout$name == "guide-box-top")
    if (length(guide_idx_top) > 0) {
        panel_rows <- g$layout$t[grepl("^panel-", g$layout$name)]
        panel_top  <- min(panel_rows)
        guide_row  <- g$layout$t[guide_idx_top[1]]
        # Zero extra padding rows between guide-box and panels
        for (r in seq(guide_row + 1L, panel_top - 1L)) {
            g$heights[r] <- unit(0, "mm")
        }
        g$heights[guide_row] <- unit(12, "mm")
        # Shift legend right: start one column right of B1 left edge
        strip_l <- min(g$layout$l[strip_idx])
        strip_r <- max(g$layout$r[strip_idx])
        g$layout$l[guide_idx_top[1]] <- strip_l + 1L
        g$layout$r[guide_idx_top[1]] <- strip_r
    }

    # Zero out the bottom guide-box row (empty now that legend is at top)
    # and all rows below the B1/B2/B3 strips to remove white space.
    strip_rows   <- g$layout$t[strip_idx]
    strip_bottom <- max(strip_rows)
    guide_idx_bot <- which(g$layout$name == "guide-box-bottom")
    n_rows <- length(g$heights)
    for (r in seq(strip_bottom + 1L, n_rows)) {
        g$heights[r] <- unit(0, "mm")
    }
    g
}


# =============================================================================
# Damage Distribution Diagnostic
# =============================================================================

#' Load p_protein_damaged for ALL proteins (not just module-assigned ones)
load_all_protein_damage <- function(kegg_dir, sample_meta) {
    samples <- intersect(discover_samples(kegg_dir), sample_meta$sample)

    dt <- rbindlist(lapply(samples, function(s) {
        dmg_path <- file.path(kegg_dir, s, "emi.protein.slim.tsv.gz")
        if (!file.exists(dmg_path)) return(NULL)
        d <- fread(dmg_path, showProgress = FALSE)[pass_mapping_filter == 1]
        d[, sample := s]
        d[, .(sample, protein_id, p_protein_damaged)]
    }), fill = TRUE)

    merge(dt, sample_meta[, .(sample, label, short_label, member_unit, site_rnk)],
          by = "sample")
}


#' 3-panel damage distribution diagnostic
#'
#' Panel 1: global density (all proteins pooled) — look for bimodality
#' Panel 2: per-sample densities overlaid — check consistency across samples
#' Panel 3: % proteins with p_damaged >= 0.5 per sample — ancient fraction
plot_damage_distribution <- function(prot_all_dt) {
    dt <- as.data.table(prot_all_dt)

    p_global <- ggplot(dt, aes(x = p_protein_damaged)) +
        geom_density(fill = "#D7E7F0", color = "#4A90C4", linewidth = 0.8) +
        geom_vline(xintercept = MODERN_THRESHOLD,  linetype = "dashed", color = "#4A90C4") +
        geom_vline(xintercept = DAMAGED_THRESHOLD, linetype = "dashed", color = "#C94A6B") +
        annotate("text", x = MODERN_THRESHOLD  - 0.02, y = Inf, label = "no signal",
                 hjust = 1, vjust = 1.5, size = 3, color = "#4A90C4") +
        annotate("text", x = DAMAGED_THRESHOLD + 0.02, y = Inf, label = "damaged",
                 hjust = 0, vjust = 1.5, size = 3, color = "#C94A6B") +
        labs(x = "p_protein_damaged", y = "Density",
             title = "Global distribution — all samples pooled") +
        theme_bw() +
        theme(strip.background = element_blank())

    p_samples <- ggplot(
        dt[, .(short_label, member_unit, p_protein_damaged)],
        aes(x = p_protein_damaged, color = member_unit, group = short_label)
    ) +
        geom_density(linewidth = 0.4, alpha = 0.8) +
        geom_vline(xintercept = MODERN_THRESHOLD,  linetype = "dashed", color = "#4A90C4") +
        geom_vline(xintercept = DAMAGED_THRESHOLD, linetype = "dashed", color = "#C94A6B") +
        annotate("text", x = MODERN_THRESHOLD  - 0.02, y = Inf, label = "no signal",
                 hjust = 1, vjust = 1.5, size = 3, color = "#4A90C4") +
        annotate("text", x = DAMAGED_THRESHOLD + 0.02, y = Inf, label = "damaged",
                 hjust = 0, vjust = 1.5, size = 3, color = "#C94A6B") +
        scale_color_brewer(palette = "Set1", name = "Site") +
        labs(x = "p_protein_damaged", y = "Density",
             title = "Per-sample distributions (one line per sample)") +
        theme_bw() +
        theme(strip.background = element_blank(), legend.position = "right")

    # Stacked proportions: damaged / ambiguous / no-signal per sample
    sample_class <- dt[, .(
        damaged  = mean(p_protein_damaged >= DAMAGED_THRESHOLD) * 100,
        ambig    = mean(p_protein_damaged >  MODERN_THRESHOLD &
                        p_protein_damaged <  DAMAGED_THRESHOLD) * 100,
        no_signal = mean(p_protein_damaged <= MODERN_THRESHOLD) * 100
    ), by = .(short_label, member_unit, site_rnk)]
    sample_class[, short_label := fct_reorder(short_label, site_rnk)]
    sample_class_long <- melt(sample_class,
        id.vars      = c("short_label", "member_unit", "site_rnk"),
        measure.vars = c("damaged", "ambig", "no_signal"),
        variable.name = "class", value.name = "pct")
    sample_class_long[, class := factor(class, levels = c("no_signal", "ambig", "damaged"))]

    p_frac <- ggplot(sample_class_long,
                     aes(x = short_label, y = pct, fill = class)) +
        geom_col() +
        scale_fill_manual(
            values = c(damaged = "#C94A6B", ambig = "#F1AF82", no_signal = "#4A90C4"),
            labels = c(damaged   = paste0("Damaged (\u2265", DAMAGED_THRESHOLD, ")"),
                       ambig     = "Ambiguous",
                       no_signal = paste0("No signal (\u2264", MODERN_THRESHOLD, ")")),
            name = NULL
        ) +
        labs(x = "", y = "% proteins",
             title = "Protein classification per sample (DART 3-class)") +
        coord_flip() +
        theme_bw() +
        theme(strip.background = element_blank(), legend.position = "bottom")

    list(global = p_global, per_sample = p_samples, frac = p_frac)
}


# =============================================================================
# Module Enrichment Analysis (cross-sample Spearman correlation)
#
# Strategy: for each KEGG module, test whether its avg_coverage across samples
# is positively correlated with the sample's global mean_damage score.
# A significant positive ρ indicates the module is more abundant in samples
# with stronger ancient DNA signal — i.e., enriched in the ancient fraction.
# This avoids the within-sample bimodal issue (n_modern=0 in all samples).
# =============================================================================

#' Load per-sample global damage score from pathway_damage_stats.tsv
#'
#' @return data.table: sample, label, mean_damage
load_sample_damage_stats <- function(kegg_dir, sample_meta) {
    samples <- intersect(discover_samples(kegg_dir), sample_meta$sample)

    dt <- rbindlist(lapply(samples, function(s) {
        path <- file.path(kegg_dir, s, "pathway_damage_stats.tsv")
        if (!file.exists(path)) { warning("Missing: ", path); return(NULL) }
        d <- fread(path, showProgress = FALSE)
        if ("mean_damage" %in% names(d)) {
            data.table(sample = s, mean_damage = d$mean_damage[1])
        } else if ("mean_score" %in% names(d)) {
            data.table(sample = s, mean_damage = d$mean_score[1])
        } else NULL
    }), fill = TRUE)

    merge(dt, sample_meta[, .(sample, label)], by = "sample")
}


#' Cross-sample Spearman correlation: module coverage vs global damage score
#'
#' For each module, tests whether avg_coverage (per sample) is positively
#' correlated with that sample's mean_damage. A significant positive ρ means
#' the module is more abundant in samples with higher ancient DNA signal.
#'
#' @param kegg_dt       data.table from load_kegg_modules()
#' @param damage_stats  data.table from load_sample_damage_stats()
#' @param module_filter optional character vector of module_names to restrict test
#' @return data.table: module, module_name, module_class, n_samples, rho, p.value, q.value
test_enrichment_spearman <- function(kegg_dt, damage_stats, module_filter = NULL) {
    dt <- merge(as.data.table(kegg_dt),
                damage_stats[, .(sample, mean_damage)],
                by = "sample")

    # Deduplicate: same module may appear in multiple module_classes (e.g. Borrell + Woodcroft).
    # Keep the row with the highest avg_coverage per (sample, module_name).
    dt <- dt[dt[, .I[which.max(avg_coverage)], by = .(sample, module_name)]$V1]

    modules <- unique(dt[, .(module, module_name, module_class)])
    if (!is.null(module_filter)) modules <- modules[module_name %in% module_filter]

    results <- rbindlist(lapply(seq_len(nrow(modules)), function(i) {
        mdt <- dt[module_name == modules$module_name[i]]
        if (nrow(mdt) < 5) return(NULL)

        st <- tryCatch(
            cor.test(mdt$avg_coverage, mdt$mean_damage, method = "spearman", exact = FALSE),
            error = function(e) NULL
        )
        if (is.null(st)) return(NULL)

        data.table(
            module       = modules$module[i],
            module_name  = modules$module_name[i],
            module_class = modules$module_class[i],
            n_samples    = nrow(mdt),
            rho          = as.numeric(st$estimate),
            p.value      = st$p.value
        )
    }), fill = TRUE)

    if (nrow(results) > 0) results[, q.value := p.adjust(p.value, method = "BH")]
    results
}


#' Load module-specific ancient read fraction per sample
#'
#' Uses emi.functional.tsv (n_ancient, n_reads per KO) joined with enzyme_hits_in_module
#' from anvio_modules.txt. For each (sample, module), computes:
#'   module_ancient_frac = sum(n_ancient for KOs in module) / sum(n_reads for KOs in module)
#'
#' This is the fraction of reads mapping to a module's KOs that are ancient.
#' It corrects for living organisms: if living methanogens contribute heavily,
#' n_reads is high but n_ancient stays low → module_ancient_frac stays low.
#' Modules truly enriched in ancient organisms will have high module_ancient_frac
#' and this should be positively correlated with global mean_damage across samples.
#'
#' @return data.table: sample, label, module, module_name, module_class,
#'         sum_n_ancient, sum_n_reads, module_ancient_frac
load_module_ancient_coverage <- function(kegg_dir, sample_meta) {
    samples <- intersect(discover_samples(kegg_dir), sample_meta$sample)

    dt <- rbindlist(lapply(samples, function(s) {
        mod_path <- file.path(kegg_dir, s, "anvio_modules.txt")
        emi_path <- file.path(kegg_dir, s, "emi.functional.tsv")
        if (!file.exists(mod_path) || !file.exists(emi_path)) return(NULL)

        mods <- fread(mod_path, showProgress = FALSE)
        mods_ne <- mods[!is.na(enzyme_hits_in_module) & enzyme_hits_in_module != ""]
        if (nrow(mods_ne) == 0) return(NULL)

        # Explode enzyme_hits_in_module → unique KOs per module (K-numbers may repeat)
        ko_map <- mods_ne[, .(
            ko = unique(trimws(unlist(strsplit(enzyme_hits_in_module, ","))))
        ), by = .(module, module_name, module_class)]
        ko_map <- ko_map[grepl("^K[0-9]{5}$", ko)]

        # Load emi.functional.tsv: KO-level read counts
        emi <- fread(emi_path, showProgress = FALSE)
        emi <- emi[level == "group", .(ko = function_id, n_ancient, n_reads)]

        # Join and aggregate: sum reads across unique KOs per module
        joined <- merge(ko_map, emi, by = "ko", all.x = TRUE)
        joined[is.na(n_ancient), n_ancient := 0]
        joined[is.na(n_reads),   n_reads   := 0]

        agg <- joined[, .(
            sum_n_ancient = sum(n_ancient),
            sum_n_reads   = sum(n_reads)
        ), by = .(module, module_name, module_class)]
        # Module-specific ancient fraction: fraction of module reads that are ancient
        agg[sum_n_reads > 0, module_ancient_frac := sum_n_ancient / sum_n_reads]
        agg[sum_n_reads == 0, module_ancient_frac := NA_real_]
        agg[, sample := s]
        agg
    }), fill = TRUE)

    merge(dt, sample_meta[, .(sample, label)], by = "sample")
}


#' Cross-sample Spearman correlation: ancient coverage fraction vs global damage score
#'
#' For each module, tests whether its proportion of sample-wide ancient reads
#' (ancient_cov_frac from load_module_ancient_coverage) correlates with mean_damage.
#' Unlike test_enrichment_spearman(), this uses only ancient reads, removing the
#' signal from living organisms that also contribute to avg_coverage.
#'
#' @param ancient_dt    data.table from load_module_ancient_coverage()
#' @param damage_stats  data.table from load_sample_damage_stats()
#' @param module_filter optional character vector of module_names to restrict test
#' @return data.table: module, module_name, module_class, n_samples, rho, p.value, q.value
test_enrichment_spearman_ancient <- function(ancient_dt, damage_stats, module_filter = NULL) {
    dt <- merge(as.data.table(ancient_dt),
                damage_stats[, .(sample, mean_damage)],
                by = "sample")

    modules <- unique(dt[, .(module, module_name, module_class)])
    if (!is.null(module_filter)) modules <- modules[module_name %in% module_filter]

    results <- rbindlist(lapply(seq_len(nrow(modules)), function(i) {
        # Use module ID (not name) to avoid duplicates across module_classes
        mdt <- dt[module == modules$module[i]]
        if (nrow(mdt) < 5) return(NULL)

        mdt <- mdt[!is.na(module_ancient_frac)]
        if (nrow(mdt) < 5) return(NULL)

        st <- tryCatch(
            cor.test(mdt$module_ancient_frac, mdt$mean_damage,
                     method = "spearman", exact = FALSE),
            error = function(e) NULL
        )
        if (is.null(st)) return(NULL)

        data.table(
            module       = modules$module[i],
            module_name  = modules$module_name[i],
            module_class = modules$module_class[i],
            n_samples    = nrow(mdt),
            rho          = as.numeric(st$estimate),
            p.value      = st$p.value
        )
    }), fill = TRUE)

    if (nrow(results) > 0) results[, q.value := p.adjust(p.value, method = "BH")]
    results
}


#' Volcano plot: Spearman ρ vs -log10(FDR q), highlighting figure modules
plot_enrichment_volcano <- function(enr_dt, figure_module_names = NULL, title = "") {
    dt <- as.data.table(enr_dt)
    dt[, sig       := q.value < 0.1 & rho > 0]
    dt[, is_figure := module_name %in% figure_module_names]
    dt[, log10q    := -log10(q.value + 1e-300)]

    label_dt <- dt[sig == TRUE | is_figure == TRUE]

    ggplot(dt, aes(x = rho, y = log10q)) +
        geom_point(aes(color = sig, shape = is_figure), alpha = 0.7, size = 2) +
        geom_text(data = label_dt, aes(label = module_name),
                  size = 2.2, hjust = -0.1, check_overlap = TRUE) +
        geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "gray50") +
        geom_vline(xintercept = 0,           linetype = "dashed", color = "gray50") +
        scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "#C94A6B"),
                           labels = c("n.s.", "FDR < 0.1")) +
        scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18),
                           labels = c("Other", "Figure module")) +
        labs(x = "Spearman \u03c1 (coverage ~ mean damage score)",
             y = "-log10(FDR q-value)", title = title,
             color = NULL, shape = NULL) +
        theme_bw() +
        theme(legend.position = "bottom", strip.background = element_blank())
}


#' Lollipop of Spearman ρ per module (top modules + figure modules)
plot_enrichment_lollipop <- function(enr_dt, figure_module_names = NULL, top_n = 40) {
    dt <- as.data.table(enr_dt)
    dt[, is_figure := module_name %in% figure_module_names]

    top_mods <- dt[order(q.value)][seq_len(min(top_n, .N)), module_name]
    dt_plot  <- dt[module_name %in% top_mods | is_figure == TRUE]
    dt_plot[, module_label := fct_reorder(module_name, rho)]

    ggplot(dt_plot, aes(x = rho, y = module_label)) +
        geom_segment(aes(xend = 0, yend = module_label), color = "gray80") +
        geom_point(aes(size  = n_samples,
                       color = -log10(q.value + 1e-10),
                       shape = is_figure)) +
        scale_color_gradientn(colors = c("gray70", "#F1AF82", "#C94A6B", "#9B2471"),
                              name   = "-log10(q)") +
        scale_size_continuous(range = c(1, 6), name = "N samples") +
        scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18),
                           labels = c("Other", "Figure module"), name = NULL) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        labs(x = "Spearman \u03c1 (coverage ~ mean damage score)", y = "") +
        theme_bw() +
        theme(strip.background = element_blank(), legend.position = "right")
}


# =============================================================================
# Key Enzyme Analysis
# =============================================================================

#' Load key enzyme (KO-level) ancient read data per sample
#'
#' @return data.table: sample, label, short_label, member_unit, site_rnk,
#'         ko, gene, pathway, pathway_order, n_reads, n_ancient,
#'         mean_posterior, ancient_pct
load_key_enzyme_data <- function(kegg_dir, sample_meta,
                                 ko_table = KEY_KOS_TABLE) {
    samples <- intersect(discover_samples(kegg_dir), sample_meta$sample)

    dt <- rbindlist(lapply(samples, function(s) {
        path <- file.path(kegg_dir, s, "emi.functional.tsv")
        if (!file.exists(path)) return(NULL)
        d <- fread(path, showProgress = FALSE)
        d <- d[level == "group" & function_id %in% ko_table$ko]
        if (nrow(d) == 0) return(NULL)
        d[, .(sample = s, ko = function_id, n_reads, n_ancient, mean_posterior)]
    }), fill = TRUE)

    dt <- merge(dt, ko_table, by = "ko")
    dt <- merge(dt, sample_meta[, .(sample, label, short_label, member_unit, site_rnk)],
                by = "sample")
    dt[n_reads > 0, ancient_pct := n_ancient / n_reads * 100]
    dt
}


#' Dot plot of key enzyme ancient signal per sample
#'
#' Mirrors the layout of plot_functional_heatmap:
#'   rows    = KOs ordered by pathway (matching DATA_ALL_ORDER) within class facets
#'   columns = samples ordered by site_rnk within member_unit (B1/B2/B3) facets
#'   color   = ancient_pct (% of reads classified as ancient by DART)
#'   size    = log10(n_reads)
plot_key_enzyme_dotplot <- function(enzyme_dt) {
    dt <- as.data.table(enzyme_dt)

    # Row factor: genes ordered bottom→top matching heatmap top-to-bottom.
    # Central carbon pathways come first (top of dotplot), Methane last (bottom).
    # Within each pathway, preserve KEY_KOS_TABLE row order (top gene = last in table).
    ko_order <- unique(dt[order(class_order, pathway_order, -ko_idx),
                          paste0(gene, " (", ko, ")")])
    dt[, ko_label := factor(paste0(gene, " (", ko, ")"), levels = ko_order)]

    # Column factor: samples ordered by site_rnk within member_unit (B1 → B2 → B3)
    dt[, sample_label := fct_reorder(short_label, site_rnk)]
    dt[, member_unit  := fct_relevel(member_unit, c("B1", "B2", "B3"))]

    dt[, log_reads := log10(n_reads + 1)]

    # Pathway factor: hardcoded to match heatmap top-to-bottom (level 1 = top facet).
    dt[, pathway := factor(pathway, levels = DOTPLOT_PATHWAY_ORDER)]

    ggplot(dt, aes(x = sample_label, y = ko_label)) +
        geom_point(aes(size = log_reads, color = ancient_pct), alpha = 0.85) +
        scale_color_gradientn(
            colors   = c("#D7E7F0", "#F7E2B5", "#F1AF82", "#C94A6B", "#9B2471"),
            limits   = c(0, 100),
            name     = "Ancient reads (%)",
            guide    = guide_colorbar(barwidth = 8, barheight = 0.5,
                                      title.position = "top")
        ) +
        scale_size_continuous(
            range  = c(0.5, 5),
            name   = "log\u2081\u2080(reads)",
            breaks = c(1, 2, 3),
            labels = c("10", "100", "1000")
        ) +
        facet_grid(pathway ~ member_unit, scales = "free", space = "free") +
        labs(x = "", y = "") +
        theme_bw(base_size = 7) +
        theme(
            axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                            size = 6),
            axis.text.y      = element_text(size = 6),
            strip.text.y     = element_text(size = 6, angle = 0, hjust = 0),
            strip.text.x     = element_text(size = 7),
            legend.text      = element_text(size = 6),
            legend.title     = element_text(size = 6),
            strip.background = element_blank(),
            legend.position  = "top",
            panel.grid.major = element_line(linewidth = 0.2),
            panel.spacing    = unit(0.2, "lines")
        )
}


# =============================================================================
# Main
# =============================================================================

main <- function(
    kegg_dir     = "./results/functional_agp/kegg",
    cazy_dir     = "./results/functional_agp/cazy",
    cdata_path   = "./data/cdata/KapK-cdata-manuscript-20221211.tsv",
    nobloom_path = "./data/cdata/KapK-cdata-manuscript-20221211-labels-nobloom.tsv",
    output_dir   = "./results/functional_agp"
) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # Metadata
    kapk_cdata    <- read_tsv(cdata_path,   show_col_types = FALSE)
    label_nobloom <- read_tsv(nobloom_path, show_col_types = FALSE)
    kapk_cdata    <- kapk_cdata %>% filter(figure_names %in% label_nobloom$label)

    # Derive short labels: site + last part of KapK ID
    # e.g. "75_B1_83_L2_KapK-205A"    → "75_205A"
    #      "119_B3_116_L1_KapK-12-1-41" → "119_41"
    kapk_cdata <- kapk_cdata %>%
        mutate(short_label = paste(site, sub(".*-", "", figure_names), sep = "_"))

    # sample_meta: maps 10-char folder name → full label (MD5) + figure metadata
    sample_meta <- as.data.table(kapk_cdata)[,
        .(sample = substr(label, 1, 10), label, short_label, member_unit, site_rnk)]

    # Sample order matching virome figures (08d): B1→B2→B3, then site_rnk within unit,
    # preserving paper_sample_mapping.tsv row order as secondary sort within same rnk.
    sample_map   <- fread("./data/cdata/paper_sample_mapping.tsv")
    kc_slim      <- sample_meta[, .(file_label = sample, member_unit, site_rnk)]
    meta_ordered <- merge(sample_map, kc_slim, by = "file_label")
    sample_order <- unique(meta_ordered[
        order(match(member_unit, c("B1", "B2", "B3")), site_rnk)
    ]$short_label)

    cat("Loading KEGG modules (all)...\n")
    kegg_all <- load_kegg_modules(kegg_dir, sample_meta)
    cat("  ", uniqueN(kegg_all$sample), "samples,",
        uniqueN(kegg_all$module_name), "modules\n")

    cat("Computing KEGG damage (mean_p_damaged per module)...\n")
    kegg_dmg_stats <- compute_kegg_damage(kegg_dir, sample_meta)

    cat("Loading CAZy abundances (all)...\n")
    cazy_all <- load_cazy_abundances(cazy_dir, sample_meta)
    cat("  ", uniqueN(cazy_all$sample), "samples,",
        uniqueN(cazy_all$group), "families\n")

    cat("Computing CAZy damage (mean_p_damaged per family)...\n")
    cazy_dmg_stats <- compute_cazy_damage(cazy_dir, sample_meta)

    # Apply completeness filter (paper: > 0.8) before building subsets
    kegg_all_filt <- as.data.table(kegg_all)[stepwise_module_completeness > 0.8]

    # Build "damaged" subsets: completeness-filtered rows where mean_p_damaged >= threshold
    kegg_dmg_keys <- kegg_dmg_stats[mean_p_damaged >= DAMAGE_THRESHOLD, .(sample, module)]
    kegg_dmg      <- merge(kegg_all_filt, kegg_dmg_keys, by = c("sample", "module"))

    cazy_dmg_keys <- cazy_dmg_stats[mean_p_damaged >= DAMAGE_THRESHOLD, .(sample, group)]
    cazy_dmg      <- merge(as.data.table(cazy_all), cazy_dmg_keys, by = c("sample", "group"))

    cat("Building functional info...\n")
    func_all <- build_functional_info(kegg_all_filt, as.data.table(cazy_all))
    func_dmg <- build_functional_info(kegg_dmg, cazy_dmg)

    cat("Processing functional data...\n")
    processed <- process_functional_data(func_all, func_dmg, kapk_cdata, sample_order)
    saveRDS(processed, file.path(output_dir, "heatmap_processed.rds"))
    cat("  Saved: heatmap_processed.rds\n")

    cat("Plotting heatmap...\n")
    p <- plot_functional_heatmap(processed$all, processed$dmg)
    ggsave(file.path(output_dir, "fig_functional_heatmap.pdf"), p,
           width = 200, height = 140, units = "mm")
    ggsave(file.path(output_dir, "fig_functional_heatmap.png"), p,
           width = 200, height = 140, units = "mm", dpi = 200, bg = "white")
    cat("  Saved: fig_functional_heatmap.pdf / .png\n")

    # ------------------------------------------------------------------
    # Damage distribution diagnostic
    # ------------------------------------------------------------------
    cat("Loading all protein damage scores for distribution diagnostic...\n")
    prot_all_dt <- load_all_protein_damage(kegg_dir, sample_meta)
    cat("  ", nrow(prot_all_dt), "proteins across",
        uniqueN(prot_all_dt$sample), "samples\n")

    cat("Plotting damage distribution...\n")
    dist_plots <- plot_damage_distribution(prot_all_dt)
    pdf(file.path(output_dir, "fig_damage_distribution.pdf"), width = 10, height = 5)
    print(dist_plots$global)
    print(dist_plots$per_sample)
    print(dist_plots$frac)
    dev.off()
    cat("  Saved: fig_damage_distribution.pdf\n")

    # ------------------------------------------------------------------
    # Module enrichment analysis (cross-sample Spearman correlation)
    # ------------------------------------------------------------------
    cat("Loading per-sample damage stats...\n")
    damage_stats <- load_sample_damage_stats(kegg_dir, sample_meta)
    cat("  ", nrow(damage_stats), "samples, mean_damage range:",
        round(min(damage_stats$mean_damage), 3), "-",
        round(max(damage_stats$mean_damage), 3), "\n")

    # ------------------------------------------------------------------
    # Approach 1: avg_coverage (all reads) vs mean_damage
    # ------------------------------------------------------------------
    cat("Testing enrichment: Spearman (avg_coverage ~ mean_damage, all modules)...\n")
    enr_all <- test_enrichment_spearman(kegg_all, damage_stats)
    fwrite(enr_all, file.path(output_dir, "enrichment_spearman_all.tsv"), sep = "\t")
    cat("  ", nrow(enr_all[q.value < 0.1 & rho > 0]),
        "modules enriched in high-damage samples (FDR < 0.1)\n")

    cat("  Figure module results (avg_coverage):\n")
    enr_fig <- test_enrichment_spearman(kegg_all, damage_stats, FIGURE_MODULES_KEGG)
    fwrite(enr_fig, file.path(output_dir, "enrichment_spearman_figure.tsv"), sep = "\t")
    print(enr_fig[order(p.value), .(module_name, n_samples, rho, p.value, q.value)])

    # ------------------------------------------------------------------
    # Approach 2: ancient read fraction per module vs mean_damage
    # Corrects for living organisms contributing to avg_coverage
    # ------------------------------------------------------------------
    cat("Loading module ancient read coverage (emi.functional.tsv)...\n")
    ancient_cov <- load_module_ancient_coverage(kegg_dir, sample_meta)
    cat("  ", uniqueN(ancient_cov$sample), "samples,",
        uniqueN(ancient_cov$module), "modules\n")

    cat("Testing enrichment: Spearman (ancient_cov_frac ~ mean_damage, all modules)...\n")
    enr_anc_all <- test_enrichment_spearman_ancient(ancient_cov, damage_stats)
    fwrite(enr_anc_all, file.path(output_dir, "enrichment_spearman_ancient_all.tsv"), sep = "\t")
    cat("  ", nrow(enr_anc_all[q.value < 0.1 & rho > 0]),
        "modules enriched in high-damage samples (FDR < 0.1)\n")

    cat("  Figure module results (ancient_cov_frac):\n")
    enr_anc_fig <- test_enrichment_spearman_ancient(ancient_cov, damage_stats, FIGURE_MODULES_KEGG)
    fwrite(enr_anc_fig, file.path(output_dir, "enrichment_spearman_ancient_figure.tsv"), sep = "\t")
    print(enr_anc_fig[order(p.value), .(module_name, n_samples, rho, p.value, q.value)])

    # ------------------------------------------------------------------
    # Figures (use ancient-based enrichment as primary)
    # ------------------------------------------------------------------
    cat("Plotting enrichment figures...\n")

    p_vol_anc <- plot_enrichment_volcano(
        enr_anc_all, FIGURE_MODULES_KEGG,
        title = "Module enrichment in ancient fraction (ancient reads, Spearman \u03c1)"
    )
    ggsave(file.path(output_dir, "fig_enrichment_volcano_ancient.pdf"),
           p_vol_anc, width = 10, height = 8)

    p_lol_anc <- plot_enrichment_lollipop(enr_anc_all, FIGURE_MODULES_KEGG)
    ggsave(file.path(output_dir, "fig_enrichment_lollipop_ancient.pdf"),
           p_lol_anc, width = 10, height = 12)

    # Keep avg_coverage-based figures for comparison
    p_vol <- plot_enrichment_volcano(
        enr_all, FIGURE_MODULES_KEGG,
        title = "Module enrichment (all reads avg_coverage, Spearman \u03c1)"
    )
    ggsave(file.path(output_dir, "fig_enrichment_volcano.pdf"),
           p_vol, width = 10, height = 8)

    p_lol <- plot_enrichment_lollipop(enr_all, FIGURE_MODULES_KEGG)
    ggsave(file.path(output_dir, "fig_enrichment_lollipop.pdf"),
           p_lol, width = 10, height = 12)

    # ------------------------------------------------------------------
    # Key enzyme ancient read analysis
    # ------------------------------------------------------------------
    cat("Loading key enzyme data...\n")
    enzyme_dt <- load_key_enzyme_data(kegg_dir, sample_meta)
    cat("  ", uniqueN(enzyme_dt$ko), "KOs across",
        uniqueN(enzyme_dt$sample), "samples\n")

    fwrite(enzyme_dt[order(pathway_order, gene, site_rnk),
                     .(ko, gene, pathway, short_label, member_unit,
                       n_reads, n_ancient, ancient_pct = round(ancient_pct, 1),
                       mean_posterior)],
           file.path(output_dir, "key_enzyme_ancient.tsv"), sep = "\t")
    cat("  Saved: key_enzyme_ancient.tsv\n")

    cat("Plotting key enzyme dotplot...\n")
    p_enz <- plot_key_enzyme_dotplot(enzyme_dt)
    ggsave(file.path(output_dir, "fig_key_enzyme_dotplot.pdf"),
           p_enz, width = 183, height = 230, units = "mm")
    cat("  Saved: fig_key_enzyme_dotplot.pdf\n")

    cat("Done.\n")
    invisible(list(
        processed   = processed,
        enr_all     = enr_all,     enr_fig     = enr_fig,
        enr_anc_all = enr_anc_all, enr_anc_fig = enr_anc_fig,
        enzyme_dt   = enzyme_dt
    ))
}

if (sys.nframe() == 0) {
    main()
}
