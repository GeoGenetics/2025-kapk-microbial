# 08d--virome-figures.R
# Publication-quality viral figures (Nature format) using DART/EMI data.
#
# Three figures:
#   fig_virome_composition   — stacked bar of phyla relative abundance per sample
#   fig_virome_phyla         — 3-panel pointrange: n_refs / divergence / damage by phylum
#   fig_virome_ecosystem     — DART damage by IMGVR ecosystem of origin
#
# Data source: results/virome_agp/functional/{sample}/viral_emi.functional.tsv
# Filtering  : rank filter (prot_coverage >= 0.20, n_genes > 1 or Cressdna)
#              domain == d__Viruses, phylum in vir_phyla

library(tidyverse)
library(data.table)
library(janitor)
library(showtext)
library(patchwork)

source("./libs/lib.R")
source("./libs/lib-vir.R")
showtext_auto()


# =============================================================================
# Constants
# =============================================================================

PHYLA_COLORS <- c(
    "p__Uroviricota"        = "#457B9D",
    "p__Cressdnaviricota"   = "#C1666B",
    "p__Hofneiviricota"     = "#5C8A5E",
    "p__Phixviricota"       = "#D4A853",
    "p__Nucleocytoviricota" = "#7B6D8D",
    "Methanogenic viruses"  = "#8B6150"
)

PHYLA_LABELS <- c(
    "p__Uroviricota"        = "Uroviricota",
    "p__Cressdnaviricota"   = "Cressdnaviricota",
    "p__Hofneiviricota"     = "Hofneiviricota",
    "p__Phixviricota"       = "Phixviricota",
    "p__Nucleocytoviricota" = "Nucleocytoviricota",
    "Methanogenic viruses"  = "Methanogenic viruses"
)

UNIT_COLORS <- c("B1" = "#4097AA", "B2" = "#C94A6B", "B3" = "#AFCCB8")
UNIT_LEVELS <- c("B1", "B2", "B3")

ECO_BROAD_COLORS <- c(
    "Freshwater"      = "#74B9D5",
    "Marine"          = "#1A5276",
    "Terrestrial"     = "#7A9E67",
    "Host-associated" = "#D4856A",
    "Saline/Thermal"  = "#E2C05A",
    "Unclassified"    = "#C8C8C8"
)
ECO_BROAD_LEVELS <- names(ECO_BROAD_COLORS)

map_ecosystem_broad <- function(eco) {
    dplyr::case_when(
        is.na(eco) | eco == ";;;"                                  ~ "Unclassified",
        grepl("Aquatic;Freshwater|Aquatic;Floodplain|Deep subsurface;Groundwater", eco) ~ "Freshwater",
        grepl("Aquatic;Marine|Aquatic;Oceanic",                    eco) ~ "Marine",
        grepl("Terrestrial;",                                      eco) ~ "Terrestrial",
        grepl("Host-associated",                                   eco) ~ "Host-associated",
        grepl("Non-marine Saline|Hypersaline|Thermal springs",     eco) ~ "Saline/Thermal",
        TRUE                                                            ~ "Unclassified"
    )
}

ECOSYSTEMS <- c(
    "Environmental;Terrestrial;Soil;Wetlands",
    "Environmental;Aquatic;Marine;Coastal",
    "Environmental;Aquatic;Freshwater;River",
    "Environmental;Aquatic;Freshwater;Lake",
    "Environmental;Aquatic;Marine;Wetlands",
    "Environmental;Aquatic;Freshwater;Ice",
    "Environmental;Aquatic;Non-marine Saline and Alkaline;Saline",
    "Environmental;Aquatic;Deep subsurface;Groundwater"
)

ECOSYSTEMS_LABELS <- c(
    "Environmental;Terrestrial;Soil;Wetlands"                      = "Wetland soil",
    "Environmental;Aquatic;Marine;Coastal"                         = "Marine coastal",
    "Environmental;Aquatic;Freshwater;River"                       = "Freshwater river",
    "Environmental;Aquatic;Freshwater;Lake"                        = "Freshwater lake",
    "Environmental;Aquatic;Marine;Wetlands"                        = "Marine wetland",
    "Environmental;Aquatic;Freshwater;Ice"                         = "Freshwater ice",
    "Environmental;Aquatic;Non-marine Saline and Alkaline;Saline"  = "Saline water",
    "Environmental;Aquatic;Deep subsurface;Groundwater"            = "Groundwater"
)


# =============================================================================
# Data loading
# =============================================================================

load_metadata <- function(path) fread(path, showProgress = FALSE)

load_taxonomy <- function(path) {
    dt <- fread(path, col.names = c("reference", "tax_string"), showProgress = FALSE)
    dt[, c("domain", "lineage", "kingdom", "phylum", "class",
           "order", "family", "genus", "species", "strain") :=
        tstrsplit(tax_string, ";", fixed = TRUE)]
    dt[, tax_string := NULL]
    dt
}

load_imgvr <- function(path) {
    dt <- fread(path, showProgress = FALSE) |> clean_names()
    dt[, n_cds := as.integer(
        sub("^[^;]+;([^;]+);.*", "\\1",
            gene_content_total_genes_cds_t_rna_ge_nomad_marker)
    )]
    dt[, ecosystem := fifelse(
        ecosystem_classification == "" | ecosystem_classification == "Unclassified",
        NA_character_, ecosystem_classification
    )]
    dt[, .(uvig, n_cds, ecosystem)]
}

load_emi_functional <- function(functional_dir, valid_samples) {
    files <- list.files(functional_dir, pattern = "^viral_emi\\.functional\\.tsv$",
                        recursive = TRUE, full.names = TRUE)
    samples <- basename(dirname(files))
    # keep only sample-level subdirs (exclude stray top-level file)
    keep <- samples %in% valid_samples
    files <- files[keep]; samples <- samples[keep]

    dt <- rbindlist(mapply(function(f, s) {
        d <- fread(f, showProgress = FALSE)[level == "group"]
        d[, sample := s]
        d
    }, files, samples, SIMPLIFY = FALSE), fill = TRUE)
    setnames(dt, "function_id", "reference")
    dt
}


# =============================================================================
# Filtering (identical to 08c)
# =============================================================================

apply_rank_filter <- function(emi_data, imgvr) {
    dt <- merge(emi_data, imgvr[, .(uvig, n_cds)],
                by.x = "reference", by.y = "uvig", all.x = TRUE)

    # IMGVR references: standard prot_coverage rank filter
    imgvr_dt <- dt[!is.na(n_cds) & n_cds > 0]
    imgvr_dt[, prot_coverage := n_genes / n_cds]
    imgvr_dt[, rank := fcase(
        prot_coverage >= 0.25 & n_genes > 1,                        "green",
        prot_coverage >= 0.2 & prot_coverage < 0.25 & n_genes > 1,  "yellow",
        n_genes == 1 & prot_coverage >= 0.2,                         "grey",
        default = "red"
    )]
    imgvr_dt <- imgvr_dt[rank %in% c("green", "yellow", "grey")]

    # Custom vir* references (not in IMGVR): keep all with n_genes >= 1
    vir_dt <- dt[is.na(n_cds) & grepl("^vir", reference)]
    vir_dt[, prot_coverage := NA_real_]
    vir_dt[, rank := fifelse(n_genes > 1, "green", "grey")]

    rbind(imgvr_dt, vir_dt, fill = TRUE)
}


# =============================================================================
# Nature theme
# =============================================================================

theme_nature <- function(base_size = 9) {
    theme_bw(base_size = base_size) +
        theme(
            axis.text.x      = element_text(size = 8),
            axis.text.y      = element_text(size = 8),
            axis.title       = element_text(size = 9),
            legend.text      = element_text(size = 8),
            legend.title     = element_text(size = 8),
            strip.text       = element_text(size = 9),
            strip.background = element_blank(),
            panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
            panel.grid.minor = element_blank(),
            panel.spacing    = unit(0.3, "lines")
        )
}


# =============================================================================
# Figure A: Composition stacked bar
# =============================================================================

plot_composition <- function(dat, sample_order, phyla_order) {
    # Use same order as bottom panels (by n_ref count)
    phyla_labels_ordered <- PHYLA_LABELS[phyla_order]

    comp <- dat[, .(n_reads = sum(n_reads)), by = .(short_label, member_unit, phylum)]
    comp[, total := sum(n_reads), by = .(short_label)]
    comp[, prop  := n_reads / total]

    comp[, phylum_label := factor(PHYLA_LABELS[phylum], levels = phyla_labels_ordered)]
    comp[, short_label  := factor(short_label, levels = sample_order)]
    comp[, member_unit  := factor(member_unit, levels = UNIT_LEVELS)]

    ggplot(comp, aes(x = short_label, y = prop, fill = phylum_label)) +
        geom_bar(stat = "identity", width = 0.88, colour = "black", linewidth = 0.1, alpha = 0.7,
                 position = position_stack(reverse = TRUE)) +
        facet_grid(. ~ member_unit, scales = "free_x", space = "free_x") +
        scale_fill_manual(values = setNames(PHYLA_COLORS, PHYLA_LABELS), name = NULL,
                          guide = guide_legend(nrow = 1)) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                           expand = expansion(mult = c(0, 0.02)),
                           breaks = c(0, 0.25, 0.5, 0.75, 1)) +
        theme_nature() +
        theme(
            axis.text.x        = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
            axis.title.y       = element_text(size = 9, angle = 90),
            legend.position    = "top",
            legend.justification = "left",
            legend.text        = element_text(size = 8, face = "italic"),
            legend.key.size    = unit(3, "mm"),
            legend.margin      = margin(0, 0, 0, 0),
            legend.box.spacing = unit(1, "mm"),
            strip.text         = element_text(size = 9, face = "bold"),
            panel.spacing.x    = unit(1.5, "mm"),
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_blank()
        ) +
        labs(x = NULL, y = "Relative read abundance")
}


# =============================================================================
# Figure C: Ecosystem composition stacked bar
# =============================================================================

plot_ecosystem_composition <- function(dat, sample_order) {
    comp <- as_tibble(dat) |>
        mutate(eco_broad = map_ecosystem_broad(ecosystem)) |>
        group_by(short_label, member_unit, eco_broad) |>
        summarise(n_reads = sum(n_reads), .groups = "drop") |>
        group_by(short_label) |>
        mutate(total = sum(n_reads), prop = n_reads / total) |>
        ungroup() |>
        mutate(
            eco_broad   = factor(eco_broad,   levels = ECO_BROAD_LEVELS),
            short_label = factor(short_label, levels = sample_order),
            member_unit = factor(member_unit, levels = UNIT_LEVELS)
        )

    ggplot(comp, aes(x = short_label, y = prop, fill = eco_broad)) +
        geom_bar(stat = "identity", width = 0.88, colour = "black", linewidth = 0.1,
                 alpha = 0.7, position = position_stack(reverse = TRUE)) +
        facet_grid(. ~ member_unit, scales = "free_x", space = "free_x") +
        scale_fill_manual(values = ECO_BROAD_COLORS, name = NULL,
                          guide = guide_legend(nrow = 1)) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                           expand = expansion(mult = c(0, 0.02)),
                           breaks = c(0, 0.25, 0.5, 0.75, 1)) +
        theme_nature() +
        theme(
            axis.text.x        = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
            axis.title.y       = element_text(size = 9, angle = 90),
            legend.text        = element_text(size = 8),
            legend.position    = "top",
            legend.justification = "left",
            legend.key.size    = unit(3, "mm"),
            legend.margin      = margin(0, 0, 0, 0),
            legend.box.spacing = unit(1, "mm"),
            strip.text         = element_blank(),
            panel.spacing.x    = unit(1.5, "mm"),
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_blank()
        ) +
        labs(x = NULL, y = "Relative read abundance")
}


# =============================================================================
# Figure B: Phyla — n_refs / AA identity / DART damage
#   n_refs: jitter + median crossbar coloured by depth unit (22 pts/phylum)
#   identity & DART: violin coloured by phylum (pooled across units)
# =============================================================================

plot_phyla_stats <- function(dat, phyla_order) {
    phyla_levels <- rev(gsub("p__", "", phyla_order))   # top-to-bottom on Y
    phyla_fills  <- setNames(PHYLA_COLORS, PHYLA_LABELS)

    # Shared Y-axis theme: italic phylum labels, no ticks clutter
    y_shared <- theme(
        axis.text.y  = element_text(size = 8, face = "italic"),
        axis.ticks.y = element_line(linewidth = 0.2)
    )
    y_hidden <- theme(
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank()
    )

    base <- dat |>
        filter(phylum %in% phyla_order) |>
        mutate(
            phylum_label = factor(gsub("p__", "", phylum), levels = phyla_levels),
            member_unit  = factor(member_unit, levels = c("B1", "B2", "B3"))
        )

    # p_nref: one point per sample × phylum — shows inter-sample variation
    nref_dat <- base |>
        group_by(sample, member_unit, phylum_label) |>
        summarise(n_ref = n_distinct(reference), .groups = "drop")

    p_nref <- ggplot(nref_dat, aes(y = phylum_label, x = n_ref, colour = member_unit)) +
        geom_point(alpha = 0.7, size = 1.1,
                   position = position_jitter(height = 0.18, width = 0, seed = 42)) +
        stat_summary(aes(group = phylum_label), colour = "grey15", fill = "grey15",
                     fun = median, geom = "crossbar",
                     fun.min = median, fun.max = median,
                     width = 0.45, linewidth = 0.25) +
        scale_colour_manual(values = UNIT_COLORS, name = NULL,
                            guide = guide_legend(override.aes = list(size = 2.5, alpha = 1),
                                                 nrow = 1)) +
        scale_x_log10(labels = scales::label_comma(), breaks = c(1, 10, 100)) +
        theme_nature() + y_shared +
        theme(legend.position        = c(0.97, 0.03),
              legend.justification   = c(1, 0),
              legend.direction       = "horizontal",
              legend.background      = element_rect(fill = "white", colour = NA, linewidth = 0),
              legend.key.size        = unit(1.2, "mm"),
              legend.key.spacing.x   = unit(0.5, "mm"),
              legend.text            = element_text(size = 6),
              legend.margin          = margin(0, 0, 0, 0)) +
        labs(y = NULL, x = "Unique viral refs.")

    # p_iden: violin per phylum, coloured by phylum — shows AA identity distribution
    p_iden <- base |>
        ggplot(aes(y = phylum_label, x = avg_identity * 100,
                   fill = phylum_label, colour = phylum_label)) +
        geom_violin(alpha = 0.4, width = 0.8, trim = TRUE, linewidth = 0.2,
                    bw = 0.35) +
        geom_jitter(alpha = 0.04, size = 0.12, height = 0.25, width = 0) +
        stat_summary(fun = median, geom = "point", shape = 21,
                     size = 1.8, fill = "white", colour = "grey20", stroke = 0.5) +
        scale_fill_manual(values   = phyla_fills, guide = "none") +
        scale_colour_manual(values = phyla_fills, guide = "none") +
        scale_x_continuous(labels = scales::label_number(suffix = "%"),
                           limits = c(85, 100)) +
        theme_nature() + y_hidden +
        labs(y = NULL, x = "AA identity to reference")

    # p_dmg: violin per phylum, coloured by phylum — DART authentication
    p_dmg <- base |>
        ggplot(aes(y = phylum_label, x = mean_posterior,
                   fill = phylum_label, colour = phylum_label)) +
        geom_violin(alpha = 0.4, width = 0.8, trim = FALSE, linewidth = 0.2,
                    bw = 0.02) +
        geom_jitter(alpha = 0.04, size = 0.12, height = 0.25, width = 0) +
        stat_summary(fun = median, geom = "point", shape = 21,
                     size = 1.8, fill = "white", colour = "grey20", stroke = 0.5) +
        scale_fill_manual(values   = phyla_fills, guide = "none") +
        scale_colour_manual(values = phyla_fills, guide = "none") +
        coord_cartesian(xlim = c(0.85, 1.02)) +
        scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                           breaks = c(0.85, 0.9, 0.95, 1.0)) +
        theme_nature() + y_hidden +
        labs(y = NULL, x = "Protein damage probability")

    list(n_ref = p_nref, identity = p_iden, damage = p_dmg)
}


# Assemble phyla stats as a patchwork row (facet_grid-like)
assemble_phyla_row <- function(phyla_list) {
    p1 <- phyla_list$n_ref    + labs(tag = "B") + theme(plot.margin = margin(2, 3, 2, 2),
                                                        plot.tag = element_text(size = 11, face = "bold"))
    p2 <- phyla_list$identity + theme(plot.margin = margin(2, 3, 2, 3))
    p3 <- phyla_list$damage   + theme(plot.margin = margin(2, 2, 2, 3))
    (p1 | p2 | p3) + plot_layout(widths = c(1, 1, 1))
}


# =============================================================================
# Combined: composition (top) + phyla stats 3 panels (bottom)
# =============================================================================

plot_combined <- function(dat, sample_order, phyla_order) {
    p_comp  <- plot_composition(dat, sample_order, phyla_order) + labs(tag = "A") +
                   theme(plot.tag = element_text(size = 11, face = "bold"))
    p_eco   <- plot_ecosystem_composition(dat, sample_order) + labs(tag = "C") +
                   theme(plot.tag = element_text(size = 11, face = "bold"))
    p_phyla <- plot_phyla_stats(as_tibble(dat), phyla_order)
    p_bot   <- assemble_phyla_row(p_phyla)
    free(p_comp) / p_bot / free(p_eco) +
        plot_layout(heights = c(1, 1, 1)) +
        plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2)))
}


# =============================================================================
# Figure C: Environmental origin — AA identity by IMGVR ecosystem
# =============================================================================

plot_ecosystem <- function(dat, ecosystems) {
    eco_dat <- dat |>
        filter(!is.na(ecosystem), ecosystem %in% ecosystems) |>
        mutate(eco_label = ECOSYSTEMS_LABELS[ecosystem])

    bubble <- eco_dat |>
        group_by(eco_label, phylum) |>
        summarise(
            med_identity = median(avg_identity * 100),
            n_ref        = n_distinct(reference),
            .groups      = "drop"
        ) |>
        mutate(
            phylum_label = PHYLA_LABELS[phylum],
            phylum_label = factor(phylum_label, levels = PHYLA_LABELS),
            eco_label    = fct_reorder(eco_label, med_identity, .fun = median)
        )

    ggplot(bubble, aes(x = med_identity, y = eco_label,
                       size = n_ref, colour = phylum_label)) +
        geom_point(alpha = 0.82) +
        scale_colour_manual(values = setNames(PHYLA_COLORS, PHYLA_LABELS), name = NULL) +
        scale_size_area(name = "Viral refs.", max_size = 9,
                        breaks = c(1, 5, 25, 100)) +
        scale_x_continuous(labels = scales::label_number(suffix = "%"),
                           limits = c(90, 100),
                           breaks = c(90, 92, 94, 96, 98, 100)) +
        theme_nature() +
        theme(
            legend.position    = "right",
            legend.key.size    = unit(3, "mm"),
            panel.grid.major.y = element_line(linewidth = 0.15, colour = "grey92"),
            panel.grid.major.x = element_blank()
        ) +
        labs(x = "Median AA identity to reference (%)", y = NULL)
}


# =============================================================================
# Option 2: Streamgraph — stacked area, samples ordered shallow → deep on X
# =============================================================================

plot_streamgraph <- function(dat, sample_order) {
    comp <- dat[, .(n_reads = sum(n_reads)), by = .(short_label, member_unit, phylum)]
    comp[, total := sum(n_reads), by = short_label]
    comp[, prop  := n_reads / total]

    comp[, phylum_label := factor(PHYLA_LABELS[phylum], levels = rev(PHYLA_LABELS))]
    comp[, short_label  := factor(short_label, levels = sample_order)]
    comp[, x_pos        := as.integer(short_label)]

    # Unit boundaries and midpoints for annotations
    samp_units <- unique(comp[, .(short_label, member_unit, x_pos)])
    b1_max  <- max(samp_units[member_unit == "B1", x_pos])
    b2_max  <- max(samp_units[member_unit == "B2", x_pos])
    unit_mids <- samp_units[, .(mid = mean(x_pos)), by = member_unit][
        order(match(member_unit, c("B1", "B2", "B3")))]

    ggplot(comp, aes(x = x_pos, y = prop, fill = phylum_label)) +
        geom_area(position = "stack", colour = "white", linewidth = 0.15,
                  alpha = 0.92) +
        geom_vline(xintercept = c(b1_max + 0.5, b2_max + 0.5),
                   colour = "grey45", linewidth = 0.35, linetype = "dashed") +
        annotate("text", x = unit_mids$mid, y = 1.05,
                 label = unit_mids$member_unit,
                 size = 6 / .pt, fontface = "bold", colour = "grey25",
                 hjust = 0.5) +
        scale_fill_manual(values = rev(setNames(PHYLA_COLORS, PHYLA_LABELS)),
                          name = NULL,
                          guide = guide_legend(reverse = TRUE, ncol = 2,
                                               byrow = FALSE)) +
        scale_x_continuous(breaks = seq_along(sample_order),
                           labels = sample_order,
                           expand = expansion(mult = c(0.01, 0.01))) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                           expand = expansion(mult = c(0, 0.08)),
                           breaks = c(0, 0.25, 0.5, 0.75, 1)) +
        coord_cartesian(clip = "off") +
        theme_nature() +
        theme(
            axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1,
                                               size = 5.0),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(linewidth = 0.15, colour = "grey90"),
            panel.border       = element_rect(colour = "grey75", fill = NA,
                                              linewidth = 0.3),
            legend.position    = "bottom",
            legend.key.size    = unit(2.5, "mm"),
            legend.text        = element_text(size = 5.5),
            legend.margin      = margin(0, 0, 0, 0)
        ) +
        labs(x = NULL, y = "Relative read abundance")
}


# =============================================================================
# Option 3: Dot/strip by phylum — abundance distribution, coloured by B unit
# =============================================================================

plot_lollipop <- function(dat) {
    samp_units <- unique(dat[, .(short_label, member_unit)])

    comp <- dat[, .(n_reads = sum(n_reads)), by = .(short_label, phylum)]
    comp[, total := sum(n_reads), by = short_label]
    comp[, prop  := n_reads / total]

    # Include zero-abundance phyla per sample
    all_combos <- CJ(short_label = samp_units$short_label,
                     phylum      = names(PHYLA_COLORS))
    comp_full  <- merge(all_combos, comp[, .(short_label, phylum, prop)],
                        by = c("short_label", "phylum"), all.x = TRUE)
    comp_full  <- merge(comp_full, samp_units, by = "short_label")
    comp_full[is.na(prop), prop := 0]

    comp_full[, phylum_label := factor(PHYLA_LABELS[phylum], levels = PHYLA_LABELS)]
    comp_full[, member_unit  := factor(member_unit, levels = c("B1", "B2", "B3"))]

    ggplot(comp_full, aes(x = phylum_label, y = prop, colour = member_unit)) +
        geom_point(alpha = 0.75, size = 1.5,
                   position = position_jitter(width = 0.2, height = 0, seed = 42)) +
        stat_summary(aes(group = phylum_label), colour = "grey15",
                     fun = median, geom = "crossbar",
                     fun.min = median, fun.max = median,
                     width = 0.55, linewidth = 0.28) +
        scale_colour_manual(values = UNIT_COLORS, name = NULL,
                            guide = guide_legend(
                                override.aes = list(size = 2.5, alpha = 1))) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                           limits = c(-0.02, 1),
                           breaks = c(0, 0.25, 0.5, 0.75, 1)) +
        theme_nature() +
        theme(
            axis.text.x        = element_text(angle = 35, hjust = 1, vjust = 1,
                                               size = 6, face = "italic"),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(linewidth = 0.15, colour = "grey90"),
            panel.border       = element_rect(colour = "grey75", fill = NA,
                                              linewidth = 0.3),
            legend.position    = c(0.88, 0.82),
            legend.background  = element_rect(fill = "white", colour = NA,
                                              linewidth = 0),
            legend.key.size    = unit(2, "mm"),
            legend.text        = element_text(size = 5.5)
        ) +
        labs(x = NULL, y = "Relative read abundance")
}


# =============================================================================
# Fig-4 panel: horizontal stacked bars (Figure-4B style)
# =============================================================================

plot_integrated <- function(dat, sample_order) {
    samp_levels <- rev(sample_order)
    unit_levels <- c("B3", "B2", "B1")

    ddt <- copy(dat)
    ddt[, short_label := factor(short_label, levels = samp_levels)]
    ddt[, member_unit := factor(member_unit, levels = unit_levels)]

    comp <- ddt[, .(n_reads = sum(n_reads)), by = .(short_label, member_unit, phylum)]
    comp[, total        := sum(n_reads), by = short_label]
    comp[, prop         := n_reads / total]
    comp[, phylum_label := factor(PHYLA_LABELS[phylum], levels = PHYLA_LABELS)]

    # Unit colour sidebar: narrow panel with one tile per sample, coloured by unit
    sidebar <- unique(ddt[, .(short_label, member_unit)])
    sidebar[, x := 1]

    p_side <- ggplot(sidebar, aes(x = x, y = short_label, fill = member_unit)) +
        geom_tile(width = 1, colour = NA) +
        facet_grid(member_unit ~ ., scales = "free_y", space = "free_y",
                   switch = "y") +
        scale_fill_manual(values = UNIT_COLORS, guide = "none") +
        scale_x_continuous(expand = c(0, 0)) +
        theme_void() +
        theme(
            strip.text.y.left = element_text(size = 7, face = "bold", angle = 0,
                                             hjust = 1, margin = margin(r = 2),
                                             colour = "grey20"),
            strip.placement   = "outside",
            strip.background  = element_blank(),
            panel.spacing.y   = unit(1, "mm"),
            plot.margin       = margin(4, 1, 4, 2)
        )

    p_bars <- ggplot(comp, aes(y = short_label, x = prop, fill = phylum_label)) +
        geom_bar(stat = "identity", width = 0.88, colour = NA) +
        facet_grid(member_unit ~ ., scales = "free_y", space = "free_y") +
        scale_fill_manual(values = setNames(PHYLA_COLORS, PHYLA_LABELS), name = NULL,
                          guide = guide_legend(ncol = 1, byrow = FALSE)) +
        scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                           expand = expansion(mult = c(0, 0.01)),
                           breaks = c(0, 0.5, 1)) +
        theme_nature() +
        theme(
            panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey88"),
            panel.grid.major.y = element_blank(),
            panel.border       = element_rect(colour = "grey70", fill = NA,
                                              linewidth = 0.35),
            strip.text         = element_blank(),
            strip.background   = element_blank(),
            panel.spacing.y    = unit(1, "mm"),
            axis.text.y        = element_text(size = 5.5),
            axis.text.x        = element_text(size = 6),
            legend.position    = "right",
            legend.key.size    = unit(3, "mm"),
            legend.text        = element_text(size = 6),
            legend.margin      = margin(0, 0, 0, 2),
            plot.margin        = margin(4, 2, 4, 0)
        ) +
        labs(y = NULL, x = "Relative read abundance")

    ggpubr::ggarrange(p_side, p_bars, ncol = 2, align = "hv",
                      widths = c(0.12, 1))
}


# =============================================================================
# Bottom-left: AA similarity by phylum (Fig-6D style)
#   Jittered points per reference × sample, coloured by B unit, median bar
# =============================================================================

plot_sim_phylum <- function(dat, phyla_order) {
    phyla_levels <- rev(gsub("p__", "", phyla_order))

    base <- as_tibble(dat) |>
        filter(phylum %in% phyla_order) |>
        mutate(
            phylum_label = factor(gsub("p__", "", phylum), levels = phyla_levels),
            member_unit  = factor(member_unit, levels = c("B1", "B2", "B3"))
        )

    ggplot(base, aes(y = phylum_label, x = avg_identity * 100,
                     colour = member_unit)) +
        geom_jitter(alpha = 0.25, size = 0.5, height = 0.22, width = 0) +
        stat_summary(aes(group = phylum_label), colour = "grey15",
                     fun = median, geom = "crossbar",
                     fun.min = median, fun.max = median,
                     width = 0.55, linewidth = 0.28) +
        scale_colour_manual(values = UNIT_COLORS, name = NULL,
                            guide = guide_legend(
                                override.aes = list(size = 2.5, alpha = 1))) +
        scale_x_continuous(labels = scales::label_number(suffix = "%"),
                           limits = c(80, 100),
                           breaks = c(80, 90, 100)) +
        theme_nature() +
        theme(
            panel.grid.major.x = element_line(linewidth = 0.15, colour = "grey90"),
            panel.grid.major.y = element_blank(),
            panel.border       = element_rect(colour = "grey70", fill = NA,
                                              linewidth = 0.35),
            axis.text.y        = element_text(size = 6, face = "italic"),
            legend.position    = c(0.2, 0.12),
            legend.background  = element_rect(fill = "white", colour = NA,
                                              linewidth = 0),
            legend.key.size    = unit(2, "mm"),
            legend.text        = element_text(size = 5.5)
        ) +
        labs(y = NULL, x = "AA similarity (%)")
}


# =============================================================================
# Bottom-right: AA similarity by IMGVR ecosystem (Fig-6E right style)
#   Points per phylum per ecosystem, coloured by phylum, median + IQR
# =============================================================================

plot_sim_ecosystem <- function(dat, ecosystems) {
    eco_dat <- as_tibble(dat) |>
        filter(!is.na(ecosystem), ecosystem %in% ecosystems) |>
        mutate(
            eco_label    = ECOSYSTEMS_LABELS[ecosystem],
            phylum_label = factor(PHYLA_LABELS[phylum], levels = PHYLA_LABELS)
        )

    # One point per unique reference × ecosystem, averaged across samples
    eco_refs <- eco_dat |>
        group_by(eco_label, phylum_label, reference) |>
        summarise(identity = mean(avg_identity * 100), .groups = "drop") |>
        mutate(eco_label = fct_reorder(eco_label, identity, .fun = median))

    eco_stats <- eco_refs |>
        group_by(eco_label, phylum_label) |>
        summarise(med = median(identity),
                  q25 = quantile(identity, 0.25),
                  q75 = quantile(identity, 0.75),
                  .groups = "drop")

    ggplot(eco_refs, aes(y = eco_label, x = identity, colour = phylum_label)) +
        geom_jitter(alpha = 0.3, size = 0.5, height = 0.2, width = 0) +
        geom_segment(data = eco_stats,
                     aes(x = q25, xend = q75, yend = eco_label),
                     linewidth = 0.5, alpha = 0.55) +
        geom_point(data = eco_stats, aes(x = med), size = 1.6) +
        scale_colour_manual(values = setNames(PHYLA_COLORS, PHYLA_LABELS),
                            name = NULL) +
        scale_x_continuous(labels = scales::label_number(suffix = "%"),
                           limits = c(80, 100),
                           breaks = c(80, 90, 100)) +
        theme_nature() +
        theme(
            panel.grid.major.x = element_line(linewidth = 0.15, colour = "grey90"),
            panel.grid.major.y = element_blank(),
            panel.border       = element_rect(colour = "grey70", fill = NA,
                                              linewidth = 0.35),
            axis.text.y        = element_text(size = 6),
            legend.position    = "right",
            legend.key.size    = unit(2.5, "mm"),
            legend.text        = element_text(size = 5.5)
        ) +
        labs(y = NULL, x = "AA similarity (%)")
}


# =============================================================================
# Combined three-panel figure
# =============================================================================

plot_virome_panel <- function(dat, sample_order, phyla_order, ecosystems) {
    p_top    <- plot_integrated(dat, sample_order)
    p_bottom <- ggpubr::ggarrange(
        plot_sim_phylum(dat, phyla_order),
        plot_sim_ecosystem(dat, ecosystems),
        ncol = 2, nrow = 1, align = "hv",
        widths = c(1, 1.4)
    )
    ggpubr::ggarrange(p_top, p_bottom, nrow = 2, heights = c(2.2, 1.8))
}


# =============================================================================
# Main
# =============================================================================

main <- function(
    functional_dir = "./results/virome_agp/functional",
    output_dir     = "./results/virome_agp"
) {
    # ---- Metadata -----------------------------------------------------------
    kapk_cdata    <- load_metadata("./data/cdata/KapK-cdata-manuscript-20221211.tsv")
    label_nobloom <- load_metadata("./data/cdata/KapK-cdata-manuscript-20221211-labels-nobloom.tsv")
    kapk_cdata    <- kapk_cdata[figure_names %in% label_nobloom$label]
    kapk_cdata[, short_label := paste(site, sub(".*-", "", figure_names), sep = "_")]

    sample_map  <- load_metadata("./data/cdata/paper_sample_mapping.tsv")
    kapk_lookup <- kapk_cdata[, .(label, figure_names, member_unit, site_rnk, site,
                                   file_label = substr(label, 1, 10))]
    sample_meta <- merge(sample_map, kapk_lookup, by = "file_label")
    setnames(sample_meta, "file_label", "sample")

    # Sample order for x-axis (member_unit then site_rnk within unit)
    sample_order <- unique(
        sample_meta[order(match(member_unit, c("B1", "B2", "B3")), site_rnk)]$short_label
    )

    tax_info <- load_taxonomy("./data/taxonomy/hires-organelles-viruses-arctic.tax.tsv")
    imgvr    <- load_imgvr("./data/cdata/IMGVR_all_Sequence_information-high_confidence.tsv")

    # ---- Load & filter EMI data ---------------------------------------------
    cat("Loading viral_emi.functional.tsv...\n")
    emi_data <- load_emi_functional(functional_dir, valid_samples = sample_meta$sample)
    cat("  Raw:", nrow(emi_data), "ref×sample pairs,",
        uniqueN(emi_data$sample), "samples\n")

    emi_filt <- apply_rank_filter(emi_data, imgvr)

    # Join taxonomy — all.x = TRUE keeps custom vir* refs (not in tax_info)
    dat <- merge(emi_filt, tax_info[, .(reference, domain, phylum)],
                 by = "reference", all.x = TRUE)
    # Keep d__Viruses + unclassified custom vir* refs
    dat <- dat[domain == "d__Viruses" | (is.na(domain) & grepl("^vir", reference))]
    dat <- merge(dat, sample_meta[, .(sample, short_label, member_unit, site_rnk)],
                 by = "sample")
    dat <- merge(dat, imgvr[, .(uvig, ecosystem)],
                 by.x = "reference", by.y = "uvig", all.x = TRUE)

    # Promote Cressdna grey → blue, drop remaining grey (IMGVR refs only)
    dat[rank == "grey" & phylum == "p__Cressdnaviricota", rank := "blue"]
    dat <- dat[rank != "grey"]

    # Classify custom vir* references (same logic as 08--virome.R):
    #   vir291/vir303 → p__Uroviricota (Caudoviricetes)
    #   all other vir* → "Methanogenic viruses"
    dat[is.na(domain), domain := "d__Viruses"]
    dat[reference %in% c("vir291", "vir303"), phylum := "p__Uroviricota"]
    dat[is.na(phylum) & grepl("^vir", reference), phylum := "Methanogenic viruses"]

    # DART authentication filter: only references where reads authenticate as ancient
    dat_filt <- dat[rank %in% c("green", "yellow", "blue") &
                    phylum %in% vir_phyla &
                    mean_posterior >= 0.7]
    cat("  After filter:", nrow(dat_filt), "ref×sample pairs,",
        uniqueN(dat_filt$reference), "unique UViGs,",
        uniqueN(dat_filt$sample), "samples\n")

    # Phyla order (by n unique references)
    phyla_order <- dat_filt[, .(n = uniqueN(reference)), by = phylum][
        order(-n), phylum]
    cat("  Phyla:", paste(gsub("p__", "", phyla_order), collapse = ", "), "\n")

    # ---- Figures ------------------------------------------------------------
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    save_fig <- function(p, name, w, h) {
        ggsave(file.path(output_dir, paste0(name, ".pdf")),
               p, width = w, height = h, units = "mm")
        ggsave(file.path(output_dir, paste0(name, ".png")),
               p, width = w, height = h, units = "mm", dpi = 250, bg = "white")
        cat("  Saved:", name, "\n")
    }

    # Individual figures (kept for supplementary)
    cat("Plotting composition...\n")
    p_comp <- plot_composition(dat_filt, sample_order, phyla_order)
    save_fig(p_comp, "fig_virome_composition", 183, 68)

    cat("Plotting phyla stats...\n")
    p_phyla <- plot_phyla_stats(as_tibble(dat_filt), phyla_order)
    p_phyla_combined <- assemble_phyla_row(p_phyla)
    save_fig(p_phyla_combined, "fig_virome_phyla", 183, 85)

    cat("Plotting combined composition+phyla figure...\n")
    p_combined <- plot_combined(dat_filt, sample_order, phyla_order)
    save_fig(p_combined, "fig_virome_combined", 183, 165)

    cat("Plotting ecosystem origin...\n")
    p_eco <- plot_ecosystem(as_tibble(dat_filt), ECOSYSTEMS)
    save_fig(p_eco, "fig_virome_ecosystem", 140, 65)

    cat("Plotting integrated figure...\n")
    p_int <- plot_integrated(dat_filt, sample_order)
    save_fig(p_int, "fig_virome_main", 95, 115)

    cat("Plotting combined panel figure...\n")
    p_panel <- plot_virome_panel(dat_filt, sample_order, phyla_order, ECOSYSTEMS)
    save_fig(p_panel, "fig_virome_panel", 140, 170)

    cat("Plotting streamgraph...\n")
    p_stream <- plot_streamgraph(dat_filt, sample_order)
    save_fig(p_stream, "fig_virome_stream", 110, 65)

    cat("Plotting lollipop...\n")
    p_lollipop <- plot_lollipop(dat_filt)
    save_fig(p_lollipop, "fig_virome_lollipop", 75, 60)

    cat("Done.\n")

    invisible(list(
        dat_filt    = dat_filt,
        phyla_order = phyla_order,
        figures     = list(comp = p_comp, phyla = p_phyla, eco = p_eco,
                           main = p_int)
    ))
}

if (sys.nframe() == 0) main()
