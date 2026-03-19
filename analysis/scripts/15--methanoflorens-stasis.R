#!/usr/bin/env Rscript
# Fig 7 — Evolutionary stasis of Methanoflorens at the Kap København Formation
#
# 3-panel figure:
#   a) Post-mortem DNA damage (ancient vs deconvolved modern fractions)
#   b) Fragment length distributions (ancient = short, modern = longer, S17 = longest)
#   c) IQ-TREE ML focal subtree with BEAST2-calibrated branch lengths and
#      direct pairwise SNP counts (nucmer alignment vs ancient_consensus)

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2); library(patchwork); library(dplyr)
  library(tidyr); library(readr); library(scales)
  library(showtext)
})
font_add("DroidSans",
  regular = "/usr/share/fonts/google-droid/DroidSans.ttf",
  bold    = "/usr/share/fonts/google-droid/DroidSans-Bold.ttf")
showtext_auto()
FONT <- "DroidSans"

# ── Paths ──────────────────────────────────────────────────────────────────────
# Data lives in data/methanoflorens/ relative to the analysis/ working directory.
# ERDA package mirrors this structure so the script runs without modification.
BASE         <- "./data/methanoflorens"
MCC_FILE     <- file.path(BASE, "beast2_nodec/tipdating/codon_tipdating_mcc.tree")
IQTREE_FILE  <- file.path(BASE, "phylogenomics_79taxa_bak/08_trees/codon12_ml_tree.treefile")
ANI_FILE     <- file.path(BASE, "fastani/deconvolved_vs_refs.tsv")
FRAG_FILE    <- file.path(BASE, "deconvolve/read_length_histogram.tsv")
DMG_FILE     <- file.path(BASE, "deconvolve/damage_model.tsv")
SMILEY_FILE  <- file.path(BASE, "deconvolve/smiley_data.tsv")
S17_FRAG     <- file.path(BASE, "s17_assembly/deconvolve/read_length_histogram.tsv")
S17_SMILEY   <- file.path(BASE, "s17_assembly/deconvolve/smiley_data.tsv")
SNP_FILE     <- file.path(BASE, "phylogenomics_79taxa_bak/11_focal_snps/focal_pairwise_snp_counts.tsv")
ANI_FILE_FOC <- file.path(BASE, "phylogenomics_79taxa_bak/11_focal_snps/focal_pairwise_ani.tsv")
OUTDIR       <- "./results/figures"
META_FILE    <- file.path(BASE, "methanoflorens_metadata.txt")

# ── BEAST2 log files ──────────────────────────────────────────────────────────
NODEC_BASE      <- file.path(BASE, "beast2_nodec")
LOG_TIPDATING   <- file.path(NODEC_BASE, "tipdating/codon_tipdating.log")
LOG_2CLOCK      <- file.path(NODEC_BASE, "2clock_stasis/codon_2clock_stasis_nodec.log")
LOG_UCLD        <- file.path(NODEC_BASE, "ucld_tipdating/codon_tipdating.log")
LOG_FREEDATE_R3 <- file.path(NODEC_BASE, "freedate_rep3/codon_freedate_nodec.log")
LOG_FREEDATE_R4 <- file.path(NODEC_BASE, "freedate_rep4/codon_freedate_nodec.log")
PS_DIR          <- NODEC_BASE   # path sampling beast_ps.out files live in ps_*/
GEO_AGE_YR      <- 2e6   # ancient_consensus tip fixed at 2 Ma

# ── Pairwise SNP counts and ANI (from Snakemake workflow: phylogenomics/11_focal_snps/)
# Regenerate: snakemake phylogenomics/11_focal_snps/focal_pairwise_snp_counts.tsv
#             snakemake phylogenomics/11_focal_snps/focal_pairwise_ani.tsv
snp_tbl <- read_tsv(SNP_FILE, show_col_types = FALSE)
SNP_MOD <- snp_tbl$nondamage_snps[snp_tbl$pair == "ancient_vs_modern"]
SNP_S17 <- snp_tbl$nondamage_snps[snp_tbl$pair == "ancient_vs_s17"]
ani_tbl <- read_tsv(ANI_FILE_FOC, show_col_types = FALSE)
ANI_MOD <- ani_tbl$ani[ani_tbl$pair == "ancient_vs_modern"]
ANI_S17 <- ani_tbl$ani[ani_tbl$pair == "ancient_vs_s17"]

# ── Figure dimensions (Nature double-column) ───────────────────────────────────
FIG_W <- 183
FIG_H <- 155

# ── Colour palette ─────────────────────────────────────────────────────────────
COL_ANCIENT  <- "#A63820"   # deep terracotta — anchor color, stays saturated
COL_MODERN   <- "#6B9AB8"   # muted steel-blue — KapK-dec (deconvolved modern)
COL_RELATIVE <- "#6B8C6E"   # dusty sage       — KapK-mod (S17 living)

# ── Shared theme ───────────────────────────────────────────────────────────────
theme_nat <- function(base_size = 9) {
  theme_bw(base_size = base_size, base_family = FONT) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(linewidth = 0.2, colour = "grey92"),
      panel.border       = element_rect(linewidth = 0.4, fill = NA, colour = "grey30"),
      axis.ticks         = element_line(linewidth = 0.3, colour = "grey40"),
      axis.ticks.length  = unit(1.5, "pt"),
      axis.text          = element_text(size = 8, colour = "grey20", family = FONT),
      axis.title         = element_text(size = 9, colour = "grey10", family = FONT),
      legend.text        = element_text(size = 6.5, family = FONT),
      legend.title       = element_blank(),
      legend.key.size    = unit(9, "pt"),
      legend.background  = element_rect(fill = alpha("white", 0.9),
                                        colour = "grey80", linewidth = 0.2),
      legend.margin      = margin(2, 4, 2, 4, "pt"),
      plot.margin        = margin(3, 4, 2, 3, "pt"),
      plot.tag           = element_text(size = 11, face = "bold", colour = "grey5", family = FONT)
    )
}

# ── Data loaders ───────────────────────────────────────────────────────────────
load_damage_data <- function(path) {
  lines <- readLines(path)
  parse_section <- function(header_pat) {
    idx   <- grep(header_pat, lines, fixed = TRUE)
    if (!length(idx)) return(NULL)
    block <- lines[(idx + 1):(idx + 16)]
    read.table(text = paste(block, collapse = "\n"),
               header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }
  bind_rows(
    parse_section("5' end C->T") %>%
      transmute(pos = position, rate, ci_lo = ci_lower, ci_hi = ci_upper,
                end = "5\u2032 C\u2192T"),
    parse_section("3' end G->A") %>%
      transmute(pos = position, rate, ci_lo = ci_lower, ci_hi = ci_upper,
                end = "3\u2032 G\u2192A")
  )
}

load_damage_params <- function(path) {
  lines <- readLines(path)
  get_v <- function(k) {
    l <- grep(paste0("^", k, "\t"), lines, value = TRUE)
    if (!length(l)) NA_real_ else as.numeric(strsplit(l[1], "\t")[[1]][2])
  }
  list(amp5  = get_v("amplitude_5p"), lam5  = get_v("lambda_5p"),
       base5 = get_v("baseline_5p"),  amp3  = get_v("amplitude_3p"),
       lam3  = get_v("lambda_3p"),    base3 = get_v("baseline_3p"))
}

load_smiley_excess <- function(path, fraction = "mod") {
  sm   <- read_tsv(path, show_col_types = FALSE)
  col5 <- paste0(fraction, "_5p")
  col3 <- paste0(fraction, "_3p")
  bg5  <- mean(sm[[col5]][10:15])
  bg3  <- mean(sm[[col3]][10:15])
  sm %>% filter(pos <= 15) %>%
    transmute(pos,
              rate_5p = (.data[[col5]] - bg5) * 100,
              rate_3p = (.data[[col3]] - bg3) * 100)
}

# ── Panel A: DNA damage — mirrored smiley plot ────────────────────────────────
# Classic layout: 5' C→T on left (pos 1→15), 3' G→A on right (pos 1→15 mirrored).
# Ancient damage (KapK-anc) is prominent; modern fractions (KapK-dec, KapK-mod)
# are shown at absolute scale as muted background lines.
plot_damage <- function(dmg_path, kapk_smiley, s17_smiley) {
  dmg <- load_damage_data(dmg_path)
  p   <- load_damage_params(dmg_path)

  # Smiley layout: pos=1 at OUTER edges, pos=15 near centre divider.
  # 5' C→T: xpos = -(16 - pos)  →  pos=1 → x=-15 (far left), pos=15 → x=-1
  # 3' G→A: xpos =  (16 - pos)  →  pos=1 → x=+15 (far right), pos=15 → x=+1
  anc_dat <- bind_rows(
    dmg %>% filter(end == "5\u2032 C\u2192T") %>%
      transmute(xpos = -(16L - pos), rate, ci_lo, ci_hi, series = "KapK-anc"),
    dmg %>% filter(end == "3\u2032 G\u2192A") %>%
      transmute(xpos =  (16L - pos), rate, ci_lo, ci_hi, series = "KapK-anc")
  )

  # Scale modern excess to match ancient at pos=1 (so shapes are visible)
  smiley_anc <- load_smiley_excess(kapk_smiley, "anc")
  anc_5p_1   <- filter(dmg, end == "5\u2032 C\u2192T", pos == 1)$rate * 100
  anc_3p_1   <- filter(dmg, end == "3\u2032 G\u2192A", pos == 1)$rate * 100
  scale_5p   <- anc_5p_1 / max(smiley_anc$rate_5p[1], 1e-3)
  scale_3p   <- anc_3p_1 / max(smiley_anc$rate_3p[1], 1e-3)

  make_mod_smiley <- function(smiley_path, label) {
    sm <- load_smiley_excess(smiley_path, "mod")
    bind_rows(
      sm %>% transmute(xpos = -(16L - pos),
                       rate = pmax(0, rate_5p * scale_5p) / 100,
                       ci_lo = NA_real_, ci_hi = NA_real_, series = label),
      sm %>% transmute(xpos =  (16L - pos),
                       rate = pmax(0, rate_3p * scale_3p) / 100,
                       ci_lo = NA_real_, ci_hi = NA_real_, series = label)
    )
  }

  mod_dec <- make_mod_smiley(kapk_smiley, "KapK-dec")
  mod_s17 <- make_mod_smiley(s17_smiley,  "KapK-mod")

  all_dat <- bind_rows(mod_dec, mod_s17, anc_dat) %>%
    mutate(series = factor(series, levels = c("KapK-anc", "KapK-dec", "KapK-mod")))

  # Fitted decay curves (5' left→right, 3' right→left)
  xs    <- seq(1, 15, 0.1)
  model <- bind_rows(
    tibble(xpos = -(16 - xs), rate = p$amp5 * exp(-p$lam5 * (xs - 1)) + p$base5, series = "KapK-anc"),
    tibble(xpos =  (16 - xs), rate = p$amp3 * exp(-p$lam3 * (xs - 1)) + p$base3, series = "KapK-anc")
  ) %>% mutate(series = factor(series, levels = levels(all_dat$series)))

  # Same palette as panel B
  ser_col <- c(
    "KapK-anc" = COL_ANCIENT,
    "KapK-dec" = COL_MODERN,
    "KapK-mod" = COL_RELATIVE
  )

  ggplot(all_dat, aes(x = xpos, y = rate, colour = series)) +
    # Subtle centre divider
    geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey85") +
    # CI ribbon: ancient only, very subtle
    geom_ribbon(data = filter(all_dat, series == "KapK-anc", !is.na(ci_lo)),
                aes(ymin = ci_lo, ymax = ci_hi, fill = series),
                alpha = 0.10, colour = NA, show.legend = FALSE) +
    # Modern reference: thin, muted, no points
    geom_line(data = filter(all_dat, series != "KapK-anc"),
              linewidth = 0.38, alpha = 0.60) +
    # Ancient: prominent line + points
    geom_line(data = filter(all_dat, series == "KapK-anc"),
              linewidth = 0.55) +
    geom_point(data = filter(all_dat, series == "KapK-anc"),
               size = 1.0, shape = 16) +
    # Structural end labels (larger) + mutation type below
    annotate("text", x = -13.5, y = 0.42, label = "5' end",
             size = 2.6, colour = "grey20", hjust = 0.5, fontface = "bold") +
    annotate("text", x =  13.5, y = 0.42, label = "3' end",
             size = 2.6, colour = "grey20", hjust = 0.5, fontface = "bold") +
    annotate("text", x = -13.5, y = 0.385, label = "C>T",
             size = 2.0, colour = "grey45", hjust = 0.5, fontface = "italic") +
    annotate("text", x =  13.5, y = 0.385, label = "G>A",
             size = 2.0, colour = "grey45", hjust = 0.5, fontface = "italic") +
    scale_colour_manual(values = ser_col, name = NULL,
                        breaks = c("KapK-anc", "KapK-dec", "KapK-mod")) +
    scale_fill_manual(values   = ser_col, guide = "none") +
    guides(colour = "none") +
    scale_x_continuous(
      breaks = c(-15, -11, -6, 6, 11, 15),
      labels = c("1", "5", "10", "10", "5", "1"),
      limits = c(-16, 16), expand = c(0, 0)
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 0.44),
      breaks = c(0, 0.1, 0.2, 0.3, 0.4)
    ) +
    coord_cartesian(clip = "off") +
    labs(x = "Position from read terminus", y = "Misincorporation frequency") +
    theme_nat() +
    theme(
      panel.grid.major   = element_blank(),
      panel.border       = element_blank(),
      axis.line          = element_line(linewidth = 0.35, colour = "grey30"),
      legend.position    = "none"
    )
}

# ── Panel B: Fragment lengths ──────────────────────────────────────────────────
plot_fragment_lengths <- function(kapk_path, s17_path) {
  load_hist <- function(path, cols) {
    read_tsv(path, show_col_types = FALSE) %>%
      mutate(pos = (length_min + length_max) / 2) %>%
      filter(pos >= 25, pos <= 600) %>%
      select(pos, all_of(cols))
  }

  # Raw read counts per fraction (for legend labels)
  kapk_raw <- read_tsv(kapk_path, show_col_types = FALSE)
  s17_raw  <- read_tsv(s17_path,  show_col_types = FALSE)
  n_anc <- sum(kapk_raw$ancient, na.rm = TRUE)
  n_mod <- sum(kapk_raw$modern,  na.rm = TRUE)
  n_s17 <- sum(s17_raw$modern,   na.rm = TRUE)

  fmt_n <- function(n) {
    if (n >= 1e6) sprintf("%.1fM", n / 1e6)
    else          sprintf("%.0fk", n / 1e3)
  }

  lbl_anc <- sprintf("KapK-anc (%s reads)", fmt_n(n_anc))
  lbl_mod <- sprintf("KapK-dec (%s reads)", fmt_n(n_mod))
  lbl_s17 <- sprintf("KapK-mod (%s reads)", fmt_n(n_s17))

  kapk <- load_hist(kapk_path, c("ancient", "modern"))
  s17  <- load_hist(s17_path, "modern") %>% rename(s17_modern = modern)

  frag <- kapk %>%
    left_join(s17, by = "pos") %>%
    rename(Ancient = ancient,
           `Kap Kbh. modern` = modern,
           `S17 modern`      = s17_modern) %>%
    pivot_longer(c(Ancient, `Kap Kbh. modern`, `S17 modern`),
                 names_to = "fraction", values_to = "count") %>%
    group_by(fraction) %>%
    mutate(freq = count / sum(count, na.rm = TRUE) * 100) %>%
    ungroup() %>%
    mutate(fraction = factor(fraction,
                             levels = c("Ancient", "Kap Kbh. modern", "S17 modern")))

  means <- frag %>%
    group_by(fraction) %>%
    summarise(mu = sum(pos * freq / 100, na.rm = TRUE), .groups = "drop") %>%
    mutate(label_y = c(9.8, 8.4, 7.0)[as.integer(fraction)])

  # Legend labels with read counts
  frac_lvls <- c(lbl_anc, lbl_mod, lbl_s17)
  frag <- frag %>%
    mutate(fraction_lbl = factor(
      case_when(fraction == "Ancient"          ~ lbl_anc,
                fraction == "Kap Kbh. modern"  ~ lbl_mod,
                fraction == "S17 modern"        ~ lbl_s17),
      levels = frac_lvls
    ))
  means <- means %>%
    mutate(fraction_lbl = factor(
      case_when(fraction == "Ancient"          ~ lbl_anc,
                fraction == "Kap Kbh. modern"  ~ lbl_mod,
                fraction == "S17 modern"        ~ lbl_s17),
      levels = frac_lvls
    ))

  frac_cols <- setNames(
    c(COL_ANCIENT, COL_MODERN, COL_RELATIVE),
    frac_lvls
  )
  # Simplified legend labels (same names as panel A, read counts in caption)
  simple_lvls <- c("KapK-anc", "KapK-dec", "S17")
  simple_cols <- setNames(c(COL_ANCIENT, COL_MODERN, COL_RELATIVE), simple_lvls)
  frag <- frag %>%
    mutate(series = factor(
      case_when(fraction == "Ancient"         ~ "KapK-anc",
                fraction == "Kap Kbh. modern" ~ "KapK-dec",
                fraction == "S17 modern"      ~ "S17"),
      levels = simple_lvls
    ))
  means <- means %>%
    mutate(series = factor(
      case_when(fraction == "Ancient"         ~ "KapK-anc",
                fraction == "Kap Kbh. modern" ~ "KapK-dec",
                fraction == "S17 modern"      ~ "S17"),
      levels = simple_lvls
    ))

  ggplot(frag, aes(x = pos, y = freq, colour = series, fill = series)) +
    # Subtle fills — legend hidden here, shown via geom_line only
    geom_area(alpha = 0.07, linewidth = 0, position = "identity", show.legend = FALSE) +
    # Clean lines on top — legend from these only
    geom_line(linewidth = 0.65, position = "identity") +
    # Mean lines: staircase segments stopping at each label's height
    geom_segment(data = means, aes(x = mu, xend = mu, y = 0, yend = label_y,
                                   colour = series),
                 linetype = "dashed", linewidth = 0.3, alpha = 0.45, show.legend = FALSE) +
    # Mean labels: staircase inside plot, one per fraction at different heights
    geom_text(data = means,
              aes(x = mu + 2, y = label_y,
                  label = sprintf("%.0f bp", mu), colour = series),
              size = 1.75, hjust = 0, vjust = 1, show.legend = FALSE) +
    scale_fill_manual(values   = simple_cols, name = NULL) +
    scale_colour_manual(values = simple_cols, name = NULL) +
    guides(
      fill   = "none",
      colour = guide_legend(override.aes = list(linewidth = 0.7, alpha = 1))
    ) +
    scale_x_continuous(breaks = seq(100, 600, 100), limits = c(25, 600),
                       expand = c(0.01, 0)) +
    scale_y_continuous(expand = c(0.01, 0), limits = c(0, 11.5)) +
    labs(x = "Fragment length (bp)", y = "Frequency (%)") +
    theme_nat() +
    theme(
      panel.grid.major   = element_blank(),
      panel.border       = element_blank(),
      axis.line          = element_line(linewidth = 0.35, colour = "grey30"),
      legend.position    = "top",
      legend.direction   = "horizontal",
      legend.background  = element_blank(),
      legend.key.size    = unit(7, "pt"),
      legend.spacing.x   = unit(3, "pt"),
      legend.text        = element_text(size = 6.5, family = FONT),
      legend.title       = element_blank()
    )
}

# ── Panel C: Focal subtree + BEAST2 tip annotations ────────────────────────────
plot_combined_tree <- function(iqtree_file, mcc_file, ani_file) {
  suppressMessages({

  # ── IQ-TREE focal subtree ──────────────────────────────────────────────────
  iqtree     <- read.tree(iqtree_file)
  focal_tips <- c("ancient_consensus", "s17_ancient_relative",
                  "NAY3300025461b7", "SOIL100000032",
                  "GCA054338395", "GCA054338475", "GCA963719035")
  focal_tips <- intersect(focal_tips, iqtree$tip.label)
  # Round-trip through newick to give treeio a clean edge matrix (avoids
  # "Invalid edge matrix for <phylo>" from keep.tip() artefacts)
  sub_iq     <- read.tree(text = write.tree(keep.tip(iqtree, focal_tips)))

  # ── BEAST2 terminal branch lengths (years → Ma) ───────────────────────────
  bt      <- read.nexus(mcc_file)
  N       <- Ntip(bt)
  te      <- which(bt$edge[, 2] <= N)
  tl      <- bt$tip.label[bt$edge[te, 2]]
  tbl     <- bt$edge.length[te] / 1e6
  get_bl  <- function(nm) { v <- tbl[tl == nm]; if (length(v)) v[1] else NA_real_ }

  anc_bl <- get_bl("ancient_consensus")
  mod_bl <- get_bl("modern_consensus")

  bl_str <- function(nm) {
    bl <- get_bl(nm)
    if (is.na(bl))    return("")
    if (bl >= 0.1)    sprintf("%.1f Ma", bl)
    else              sprintf("%.0f kyr", bl * 1e3)
  }

  # ── Tip metadata ──────────────────────────────────────────────────────────
  name_map <- c(
    ancient_consensus    = "KapK-anc",
    s17_ancient_relative = "KapK-mod",
    NAY3300025461b7      = "NAY3300025461",
    SOIL100000032        = "SOIL100000032",
    GCA054338395         = "GCA054338395",
    GCA054338475         = "GCA054338475",
    GCA963719035         = "GCA963719035"
  )

  # ── Biome metadata ────────────────────────────────────────────────────────
  meta <- tryCatch({
    raw <- read.table(META_FILE, sep = "\t", header = FALSE, quote = "",
                      comment.char = "#", stringsAsFactors = FALSE,
                      col.names = c("id", "completeness", "contamination",
                                    "strain_hetero", "classification",
                                    "ncbi_name", "sample_description"))
    raw %>% mutate(id = trimws(id))
  }, error = function(e) data.frame())

  # Classify biome from sample_description
  classify_biome <- function(desc) {
    desc <- tolower(desc)
    if (grepl("permafrost", desc))               "Permafrost"
    else if (grepl("peatland|peat|bog|fen|sphagnum|spruce", desc)) "Peatland"
    else if (grepl("tundra", desc))              "Tundra"
    else if (grepl("arctic|barrow|alaska", desc)) "Arctic peat"
    else if (grepl("freshwater|sediment.*fresh|fresh.*sediment|urban freshwater|natural freshwater|aquatic", desc)) "Freshwater"
    else if (grepl("sediment", desc))            "Sediment"
    else if (grepl("rice|wastewater|bioreactor", desc)) "Other"
    else                                          "Soil"
  }

  biome_from_id <- function(nm) {
    row <- meta[meta$id == nm, ]
    if (nrow(row) == 0) return(NA_character_)
    classify_biome(row$sample_description[1])
  }

  # Hardcoded focal biomes
  focal_biomes <- c(
    ancient_consensus    = "Kap K\u00f8benhavn",
    s17_ancient_relative = "Kap K\u00f8benhavn"
  )

  biome_cols <- c(
    "Kap K\u00f8benhavn" = COL_ANCIENT,
    "Permafrost"      = "#8AAFC2",
    "Arctic peat"     = "#6B8C6B",
    "Peatland"        = "#8BAD6A",
    "Tundra"          = "#9B8FB8",
    "Freshwater"      = "#6AA898",
    "Sediment"        = "#B09070",
    "Soil"            = "#C4A87A",
    "Other"           = "grey75"
  )

  # ANI from file (for non-focal tips)
  ani_df <- tryCatch(
    read_tsv(ani_file,
             col_names      = c("query", "ref", "ani", "frag_mapped", "frag_total"),
             show_col_types = FALSE) %>%
      filter(grepl("ancient_consensus", query)) %>%
      mutate(tip = sub("\\.fna$", "", basename(ref))),
    error = function(e) tibble()
  )
  get_ani <- function(nm) {
    v <- ani_df$ani[ani_df$tip == nm]
    if (length(v)) v[1] else NA_real_
  }

  # Build display label per tip
  # NOTE: Removed kyr/Ma terminal branch labels - confusing for readers

  # Now showing: SNP counts + ANI for focal tips, ANI only for references
  make_label <- function(nm) {
    disp <- name_map[nm]
    if (is.na(disp)) disp <- nm
    if (nm == "ancient_consensus") {
      # Ancient tip: just the name (dated to ~2 Ma via tip constraint)
      disp
    } else if (nm == "s17_ancient_relative") {
      sprintf("%s | %s SNPs | %.1f%% ANI",
              disp, format(SNP_S17, big.mark = ","), ANI_S17)
    } else {
      ani <- get_ani(nm)
      if (!is.na(ani)) sprintf("%s | %.1f%% ANI", disp, ani) else disp
    }
  }

  tip_meta <- tibble(
    label      = sub_iq$tip.label,
    group      = case_when(
      label == "ancient_consensus"    ~ "ancient",
      label == "s17_ancient_relative" ~ "relative",
      TRUE                            ~ "other"
    ),
    disp_label = sapply(label, make_label),
    tip_face   = ifelse(label %in% c("ancient_consensus", "s17_ancient_relative"),
                        "bold", "plain"),
    biome      = factor(
      ifelse(label %in% names(focal_biomes),
             focal_biomes[label],
             sapply(label, biome_from_id)),
      levels = c("Kap K\u00f8benhavn", "Tundra", "Permafrost", "Peatland", "Freshwater")
    )
  )

  # tip_cols/shapes/sizes no longer used — all tips rendered via biome geom_point

  # ── Branch support (internal nodes only, exclude root) ────────────────────
  # Convert to treedata to avoid ggtree 3.6/treeio "Invalid edge matrix" warnings
  base_p  <- ggtree(treeio::as.treedata(sub_iq), color = "grey35", linewidth = 0.35)
  td      <- base_p$data
  n_tips  <- Ntip(sub_iq)
  max_x   <- max(td$x, na.rm = TRUE)

  support_dat <- td %>%
    filter(!isTip, !is.na(label), nzchar(label),
           node != n_tips + 1, grepl("/", label)) %>%
    mutate(
      # Parse SH-aLRT / UFBoot
      sh   = as.numeric(sub("/.*", "", label)),
      uf   = as.numeric(sub(".*/", "", label)),
      # Only show well-supported nodes
      show = (sh >= 80 | uf >= 95)
    ) %>%
    filter(show)

  # ── Stasis clade highlight ────────────────────────────────────────────────
  stasis_tips <- intersect(c("ancient_consensus", "s17_ancient_relative"),
                           sub_iq$tip.label)
  mrca_stasis <- getMRCA(sub_iq, stasis_tips)
  stasis_node <- td %>% filter(node == mrca_stasis)
  stasis_tips_y <- td %>% filter(isTip, label %in% stasis_tips)

  y_lo  <- min(stasis_tips_y$y) - 0.45
  y_hi  <- max(stasis_tips_y$y) + 0.45

  # Subtle amber shading behind KapK clade
  shade_dat <- data.frame(
    xmin = pmax(0, stasis_node$x - max_x * 0.02),
    xmax = max_x * 1.01,
    ymin = y_lo,
    ymax = y_hi
  )

  # "~2 Ma" label: above clade, centred on MRCA x
  mrca_label_dat <- data.frame(
    x     = stasis_node$x,
    y     = max(stasis_tips_y$y) + 0.55,
    label = "~2 Ma"
  )

  # Arrow from "~2 Ma" label down to MRCA node
  arrow_dat <- data.frame(
    x    = stasis_node$x,
    xend = stasis_node$x,
    y    = max(stasis_tips_y$y) + 0.40,
    yend = stasis_node$y + 0.13
  )

  # Rate suppression annotation — 82× strict lower bound, 244× 2-clock conditional
  rate_ann_dat <- data.frame(
    x     = stasis_node$x,
    y     = max(stasis_tips_y$y) + 0.98,
    label = "82\u2013244\u00d7 rate suppression"
  )

  biome_x_offset <- max_x * 0.02

  # ── Assemble ──────────────────────────────────────────────────────────────
  base_p %<+% tip_meta +
    # All tips as biome-coloured squares — consistent shape across KapK and references
    geom_point(aes(x = x + biome_x_offset, fill = biome),
               data = function(d) d[d$isTip & !is.na(d$isTip), ],
               shape = 22, size = 2.0, colour = "black", stroke = 0.3,
               show.legend = TRUE) +
    geom_tiplab(aes(label = disp_label, fontface = tip_face),
                colour = "grey15", size = 2.5, offset = max_x * 0.05, hjust = 0,
                family = FONT) +
    # Branch support: darker for legibility at journal scale
    geom_text(data = support_dat,
              aes(x = x, y = y, label = label),
              size = 1.8, hjust = 1.1, vjust = -0.5,
              colour = "grey35", family = FONT, inherit.aes = FALSE) +
    scale_fill_manual(values   = biome_cols, name = NULL,
                      na.value = "grey80",
                      guide = guide_legend(
                        override.aes = list(size = 2.5, shape = 22,
                                            colour = "black", stroke = 0.3)
                      )) +
    scale_x_continuous(expand  = expansion(mult = c(0.02, 1.2))) +
    coord_cartesian(clip = "off") +
    labs(x = "Nucleotide substitutions per site", y = NULL) +
    theme_tree2(base_size = 7) +
    theme(
      text         = element_text(family = FONT),
      axis.text.x  = element_text(size = 6, colour = "grey20"),
      axis.title.x = element_text(size = 7, colour = "grey10"),
      axis.line.x  = element_line(linewidth = 0.3, colour = "grey35"),
      axis.ticks.x = element_line(linewidth = 0.3, colour = "grey35"),
      plot.margin  = margin(3, 8, 2, 10, "pt"),
      legend.position      = "top",
      legend.direction     = "horizontal",
      legend.background    = element_blank(),
      legend.text          = element_text(size = 6.5, family = FONT, margin = margin(l = 1, r = 3, unit = "pt")),
      legend.title         = element_blank(),
      legend.key.size      = unit(5, "pt"),
      legend.key.width     = unit(6, "pt"),
      legend.spacing.x     = unit(0.5, "pt"),
      plot.tag             = element_text(size = 11, face = "bold", colour = "grey5")
    )

  }) # end suppressMessages
}

# ── BEAST2 log helpers ─────────────────────────────────────────────────────────
read_beast_log <- function(path, burnin = 0.10) {
  lines <- readLines(path)
  hdr_i <- grep("^Sample\t", lines)
  if (!length(hdr_i)) stop("No 'Sample' header in ", path)
  hdr   <- strsplit(lines[hdr_i[1]], "\t")[[1]]
  data_lines <- lines[(hdr_i[1] + 1):length(lines)]
  data_lines <- data_lines[!grepl("^#", data_lines) & nzchar(data_lines)]
  dat   <- read.table(text = paste(data_lines, collapse = "\n"),
                      header = FALSE, col.names = hdr, sep = "\t",
                      stringsAsFactors = FALSE)
  n_drop <- ceiling(nrow(dat) * burnin)
  dat[-(1:n_drop), , drop = FALSE]
}

hpd95 <- function(x) {
  x  <- sort(x[!is.na(x)])
  n  <- length(x)
  w  <- floor(0.95 * n)
  gap <- x[(w + 1):n] - x[1:(n - w)]
  i  <- which.min(gap)
  c(lo = x[i], hi = x[i + w])
}

parse_ps_marginal_L <- function(ps_run_dir) {
  path <- file.path(PS_DIR, ps_run_dir, "beast_ps.out")
  if (!file.exists(path)) return(NA_real_)
  lines <- readLines(path)
  hit   <- grep("marginal L estimate =", lines, value = TRUE)
  if (!length(hit)) return(NA_real_)
  as.numeric(sub(".*marginal L estimate = ", "", tail(hit, 1)))
}

compute_bayes_factors <- function() {
  # ── PRIMARY: model comparison at 2 Ma fixed, pos12 only ──────────────────
  # Reference = strict clock; alternative = UCLD relaxed clock.
  # Both use pos12-only data (filter 1-4406\3,2-4406\3) for valid comparison.
  #
  # NOTE: ps_2clock_stasis (FlexibleLocalClockModel) is NOT included here.
  # That run has catastrophic MCMC mixing failure (treeLikelihood ESS=8),
  # with posterior treeLik=-106,442 vs strict -61,118 on identical data — a
  # degenerate region caused by the near-zero KapK rate inflating clockRateBg
  # 3×, trapping the tree topology. The 244× rate ratio from the 2clock MCMC
  # is reported separately as a conditional posterior estimate (below), not as
  # a formal Bayes factor.
  primary_models <- c("strict_2Ma_pos12", "ucld_2Ma_pos12")
  primary_dirs   <- paste0("ps_", primary_models)
  primary_lml    <- setNames(vapply(primary_dirs, parse_ps_marginal_L, numeric(1)),
                              primary_models)

  ps_tbl <- function(labels, dirs, lml, ref_label) {
    ref <- lml[ref_label]
    data.frame(
      model      = labels,
      log_ML     = lml,
      two_lnBF   = 2 * (lml - ref),   # positive = better than reference
      steps_done = vapply(dirs, function(d) {
        f <- file.path(PS_DIR, d, "beast_ps.out")
        if (!file.exists(f)) return(0L)
        sum(grepl("Finished step", readLines(f)))
      }, integer(1)),
      row.names  = NULL
    )
  }

  primary_tbl <- ps_tbl(primary_models, primary_dirs, primary_lml, "strict_2Ma_pos12")

  message("\n── Path sampling: PRIMARY — model comparison at 2 Ma (pos12 only) ──────────")
  message("  Reference: strict clock. Positive 2×ln BF = UCLD better than strict.")
  message("  (2clock FlexibleLocalClock excluded — MCMC ESS=8, degenerate posterior)")
  message(sprintf("  %-20s  %14s  %14s  %7s", "Model", "ln ML", "2×ln BF vs strict", "Steps"))
  for (i in seq_len(nrow(primary_tbl))) {
    r     <- primary_tbl[i, ]
    ml_s  <- if (is.na(r$log_ML))  sprintf("%14s", "–") else sprintf("%14.2f", r$log_ML)
    bf_s  <- if (r$model == "strict_2Ma_pos12" || is.na(r$two_lnBF))
               sprintf("%14s", "0 (ref)")
             else
               sprintf("%14.2f", r$two_lnBF)
    message(sprintf("  %-20s  %s  %s  %4d/64", r$model, ml_s, bf_s, r$steps_done))
  }
  message("  2×ln BF > 10 = decisive; > 6 = strong; > 2 = substantial\n")

  # ── SECONDARY: fixed-age sensitivity under strict clock ──────────────────
  # Interpretation: NOT "which age is correct" — younger ages fit strict clock
  # better because the anomalous ancient branch is a smaller outlier. This
  # confirms rate heterogeneity, not a preferred geological age.
  sec_ages <- c("0Ma", "0.5Ma", "1Ma", "2Ma")
  sec_dirs <- paste0("ps_", sec_ages)
  sec_lml  <- setNames(vapply(sec_dirs, parse_ps_marginal_L, numeric(1)), sec_ages)
  sec_tbl  <- ps_tbl(sec_ages, sec_dirs, sec_lml, "2Ma")

  message("── Path sampling: SECONDARY — fixed-age sensitivity (strict clock only) ──────")
  message("  Reference: 2 Ma. Positive 2×ln BF = that age better than 2 Ma under strict clock.")
  message("  (Expected: younger ages score better — confirms rate outlier, not age preference)")
  message(sprintf("  %-8s  %14s  %14s  %7s", "Age", "ln ML", "2×ln BF vs 2Ma", "Steps"))
  for (i in seq_len(nrow(sec_tbl))) {
    r    <- sec_tbl[i, ]
    ml_s <- if (is.na(r$log_ML)) sprintf("%14s", "–") else sprintf("%14.2f", r$log_ML)
    bf_s <- if (r$model == "2Ma" || is.na(r$two_lnBF))
              sprintf("%14s", "0 (ref)")
            else
              sprintf("%14.2f", r$two_lnBF)
    message(sprintf("  %-8s  %s  %s  %4d/64", r$model, ml_s, bf_s, r$steps_done))
  }
  message("")

  invisible(list(primary = primary_tbl, secondary = sec_tbl))
}

# ── Evolutionary stasis statistics ────────────────────────────────────────────
compute_stasis_stats <- function() {
  message("\n── Evolutionary stasis statistics ──────────────────────────────────────────")
  bg_med <- NA_real_

  # 1. Background clock rate — strict clock tipdating (2 Ma fixed)
  if (file.exists(LOG_TIPDATING)) {
    td     <- read_beast_log(LOG_TIPDATING)
    bg_med <- median(td$clockRate)
    bg_hpd <- hpd95(td$clockRate)
    message(sprintf("Background rate (tipdating):  %.3g [%.3g\u2013%.3g] sub/site/yr  (n=%d post-burnin)",
                    bg_med, bg_hpd["lo"], bg_hpd["hi"], nrow(td)))
  }

  # 2. Ancient branch molecular duration from tipdating MCC tree
  if (file.exists(MCC_FILE)) {
    bt     <- read.nexus(MCC_FILE)
    N      <- Ntip(bt)
    te     <- which(bt$edge[, 2] <= N)
    tl     <- bt$tip.label[bt$edge[te, 2]]
    tbl_yr <- bt$edge.length[te]
    anc_yr <- tbl_yr[tl == "ancient_consensus"][1]
    fold   <- GEO_AGE_YR / anc_yr
    eff    <- if (!is.na(bg_med)) bg_med * (anc_yr / GEO_AGE_YR) else NA_real_
    message(sprintf("Ancient branch (MCC):         %.0f yr  (%.4f Ma)", anc_yr, anc_yr / 1e6))
    message(sprintf("Geological age (tip prior):   %.0f yr  (%.1f Ma)", GEO_AGE_YR, GEO_AGE_YR / 1e6))
    message(sprintf("Fold stasis (geo/branch):     %.0f\u00d7", fold))
    if (!is.na(eff))
      message(sprintf("Effective ancient rate:       %.3g sub/site/yr", eff))
  }

  # 3. Freedate — tip date as free parameter (two replicate runs)
  message("")
  for (rep_info in list(list(LOG_FREEDATE_R3, "rep3"), list(LOG_FREEDATE_R4, "rep4"))) {
    log_path <- rep_info[[1]]
    rep_lbl  <- rep_info[[2]]
    if (!file.exists(log_path)) next
    fd <- read_beast_log(log_path)
    # Column: height(ancient_consensus) — use make.names-safe accessor
    h_col <- "height.ancient_consensus."   # read.table converts () to .
    if (!h_col %in% names(fd)) {
      # Try the literal name (data read via scan/col.names preserves it)
      h_col <- grep("ancient_consensus", names(fd), value = TRUE)[1]
    }
    if (!is.na(h_col) && h_col %in% names(fd)) {
      h_vec <- fd[[h_col]]
      h_med <- median(h_vec)
      h_hpd <- hpd95(h_vec)
      message(sprintf("Freedate %s tip age:          %.4f [%.4f\u2013%.4f] Ma  (n=%d post-burnin)",
                      rep_lbl, h_med / 1e6, h_hpd["lo"] / 1e6, h_hpd["hi"] / 1e6, nrow(fd)))
    }
  }

  # 4. Two-clock stasis model — separate background and stasis rates
  message("")
  if (file.exists(LOG_2CLOCK)) {
    cl      <- read_beast_log(LOG_2CLOCK)
    bg2_med <- median(cl$clockRateBg); bg2_hpd <- hpd95(cl$clockRateBg)
    mod_med <- median(cl$clockRateMod); mod_hpd <- hpd95(cl$clockRateMod)
    ratio_vec  <- cl$clockRateBg / cl$clockRateMod
    fold2_med  <- median(ratio_vec)
    fold2_hpd  <- hpd95(ratio_vec)
    p_gt10     <- mean(ratio_vec > 10)
    p_gt100    <- mean(ratio_vec > 100)
    message(sprintf("2clock clockRateBg:           %.3g [%.3g\u2013%.3g] sub/site/yr  (n=%d post-burnin)",
                    bg2_med, bg2_hpd["lo"], bg2_hpd["hi"], nrow(cl)))
    message(sprintf("2clock clockRateMod (stasis): %.3g [%.3g\u2013%.3g] sub/site/yr",
                    mod_med, mod_hpd["lo"], mod_hpd["hi"]))
    message(sprintf("Rate ratio posterior:         median=%.0f\u00d7 [%.0f\u2013%.0f\u00d7] (95%% HPD)",
                    fold2_med, fold2_hpd["lo"], fold2_hpd["hi"]))
    message(sprintf("  P(ratio > 10\u00d7)  = %.3f", p_gt10))
    message(sprintf("  P(ratio > 100\u00d7) = %.3f", p_gt100))
  }

  # 5. UCLD relaxed clock — rate heterogeneity diagnostic
  message("")
  if (file.exists(LOG_UCLD)) {
    uc      <- read_beast_log(LOG_UCLD)
    um_med  <- median(uc$ucldMean);  um_hpd  <- hpd95(uc$ucldMean)
    usd_med <- median(uc$ucldStdev); usd_hpd <- hpd95(uc$ucldStdev)
    message(sprintf("UCLD mean rate:               %.3g [%.3g\u2013%.3g] sub/site/yr  (n=%d post-burnin)",
                    um_med, um_hpd["lo"], um_hpd["hi"], nrow(uc)))
    message(sprintf("UCLD SD:                      %.3f [%.3f\u2013%.3f]  (>>1 = strong rate heterogeneity)",
                    usd_med, usd_hpd["lo"], usd_hpd["hi"]))
  }

  # 6. Path sampling Bayes factors (fixed-age sensitivity)
  compute_bayes_factors()

  # 7. SNP context
  message(sprintf("SNPs ancient vs KapK-dec:     %d", SNP_MOD))
  message(sprintf("SNPs ancient vs KapK-mod:     %d", SNP_S17))
  message("────────────────────────────────────────────────────────────────────────────\n")
}

# ── Main ───────────────────────────────────────────────────────────────────────
main <- function() {
  dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

  compute_stasis_stats()

  message("Building panels...")
  pA <- plot_damage(DMG_FILE, SMILEY_FILE, S17_SMILEY)
  pB <- plot_fragment_lengths(FRAG_FILE, S17_FRAG)
  pC <- plot_combined_tree(IQTREE_FILE, MCC_FILE, ANI_FILE)

  message("Assembling figure...")
  fig <- ((pA / pB) | pC) +
    plot_layout(widths = c(1, 1.5)) +
    plot_annotation(
      tag_levels = "A",
      theme = theme(
        plot.tag = element_text(size = 11, face = "bold", colour = "grey5")
      )
    )

  out_pdf <- file.path(OUTDIR, "fig7-methanoflorens-stasis.pdf")
  out_png <- file.path(OUTDIR, "fig7-methanoflorens-stasis.png")

  ggsave(out_pdf, fig, width = FIG_W, height = FIG_H,
         units = "mm", device = cairo_pdf)
  message("Saved: ", out_pdf)
  ggsave(out_png, fig, width = FIG_W, height = FIG_H,
         units = "mm", dpi = 300, bg = "white", device = ragg::agg_png)
  message("Saved: ", out_png)
}

main()
