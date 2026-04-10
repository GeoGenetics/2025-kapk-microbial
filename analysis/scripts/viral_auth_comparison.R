suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(patchwork)
})

setwd("/maps/projects/fernandezguerra/apps/repos/2025-kapk-microbial/analysis")

# Load viral-auth per-sample results
ps <- fread("results/virome_agp/viral_auth.per_sample.tsv", showProgress = FALSE) |>
    as_tibble()

# Old rank filter
rank_filter <- ps |>
    mutate(rank = case_when(
        prot_coverage >= 0.25 & n_genes > 1                          ~ "green",
        prot_coverage >= 0.20 & prot_coverage < 0.25 & n_genes > 1  ~ "yellow",
        n_genes == 1 & prot_coverage >= 0.20                         ~ "grey",
        TRUE                                                          ~ "red"
    ))

old_pass <- rank_filter |> filter(rank != "red")
new_pass <- ps |> filter(qvalue <= 0.05)

cat("Old rank filter (ref×sample):", nrow(old_pass), "\n")
cat("Viral-auth significant:       ", nrow(new_pass), "\n")
cat("Old unique refs:", n_distinct(old_pass$reference), "\n")
cat("New unique refs:", n_distinct(new_pass$reference), "\n")

# Label each (ref, sample) pair by which filter it passes
ps_labeled <- ps |>
    left_join(rank_filter |> select(sample, reference, rank), by = c("sample", "reference")) |>
    mutate(
        pass_old = rank %in% c("green", "yellow", "grey"),
        pass_new = qvalue <= 0.05,
        category = case_when(
            pass_old & pass_new  ~ "Both",
            pass_old & !pass_new ~ "Rank filter only\n(removed by viral-auth)",
            !pass_old & pass_new ~ "Viral-auth only\n(missed by rank filter)",
            TRUE                 ~ "Neither"
        )
    )

cat_colors <- c(
    "Both"                             = "#2166AC",
    "Rank filter only\n(removed by viral-auth)" = "#D1E5F0",
    "Viral-auth only\n(missed by rank filter)"  = "#B2182B",
    "Neither"                          = "grey88"
)

# ── Panel A: prot_coverage vs n_cds, coloured by category ─────────────────
p_scatter <- ps_labeled |>
    filter(category != "Neither") |>
    mutate(category = factor(category, levels = names(cat_colors))) |>
    ggplot(aes(x = n_cds, y = prot_coverage, colour = category)) +
    geom_point(alpha = 0.3, size = 0.6) +
    geom_hline(yintercept = 0.20, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    annotate("text", x = 2500, y = 0.21, label = "20% threshold",
             size = 2.5, colour = "grey40", hjust = 1) +
    scale_colour_manual(values = cat_colors, name = NULL) +
    scale_x_log10(labels = scales::label_comma()) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1)) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          legend.position  = "bottom",
          legend.text      = element_text(size = 7)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
    labs(x = "Reference genome size (CDS)", y = "Proteome coverage",
         title = "A  Proteome coverage vs genome size")

# ── Panel B: n_genes distribution by category ──────────────────────────────
p_ngenes <- ps_labeled |>
    filter(category != "Neither", n_genes <= 20) |>
    mutate(category = factor(category, levels = names(cat_colors))) |>
    count(category, n_genes) |>
    group_by(category) |>
    mutate(prop = n / sum(n)) |>
    ggplot(aes(x = n_genes, y = prop, fill = category)) +
    geom_col(width = 0.8, show.legend = FALSE) +
    facet_wrap(~category, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = cat_colors) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(breaks = 1:20) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor  = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background   = element_blank(),
          strip.text         = element_text(size = 7)) +
    labs(x = "Detected genes (n_genes)", y = "Fraction of detections",
         title = "B  Gene detection count distribution")

# ── Panel C: Poisson λ vs n_genes, showing significance boundary ───────────
p_poisson <- ps_labeled |>
    filter(category != "Neither") |>
    mutate(category = factor(category, levels = names(cat_colors))) |>
    ggplot(aes(x = lambda, y = n_genes, colour = category)) +
    geom_point(alpha = 0.25, size = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey30", linewidth = 0.4) +
    annotate("text", x = 3, y = 3.4, label = "n_genes = lambda (null expectation)",
             size = 2.5, colour = "grey30", angle = 30) +
    scale_colour_manual(values = cat_colors, name = NULL) +
    scale_x_log10() + scale_y_log10() +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          legend.position  = "none") +
    labs(x = "Poisson null (lambda = n_cds x p_f)",
         y = "Observed n_genes",
         title = "C  Observed vs expected genes under EM noise floor")

# ── Panel D: summary counts ────────────────────────────────────────────────
summary_df <- tibble(
    filter   = c("Rank filter\n(20% threshold)", "Viral-auth\n(Poisson FDR 5%)"),
    n_unique = c(n_distinct(old_pass$reference), n_distinct(new_pass$reference))
)

p_summary <- ggplot(summary_df, aes(x = filter, y = n_unique, fill = filter)) +
    geom_col(width = 0.55, colour = "black", linewidth = 0.3) +
    geom_text(aes(label = scales::comma(n_unique)), vjust = -0.4, size = 3) +
    scale_fill_manual(values = c("#D1E5F0", "#2166AC")) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.12))) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor  = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position    = "none") +
    labs(x = NULL, y = "Unique viral references",
         title = "D  Total detected references")

# ── Assemble ───────────────────────────────────────────────────────────────
fig <- (p_scatter | p_ngenes) / (p_poisson | p_summary) +
    plot_layout(heights = c(1.4, 1))

out <- "results/virome_agp/viral_auth_comparison.pdf"
ggsave(out, fig, width = 10, height = 8)
cat("Saved:", out, "\n")

# Also PNG for quick view
ggsave(sub(".pdf", ".png", out), fig, width = 10, height = 8, dpi = 150)
cat("Saved:", sub(".pdf", ".png", out), "\n")

# Print breakdown
cat("\n--- Category breakdown (ref×sample pairs) ---\n")
ps_labeled |>
    filter(category != "Neither") |>
    count(category, sort = TRUE) |>
    print()

cat("\n--- Singleton (n_genes=1) fraction by category ---\n")
ps_labeled |>
    filter(category != "Neither") |>
    group_by(category) |>
    summarise(singleton_frac = mean(n_genes == 1), n = n()) |>
    print()
