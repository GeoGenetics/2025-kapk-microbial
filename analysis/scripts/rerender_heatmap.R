suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(janitor)
  library(showtext)
  library(ggh4x)
})

setwd("/maps/projects/fernandezguerra/apps/repos/2025-kapk-microbial/analysis")
source("./libs/lib.R")
showtext_auto()
source("./scripts/07a--agp-functional.R")

# Replicate the metadata construction from main()
cdata_path   <- "./data/cdata/KapK-cdata-manuscript-20221211.tsv"
nobloom_path <- "./data/cdata/KapK-cdata-manuscript-20221211-labels-nobloom.tsv"

kapk_cdata    <- read_tsv(cdata_path,   show_col_types = FALSE)
label_nobloom <- read_tsv(nobloom_path, show_col_types = FALSE)
kapk_cdata    <- kapk_cdata %>% filter(figure_names %in% label_nobloom$label)
kapk_cdata    <- kapk_cdata %>%
    mutate(short_label = paste(site, sub(".*-", "", figure_names), sep = "_"))

sample_meta <- as.data.table(kapk_cdata)[,
    .(sample = substr(label, 1, 10), label, short_label, member_unit, site_rnk)]

sample_map   <- fread("./data/cdata/paper_sample_mapping.tsv")
kc_slim      <- sample_meta[, .(file_label = sample, member_unit, site_rnk)]
meta_ordered <- merge(sample_map, kc_slim, by = "file_label")
sample_order <- unique(meta_ordered[
    order(match(member_unit, c("B1", "B2", "B3")), site_rnk)
]$short_label)
cat("sample_order:", paste(sample_order, collapse = ", "), "\n")

processed <- readRDS("./results/functional_agp/heatmap_processed.rds")
processed$all$label <- factor(as.character(processed$all$label), levels = sample_order)
processed$dmg$label <- factor(as.character(processed$dmg$label), levels = sample_order)

p <- plot_functional_heatmap(processed$all, processed$dmg)

out_dir <- "./results/functional_agp"
ggsave(file.path(out_dir, "fig_functional_heatmap.pdf"), p,
       width = 200, height = 140, units = "mm")
ggsave(file.path(out_dir, "fig_functional_heatmap.png"), p,
       width = 200, height = 140, units = "mm", dpi = 200, bg = "white")
cat("Saved: fig_functional_heatmap.pdf / .png\n")
