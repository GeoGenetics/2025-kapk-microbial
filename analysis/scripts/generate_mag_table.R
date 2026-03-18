#!/usr/bin/env Rscript
# Generate comprehensive MAG quality table for the KapK paper.
# Sources: AMBER paper_run/ outputs (CheckM2, GUNC, GTDB-Tk, damage-per-bin)
#
# MIMAG quality thresholds (AMBER workflow / Parks et al. 2022):
#   HQ: completeness ≥ 90%,  contamination < 5%
#   MQ: completeness ≥ 50%,  contamination < 10%
#   LQ: everything else
# Note: strict MIMAG HQ also requires 5S/16S/23S rRNA + ≥18 tRNAs; this is
# rarely met by aDNA MAGs due to fragmentation and is not applied here.

suppressPackageStartupMessages({
    library(tidyverse)
    library(openxlsx)
    library(janitor)
})

PAPER_RUN <- "/projects/caeg/people/kbd606/scratch/kapk-assm/amber/paper_run"

# ── 1. CheckM2 ────────────────────────────────────────────────────────────────

checkm <- read_tsv(
    file.path(PAPER_RUN, "checkm2_final", "quality_report.tsv"),
    show_col_types = FALSE
) |>
    select(
        mag               = Name,
        completeness      = Completeness,
        contamination     = Contamination,
        completeness_model = Completeness_Model_Used,
        genome_size       = Genome_Size,
        gc_content        = GC_Content,
        contig_n50        = Contig_N50,
        max_contig_length = Max_Contig_Length,
        total_contigs     = Total_Contigs,
        total_cds         = Total_Coding_Sequences,
        coding_density    = Coding_Density
    )

# ── 2. GUNC ───────────────────────────────────────────────────────────────────

gunc <- read_tsv(
    file.path(PAPER_RUN, "gunc", "GUNC.progenomes_2.1.maxCSS_level.tsv"),
    show_col_types = FALSE
) |>
    select(
        mag                         = genome,
        gunc_taxonomic_level        = taxonomic_level,
        gunc_css                    = clade_separation_score,
        gunc_contamination_portion  = contamination_portion,
        gunc_reference_rep_score    = reference_representation_score,
        gunc_pass                   = pass.GUNC
    )

# ── 3. GTDB-Tk ────────────────────────────────────────────────────────────────

parse_gtdbtk <- function(path) {
    if (!file.exists(path)) return(tibble())
    df <- read_tsv(path, show_col_types = FALSE)
    if (nrow(df) == 0) return(tibble())
    df |>
        select(mag = user_genome, gtdb_classification = classification) |>
        separate(gtdb_classification,
                 into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
                 sep  = ";",
                 fill = "right",
                 remove = FALSE) |>
        mutate(across(domain:species, ~ sub("^[a-z]__", "", .x)))
}

gtdbtk <- bind_rows(
    parse_gtdbtk(file.path(PAPER_RUN, "gtdbtk", "gtdbtk.bac120.summary.tsv")),
    parse_gtdbtk(file.path(PAPER_RUN, "gtdbtk", "gtdbtk.ar53.summary.tsv"))
)

# ── 4. Damage ─────────────────────────────────────────────────────────────────

damage <- read_tsv(
    file.path(PAPER_RUN, "publication", "resolve_damage_per_bin.tsv"),
    show_col_types = FALSE
) |>
    select(
        mag          = bin,
        p_ancient,
        d_max_5prime = lambda5,
        d_max_3prime = lambda3,
        ct_1p,
        ga_1p,
        frag_mean,
        damage_class
    )

# ── 5. Merge + quality classification ─────────────────────────────────────────

mag_table <- checkm |>
    left_join(gunc,   by = "mag") |>
    left_join(gtdbtk, by = "mag") |>
    left_join(damage, by = "mag") |>
    mutate(
        quality = case_when(
            completeness >= 90 & contamination < 5  ~ "HQ",
            completeness >= 50 & contamination < 10 ~ "MQ",
            TRUE                                     ~ "LQ"
        ),
        # Flag MQ/HQ bins that fail GUNC (contamination may be overestimated)
        gunc_fail_note = if_else(
            quality %in% c("HQ", "MQ") & !is.na(gunc_pass) & !gunc_pass,
            "GUNC fail", NA_character_
        )
    ) |>
    # Round numeric columns
    mutate(
        across(c(completeness, contamination, gc_content, coding_density,
                 gunc_css, gunc_contamination_portion, gunc_reference_rep_score,
                 p_ancient, d_max_5prime, d_max_3prime, ct_1p, ga_1p, frag_mean),
               ~ round(.x, 4))
    ) |>
    arrange(factor(quality, levels = c("HQ", "MQ", "LQ")), desc(completeness))

# ── 6. Print summary ──────────────────────────────────────────────────────────

cat("\n=== MAG Quality Summary ===\n")
cat("Total MAGs:", nrow(mag_table), "\n")
print(table(mag_table$quality))
cat("\nDamage class breakdown (HQ+MQ only):\n")
print(table(mag_table |> filter(quality %in% c("HQ","MQ")) |> pull(damage_class)))

cat("\nHQ bins:\n")
mag_table |>
    filter(quality == "HQ") |>
    select(mag, completeness, contamination, genus, species, damage_class, gunc_pass) |>
    print(n = Inf)

# ── 7. Save outputs ───────────────────────────────────────────────────────────

out_dir <- file.path(
    "/maps/projects/fernandezguerra/apps/repos/2025-kapk-microbial/analysis",
    "results", "mag_quality"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tsv_path <- file.path(out_dir, "mag_quality_table.tsv")
write_tsv(mag_table, tsv_path)
cat("\nSaved TSV:", tsv_path, "\n")

# Supplementary Excel
wb <- createWorkbook()
addWorksheet(wb, "MAG quality")

# Header style
hs <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2",
                  border = "Bottom", wrapText = FALSE)

writeData(wb, "MAG quality", mag_table |> clean_names(case = "sentence"),
          headerStyle = hs)

setColWidths(wb, "MAG quality", cols = 1:ncol(mag_table), widths = "auto")

saveWorkbook(wb, "../supp-tab-v2/sup_table_9.xlsx", overwrite = TRUE)
cat("Saved sup_table_9.xlsx\n")
cat("  HQ:", sum(mag_table$quality == "HQ"),
    "| MQ:", sum(mag_table$quality == "MQ"),
    "| LQ:", sum(mag_table$quality == "LQ"), "\n")
