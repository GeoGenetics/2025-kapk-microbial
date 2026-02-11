library(tidyverse)
library(purrr)
library(data.table)
library(pbmcapply)
library(unixtools)
library(showtext)
showtext_auto()

set.tempdir("/maps/projects/fernandezguerra/scratch/sandbox/tmp")
data_filt_folder <- "/maps/projects/lundbeck/people/kbd606/projects/kapk/results/taxonomic-profiling-filtering-benchmark-20221202/standard"
data_dmg_folder <- "/maps/projects/lundbeck/people/kbd606/projects/kapk/results/taxonomic-profiling-dmg-benchmark-20221202/standard/"
# Read in the data
k <- c(50, 100, 250, 500, 750, 1000)
ani <- c(90, 91, 92, 93, 94, 95, 96)

k_colors <- c(
    "k50" = "#1f77b4",
    "k100" = "#ff7f0e",
    "k250" = "#2ca02c",
    "k500" = "#d62728",
    "k1000" = "#9467bd",
    "k750" = "#8c564b"
)

# read taxonomic annotations
tax_info <- read_tsv("./data/taxonomy/db/hires-organelles-viruses.tax.tsv",
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
            "species"
        )
    )

rhc_files <- list.files(
    path = data_filt_folder,
    pattern = "*dedup_read-hits-count.tsv.gz",
    full.names = TRUE,
    recursive = TRUE
)
sf_files <- list.files(
    path = data_filt_folder,
    pattern = "*stats-filtered.tsv.gz",
    full.names = TRUE,
    recursive = TRUE
)

dmg_files <- list.files(
    path = data_dmg_folder,
    pattern = "*weight-1.csv.gz",
    full.names = TRUE,
    recursive = TRUE
)

dmg_files_lca <- dmg_files[grepl("lca", dmg_files)]
dmg_files_global <- dmg_files[grepl("global", dmg_files)]
dmg_files_local <- dmg_files[grepl("local", dmg_files)]

file <- dmg_files[1]
# Read in the files
read_files <- function(file, file_type = "rc", threads = 4) {
    if (file_type %in% c("sf", "rc")) {
        k <- str_match(
            string = file,
            pattern = ".*/k(\\d.+)/(\\d+)/(\\S{10})\\..*"
        )[1, 2][1]
        ani <- str_match(
            string = file,
            pattern = ".*/k(\\d.+)/(\\d+)/(\\S{10})\\..*"
        )[1, 3][1]
        label <- str_match(
            string = file,
            pattern = ".*/k(\\d.+)/(\\d+)/(\\S{10})\\..*"
        )[1, 4][1]
    } else if (file_type == "dmg_global") {
        k <- str_match(
            string = file,
            pattern = ".*/k(\\d.+)/(\\d+)/global/(\\S{10})\\..*"
        )[1, 2][1]
        ani <- str_match(
            string = file,
            pattern = ".*/k(\\d.+)/(\\d+)/global/(\\S{10})\\..*"
        )[1, 3][1]
        label <- str_match(
            string = file,
            pattern = ".*/k(\\d.+)/(\\d+)/global/(\\S{10})\\..*"
        )[1, 4][1]
    } else if (file_type == "dmg_lca") {
        k <- str_match(
            string = file,
            pattern = ".*/k(\\d.+)/(\\d+)/lca/(\\S{10})\\..*"
        )[1, 2][1]
        ani <- str_match(
            string = file,
            pattern = ".*/k(\\d.+)/(\\d+)/lca/(\\S{10})\\..*"
        )[1, 3][1]
        label <- str_match(
            string = file,
            pattern = ".*/k(\\d.+)/(\\d+)/lca/(\\S{10})\\..*"
        )[1, 4][1]
    } else if (file_type == "dmg_local") {
        k <- str_match(
            string = file,
            pattern = ".*/k(\\d.+)/(\\d+)/local/(\\S{10})\\..*"
        )[1, 2][1]
        ani <- str_match(
            string = file,
            pattern = ".*/k(\\d.+)/(\\d+)/local/(\\S{10})\\..*"
        )[1, 3][1]
        label <- str_match(
            string = file,
            pattern = ".*/k(\\d.+)/(\\d+)/local/(\\S{10})\\..*"
        )[1, 4][1]
    }

    data <- fread(file, nThread = threads)


    if (file_type == "rc") {
        n_reads <- data |> nrow()
        n_reads_saturated <- data |>
            filter(count >= as.numeric(k)) |>
            nrow()
        data <- tibble(
            label = label,
            k = k,
            ani = ani,
            n_reads = n_reads,
            n_reads_saturated = n_reads_saturated,
            prop = n_reads_saturated / n_reads,
        )
    } else if (file_type == "sf") {
        data <- data |>
            mutate(
                label = label,
                k = k,
                ani = ani
            )
    } else if (file_type %in% c("dmg_lca", "dmg_global", "dmg_local")) {
        data <- data |>
            mutate(
                label = label,
                k = k,
                ani = ani
            )
    }

    return(data)
}

# data_prop <- pbmclapply(rhc_files, read_files, mc.cores = 20)
# data_prop <- data_prop |> bind_rows()
# write_tsv(data_prop, file = "./data/benchmarks/saturation.tsv.gz")
data_prop <- read_tsv("./data/benchmarks/saturation.tsv")
data_prop |>
    group_by(k, ani) |>
    count() |>
    View()

data_prop |>
    filter(k == 1000, ani == 92) |>
    mutate(prop = prop * 100) |>
    View()

data_prop |>
    mutate(
        k = paste0("k", k),
        ani = as.character(ani)
    ) |>
    group_by(k, ani) |>
    mutate(
        k = as.character(k),
        ani = as.character(ani)
    ) |>
    summarise(
        prop_mean = mean(prop),
        prop_sd = sd(prop),
        prop_median = median(prop),
        prop_mad = mad(prop),
    ) |>
    ungroup() |>
    mutate(k = fct_relevel(
        k,
        rev(c("k1000", "k750", "k500", "k250", "k100", "k50"))
    )) |>
    View()
ggplot(aes(x = ani, y = prop_median, color = k, group = k)) +
    geom_line() +
    geom_point(aes(x = ani, y = prop_median),
        shape = 21, color = "black", fill = "grey", size = 2
    ) +
    scale_y_continuous(labels = scales::percent) +
    xlab("Read ANI filter") +
    ylab("Proportion of reads saturated") +
    theme() +
    # guides(colour = guide_legend(nrow = 1), size = guide_legend(nrow = 1)) +
    theme_bw() +
    theme(
        text = element_text(size = 12),
        legend.position = "bottom",
        legend.title = element_blank(),
    ) +
    scale_color_manual(values = k_colors)

ggsave(
    filename = "./manuscript/figures/saturation.pdf",
    width = 12, height = 6
)

data_prop |>
    janitor::clean_names(case = "sentence") |>
    write_tsv("./manuscript/tables/saturation.tsv")

# data_prop |>
#     mutate(
#         k = as.character(k),
#         ani = as.character(ani)
#     ) |>
#     select(label, ani, n_reads) |>
#     group_by(label, ani) |>
#     summarise(
#         n_reads = round(mean(n_reads)),
#     ) |>
#     group_by(label) |>
#     mutate(rdiff = abs((n_reads - lag(n_reads))) / first(n_reads))
# group_by(label) |>
#     mutate(abs_diff = abs(rdiff - lag(rdiff))) |>
#     ungroup() |>
#     ggplot(aes(x = ani, y = abs_diff, group = label, color = label)) +
#     geom_line() +
#     geom_point(
#         shape = 21, color = "black", fill = "grey"
#     )



data_prop |>
    mutate(
        k = as.character(k),
        ani = as.character(ani)
    ) |>
    select(label, ani, n_reads) |>
    group_by(label, ani) |>
    summarise(
        n_reads = round(mean(n_reads)),
    ) |>
    ungroup() |>
    ggplot(aes(x = ani, y = n_reads, group = label, color = label)) +
    geom_line() +
    geom_point(aes(x = ani, y = n_reads),
        shape = 21, color = "black", fill = "grey"
    ) +
    xlab("Read ANI filter") +
    ylab("Number of reads") +
    theme_bw() +
    guides(size = guide_legend(nrow = 1))

map_dfr_progress <- function(.x, .f, ..., .id = NULL) {
    .f <- purrr::as_mapper(.f, ...)
    pb <- progress::progress_bar$new(
        total = length(.x),
        format = " [:bar] :current/:total (:percent) eta: :eta",
        force = TRUE
    )

    f <- function(...) {
        pb$tick()
        .f(...)
    }
    purrr::map_dfr(.x, f, ..., .id = .id)
}

stats <- map_dfr_progress(sf_files, read_files,
    file_type = "sf", threads = 4
)

# stats |> filter(label == "7fe94f44b2")

# write_tsv(stats, file = "./data/benchmarks/reference-filtering-k.tsv.gz")
# stats <- read_tsv("./data/benchmarks/reference-filtering-k.tsv.gz")

nrefs_euk <- tax_info |>
    filter(domain == "d__Eukaryota") |>
    nrow()
nrefs_bact <- tax_info |>
    filter(domain == "d__Bacteria") |>
    nrow()
nrefs_arc <- tax_info |>
    filter(domain == "d__Archaea") |>
    nrow()
nrefs_vir <- tax_info |>
    filter(domain == "d__Viruses") |>
    nrow()

nrefs <- tibble(
    domain = c("d__Eukaryota", "d__Bacteria", "d__Archaea", "d__Viruses"),
    nrefs = c(nrefs_euk, nrefs_bact, nrefs_arc, nrefs_vir)
)

stats |>
    filter(breadth >= 0.01) |>
    inner_join(tax_info) |>
    filter(domain != "d__Eukaryota") |>
    group_by(k, label, ani, domain) |>
    count() |>
    ungroup() |>
    inner_join(nrefs) |>
    mutate(prop = n / nrefs) |>
    group_by(k, ani, domain) |>
    summarise(
        nrefs_median = median(prop),
        nrefs_mad = mad(prop),
        nrefs_mean = mean(prop),
        nrefs_sd = sd(prop)
    ) |>
    ungroup() |>
    mutate(k = fct_relevel(
        k,
        rev(c("1000", "750", "500", "250", "100", "50"))
    )) |>
    ggplot(aes(x = ani, y = nrefs_mean, color = k, group = k)) +
    geom_line() +
    # geom_errorbar(aes(ymin = prop_mean - prop_sd,
    # ymax = prop_mean + prop_sd),
    # width = 0) +
    geom_point(aes(x = ani, y = nrefs_mean, size = nrefs_sd),
        shape = 21, color = "black", fill = "grey"
    ) +
    facet_wrap(~domain, scales = "free_y") +
    guides(colour = guide_legend(nrow = 1)) +
    xlab("Read ANI filter") +
    ylab("Number of references") +
    theme_bw() +
    theme(
        text = element_text(size = 16),
        legend.position = "bottom",
        legend.title = element_blank(),
    )


stats |>
    filter(breadth >= 0.01) |>
    inner_join(tax_info) |>
    filter(domain != "d__Eukaryota") |>
    filter(
        n_reads >= 100
    ) |>
    group_by(k, ani, label) |>
    count() |>
    ungroup() |>
    inner_join(nrefs) |>
    mutate(prop = n / nrefs) |>
    group_by(k, ani) |>
    summarise(
        prop_mean = mean(prop),
        prop_sd = sd(prop)
    ) |>
    ungroup() |>
    mutate(k = fct_relevel(
        k,
        rev(c("1000", "750", "500", "250", "100", "50"))
    )) |>
    ggplot(aes(x = ani, y = prop_mean, color = k, group = k)) +
    geom_line() +
    # geom_errorbar(aes(ymin = prop_mean - prop_sd,
    # ymax = prop_mean + prop_sd),
    # width = 0) +
    geom_point(aes(x = ani, y = prop_mean, size = prop_sd),
        shape = 21, color = "black", fill = "grey"
    ) +
    guides(colour = guide_legend(nrow = 1)) +
    xlab("Read ANI filter") +
    ylab("Number of reference") +
    theme_bw() +
    theme(legend.position = "bottom")

dmg_global <- map_dfr_progress(dmg_files_global, read_files,
    file_type = "dmg_global", threads = 4
)

dmg_global |>
    filter(
        MAP_significance > 2
    ) |>
    mutate(k = fct_relevel(
        k,
        rev(c("1000", "750", "500", "250", "100", "50"))
    )) |>
    ggplot(aes(x = ani, y = MAP_damage, color = k, group = k)) +
    geom_line() +
    geom_point(shape = 21, color = "black", fill = "grey", size = 1.5) +
    facet_wrap(~label, scales = "free", nrow = 2) +
    guides(colour = guide_legend(nrow = 1)) +
    xlab("Read ANI filter") +
    ylab("MAP damage") +
    theme_bw() +
    theme(legend.position = "bottom")


dmg_global |>
    filter(
        MAP_significance > 2
    ) |>
    group_by(k, ani) |>
    summarise(
        MAP_damage_mean = median(MAP_damage),
        MAP_damage_sd = mad(MAP_damage)
    ) |>
    ungroup() |>
    mutate(k = fct_relevel(
        k,
        rev(c("1000", "750", "500", "250", "100", "50"))
    )) |>
    ggplot(aes(x = ani, y = MAP_damage_mean, color = k, group = k)) +
    geom_line() +
    geom_point(aes(x = ani, y = MAP_damage_mean, size = MAP_damage_sd),
        shape = 21, color = "black", fill = "grey"
    ) +
    guides(colour = guide_legend(nrow = 1)) +
    xlab("Read ANI filter") +
    ylab("MAP damage") +
    theme_bw() +
    theme(legend.position = "bottom")


dmg_local <- map_dfr_progress(dmg_files_local, read_files,
    file_type = "dmg_local", threads = 4
)

dmg_local <- dmg_local |>
    as_tibble() |>
    rename(reference = tax_id) |>
    inner_join(tax_info) |>
    filter(domain != "d__Eukaryota")

# dmg_lca |>
#     filter(
#         MAP_significance > 2,
#         tax_rank == "species",
#         N_reads >= 100
#     ) |>
#     group_by(k, ani) |>
#     summarise(
#         MAP_damage_mean = mean(MAP_damage),
#         MAP_damage_sd = sd(MAP_damage)
#     ) |>
#     ungroup() |>
#     mutate(k = fct_relevel(
#         k,
#         rev(c("1000", "750", "500", "250", "100", "50"))
#     )) |>
#     ggplot(aes(x = ani, y = MAP_damage_mean, color = k, group = k)) +
#     geom_line() +
#     geom_point(aes(x = ani, y = MAP_damage_mean, size = MAP_damage_sd),
#         shape = 21, color = "black", fill = "grey"
#     ) +
#     guides(colour = guide_legend(nrow = 1)) +
#     xlab("Read ANI filter") +
#     ylab("MAP damage") +
#     theme_bw() +
#     theme(legend.position = "bottom")

dmg_local |>
    filter(
        MAP_significance > 2,
        MAP_damage > 0.1,
        N_reads >= 100
    ) |>
    group_by(k, ani, domain) |>
    summarise(
        MAP_damage_mean = mean(MAP_damage),
        MAP_damage_sd = sd(MAP_damage),
        MAP_damage_median = median(MAP_damage),
        MAP_damage_mad = mad(MAP_damage)
    ) |>
    ungroup() |>
    mutate(k = fct_relevel(
        k,
        rev(c("1000", "750", "500", "250", "100", "50"))
    )) |>
    ggplot(aes(x = ani, y = MAP_damage_median, color = k, group = k)) +
    geom_line() +
    geom_point(aes(x = ani, y = MAP_damage_median, size = MAP_damage_mad),
        shape = 21, color = "black", fill = "grey"
    ) +
    facet_wrap(~domain, scales = "free_y") +
    guides(colour = guide_legend(nrow = 1)) +
    xlab("Read ANI filter") +
    ylab("MAP damage") +
    theme_bw() +
    theme(
        text = element_text(size = 16),
        legend.position = "bottom",
        legend.title = element_blank(),
    )

dmg_local |>
    filter(
        MAP_significance > 2,
        MAP_damage > 0.1,
        N_reads >= 100
    ) |>
    group_by(k, label, ani, domain) |>
    count() |>
    ungroup() |>
    inner_join(nrefs) |>
    mutate(prop = n / nrefs) |>
    group_by(k, ani, domain) |>
    summarise(
        nrefs_median = median(prop),
        nrefs_mad = mad(prop),
        nrefs_mean = mean(prop),
        nrefs_sd = sd(prop)
    ) |>
    ungroup() |>
    mutate(k = fct_relevel(
        k,
        rev(c("1000", "750", "500", "250", "100", "50"))
    )) |>
    ggplot(aes(x = ani, y = nrefs_median, color = k, group = k)) +
    geom_line() +
    # geom_errorbar(aes(ymin = prop_mean - prop_sd,
    # ymax = prop_mean + prop_sd),
    # width = 0) +
    geom_point(aes(x = ani, y = nrefs_median, size = nrefs_mad),
        shape = 21, color = "black", fill = "grey"
    ) +
    facet_wrap(~domain, scales = "free_y") +
    guides(colour = guide_legend(nrow = 1)) +
    xlab("Read ANI filter") +
    ylab("Number of references") +
    theme_bw() +
    theme(
        text = element_text(size = 16),
        legend.position = "bottom",
        legend.title = element_blank(),
    )
