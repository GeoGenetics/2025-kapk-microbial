library(tidyverse)

get_sim_data <- function(X) {
    wf <- strsplit(X, "/")[[1]][5]
    k <- strsplit(X, "/")[[1]][6]
    id <- strsplit(X, "/")[[1]][7]
    read_tsv(X, show_col_types = FALSE) |>
        mutate(k = k, id = id, wf = wf)
}

get_dmg_sim_data <- function(X) {
    wf <- strsplit(X, "/")[[1]][5]
    k <- strsplit(X, "/")[[1]][6]
    id <- strsplit(X, "/")[[1]][7]
    read_csv(X, show_col_types = FALSE) |>
        select(label = sample, reference = tax_id, N_reads, damage, significance) |>
        mutate(k = k, id = id, wf = wf)
}

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
