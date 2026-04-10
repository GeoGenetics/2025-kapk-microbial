library(tidyverse)
library(showtext)
showtext_auto()
# Create a dataframe with the data


data <- tibble(
    sample = rep("abfe73874c", 5),
    params = c("sens", "sens-N1", "sens-N1-L20", "vsens", "vsens-N1"),
    unique_aligned = c(1652691, 2037105, 2079449, 1856161, 2119162),
    multiple_aligned = c(2278727, 3390574, 3496992, 2792593, 3550463),
    total_reads = rep(36744399, 5),
    time = c(7, 72, 228, 27, 428),
) |>
    mutate(
        total_aligned = unique_aligned + multiple_aligned
    )

data1 <- tibble(
    sample = rep("c74844188e", 5),
    params = c("sens", "sens-N1", "sens-N1-L20", "vsens", "vsens-N1"),
    unique_aligned = c(11156018, 11972421, 12138116, 11590992, 12310476),
    multiple_aligned = c(14496115, 20012420, 20674677, 17010564, 21011143),
    total_reads = rep(132383467, 5),
    time = c(30, 225, 712, 93, 1098),
) |>
    mutate(
        total_aligned = unique_aligned + multiple_aligned
    )

data2 <- tibble(
    sample = rep("d0bc152641", 5),
    params = c("sens", "sens-N1", "sens-N1-L20", "vsens", "vsens-N1"),
    unique_aligned = c(6745873, 8206319, 8424734, 7495044, 8637158),
    multiple_aligned = c(16161000, 20579716, 20988104, 18260964, 21177065),
    total_reads = rep(204107985, 5),
    time = c(30, 350, 1202, 127, 1791),
) |>
    mutate(
        total_aligned = unique_aligned + multiple_aligned
    )

plot_data <- data |>
    mutate(
        base_aligned_reads = data[data["params"] == "sens", ]$total_aligned,
        base_time = data[data["params"] == "sens", ]$time,
    ) |>
    mutate(
        total_aligned_perc = total_aligned / total_reads,
        total_aligned_x = total_aligned / base_aligned_reads,
        time_x = time / base_time
    ) |>
    filter(params != "sens") |>
    mutate(
        params = fct_relevel(params, c("vsens", "sens-N1", "sens-N1-L20", "vsens-N1"))
    )

plot_data_base <- data |>
    mutate(
        base_aligned_reads = data[data["params"] == "sens", ]$total_aligned,
        base_time = data[data["params"] == "sens", ]$time,
    ) |>
    mutate(
        total_aligned_perc = total_aligned / total_reads,
        total_aligned_x = total_aligned / base_aligned_reads,
        time_x = time / base_time
    ) |>
    filter(params == "sens") |>
    mutate(params = "baseline")

plot_data1 <- data1 |>
    mutate(
        base_aligned_reads = data1[data1["params"] == "sens", ]$total_aligned,
        base_time = data1[data1["params"] == "sens", ]$time,
    ) |>
    mutate(
        total_aligned_perc = total_aligned / total_reads,
        total_aligned_x = total_aligned / base_aligned_reads,
        time_x = time / base_time
    ) |>
    filter(params != "sens") |>
    mutate(
        params = fct_relevel(params, c("vsens", "sens-N1", "sens-N1-L20", "vsens-N1"))
    )

plot_data1_base <- data1 |>
    mutate(
        base_aligned_reads = data1[data1["params"] == "sens", ]$total_aligned,
        base_time = data1[data1["params"] == "sens", ]$time,
    ) |>
    mutate(
        total_aligned_perc = total_aligned / total_reads,
        total_aligned_x = total_aligned / base_aligned_reads,
        time_x = time / base_time
    ) |>
    filter(params == "sens") |>
    mutate(params = "baseline")

plot_data2 <- data2 |>
    mutate(
        base_aligned_reads = data2[data2["params"] == "sens", ]$total_aligned,
        base_time = data2[data2["params"] == "sens", ]$time,
    ) |>
    mutate(
        total_aligned_perc = total_aligned / total_reads,
        total_aligned_x = total_aligned / base_aligned_reads,
        time_x = time / base_time
    ) |>
    filter(params != "sens") |>
    mutate(
        params = fct_relevel(params, c("vsens", "sens-N1", "sens-N1-L20", "vsens-N1"))
    )

plot_data2_base <- data2 |>
    mutate(
        base_aligned_reads = data2[data2["params"] == "sens", ]$total_aligned,
        base_time = data2[data2["params"] == "sens", ]$time,
    ) |>
    mutate(
        total_aligned_perc = total_aligned / total_reads,
        total_aligned_x = total_aligned / base_aligned_reads,
        time_x = time / base_time
    ) |>
    filter(params == "sens") |>
    mutate(params = "baseline")



p1 <- ggplot(plot_data, aes(x = params, y = total_aligned_x, group = 1)) +
    geom_line() +
    geom_point(size = 3, shape = 21, fill = "grey", color = "black") +
    labs(
        x = "Bpwtie2 parameters",
        y = "Reads fold increase from default parameters"
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

p2 <- ggplot(plot_data1, aes(x = params, y = time_x, group = 1)) +
    geom_line() +
    geom_point(size = 3, shape = 21, fill = "grey", color = "black") +
    labs(
        x = "Bpwtie2 parameters",
        y = "Time (min) fold increase from default parameters"
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

p3 <- ggplot(plot_data2, aes(x = params, y = time_x, group = 1)) +
    geom_line() +
    geom_point(size = 3, shape = 21, fill = "grey", color = "black") +
    labs(
        x = "Bpwtie2 parameters",
        y = "Time (min) fold increase from default parameters"
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )



plot_data |>
    bind_rows(plot_data_base) |>
    bind_rows(plot_data1) |>
    bind_rows(plot_data1_base) |>
    bind_rows(plot_data2) |>
    bind_rows(plot_data2_base) |>
    janitor::clean_names(case = "sentence") |>
    write_tsv("./manuscript/tables/bowtie2-benchmarks.tsv")


coeff <- 55

ggplot(plot_data, aes(x = params, group = 1)) +
    geom_col(aes(y = total_aligned_x)) +
    geom_line(aes(y = time_x / coeff), size = 2, color = "red") + # Divide by 10 to get the same range than the temperature
    geom_point(aes(y = time_x / coeff), size = 4, shape = 21, fill = "grey", color = "black") + # Divide by 10 to get the same range than the temperature
    scale_y_continuous(

        # Features of the first axis
        name = "Reads fold increase from default parameters",

        # Add a second axis and specify its features
        sec.axis = sec_axis(~ . * coeff, name = "Time (min) fold increase from default parameters")
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 24),
    )


coeff <- 31
ggplot(plot_data1, aes(x = params, group = 1)) +
    geom_col(aes(y = total_aligned_x)) +
    geom_line(aes(y = time_x / coeff), size = 2, color = "red") + # Divide by 10 to get the same range than the temperature
    geom_point(aes(y = time_x / coeff), size = 4, shape = 21, fill = "grey", color = "black") + # Divide by 10 to get the same range than the temperature

    scale_y_continuous(

        # Features of the first axis
        name = "Reads fold increase from default parameters",

        # Add a second axis and specify its features
        sec.axis = sec_axis(~ . * coeff, name = "Time (min) fold increase from default parameters")
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 24),
    )


coeff <- 50
ggplot(plot_data2, aes(x = params, group = 1)) +
    geom_col(aes(y = total_aligned_x)) +
    geom_line(aes(y = time_x / coeff), size = 2, color = "red") + # Divide by 10 to get the same range than the temperature
    geom_point(aes(y = time_x / coeff), size = 4, shape = 21, fill = "grey", color = "black") + # Divide by 10 to get the same range than the temperature

    scale_y_continuous(

        # Features of the first axis
        name = "Reads fold increase from default parameters",

        # Add a second axis and specify its features
        sec.axis = sec_axis(~ . * coeff, name = "Time (min) fold increase from default parameters")
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 12),
    )

ggsave("./manuscript/figures/bowtie2-benchmarks.pdf", width = 12, height = 8)

# find read overlapping between the different parameters ----------------------

# read files for each set of parameters
files <- list.files(path = "data/benchmarks", pattern = "ids", full.names = TRUE)

# read files with read_tsv and parse parameters from the file name

# create a list of dataframes
dfs <- map_dfr(files, function(x) {
    read_tsv(x, col_names = "read_id") |>
        mutate(
            params = gsub("abfe73874c-k50-90-|.ids", "", basename(x))
        )
})


dfs |>
    filter(params == "sens-N1") |>
    inner_join(dfs |>
        filter(params == "sens-N1-L20") |>
        select(-params)) |>
    nrow()
# 5354814

dfs |>
    filter(params == "sens-N1") |>
    anti_join(dfs |>
        filter(params == "sens-N1-L20") |>
        select(-params)) |>
    nrow()
# 72,865


dfs |>
    filter(params == "sens-N1-L20") |>
    anti_join(dfs |>
        filter(params == "sens-N1") |>
        select(-params)) |>
    nrow()
# 221627

dfs_long <- dfs |>
    mutate(count = 1) |>
    sample_n(10000000) |>
    pivot_wider(names_from = params, values_from = count, values_fill = 0) |>
    column_to_rownames("read_id")

library(ComplexUpset)
params <- c("sens", "vsens", "sens-N1", "sens-N1-L20", "vsens-N1")
names(params) <- params
upset(
    data = dfs_long, intersect = params,
    name = "Number of reads shared between different bowtie2 parameters",
    min_size = 0,
    width_ratio = 0.125
)
