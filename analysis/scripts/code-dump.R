# d1 <- briggs_1 |>
#     filter(prob_dmg>0.9)

# eps <- c(1, 2, 3, 4, 5)

# briggs_0 <- map_dfr(eps, function(ep) {
#     read_tsv(paste0("/maps/projects/fernandezguerra/scratch/kapk-briggs/c_000000000001-nonidentical.nc.",rp,".sorted.briggs.bam"), col_names = c("read_name", "prob_ancient", "prob_dmg")) |>
#         mutate(ep = ep)
# })



# Let's check eps 10 and different isrecal values
# data_b0_all <- briggs_0 |>
#     inner_join(read_ids)



# data_b1_all |>
# mutate(class= ifelse(percid == 1, "identical", "not-identical")) |>
# ggplot(aes(prob_ancient, fill = class)) +
# geom_density()



# data_b0 <- briggs_0 |>
#     filter(ep == 10) |>
#     inner_join(read_ids) |>
#     select(-ep)

# data_b0_40 <- briggs_0 |>
#     filter(ep == 40) |>
#     inner_join(read_ids) |>
#     select(-ep)

# data_b1_10 <- briggs_1_10 |>
#     inner_join(read_ids)

# data_b1_30 <- briggs_1_30 |>
#     inner_join(read_ids)

# p1 <- ggplot(data_b0, aes(x = prob_ancient, y = prob_dmg)) +
#     # geom_point() +
#     # geom_density_2d_filled(contour_var = "ndensity") +
#     # geom_density_2d(linewidth = 0.25, colour = "black") +
#     geom_hdr() +
#     # geom_point(shape = 21, color = "grey", alpha = 0.3, size = 0.5) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # scale_fill_viridis_c(direction = -1, option = "inferno") +
#     theme_bw() +
#     xlab("AN (-isrecal 0 -eps 0.1)") +
#     ylab("PD (-isrecal 0 -eps 0.1)") +
#     theme(
#         legend.position = "top",
#         legend.title = element_blank(),
#         text = element_text(size = 12),
#         strip.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.text.x = element_blank(),
#     )


# p2 <- ggplot(data_b1, aes(x = prob_ancient, y = prob_dmg)) +
#     # geom_point() +
#     # geom_density_2d_filled(contour_var = "ndensity") +
#     # geom_density_2d(linewidth = 0.25, colour = "black") +
#     geom_hdr() +
#     # geom_point(shape = 21, color = "grey", alpha = 0.3, size = 0.5) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # scale_fill_viridis_c(direction = -1, option = "inferno") +
#     theme_bw() +
#     xlab("AN (-isrecal 1 -eps 0.1)") +
#     ylab("PD (-isrecal 1 -eps 0.1)") +
#     theme(
#         legend.position = "top",
#         legend.title = element_blank(),
#         text = element_text(size = 12),
#         strip.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.text.x = element_blank(),
#     )



# p3 <- data_b0 |>
#     select(read_name, prob_ancient_0 = prob_ancient) |>
#     inner_join(data_b1 |> select(read_name, prob_ancient_1 = prob_ancient), by = "read_name") |>
#     ggplot(aes(x = prob_ancient_0, y = prob_ancient_1)) +
#     # geom_point() +
#     # geom_density_2d_filled(contour_var = "ndensity") +
#     # geom_density_2d(linewidth = 0.25, colour = "black") +
#     geom_hdr() +
#     # geom_point(shape = 21, color = "grey", alpha = 0.3, size = 0.5) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # scale_fill_viridis_c(direction = -1, option = "inferno") +
#     theme_bw() +
#     xlab("AN -isrecal 0") +
#     ylab("AN -isrecal 1") +
#     theme(
#         legend.position = "top",
#         legend.title = element_blank(),
#         text = element_text(size = 12),
#         strip.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.text.x = element_blank(),
#     )


# p4 <- data_b0 |>
#     select(read_name, prob_dmg_0 = prob_dmg) |>
#     inner_join(data_b1 |> select(read_name, prob_dmg_1 = prob_dmg), by = "read_name") |>
#     ggplot(aes(x = prob_dmg_0, y = prob_dmg_1)) +
#     # geom_point() +
#     # geom_density_2d_filled(contour_var = "ndensity") +
#     # geom_density_2d(linewidth = 0.25, colour = "black") +
#     geom_hdr() +
#     # geom_point(shape = 21, color = "grey", alpha = 0.3, size = 0.5) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # scale_fill_viridis_c(direction = -1, option = "inferno") +
#     theme_bw() +
#     xlab("PD -isrecal 0") +
#     ylab("PD -isrecal 1") +
#     theme(
#         legend.position = "top",
#         legend.title = element_blank(),
#         text = element_text(size = 12),
#         strip.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.text.x = element_blank(),
#     )




# p5 <- data_b0 |>
#     inner_join(read_ids) |>
#     ggplot(aes(x = percid, y = prob_ancient)) +
#     # geom_point() +
#     # geom_density_2d_filled(contour_var = "ndensity") +
#     # geom_density_2d(linewidth = 0.25, colour = "black") +
#     geom_hdr() +
#     # geom_vline(xintercept = 0.995, linetype = "dashed", color = "red") +
#     # geom_hline(yintercept = 0.975, linetype = "dashed", color = "red") +
#     # geom_point(shape = 21, color = "grey", alpha = 0.3, size = 0.5) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # scale_fill_viridis_c(direction = -1, option = "inferno") +
#     theme_bw() +
#     xlab("ANIr") +
#     ylab("AN (-isrecal 0; -eps 0.1)") +
#     theme(
#         legend.position = "top",
#         legend.title = element_blank(),
#         text = element_text(size = 12),
#         strip.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.text.x = element_blank(),
#     )


# p6 <- data_b1 |>
#     inner_join(read_ids) |>
#     ggplot(aes(x = percid, y = prob_ancient)) +
#     # geom_point() +
#     # geom_density_2d_filled(contour_var = "ndensity") +
#     # geom_density_2d(linewidth = 0.25, colour = "black") +
#     geom_hdr() +
#     # geom_vline(xintercept = 0.995, linetype = "dashed", color = "red") +
#     # geom_hline(yintercept = 0.975, linetype = "dashed", color = "red") +
#     # geom_point(shape = 21, color = "grey", alpha = 0.3, size = 0.5) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # scale_fill_viridis_c(direction = -1, option = "inferno") +
#     theme_bw() +
#     xlab("ANIr") +
#     ylab("AN (-isrecal 1; -eps 0.1)") +
#     theme(
#         legend.position = "top",
#         legend.title = element_blank(),
#         text = element_text(size = 12),
#         strip.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.text.x = element_blank(),
#     )


# p7 <- data_b0 |>
#     inner_join(read_ids) |>
#     ggplot(aes(x = percid, y = prob_ancient)) +
#     # geom_point() +
#     # geom_density_2d_filled(contour_var = "ndensity") +
#     # geom_density_2d(linewidth = 0.25, colour = "black") +
#     geom_hdr() +
#     # geom_vline(xintercept = 0.995, linetype = "dashed", color = "red") +
#     # geom_hline(yintercept = 0.975, linetype = "dashed", color = "red") +
#     # geom_point(shape = 21, color = "grey", alpha = 0.3, size = 0.5) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # scale_fill_viridis_c(direction = -1, option = "inferno") +
#     theme_bw() +
#     xlab("ANIr") +
#     ylab("AN (-isrecal 0; -eps 0.1)") +
#     theme(
#         legend.position = "top",
#         legend.title = element_blank(),
#         text = element_text(size = 12),
#         strip.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.text.x = element_blank(),
#     )
# p8 <- data_b1 |>
#     inner_join(read_ids) |>
#     ggplot(aes(x = percid, y = prob_ancient)) +
#     # geom_point() +
#     # geom_density_2d_filled(contour_var = "ndensity") +
#     # geom_density_2d(linewidth = 0.25, colour = "black") +
#     geom_hdr() +
#     # geom_vline(xintercept = 0.995, linetype = "dashed", color = "red") +
#     # geom_hline(yintercept = 0.975, linetype = "dashed", color = "red") +
#     # geom_point(shape = 21, color = "grey", alpha = 0.3, size = 0.5) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # scale_fill_viridis_c(direction = -1, option = "inferno") +
#     theme_bw() +
#     xlab("ANIr") +
#     ylab("AN (-isrecal 1; -eps 0.1)") +
#     theme(
#         legend.position = "top",
#         legend.title = element_blank(),
#         text = element_text(size = 12),
#         strip.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.text.x = element_blank(),
#     )



# p9 <- data_b0_all |>
#     mutate(ep = paste0(ep, "%")) |>
#     ggplot(aes(x = percid, y = prob_ancient)) +
#     # geom_point() +
#     # geom_density_2d_filled(contour_var = "ndensity") +
#     # geom_density_2d(linewidth = 0.25, colour = "black") +
#     geom_hdr() +
#     # geom_vline(xintercept = 0.995, linetype = "dashed", color = "red") +
#     # geom_hline(yintercept = 0.975, linetype = "dashed", color = "red") +
#     # geom_point(shape = 21, color = "grey", alpha = 0.3, size = 0.5) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # scale_fill_viridis_c(direction = -1, option = "inferno") +
#     theme_bw() +
#     xlab("ANIr") +
#     ylab("AN (-isrecal 0)") +
#     facet_grid(~ep) +
#     theme(
#         legend.position = "top",
#         legend.title = element_blank(),
#         text = element_text(size = 12),
#         strip.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#     )


# p9a <- data_b1_all |>
#     mutate(ep = paste0(ep, "%")) |>
#     ggplot(aes(x = percid, y = prob_ancient)) +
#     # geom_point() +
#     # geom_density_2d_filled(contour_var = "ndensity") +
#     # geom_density_2d(linewidth = 0.25, colour = "black") +
#     geom_hdr() +
#     # geom_vline(xintercept = 0.995, linetype = "dashed", color = "red") +
#     # geom_hline(yintercept = 0.975, linetype = "dashed", color = "red") +
#     # geom_point(shape = 21, color = "grey", alpha = 0.3, size = 0.5) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # scale_fill_viridis_c(direction = -1, option = "inferno") +
#     theme_bw() +
#     xlab("ANIr") +
#     ylab("AN (-isrecal 1)") +
#     facet_grid(~ep) +
#     theme(
#         legend.position = "top",
#         legend.title = element_blank(),
#         text = element_text(size = 12),
#         strip.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#     )

# p9a





# p9
# ggsave("manuscript/figures/b-fig4.pdf", width = 12, height = 3, dpi = 300)


# ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "top")
# ggsave("manuscript/figures/b-fig1.pdf", width = 8, height = 4, dpi = 300)

# ggpubr::ggarrange(p3, p4, ncol = 2, nrow = 1, common.legend = TRUE, legend = "top")
# ggsave("manuscript/figures/b-fig2.pdf", width = 8, height = 4, dpi = 300)

# ggpubr::ggarrange(p5, p6, ncol = 2, nrow = 1, common.legend = TRUE, legend = "top")
# ggsave("manuscript/figures/b-fig3.pdf", width = 8, height = 4, dpi = 300)

# ggpubr::ggarrange(p7, p8, ncol = 2, nrow = 1, common.legend = TRUE, legend = "top")
# ggsave("manuscript/figures/b-fig3.pdf", width = 8, height = 4, dpi = 300)



# data_b0_all |>
#     filter(ep == 30) |>
#     inner_join(read_ids) |>
#     ggplot(aes(x = percid, y = prob_ancient)) +
#     # geom_point() +
#     # geom_density_2d_filled(contour_var = "ndensity") +
#     # geom_density_2d(linewidth = 0.25, colour = "black") +
#     geom_hdr() +
#     geom_vline(xintercept = id_threshold, linetype = "dashed", color = "red") +
#     geom_hline(yintercept = an_threshold, linetype = "dashed", color = "red") +
#     # geom_point(shape = 21, color = "grey", alpha = 0.3, size = 0.5) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # scale_fill_viridis_c(direction = -1, option = "inferno") +
#     theme_bw() +
#     xlab("ANIr") +
#     ylab("AN (-isrecal 1; -eps 0.1)") +
#     theme(
#         legend.position = "top",
#         legend.title = element_blank(),
#         text = element_text(size = 12),
#         strip.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.text.x = element_blank(),
#     )

# p10 <- data_b1_10 |>
#     inner_join(read_ids) |>
#     ggplot(aes(x = percid, y = prob_ancient)) +
#     # geom_point() +
#     # geom_density_2d_filled(contour_var = "ndensity") +
#     # geom_density_2d(linewidth = 0.25, colour = "black") +
#     geom_hdr() +
#     # geom_vline(xintercept = id_threshold, linetype = "dashed", color = "red") +
#     # geom_hline(yintercept = an_threshold, linetype = "dashed", color = "red") +
#     # geom_point(shape = 21, color = "grey", alpha = 0.3, size = 0.5) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # scale_fill_viridis_c(direction = -1, option = "inferno") +
#     theme_bw() +
#     xlab("ANIr") +
#     ylab("AN (-isrecal 1; -eps 0.1)") +
#     theme(
#         legend.position = "top",
#         legend.title = element_blank(),
#         text = element_text(size = 12),
#         strip.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.text.x = element_blank(),
#     )

# p11 <- data_b1_30 |>
#     inner_join(read_ids) |>
#     ggplot(aes(x = percid, y = prob_ancient)) +
#     # geom_point() +
#     # geom_density_2d_filled(contour_var = "ndensity") +
#     # geom_density_2d(linewidth = 0.25, colour = "black") +
#     geom_hdr() +
#     geom_vline(xintercept = id_threshold, linetype = "dashed", color = "red") +
#     geom_hline(yintercept = an_threshold, linetype = "dashed", color = "red") +
#     # geom_point(shape = 21, color = "grey", alpha = 0.3, size = 0.5) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#     # scale_fill_viridis_c(direction = -1, option = "inferno") +
#     theme_bw() +
#     xlab("ANIr") +
#     ylab("AN (-isrecal 1; -eps 0.3)") +
#     theme(
#         legend.position = "top",
#         legend.title = element_blank(),
#         text = element_text(size = 12),
#         strip.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.text.x = element_blank(),
#     )

# p11

# ggpubr::ggarrange(p10, p11, ncol = 2, nrow = 1, common.legend = TRUE, legend = "top")



# data_avg |>
#     inner_join(read_ids) |>
#     mutate(
#         class = case_when(
#             prob_ancient >= an_threshold ~ "group0",
#             prob_ancient < an_threshold & percid >= id_threshold ~ "group2",
#             prob_ancient < an_threshold & percid < id_threshold ~ "group1"
#         )
#     ) |>
#     ggplot(aes(x = read_length)) +
#     geom_density() +
#     facet_wrap(~class, scales = "free") +
#     theme_bw()




# For the paper
data <- data_b1_all
fname <- "isrecal1"


# estimate probability threshold

probs <- seq(0, 0.98, 0.01)


find_prob <- map_dfr(probs, function(X) {
    n_reads <- data |>
        filter(prob_ancient >= X) |>
        nrow()
    tibble(prob = X, n_reads = n_reads)
})


fr <- phages |>
    as_tibble() |>
    unite(col = phage, genome, phage) |>
    select(phage, cl_name) |>
    distinct() |>
    mutate(n_phages = length(unique(phage))) |>
    group_by(cl_name) |>
    summarise(n = n(), perc = n / n_phages) |>
    ungroup() |>
    distinct() |>
    arrange(desc(perc)) |>
    mutate(pos = row_number())


df <- df |> mutate(second_derivative = lead(n_reads, 2) - 2 * lead(n_reads) + n_reads)



d2 <- with(find_prob, diff(n_reads) / diff(prob))

d1 <- diff(find_prob$n_reads)
k <- which.max(abs(diff(d1) / diff(find_prob$n_reads[-1])))
p <- find_prob |> filter(rank == k)

p <- tibble(prob = probs, n_reads = c(0, d1)) |>
    ggplot(aes(probs, n_reads)) +
    geom_line()




drv <- function(x, y) c(NA, (diff(y) / diff(x)))
middle_pts <- function(x) c(NA, (x[-1] - diff(x) / 2))

find_prob |>
    mutate(second_d = drv(middle_pts(prob), drv(n_reads, prob)))

find_prob <- find_prob |>
    mutate(rank = row_number())

get_inflection <- function(df) {
    log_total <- log_rank <- total <- NULL
    df_fit <- df |>
        transmute(
            log_total = log10(n_reads),
            log_rank = log10(rank)
        )
    d1n <- diff(df_fit$log_total) / diff(df_fit$log_rank)
    right.edge <- which.min(d1n)
    10^(df_fit$log_total[right.edge])
}
inflection <- get_inflection(find_prob)

df <- find_prob |> filter(rank < 101)
annot <- tibble(
    inflection = inflection,
    rank_cutoff = max(df$rank[df$n_reads > inflection])
)
ggplot(df, aes(n_reads, rank)) +
    geom_path() +
    geom_vline(aes(xintercept = inflection),
        data = annot, linetype = 2,
        color = "gray40"
    ) +
    geom_hline(aes(yintercept = rank_cutoff),
        data = annot, linetype = 2,
        color = "gray40"
    ) +
    geom_text(
        aes(inflection, rank_cutoff,
            label = paste(rank_cutoff, "'cells'")
        ),
        data = annot, vjust = 1
    ) +
    scale_x_log10() +
    scale_y_log10() +
    labs(y = "Rank", x = "Total UMIs") +
    annotation_logticks()

library(PCAtools)
# identify the elbow point
elbow_point <- find_prob |> filter(n_reads == find_prob$n_reads[findElbowPoint(find_prob$n_reads)])

ggthemr::ggthemr(layout = "scientific", palette = "fresh")

find_prob |>
    ggplot(aes(prob, n_reads), group = 1) +
    geom_line() +
    geom_point(size = 3, shape = 21, color = "black", fill = "red") +
    gghighlight::gghighlight(prob == 0.95, unhighlighted_params = list(fill = "grey")) +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    scale_y_continuous(labels = scales::comma, trans = "log10") +
    xlab("") +
    ylab("Number of reads")
p
plotly::ggplotly(p)


# Thresholds
# isrecal 0:
id_threshold <- 1
an_threshold <- 0.95

# # isrecal 1:
# id_threshold <- 1
# an_threshold <- 0.95

# Get average
data_avg <- data |>
    select(-ep, -percid) |>
    group_by(read_name) |>
    summarize(
        prob_ancient = mean(prob_ancient),
        prob_dmg = mean(prob_dmg)
    ) |>
    ungroup() |>
    inner_join(read_ids)



# g1_min_x <- df |> filter(subgroup == 1) |> pull(x) |> min()
# g1_max_x <- df |> filter(subgroup == 1) |> pull(x) |> max()
# g2_min_x <- df |> filter(subgroup == 2) |> pull(x) |> min()
# g2_max_x <- df |> filter(subgroup == 2) |> pull(x) |> max()
# g3_min_x <- df |> filter(subgroup == 3) |> pull(x) |> min()
# g3_max_x <- df |> filter(subgroup == 3) |> pull(x) |> max()

# g1_min_y <- df |> filter(subgroup == 1) |> pull(y) |> min()
# g1_max_y <- df |> filter(subgroup == 1) |> pull(y) |> max()
# g2_min_y <- df |> filter(subgroup == 2) |> pull(y) |> min()
# g2_max_y <- df |> filter(subgroup == 2) |> pull(y) |> max()
# g3_min_y <- df |> filter(subgroup == 3) |> pull(y) |> min()
# g3_max_y <- df |> filter(subgroup == 3) |> pull(y) |> max()

# maxmin_df <- tibble(group = c("g1", "g2", "g3"),
#                     prob_dmg_min = c(g1_min_x, g2_min_x, g3_min_x),
#                     prob_dmg_max = c(g1_max_x, g2_max_x, g3_max_x),
#                     percid_min = c(g1_min_y, g2_min_y, g3_min_y),
#                     percid_max = c(g1_max_y, g2_max_y, g3_max_y))

# df$colour[df$x>0.75] <- "red"
# df$subgroup[df$x>0.75] <- 2



# g1 <- data |>
#     filter(prob_dmg >= (maxmin_df |> filter(group == "g1") |> pull(prob_dmg_min)) &
#             prob_dmg <= (maxmin_df |> filter(group == "g1") |> pull(prob_dmg_max)) &
#             percid >= (maxmin_df |> filter(group == "g1") |> pull(percid_min)) &
#             percid <= (maxmin_df |> filter(group == "g1") |> pull(percid_max))) |>
#             mutate(group = "g1")

# g2 <- data |>
#     filter(prob_dmg >= (maxmin_df |> filter(group == "g2") |> pull(prob_dmg_min)) &
#             prob_dmg <= (maxmin_df |> filter(group == "g2") |> pull(prob_dmg_max)) &
#             percid >= (maxmin_df |> filter(group == "g2") |> pull(percid_min)) &
#             percid <= (maxmin_df |> filter(group == "g2") |> pull(percid_max))) |>
#             mutate(group = "g2")

# g3 <- data |>
#     filter(prob_dmg >= (maxmin_df |> filter(group == "g3") |> pull(prob_dmg_min)) &
#             prob_dmg <= (maxmin_df |> filter(group == "g3") |> pull(prob_dmg_max)) &
#             percid >= (maxmin_df |> filter(group == "g3") |> pull(percid_min)) &
#             percid <= (maxmin_df |> filter(group == "g3") |> pull(percid_max))) |>
#             mutate(group = "g3")


# bind_rows(g1,g2,g3) |>
# right_join(data) |>
# filter(is.na(group)) |>
# ggplot(aes(prob_dmg, prob_ancient)) +
# geom_point()

# g2 |>
# ggplot(aes(prob_dmg, prob_ancient)) +
# geom_hdr()

data |>
    inner_join(read_ids) |>
    # mutate(
    #     class = case_when(
    #         prob_ancient >= an_threshold & percid < id_threshold ~ "group0",
    #         prob_ancient < an_threshold & percid < id_threshold ~ "group1",
    #         percid >= id_threshold ~ "group2",
    #     )
    # ) |>
    mutate(class = case_when(
        prob_ancient >= an_threshold ~ "group0",
        prob_ancient < an_threshold ~ "group1",
    )) |>
    group_by(class) |>
    count()

data |>
    inner_join(read_ids) |>
    # mutate(
    #     class = case_when(
    #         prob_ancient >= an_threshold & percid < id_threshold ~ "group0",
    #         prob_ancient < an_threshold & percid < id_threshold ~ "group1",
    #         percid >= id_threshold ~ "group2",
    #     )
    # ) |>
    mutate(class = case_when(
        prob_ancient >= an_threshold ~ "group0",
        prob_ancient < an_threshold ~ "group1",
    )) |>
    filter(percid == 1, class == "group1") |>
    ggplot(aes(read_length, fill = class)) +

    # geom_hdr()
    geom_density(alpha = 0.5)


data |>
    inner_join(read_ids) |>
    # mutate(
    #     class = case_when(
    #         prob_ancient >= an_threshold & percid < id_threshold ~ "group0",
    #         prob_ancient < an_threshold & percid < id_threshold ~ "group1",
    #         percid >= id_threshold ~ "group2",
    #     )
    # ) |>
    mutate(class = case_when(
        prob_ancient >= an_threshold ~ "group0",
        prob_ancient < an_threshold ~ "group1",
    )) |>
    ggplot(aes(read_length, fill = class)) +
    geom_density(alpha = 0.5)

read_ids |>
    ggplot(aes(percid, read_length)) +
    geom_hdr()



# Create groups and export them
data_avg |>
    inner_join(read_ids) |>
    # mutate(
    #     class = case_when(
    #         prob_ancient >= an_threshold & percid < id_threshold ~ "group0",
    #         prob_ancient < an_threshold & percid < id_threshold ~ "group1",
    #         percid >= id_threshold ~ "group2",
    #     )
    # ) |>
    mutate(class = case_when(
        prob_ancient >= an_threshold ~ "group0",
        prob_ancient < an_threshold ~ "group1",
    )) |>
    select(read_name, class) |>
    write_tsv(paste0("./results/taxonomy/best-match-profiling/ancient-modern-briggs-classification-", fname, ".tsv.gz"), col_names = FALSE)

data |>
    # mutate(
    #     class = case_when(
    #         prob_ancient >= an_threshold & percid < id_threshold ~ "group0",
    #         prob_ancient < an_threshold & percid < id_threshold ~ "group1",
    #         percid >= id_threshold ~ "group2",
    #     )
    # ) |>
    mutate(class = case_when(
        prob_ancient >= an_threshold ~ "group0",
        prob_ancient < an_threshold ~ "group1",
    )) |>
    group_by(class) |>
    count()
