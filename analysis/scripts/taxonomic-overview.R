library(tidyverse)
library(janitor)
library(ggthemr)
library(ggpubr)
library(ggsankey)
source("libs/lib.R")
# Get sample data
kapk_cdata <- readxl::read_xlsx("data/cdata/KapK_samples-20210702.xlsx") |>
  clean_names() |> View()
  select(member_unit, collapsed_read_files_path_29bp_and_q_29_only_collapsed_reads, collapsed_read_files_path_29bp_and_q_29_only_collapsed_reads_md5sum, figure_names, site)

kapk_cdata <- kapk_cdata |>
  mutate(label = basename(collapsed_read_files_path_29bp_and_q_29_only_collapsed_reads),
         label = gsub("\\..*","",label)) |>
  distinct() |>
  select(label, collapsed_read_files_path_29bp_and_q_29_only_collapsed_reads_md5sum, member_unit, figure_names, site) |>
  rename(label_orig = label, label = collapsed_read_files_paths_29bp_and_q_29_only_collapsed_reads_md5sum) |>
  mutate(site = as.character(site),
         site = ifelse(site == "NA", "nosite", site))

site_order <- tibble(site = c("119", "50", "69", "75", "74", "nosite"),
                     site_rnk = 1:6)
# Initial stats -----------------------------------------------------------

# Read different steps stats
# Initial number of reads
initial_stats <- read_tsv("data/stats/all.stats-initial-summary.tsv.gz") |>
  clean_names()

# After extension and dereplication (used for profiling)
extended_stats <- read_tsv("data/stats/all.stats-extension-summary.tsv.gz") |>
  clean_names()

# After extension and dereplication (used for profiling)
derep_stats <- read_tsv("data/stats/all.stats-derep-summary.tsv.gz") |>
  clean_names()

# The results after profiling with the hires BD, only bac and arc
hires2motus_stats <- read_tsv("data/stats/all.stats-hires2motus-summary.tsv.gz") |>
  clean_names()

initial_rl_mode <- estimate_mode(initial_stats$avg_len)
initial_rl_avg <- mean(initial_stats$avg_len)
derep_rl_mode <- estimate_mode(derep_stats$avg_len)
derep_rl_avg <- mean(derep_stats$avg_len)

ggthemr(palette = "dust")
initial_stats |>
  mutate(class = "initial") |>
  bind_rows(derep_stats |>
              mutate(class = "derep")) |>
  bind_rows(extended_stats |>
              mutate(class = "extended")) |>
  mutate(class = fct_rev(class)) |>
  ggplot(aes(avg_len, fill = class)) +
  geom_density(alpha = 0.5, color = "black") +
  theme_light() +
  theme(legend.position = "top") +
  xlab("Average read length (nt)") +
  ylab("Density") +
  scale_color_discrete(name="")



# Number of reads ---------------------------------------------------------
# |class    |    num_seqs|
# |:--------|-----------:|
# |initial  | 10631904579|
# |extended | 10631904579|
# |derep    |  4867068849|
# 45% of reads used

initial_stats |>
  mutate(class = "initial") |>
  bind_rows(derep_stats |>
              mutate(class = "derep")) |>
  bind_rows(extended_stats |>
              mutate(class = "extended")) |>
  mutate(class = fct_rev(class)) |>
  group_by(class) |>
  summarise(num_seqs = sum(num_seqs)) |>
  ggplot(aes(class, num_seqs, fill = class)) +
  geom_col(alpha = 0.5, color = "black") +
  scale_y_continuous(labels = scales::comma) +
  xlab("") +
  ylab("Number of reads") +
  scale_color_discrete(name="") +
  theme(legend.position = "top")


# Coarse taxonomic classification -----------------------------------------
# Here we label reads at the domain level, using a coarse DB

k2_coarse <- read_tsv(file = "data/kraken2/kraken2-coarse-summary.tsv.gz")

k2_coarse_class <- k2_coarse |>
  group_by(label) |>
  dplyr::summarise(rclass = sum(counts)) |>
  ungroup() |>
  inner_join(derep_stats |> select(label, num_seqs)) |>
  mutate(unclassified = num_seqs - rclass)


k2_coarse_class_sum <- k2_coarse |>
  mutate(class = ifelse(domain == "Eukaryota", "Eukaryota", "hires")) |>
  group_by(label, class) |>
  dplyr::summarise(counts = sum(counts)) |>
  ungroup() |>
  pivot_wider(names_from = "class", values_from = counts, values_fill = 0) |>
  inner_join(derep_stats |> select(label, num_seqs)) |>
  mutate(hires = num_seqs - Eukaryota)

dom_h <- c("unclassified", "root", "Eukaryota", "Bacteria", "Archaea", "Viruses")
ggthemr::ggthemr(palette = "fresh")
k2_coarse |>
  bind_rows(k2_coarse_class |>
              select(label, unclassified) |>
              dplyr::rename(counts = unclassified) |>
              mutate(domain = "unclassified")) |>
  group_by(label) |>
  mutate(prop = counts/sum(counts)) |>
  ungroup() |>
  inner_join(kapk_cdata) |>
  inner_join(site_order) |>
  mutate(domain = fct_relevel(domain, rev(domains)),
         label = fct_reorder(label, site_rnk)) |>
  #filter(sample_type == "Lake sediment") |>
  ggplot(aes(label, prop, fill = domain)) +
  geom_col(width = 1, size = 0) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "top") +
  xlab("Samples") +
  ylab("Proportion") +
  facet_wrap(~member_unit, scales = "free_x") +
  scale_y_continuous(labels = scales::percent) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

dom_h <- c("unclassified", "root", "Eukaryota", "Bacteria", "Archaea", "Viruses")
k2_coarse |>
  bind_rows(k2_coarse_class |>
              select(label, unclassified) |>
              dplyr::rename(counts = unclassified) |>
              mutate(domain = "unclassified")) |>
  inner_join(kapk_cdata) |>
  group_by(domain, member_unit) |>
  summarise(counts = sum(counts)) |>
  ungroup() |>
  group_by(member_unit) |>
  mutate(prop = counts/sum(counts),
         domain = fct_relevel(domain, rev(dom_h))) |>
  #filter(sample_type == "Lake sediment") |>
  ggplot(aes(domain, prop, fill = domain)) +
  geom_col(width = 1, size = 0) +
  theme_bw() +
  theme(legend.position = "top") +
  xlab("Domain") +
  ylab("Proportion") +
  facet_wrap(~member_unit) +
  scale_y_continuous(labels = scales::percent) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  ggpubr::rotate()


# HIGH-RES TAXONOMY -------------------------------------------------------
k2_hires <- read_tsv("data/kraken2/kraken2-hires-summary.tsv.gz")

k2_hires_class <- k2_hires |>
  group_by(label) |>
  dplyr::summarise(rclass = sum(counts)) |>
  ungroup() |>
  inner_join(k2_coarse_class_sum |> select(label, hires)) |>
  mutate(unclassified = hires - rclass)

ggthemr::ggthemr(palette = "fresh")
k2_hires |>
  bind_rows(k2_hires_class |>
              select(label, unclassified) |>
              dplyr::rename(counts = unclassified) |>
              mutate(domain = "unclassified")) |>
  group_by(label) |>
  mutate(prop = counts/sum(counts)) |>
  ungroup() |>
  inner_join(kapk_cdata) |>
  inner_join(site_order) |>
  mutate(domain = fct_relevel(domain, rev(dom_h)),
         label = fct_reorder(label, site_rnk)) |>
  #filter(sample_type == "Lake sediment") |>
  ggplot(aes(label, prop, fill = domain)) +
  geom_col(width = 1, size = 0) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "top") +
  xlab("Samples") +
  ylab("Proportion") +
  facet_wrap(~member_unit, scales = "free_x") +
  scale_y_continuous(labels = scales::percent) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

dom_h <- c("Eukaryota", "unclassified", "root", "Bacteria", "Archaea", "Viruses")
k2_hires |>
  bind_rows(k2_hires_class |>
              select(label, unclassified) |>
              dplyr::rename(counts = unclassified) |>
              mutate(domain = "unclassified")) |>
  inner_join(kapk_cdata) |>
  group_by(domain, member_unit) |>
  summarise(counts = sum(counts)) |>
  ungroup() |>
  group_by(member_unit) |>
  mutate(prop = counts/sum(counts),
         domain = fct_relevel(domain, rev(dom_h))) |>
  #filter(sample_type == "Lake sediment") |>
  ggplot(aes(domain, prop, fill = domain)) +
  geom_col(width = 1, size = 0) +
  theme_bw() +
  theme(legend.position = "top") +
  xlab("Domain") +
  ylab("Proportion") +
  facet_wrap(~member_unit) +
  scale_y_continuous(labels = scales::percent) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  ggpubr::rotate()


k2_coarse |>
  filter(domain == "Bacteria" | domain == "Archaea" | domain == "Viruses") |>
  rename(counts_coarse = counts) |>
  inner_join(k2_hires) |>
  mutate(diff = counts - counts_coarse) |>
  inner_join(kapk_cdata) |>
  inner_join(site_order) |>
  mutate(domain = fct_relevel(domain, (domains)),
         label = fct_reorder(label, site_rnk)) |>
  ggplot(aes(label, diff)) +
  geom_col(width = 1, size = 0.1, color = "black") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "top") +
  facet_grid(domain~member_unit, scales = "free") +
  xlab("Samples") +
  ylab("Read difference between hires and coarse") +
  scale_y_continuous(labels = scales::comma) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))



# Figure1 -----------------------------------------------------------------
# We want a figure that summarizes the process of taxonomically assigning
# a label to the reads. It will go from:
# initial -> coarse -> hires -> OGUs
# HOw reads are flowing from one step to another

w_stats <- read_tsv("data/woltka/woltka-summary-mapping.tsv.gz") |>
  clean_names() |>
  select(label, n_reads)




k2_coarse_class_sum <- k2_coarse |>
  filter(domain != "root") |>
  mutate(class = ifelse(domain == "Eukaryota", "Eukaryota", "hires")) |>
  group_by(label, class) |>
  dplyr::summarise(counts = sum(counts)) |>
  ungroup() |>
  pivot_wider(names_from = "class", values_from = counts, values_fill = 0) |>
  inner_join(derep_stats |> select(label, num_seqs)) |>
  mutate(hires = num_seqs - Eukaryota)

k2_coarse_abv <- k2_coarse |>
  filter(domain != "root", domain != "Eukaryota") |>
  group_by(label) |>
  dplyr::summarise(rclass = sum(counts)) |>
  ungroup() |>
  inner_join(derep_stats |> select(label, num_seqs)) |>
  mutate(unclassified = num_seqs - rclass)

k2_coarse_euk  <- k2_coarse_class_sum$Eukaryota |> sum()

k2_coarse_u <- k2_coarse_class$unclassified |> sum()


k2_hires_ab <- k2_hires |> filter(domain == "Archaea" | domain == "Bacteria") |> .$counts |> sum()
k2_hires_v <- k2_hires |> filter(domain == "Viruses") |> .$counts |> sum()
k2_hires_u <- (sum(k2_coarse_abv$unclassified) + sum(k2_coarse_abv$rclass)) - (k2_hires_ab + k2_hires_v)



w_hires_class <- w_stats$n_reads |> sum()
w_hires_uclass <- (k2_hires_ab + k2_hires_u) - w_hires_class


data_alluvial <- tribble(~coarse,     ~hires,  ~ogus,    ~counts,
        "Eukaryota", "NA",    "NA",     k2_coarse_euk,
        "Hires",     "a/b", "a/b",    w_hires_class,
        "Hires",     "u", "u",      w_hires_uclass,
        "Hires",     "Vir",   "NA",     k2_hires_v,
) |>
  write_tsv("results/tax_labels_alluvial.tsv")


# How many eukaryotic reads we classify
k2_coarse_euk <- k2_coarse |> filter(domain == "Eukaryota") |> .$counts |> sum()
k2_coarse_abv <- k2_coarse |> filter(domain == "Archaea" | domain == "Bacteria" | domain == "Viruses") |> .$counts |> sum()
k2_coarse_u <- (derep_stats$num_seqs |> sum()) - (k2_coarse_euk + k2_coarse_abv)

k2_hires_v <- k2_hires |> filter(domain == "Viruses") |> .$counts |> sum()
k2_hires_ab <- k2_hires |> filter(domain == "Archaea" | domain == "Bacteria") |> .$counts |> sum()
k2_hires_u <- (k2_coarse_abv + k2_coarse_u) - (k2_coarse_ab + k2_hires_v)

w_hires_class <- w_stats$n_reads |> sum()
w_hires_uclass <- (k2_hires_ab + k2_hires_u) - w_hires_class


data_alluvial <- tribble(~from,        ~to,    ~counts,
                         "coarse-euk", "hires-euk",   k2_coarse_euk,
                         "coarse-abv",  "hires", k2_coarse_abv,
                         "coarse-ucl", "hires", k2_coarse_u,
                         "hires", "hires-class", k2_hires_ab,
                         "hires", "hires-vir",   k2_hires_v,
                         "hires", "hires-uclass", k2_hires_u,
                         "hires-class", "ogu", k2_hires_ab,
                         "hires-uclass", "ogu", k2_hires_u,
                         "ogu", "ogu-class", w_hires_class,
                         "ogu", "ogu-uclass", w_hires_uclass,
) |>
  write_tsv("results/tax_labels_alluvial.tsv")


# Get WOLTKA profiling results --------------------------------------------

w_none <- read_tsv("data/woltka/woltka-summary-none.tsv.gz")
w_free <- read_tsv("data/woltka/woltka-summary-free.tsv.gz")
w_species <- read_tsv("data/woltka/woltka-summary-species.tsv.gz")
w_genus <- read_tsv("data/woltka/woltka-summary-genus.tsv.gz")
w_family <- read_tsv("data/woltka/woltka-summary-family.tsv.gz")
w_order <- read_tsv("data/woltka/woltka-summary-order.tsv.gz")
w_class <- read_tsv("data/woltka/woltka-summary-class.tsv.gz")
w_phylum <- read_tsv("data/woltka/woltka-summary-phylum.tsv.gz")



w_phylum_prop <- w_phylum |>
  select(label, counts, name) |>
  group_by(label) |>
  mutate(prop = counts/sum(counts)) |>
  arrange(desc(prop))

top_phylum <- w_phylum_prop |>
  group_by(name) |>
  summarise(mean_prop = mean(prop)) |>
  arrange(desc(mean_prop)) |>
  filter(mean_prop > 0.0001)
  head(n = 50)

tol21rainbow <- c("#88CCAA", "#44AAAA", "#771155", "#DDDD77", "#DD7788", "#777711",
                  "#44AA77", "#774411", "#77CCCC", "#4477AA", "#77AADD", "#AAAA44",
                  "#CC99BB", "#AA7744", "#117744", "#DDAA77", "#771122", "#AA4455",
                  "#AA4488", "#114477", "#117777")
names(tol21rainbow) <- c(top_phylum$name, "Other")

library(wesanderson)
pal <- wes_palette("Zissou1", 100, type = "continuous")

w_phylum_prop |>
  filter(name %in% top_phylum$name) |>
  #mutate(name = ifelse(name %in% top_phylum$name, name, "Other")) |>
  droplevels() |>
  group_by(name, label) |>
  summarise(prop = sum(prop)) |>
  ungroup() |>
  inner_join(kapk_cdata) |>
  inner_join(site_order) |>
  mutate(name = fct_relevel(name, rev(top_phylum$name)),
         figure_names = fct_reorder(figure_names, site_rnk)) |>
  ggplot(aes(figure_names, name,  fill = prop)) +
  geom_tile()+
  #redrawing tiles to remove cross lines from legend
  geom_tile(colour="gray30",size=0.3, show.legend = FALSE, na.rm = TRUE)+
  #remove axis labels, add title
  labs(x="",y="",title="")+
  #remove extra space
  scale_y_discrete(expand=c(0,0), position = "left")+
  #custom breaks on x-axis
  scale_x_discrete(expand=c(0,0))+
  #custom colours for cut levels and na values
  #scale_fill_gradientn(colours =rev(c("#d53e4f","#f46d43","#fdae61",
  # "#fee08b","#e6f598","#abdda4","#ddf1da")), na.value="grey90", trans = "log10", labels = scales::percent) +
  #viridis::scale_fill_viridis(option = "D", direction = 1, labels = scales::percent, trans = scales::log10_trans()) + #, trans = scales::log10_trans(), breaks = base_breaks()) +
  scale_fill_gradientn(colours = pal, labels = scales::percent, trans = "log10") + #, name="Number of restaurant", guide = guide_legend( keyheight = unit(3, units= "mm"), keywidth=unit(12, units = "mm"), label.position = "bottom", title.position = 'top', nrow=1) ) +

  #labels=c("0k", "1k", "2k", "3k", "4k")) +
  #scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(9,"YlGnBu")),na.value="grey90", trans = scales::log_trans(), breaks = base_breaks(), labels = scales::percent) +
  #mark year of vaccination
  #geom_vline(aes(xintercept = 36),size=3.4,alpha=0.24)+
  #equal aspect ratio x and y axis
  #coord_fixed() +
  #set base size for all font elements
  theme_grey(base_size=10)+
  #theme options
  theme(
    legend.position = "right",
    #remove legend title
    legend.title=element_blank(),
    #remove legend margin
    legend.spacing = grid::unit(0,"cm"),
    #change legend text properties
    legend.text=element_text(size=7,face="bold"),
    #change legend key height
    #legend.key.height=grid::unit(0.8,"cm"),
    #set a slim legend
    #legend.key.width=grid::unit(0.2,"cm"),
    #set x axis text size and colour
    axis.text.x=element_text(hjust = 1, vjust = 0.5, angle = 90, size = 6),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.ticks.y = element_blank(),
    #axis.text.y = element_blank(),
    #set y axis text colour and adjust vertical justification
    #axis.text.y=element_text(vjust = 0.2,colour=textcol, size = 6),
    #change axis ticks thickness
    axis.ticks=element_line(size=0.4),
    #change title font, size, colour and justification
    plot.title=element_blank(),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank()) +
  facet_wrap(~member_unit, scales = "free_x")




base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}



# PhilR -------------------------------------------------------------------
library(philr)
library(phyloseq)
source("libs/phil_var.R")
library(ape)
library(phytools)
# Create phyloseq object
ogu_df <- w_none |>
  select(feature, label, counts) |>
  pivot_wider(names_from = "label", values_from = counts, values_fill = 0) |>
  as.data.frame() |>
  column_to_rownames("feature")
kapk_cdata_df <- kapk_cdata |>
  mutate(label_cl = label) |>
  as.data.frame() |>
  column_to_rownames("label_cl")

gtdb_tax_bac <- read_tsv("data/woltka/bac120_taxonomy_r202.tsv", col_names = c("genome", "tax_string")) |>
  tidyr::separate(col = tax_string, sep = ";", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
gtdb_tax_arc <- read_tsv("data/woltka/ar122_taxonomy_r202.tsv", col_names = c("genome", "tax_string")) |>
  tidyr::separate(col = tax_string, sep = ";", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

gtdb_taxonomy <- gtdb_tax_bac |>
  bind_rows(gtdb_tax_arc) |>
  distinct() |>
  as.data.frame() |>
  column_to_rownames("genome")


gtdb_tree_bac <- read.newick("data/woltka/bac120_r202.tree")
gtdb_tree_arc <- read.newick("data/woltka/ar122_r202.tree")
gtdb_tree <- bind.tree(gtdb_tree_bac, gtdb_tree_arc, where = "root")

kapk_ps <- phyloseq(otu_table(ogu_df, taxa_are_rows = TRUE),
         tax_table(as.matrix(gtdb_taxonomy)),
         sample_data(kapk_cdata_df),
         phy_tree(gtdb_tree))

kapk_philr <- get_philr(kapk_ps)
kapk_philr_g <- get_philr(kapk_ps, glom="Genus")

df1 <- df |> as_tibble()
df1$label_orig <- NULL
hc <- hclust(kapk_dist, "ward.D2")
ggtree(hc, linetype='dashed', aes(color = member_unit)) %<+% df1 +
  layout_dendrogram() +
  geom_tippoint(size=5, shape=21, aes(fill=site, x=x+.5), color='black') +
  scale_fill_manual(values = tol21rainbow)


df <- as(sample_data(kapk_ps_filt), "data.frame")

perm_tina_mp <- adonis(as.dist(as.matrix(kapk_dist)[rownames(df), rownames(df)]) ~ site, data = df)
perm_tina_mp <- adonis(kapk_dist ~ member_unit, data = as(sample_data(kapk_ps_filt), "data.frame"))

perm_tina_wt <- adonis(as.dist(osd2014_pina_tina_results[[2]]$cs[selected_samples,selected_samples]) ~ water_temperature, data = osd2014_cdata_df)


library(glmnet);

sample_data(kapk_ps_filt)$is_B1 <- factor(get_variable(kapk_ps_filt, "member_unit") %in% c("B1"))

glmmod <- glmnet(kapk_philr, sample_data(kapk_ps_filt)$is_B1, family="binomial")

coefficients(glmmod)
plot(glmmod, label = TRUE)
print(glmmod)

top.coords <- as.matrix(coef(glmmod, s=glmmod$lambda[15]))
top.coords <- rownames(top.coords)[which(top.coords != 0)]
(top.coords <- top.coords[2:length(top.coords)]) # remove the intercept as a coordinate

tc.names <- sapply(top.coords, function(x) name.balance(tree, tax, x))
tc.names

kapk_philr_long <- convert_to_long(kapk_philr, get_variable(kapk_ps_filt, 'is_B1')) |>
  filter(coord %in% top.coords)

ggplot(kapk_philr_long, aes(x=labels, y=value)) +
  geom_boxplot(fill='lightgrey') +
  facet_wrap(.~coord, scales='free') +
  xlab('Human') + ylab('Balance Value') +
  theme_bw()

kapk_philr_long |>
  rename(is_B3=labels) |>
  filter(coord %in% c('n194', 'n214')) |>
  spread(coord, value) |>
  ggplot(aes(x=n194, y=n214, color=is_B3)) +
  geom_point(size=4) +
  xlab(tc.names['n194']) + ylab(tc.names['n214']) +
  theme_bw()

library(ggtree); packageVersion("ggtree")

tc.nn <- name.to.nn(tree, top.coords)
tc.colors <- tol21rainbow
p <- ggtree(tree, layout='fan') +
  geom_balance(node=tc.nn[1], fill=tc.colors[1], alpha=0.6) +
  geom_balance(node=tc.nn[2], fill=tc.colors[2], alpha=0.6) +
  geom_balance(node=tc.nn[3], fill=tc.colors[3], alpha=0.6) +
  geom_balance(node=tc.nn[4], fill=tc.colors[4], alpha=0.6) +
  geom_balance(node=tc.nn[5], fill=tc.colors[5], alpha=0.6)
p <- annotate_balance(tree, 'n950', p=p, labels = c('n220+', 'n220-'),
                      offset.text=0.15, bar=FALSE)
annotate_balance(tree, 'n602', p=p, labels = c('n602+', 'n602-'),
                 offset.text=0.15, bar=FALSE)


sample_data(kapk_ps_filt)$is_B2 <- factor(get_variable(kapk_ps_filt, "member_unit") %in% c("B2"))

glmmod <- glmnet(kapk_philr, sample_data(kapk_ps_filt)$is_B2, family="binomial")

plot(glmmod, label = TRUE)
print(glmmod)

top.coords <- as.matrix(coef(glmmod, s=0.037040))
top.coords <- rownames(top.coords)[which(top.coords != 0)]
(top.coords <- top.coords[2:length(top.coords)]) # remove the intercept as a coordinate

tc.names <- sapply(top.coords, function(x) name.balance(tree, tax, x))
tc.names

votes <- name.balance(tree, tax, 'n883', return.votes = c('up', 'down'))
votes[[c('up.votes', 'Family')]]

library(ggtree); packageVersion("ggtree")

tc.nn <- name.to.nn(tree, top.coords)
tc.colors <- tol21rainbow
p <- ggtree(tree, layout='fan') +
  geom_balance(node=tc.nn[1], fill=tc.colors[1], alpha=0.6) +
  geom_balance(node=tc.nn[2], fill=tc.colors[2], alpha=0.6) +
  geom_balance(node=tc.nn[3], fill=tc.colors[3], alpha=0.6) +
  geom_balance(node=tc.nn[4], fill=tc.colors[4], alpha=0.6)
p <- annotate_balance(tree, 'n829', p=p, labels = c('n829+', 'n829-'),
                      offset.text=0.15, bar=FALSE)
annotate_balance(tree, 'n883', p=p, labels = c('n883+', 'n883-'),
                 offset.text=0.15, bar=FALSE)



sample_data(kapk_ps_filt)$is_B3 <- factor(get_variable(kapk_ps_filt, "member_unit") %in% c("B3"))

glmmod <- glmnet(kapk_philr, sample_data(kapk_ps_filt)$is_B3, family="binomial")

plot(glmmod, label = TRUE)
print(glmmod)

top.coords <- as.matrix(coef(glmmod, s=0.21830))
top.coords <- rownames(top.coords)[which(top.coords != 0)]
(top.coords <- top.coords[2:length(top.coords)]) # remove the intercept as a coordinate

tc.names <- sapply(top.coords, function(x) name.balance(tree, tax, x))
tc.names

votes <- name.balance(tree, tax, 'n214 ', return.votes = c('up', 'down'))
votes[[c('up.votes', 'Family')]]

library(ggtree); packageVersion("ggtree")

tc.nn <- name.to.nn(tree, top.coords)
tc.colors <- tol21rainbow
p <- ggtree(tree, layout='fan') +
  geom_balance(node=tc.nn[1], fill=tc.colors[1], alpha=0.6) +
  geom_balance(node=tc.nn[2], fill=tc.colors[2], alpha=0.6) +
  geom_balance(node=tc.nn[3], fill=tc.colors[3], alpha=0.6) +
  geom_balance(node=tc.nn[4], fill=tc.colors[4], alpha=0.6)
p <- annotate_balance(tree, 'n214', p=p, labels = c('n214+', 'n214-'),
                      offset.text=0.15, bar=FALSE)
annotate_balance(tree, 'n194', p=p, labels = c('n194+', 'n194-'),
                 offset.text=0.15, bar=FALSE)

kapk_philr_long <- convert_to_long(kapk_philr, get_variable(kapk_ps_filt, 'is_B3')) |>
  filter(coord %in% top.coords)

ggplot(kapk_philr_long, aes(x=labels, y=value)) +
  geom_boxplot(fill='lightgrey') +
  facet_grid(.~coord, scales='free_x') +
  xlab('Human') + ylab('Balance Value') +
  theme_bw()

kapk_philr_long |>
  rename(is_B3=labels) |>
  filter(coord %in% c('n194', 'n214')) |>
  spread(coord, value) |>
  ggplot(aes(x=n194, y=n214, color=is_B3)) +
  geom_point(size=4) +
  xlab(tc.names['n194']) + ylab(tc.names['n214']) +
  theme_bw()


# Sourcetracker -----------------------------------------------------------
library(phyloseq)
w_none_src <- read_tsv("data/mt-st/woltka-summary-none.tsv.gz")
src_data <- read_tsv("data/mt-st/mt-samples.txt") |>
  clean_names() |>
  rename(Env = biome) |>
  mutate(SourceSink = "source") |>
  select(label, SourceSink, Env)

ogu_src_df <- w_none_src |>
  select(feature, label, counts) |>
  bind_rows(w_none |> select(feature, label, counts)) |>
  pivot_wider(names_from = "label", values_from = counts, values_fill = 0) |>
  as.data.frame() |>
  column_to_rownames("feature")

ogu_src_df_cdata <- kapk_cdata |>
  select(label, member_unit) |>
  mutate(SourceSink = "sink") |>
  rename(Env = member_unit) |>
  select(label, SourceSink, Env) |>
  bind_rows(src_data) |>
  as.data.frame() |>
  column_to_rownames("label")

st_ps <- phyloseq(otu_table(ogu_src_df, taxa_are_rows = TRUE),
                    tax_table(as.matrix(gtdb_taxonomy)),
                    sample_data(ogu_src_df_cdata),
                    phy_tree(gtdb_tree))

st_ps_g <- speedyseq::tax_glom(st_ps, "Genus")

tax_filt <- speedyseq::psmelt(st_ps_g) |>
  dplyr::select(Sample, OTU, Abundance) |>
  group_by(Sample) |>
  mutate(prop = Abundance/sum(Abundance)) |>
  ungroup() |>
  group_by(OTU) |>
  summarise(mean_prop = sum(prop)/nsamples(st_ps_g)) |>
  ungroup() |>
  filter(mean_prop >= 1e-5) |> .$OTU

st_ps_filt <- prune_taxa(tax_filt, st_ps_g)
st_ps_filt

library(biomformat)
ogu_st_biom <- biomformat::make_biom(
  data = as((otu_table(ogu_src_df, taxa_are_rows = FALSE)), "matrix"),
  matrix_element_type = "int"
)
write_biom(ogu_st_biom, biom_file = "results/ogu-st.biom")

kapk_cdata |>
  select(label, member_unit) |>
  mutate(SourceSink = "sink") |>
  rename(Env = member_unit) |>
  select(label, SourceSink, Env) |>
  bind_rows(src_data) |>
  rename("#SampleID" = label) |>
  write_tsv("results/ogu-st.map")
