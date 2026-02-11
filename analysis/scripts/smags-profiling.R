library(tidyverse)





smags_metadata <- read_tsv("data/map-smags/SMAGS-anvio/metadata.txt") |>
  select(Genome_Id, contains("Best_taxonomy"), Arctic_119_Stations, Cosmopolitan_Score_119_Stations, Taxa_Super_Groups)
# janitor::clean_names()
smags_stats <- read_tsv(file = "data/map-smags/smags-mapping-filtered.summary.tsv.gz") |>
  filter(breadth_exp_ratio >= 0.5, n_reads > 100) |>
  select(reference, n_reads, label, breadth_exp_ratio) |>
  rename(Genome_Id = reference) |>
  inner_join(kapk_cdata)
smags_tree <- ape::read.tree("data/map-smags/SMAGS-anvio/Conca_allMETDB_allSMAGs_RNAP_200aa_no_duplicate_NUCLEAR_noRhodophyta.fa.treefile")
library(ape)

# Get taxa
smags_tax <- read_tsv("https://raw.githubusercontent.com/chassenr/taxdb-integration/master/assets/TARA_SMAGs.txt") |>
  select(accession, taxon) |>
  tidyr::separate(col = taxon, sep = ";", into = c("domain", "phylum", "class", "order", "family", "genus", "species")) |>
  rename(Genome_Id = accession)

# Let;s do a quick fix of the taxonomies
smags_tax |>
  left_join(smags_metadata |> select(Genome_Id, Taxa_Super_Groups)) |>
  mutate(genus = ifelse(genus == "g__Eukaryota" & !is.na(Taxa_Super_Groups), paste0(genus, "_", Taxa_Super_Groups), genus)) |>
  mutate(family = ifelse(family == "f__Eukaryota" & !is.na(Taxa_Super_Groups), paste0(family, "_", Taxa_Super_Groups), family)) |>
  mutate(order = ifelse(order == "o__Eukaryota" & !is.na(Taxa_Super_Groups), paste0(order, "_", Taxa_Super_Groups), order)) |>
  mutate(class = ifelse(class == "c__Eukaryota" & !is.na(Taxa_Super_Groups), paste0(class, "_", Taxa_Super_Groups), class)) |>
  unite(domain, phylum, class, order, family, genus, species, col = "tax_path", sep = ";") |>
  select(-Taxa_Super_Groups) |>
  write_tsv(file = "data/map-smags/SMAGS-anvio/smags-tax-fixed.tsv", col_names = F)

# Not all SMAGS are present in the tree because they didn't have RNApol, let's do a bit of a hack
# We will use a anchor the taxonomic information from Tom's best taxonomy. We will first will
# try to get as many in the lower ranks as we can. Starting from Genus to Kingdom
# For example:
# we will identify those smags not present in the tree, get their genera and find out where in the tree they are
# then create a subtree using those tips and get the node that it closer to the root.
# We will place in this node the SMAGS that are missing

# Which genera have a taxonomy
miss_gen <- smags_stats |>
  select(Genome_Id) |>
  distinct() |>
  inner_join(smags_metadata) |>
  distinct() |>
  filter(!(Genome_Id %in% smags_tree$tip.label)) |>
  select(Best_taxonomy_GENRE, Genome_Id) |>
  filter(!is.na(Best_taxonomy_GENRE)) |>
  distinct() |>
  rename(tax = Best_taxonomy_GENRE)

miss_fam <- smags_stats |>
  select(Genome_Id) |>
  distinct() |>
  inner_join(smags_metadata) |>
  distinct() |>
  filter(
    !(Genome_Id %in% smags_tree$tip.label),
    !(Genome_Id %in% miss_gen$Genome_Id)
  ) |>
  select(Best_taxonomy_FAMILY, Genome_Id) |>
  filter(!is.na(Best_taxonomy_FAMILY)) |>
  distinct() |>
  rename(tax = Best_taxonomy_FAMILY)


miss_class <- smags_stats |>
  select(Genome_Id) |>
  distinct() |>
  inner_join(smags_metadata) |>
  distinct() |>
  filter(
    !(Genome_Id %in% smags_tree$tip.label),
    !(Genome_Id %in% miss_gen$Genome_Id),
    !(Genome_Id %in% miss_fam$Genome_Id)
  ) |>
  select(Best_taxonomy_CLASS, Genome_Id) |>
  filter(!is.na(Best_taxonomy_CLASS)) |>
  distinct() |>
  rename(tax = Best_taxonomy_CLASS)


miss_phyl <- smags_stats |>
  select(Genome_Id) |>
  distinct() |>
  inner_join(smags_metadata) |>
  distinct() |>
  filter(
    !(Genome_Id %in% smags_tree$tip.label),
    !(Genome_Id %in% miss_gen$Genome_Id),
    !(Genome_Id %in% miss_fam$Genome_Id),
    !(Genome_Id %in% miss_class$Genome_Id)
  ) |>
  select(Best_taxonomy_PHYLUM, Genome_Id) |>
  filter(!is.na(Best_taxonomy_PHYLUM)) |>
  distinct() |>
  rename(tax = Best_taxonomy_PHYLUM)


miss_king <- smags_stats |>
  select(Genome_Id) |>
  distinct() |>
  inner_join(smags_metadata) |>
  distinct() |>
  filter(
    !(Genome_Id %in% smags_tree$tip.label),
    !(Genome_Id %in% miss_gen$Genome_Id),
    !(Genome_Id %in% miss_fam$Genome_Id),
    !(Genome_Id %in% miss_class$Genome_Id),
    !(Genome_Id %in% miss_phyl$Genome_Id)
  ) |>
  select(Best_taxonomy_KINGDON, Genome_Id) |>
  filter(!is.na(Best_taxonomy_KINGDON)) |>
  distinct() |>
  rename(tax = Best_taxonomy_KINGDON)


in_tree <- smags_stats |>
  select(Genome_Id) |>
  distinct() |>
  inner_join(smags_metadata) |>
  distinct() |>
  filter((Genome_Id %in% smags_tree$tip.label))


to_drop <- setdiff(smags_tree$tip.label, unique(smags_stats$Genome_Id))

length(to_drop)

tree_filt <- drop.tip(smags_tree, tip = to_drop)
tree_filt <- smags_tree
get_tip_info <- function(X, df, phy_tree) {
  tips <- df |>
    filter(tax == X)
  to_drop <- setdiff(phy_tree$tip.label, unique(tips$Genome_Id))
  phy_tree_filt <- drop.tip(phy_tree, tip = to_drop)
  nodes <- distRoot(phy_tree_filt, tips = "all") |>
    sort() |>
    head(1) |>
    enframe(name = "tip", value = "dist") |>
    mutate(tax = X)
  node_n <- which(phy_tree$tip.label == nodes$tip)
  nodes |> mutate(node = node_n)
}
df <- in_tree |>
  select(Genome_Id, Best_taxonomy_GENRE) |>
  rename(tax = Best_taxonomy_GENRE)

gen_tips <- map_dfr(miss_gen |>
  select(tax) |>
  distinct() |> pull(tax), get_tip_info, df = df, phy_tree = tree_filt) |>
  inner_join(miss_gen)

df <- in_tree |>
  select(Genome_Id, Best_taxonomy_FAMILY) |>
  rename(tax = Best_taxonomy_FAMILY)

fam_tips <- map_dfr(miss_fam |>
  select(tax) |>
  distinct() |> pull(tax), get_tip_info, df = df, phy_tree = tree_filt) |>
  inner_join(miss_fam)


df <- in_tree |>
  select(Genome_Id, Best_taxonomy_CLASS) |>
  rename(tax = Best_taxonomy_CLASS)

class_tips <- map_dfr(miss_class |>
  select(tax) |>
  distinct() |> pull(tax), get_tip_info, df = df, phy_tree = tree_filt) |>
  inner_join(miss_class)


df <- in_tree |>
  select(Genome_Id, Best_taxonomy_PHYLUM) |>
  rename(tax = Best_taxonomy_PHYLUM)

phyl_tips <- map_dfr(miss_phyl |>
  select(tax) |>
  distinct() |> pull(tax), get_tip_info, df = df, phy_tree = tree_filt) |>
  inner_join(miss_phyl)

df <- in_tree |>
  select(Genome_Id, Best_taxonomy_KINGDON) |>
  rename(tax = Best_taxonomy_KINGDON)

king_tips <- map_dfr(miss_king |>
  select(tax) |>
  distinct() |> pull(tax), get_tip_info, df = df, phy_tree = tree_filt) |>
  inner_join(miss_king)


df_all <- bind_rows(gen_tips, fam_tips, class_tips, phyl_tips, king_tips)

new_tree <- smags_tree
for (i in 1:nrow(df_all)) {
  new_tree <- bind.tip(new_tree,
    tip.label = df_all[i, ]$Genome_Id,
    where = df_all[i, ]$node,
    edge.length = df_all[i, ]$dist
  )
}



to_drop <- setdiff(new_tree$tip.label, unique(smags_stats$Genome_Id))

smags_tree_filt <- drop.tip(new_tree, tip = to_drop)
write.tree(smags_tree_filt, file = "data/map-smags/SMAGS-anvio/smags-tree.tree")


# Read damage estimates
mtdmg_w0 <- read_tsv(file = "data/map-smags/smags-mdmg.summary.0.tsv.gz")
mtdmg_w1 <- read_tsv(file = "data/map-smags/smags-mdmg.summary.1.tsv.gz")

# Keep only those with a sigma and model
mtdmg_w1_filt <- mtdmg_w1 |>
  filter(N_reads > 100) |>
  filter(Bayesian_z >= 1.5, Bayesian_D_max >= 0.25) |>
  inner_join(kapk_cdata |> select(label, member_unit)) |>
  filter(tax_rank == "species") |>
  dplyr::select(tax_name, Bayesian_D_max, Bayesian_z, label) |>
  separate(col = tax_name, into = c("species", "Genome_Id"), sep = " ") |>
  select(-species)

dmg <- smags_stats |>
  # inner_join(smags_tax |> select(Genome_Id, genus) |> mutate(genus = gsub("g__", "",genus))) |>
  left_join(mtdmg_w1_filt) |>
  filter(Genome_Id %in% smags_tree_filt$tip.label) |>
  mutate(Bayesian_D_max = ifelse(is.na(Bayesian_D_max), 0, Bayesian_D_max)) |>
  mutate(member_unit = paste0(member_unit, "_Bayesian_D_max")) |>
  select(Genome_Id, member_unit, Bayesian_D_max) |>
  group_by(Genome_Id, member_unit) |>
  summarise(avg_D_max = max(Bayesian_D_max)) |>
  # mutate(avg_D_max = ifelse(avg_D_max < 0.1, 0, avg_D_max)) |>
  pivot_wider(names_from = member_unit, values_from = avg_D_max, values_fill = 0)


n_reads <- smags_stats |>
  select(Genome_Id, member_unit, n_reads) |>
  group_by(Genome_Id, member_unit) |>
  summarise(avg_mapped_reads = sum(n_reads)) |>
  mutate(member_unit = paste0(member_unit, "_nreads")) |>
  # select(Genome_Id, coverage_mean, ) |>
  filter(Genome_Id %in% smags_tree_filt$tip.label) |>
  pivot_wider(names_from = member_unit, values_from = avg_mapped_reads, values_fill = 0)

exp_br <- smags_stats |>
  select(Genome_Id, member_unit, breadth_exp_ratio) |>
  group_by(Genome_Id, member_unit) |>
  summarise(avg_breadth_exp_ratio = mean(breadth_exp_ratio)) |>
  mutate(member_unit = paste0(member_unit, "_breadth_exp_ratio")) |>
  # select(Genome_Id, coverage_mean, ) |>
  filter(Genome_Id %in% smags_tree_filt$tip.label) |>
  pivot_wider(names_from = member_unit, values_from = avg_breadth_exp_ratio, values_fill = 0)


n_reads |>
  inner_join(exp_br) |>
  inner_join(dmg) |>
  inner_join(smags_metadata) |>
  write_tsv("data/map-smags/SMAGS-anvio/smags-cov-dmg-ebr.tsv")


#
# mtdmg_w1_filt <- mtdmg_w1 |>
#   filter(N_reads > 100) |>
#   filter(Bayesian_z >= 1.5) |>
#   inner_join(kapk_cdata |> select(label, member_unit)) |>
#   filter(tax_rank == "species") |>
#   select(tax_name, Bayesian_D_max, Bayesian_z, label) |>
#   separate(col = tax_name, into = c("species", "Genome_Id"), sep = " ") |>
#   select(-species) |>
#   inner_join(pydamage |> filter(damage_model_p <= 0.6))
#
# pal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")
# mtdmg_w1_filt_1 <- mtdmg_w1 |>
#   filter(N_reads > 100) |>
#   filter(Bayesian_z >= 1.5, Bayesian_q > 0.5) |>
#   inner_join(kapk_cdata |> select(label, member_unit)) |>
#   filter(tax_rank == "genus") |>
#   arrange(desc(Bayesian_D_max))
#
# mtdmg_w0_filt |>
#   group_by(member_unit) |>
#   count(tax_rank) |> View()
#
# mtdmg_w0_filt |> pull(tax_name) |> unique()
#
#
get_dmg_decay <- function(df, tax, orient = "fwd") {
  # tax <- gsub("g__", "", tax)
  # df <- df |>
  #   filter(grepl(tax, tax_path))

  df_k_fwd <- df |>
    select(tax_name, label, starts_with("k+")) |>
    pivot_longer(names_to = "type", values_to = "k_fwd", c(-tax_name, -label)) |>
    mutate(x = gsub("k\\+", "", type)) |>
    select(-type)

  df_k_rev <- df |>
    select(tax_name, label, starts_with("k-")) |>
    pivot_longer(names_to = "type", values_to = "k_rev", c(-tax_name, -label)) |>
    mutate(x = gsub("k\\-", "", type)) |>
    select(-type)

  df_N_fwd <- df |>
    select(tax_name, label, starts_with("N+")) |>
    pivot_longer(names_to = "type", values_to = "N_fwd", c(-tax_name, -label, )) |>
    mutate(x = gsub("N\\+", "", type)) |>
    select(-type)

  df_N_rev <- df |>
    select(tax_name, label, starts_with("N-")) |>
    pivot_longer(names_to = "type", values_to = "N_rev", c(-tax_name, -label)) |>
    mutate(x = gsub("N\\-", "", type)) |>
    select(-type)

  df_fit_fwd <- df |>
    select(tax_name, label, starts_with("f+")) |>
    pivot_longer(names_to = "type", values_to = "f_fwd", c(-tax_name, -label)) |>
    mutate(x = gsub("f\\+", "", type)) |>
    select(-type)

  df_fit_rev <- df |>
    select(tax_name, label, starts_with("f-")) |>
    pivot_longer(names_to = "type", values_to = "f_rev", c(-tax_name, -label)) |>
    mutate(x = gsub("f\\-", "", type)) |>
    select(-type)

  dat <- df_k_fwd |>
    inner_join(df_k_rev) |>
    inner_join(df_N_fwd) |>
    inner_join(df_N_rev) |>
    inner_join(df_fit_fwd) |>
    inner_join(df_fit_rev) |>
    inner_join(df |> select(label, member_unit) |> distinct()) |>
    mutate(x = as.numeric(x)) |>
    filter(x <= 25)

  fwd_max <- dat |>
    group_by(as.character(x)) |>
    summarise(val = mean(f_fwd) + sd(f_fwd)) |>
    pull(val) |>
    max()
  fwd_min <- dat |>
    group_by(as.character(x)) |>
    summarise(val = mean(f_fwd) - sd(f_fwd)) |>
    pull(val) |>
    min()
  rev_max <- dat |>
    group_by(as.character(x)) |>
    summarise(val = mean(f_rev) + sd(f_rev)) |>
    pull(val) |>
    max()
  rev_min <- dat |>
    group_by(as.character(x)) |>
    summarise(val = mean(f_rev) - sd(f_rev)) |>
    pull(val) |>
    min()


  y_max <- ifelse(fwd_max >= rev_max, fwd_max, rev_max)
  y_min <- ifelse(fwd_min <= rev_min, fwd_min, rev_min)
  if (orient == "fwd") {
    # model <- nls(f_fwd ~ SSasymp(x,Asym, R0, lrc), data= dat )
    # new.data <- data.frame(x=seq(1, 25, by = 1))
    # interval <- as_tibble(predFit(model, newdata = new.data, interval = "confidence", level= 0.9)) |>
    #   mutate(x = new.data$x)
    #

    ggplot() +
      # geom_line(data = dat, aes(x, f_fwd, group=interaction(tax_name, label)), alpha = 0.3, color = "black", size = 0.2) +
      # geom_ribbon(data=interval, aes(x=x, ymin=lwr, ymax=upr), alpha=0.5, inherit.aes=F, fill="grey")+
      # geom_line(data = dat, aes(CtoT, prop, group=interaction(Genome_Id, label)),alpha = 0.3) +
      stat_summary(data = dat, aes(x, f_fwd), fun.data = mean_sd, geom = "ribbon", alpha = .3, size = 0) +
      stat_summary(data = dat, aes(x, f_fwd), fun = mean, geom = "line", size = 0.8) +
      # stat_smooth(data = dat, aes(x, f_rev), method = "nls", formula = y ~ SSasymp(x, Asym, R0, lrc), se = FALSE, color = "#E84635") +
      xlab("Position") +
      ylab("Frequency") +
      # expand_limits(y=y_max) +
      scale_y_continuous(limits = c(y_min, y_max)) +
      # scale_x_continuous(trans = "reverse") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  } else {
    # model <- nls(f_rev ~ SSasymp(x,Asym, R0, lrc), data= dat )
    # new.data <- data.frame(x=seq(1, 25, by = 1))
    # interval <- as_tibble(predFit(model, newdata = new.data, interval = "confidence", level= 0.9)) |>
    #   mutate(x = new.data$x)
    #

    ggplot() +
      # geom_line(data = dat, aes(x, f_rev, group=interaction(tax_name, label)), alpha = 0.3, color = "black", size = 0.2) +
      # geom_ribbon(data=interval, aes(x=x, ymin=lwr, ymax=upr), alpha=0.5, inherit.aes=F, fill="grey")+
      # geom_line(data = dat, aes(CtoT, prop, group=interaction(Genome_Id, label)),alpha = 0.3) +
      stat_summary(data = dat, aes(x, f_rev), fun.data = mean_sd, geom = "ribbon", alpha = .3, size = 0) +
      stat_summary(data = dat, aes(x, f_rev), fun = mean, geom = "line", size = 0.8) +
      # stat_smooth(data = dat, aes(x, f_rev), method = "nls", formula = y ~ SSasymp(x, Asym, R0, lrc), se = FALSE, color = "#E84635") +
      xlab("Position") +
      ylab("Frequency") +
      scale_x_continuous(trans = "reverse") +
      scale_y_continuous(limits = c(y_min, y_max), position = "right") +
      # scale_x_continuous(trans = "reverse") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }

  # ggplot(aes(as.numeric(x), f_fwd)) +
  #   geom_point(fill = "#E84635", shape = 21, alpha = 0.6, color = "black", size = 3) +
  #   # geom_smooth(method = "nls", colour = "red", formula=y ~ exp(a + b * x),
  #   #             method.args = list(start = c(a = 2, b = -1)), se = FALSE) +
  #   stat_smooth(method = "nls", formula = y ~ SSasymp(x, Asym, R0, lrc), se = FALSE, color = "#E84635")
  #
  # p1 +
  #   geom_point(aes(as.numeric(x), f_rev), fill = "#1C3C52", shape = 21, alpha = 0.6, color = "black", size = 3) +
  #   stat_smooth(method = "nls", formula = y ~ SSasymp(x, Asym, R0, lrc), se = FALSE, color = "#E84635") +
  #   ylab("k/N") +
  #   xlab("|x|")
}


tax_superg <- smags_tax |>
  inner_join(smags_metadata) |>
  select(Genome_Id, Taxa_Super_Groups, species) |>
  mutate(species = gsub("s__", "", species))

tax_g_list <- tax_superg |>
  pull(Taxa_Super_Groups) |>
  unique()


# Get Genera estimates
mtdmg_w1_filt_genus <- mtdmg_w1 |>
  filter(N_reads > 100) |>
  filter(Bayesian_z >= 1.5, Bayesian_D_max >= 0.25) |>
  inner_join(kapk_cdata |> select(label, member_unit)) |>
  filter(tax_rank == "species") |>
  arrange(desc(Bayesian_D_max))

nrank <- "genus"


X <- "Opisthokonta"

purrr::map(tax_g_list, function(X) {
  cat(X)
  nrank <- "species"
  sel_tax <- mtdmg_w1_filt_genus |>
    filter(tax_name %in% (tax_superg |> filter(Taxa_Super_Groups == X) |> pull(species))) |>
    filter(tax_rank == nrank) |>
    filter(`f+1` >= 0.1, `f-1` >= 0.1) |>
    select(tax_name, label) |>
    distinct() |>
    arrange(tax_name)
  if (nrow(sel_tax) > 0) {
    n_readsa <- mtdmg_w1_filt_genus |>
      inner_join(sel_tax) |>
      filter(tax_rank == nrank) |>
      pull(N_reads) |>
      sum()
    ggpubr::ggarrange(plotlist = list(
      get_dmg_decay(df = mtdmg_w1_filt_genus |> inner_join(sel_tax) |> filter(tax_rank == nrank), tax = tax, orient = "fwd") +
        ggtitle(paste0(X, " ", n_readsa, " FWD")),
      get_dmg_decay(df = mtdmg_w1_filt_genus |> inner_join(sel_tax) |> filter(tax_rank == nrank), tax = tax, orient = "rev") +
        ggtitle(paste0(X, " ", n_readsa, " REV"))
    ), align = "hv")
    ggsave(paste0("figures/", X, "-dmg.pdf"), plot = last_plot(), width = 8, height = 4)
  }
})

tax <- "Ascomycota"
sel_tax <- mtdmg_w1_filt_genus |>
  filter(grepl(tax, tax_path)) |>
  filter(tax_rank == nrank) |>
  filter(`f+1` >= 0.1, `f-1` >= 0.1) |>
  select(tax_name, label) |>
  distinct() |>
  arrange(tax_name)
n_readsa <- mtdmg_w1_filt_genus |>
  inner_join(sel_tax) |>
  filter(grepl(tax, tax_path)) |>
  pull(N_reads) |>
  sum()
ggpubr::ggarrange(plotlist = list(
  get_dmg_decay(df = mtdmg_w1_filt_genus |> inner_join(sel_tax), tax = tax, orient = "fwd") +
    ggtitle(paste0(tax, " ", n_readsa, " FWD")),
  get_dmg_decay(df = mtdmg_w1_filt_genus |> inner_join(sel_tax), tax = tax, orient = "rev") +
    ggtitle(paste0(tax, " ", n_readsa, " REV"))
))


tax <- "Chloropicophyceae"
sel_tax <- mtdmg_w1_filt_genus |>
  filter(grepl(tax, tax_path)) |>
  filter(tax_rank == nrank) |>
  filter(`f+1` >= 0.1, `f-1` >= 0.1) |>
  select(tax_name, label) |>
  distinct() |>
  arrange(tax_name)
n_readsa <- mtdmg_w1_filt_genus |>
  inner_join(sel_tax) |>
  filter(grepl(tax, tax_path)) |>
  pull(N_reads) |>
  sum()
ggpubr::ggarrange(plotlist = list(
  get_dmg_decay(df = mtdmg_w1_filt_genus |> inner_join(sel_tax), tax = tax, orient = "fwd") +
    ggtitle(paste0(tax, " ", n_readsa, " FWD")),
  get_dmg_decay(df = mtdmg_w1_filt_genus |> inner_join(sel_tax), tax = tax, orient = "rev") +
    ggtitle(paste0(tax, " ", n_readsa, " REV"))
))

tax <- "Haptophyta"
sel_tax <- mtdmg_w1_filt_genus |>
  filter(grepl(tax, tax_path)) |>
  filter(tax_rank == nrank) |>
  filter(`f+1` >= 0.1, `f-1` >= 0.1) |>
  select(tax_name, label) |>
  distinct() |>
  arrange(tax_name)
n_readsa <- mtdmg_w1_filt_genus |>
  inner_join(sel_tax) |>
  filter(grepl(tax, tax_path)) |>
  pull(N_reads) |>
  sum()
ggpubr::ggarrange(plotlist = list(
  get_dmg_decay(df = mtdmg_w1_filt_genus |> filter(tax_rank == nrank) |> inner_join(sel_tax), tax = tax, orient = "fwd") +
    ggtitle(paste0(tax, " ", n_readsa, " FWD")),
  get_dmg_decay(df = mtdmg_w1_filt_genus |> filter(tax_rank == nrank) |> inner_join(sel_tax), tax = tax, orient = "rev") +
    ggtitle(paste0(tax, " ", n_readsa, " REV"))
))


tax <- "MAST"
sel_tax <- mtdmg_w1_filt_genus |>
  filter(grepl(tax, tax_path)) |>
  filter(tax_rank == nrank) |>
  filter(`f+1` >= 0.1, `f-1` >= 0.1) |>
  select(tax_name, label) |>
  distinct() |>
  arrange(tax_name)
n_readsa <- mtdmg_w1_filt_genus |>
  inner_join(sel_tax) |>
  filter(grepl(tax, tax_path)) |>
  pull(N_reads) |>
  sum()
ggpubr::ggarrange(plotlist = list(
  get_dmg_decay(df = mtdmg_w1_filt_genus |> inner_join(sel_tax), tax = tax, orient = "fwd") +
    ggtitle(paste0(tax, " ", n_readsa, " FWD")),
  get_dmg_decay(df = mtdmg_w1_filt_genus |> inner_join(sel_tax), tax = tax, orient = "rev") +
    ggtitle(paste0(tax, " ", n_readsa, " REV"))
))

# smags_metadata |>
#   select(Genome_Id, contains("Best_taxo")) |>
#   filter(!grepl("METDB", Genome_Id)) |>
#   replace(is.na(.), "NA") |>
#   #filter(Best_taxonomy_PHYLUM ==  "NA") |>
#   mutate(Best_taxonomy_KINGDON = "d__Eukaryota",
#          Best_taxonomy_PHYLUM = ifelse(Best_taxonomy_PHYLUM == "NA", "p__", paste0("p__", Best_taxonomy_PHYLUM)),
#          Best_taxonomy_CLASS = ifelse(Best_taxonomy_CLASS == "NA", "c__", paste0("c__", Best_taxonomy_CLASS)),
#          Best_taxonomy_ORDER = ifelse(Best_taxonomy_ORDER == "NA", "o__", paste0("o__", Best_taxonomy_ORDER)),
#          Best_taxonomy_FAMILY = ifelse(Best_taxonomy_FAMILY == "NA", "f__", paste0("f__", Best_taxonomy_FAMILY)),
#          Best_taxonomy_GENRE = ifelse(Best_taxonomy_GENRE == "NA", "g__", paste0("g__", Best_taxonomy_GENRE)),
#          Best_taxonomy_SPECIES =paste0("s__", Genome_Id),
#          ) |>
#   unite(contains("Best_taxo"), sep = ";", col = "tax_ranks") |>
#   write_tsv(file = "data/map-smags/SMAGS-anvio/smags-tax.tsv", col_names = FALSE)
