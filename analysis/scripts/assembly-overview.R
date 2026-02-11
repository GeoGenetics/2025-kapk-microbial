library(tidyverse)
library(janitor)
library(ggthemr)
library(ggpubr)
library(lvplot)
source("libs/lib.R")
library(ggbeeswarm)
# Get sample data
kapk_cdata <- readxl::read_xlsx("data/cdata/KapK_samples-20210702.xlsx") |>
  clean_names() |>
  select(member_unit, collapsed_read_files_path_29bp_and_q_29_only_collapsed_reads, collapsed_read_files_path_29bp_and_q_29_only_collapsed_reads_md5sum, figure_names, site)

kapk_cdata <- kapk_cdata |>
  mutate(label = basename(collapsed_read_files_path_29bp_and_q_29_only_collapsed_reads),
         label = gsub("\\..*","",label)) |>
  distinct() |>
  select(label, collapsed_read_files_path_29bp_and_q_29_only_collapsed_reads_md5sum, member_unit, figure_names, site) |>
  rename(label_orig = label, label = collapsed_read_files_path_29bp_and_q_29_only_collapsed_reads_md5sum) |>
  mutate(site = as.character(site),
         site = ifelse(site == "NA", "nosite", site),
         figure_names = paste0(figure_names, "-" ,substr(label, 1, 5))) |>
  rowwise() |>
  mutate(wf_label=substr(digest::digest(label, algo="md5",  serialize = FALSE), 1, 10))

site_order <- tibble(site = c("119", "50", "69", "75", "74", "nosite"),
                     site_rnk = 1:6)


# Read assembly refined stats ---------------------------------------------

assm_ref_stats <- read_tsv("data/stats/all.stats-assm-refined-summary.tsv.gz") |>
  mutate(class = "Refined")
assm_comb_stats <- read_tsv("data/stats/all.stats-assm-combined-summary.tsv.gz") |>
  mutate(class = "Combined")

# combined stats
# n contigs: 1,970,953
assm_comb_stats |> nrow()
summary(assm_comb_stats$length) |> enframe() |> knitr::kable()

# |name    |  value|
# |:-------|------:|
# |Min.    |   1000|
# |1st Qu. |   1239|
# |Median  |   1641|
# |Mean    |   2389|
# |3rd Qu. |   2512|
# |Max.    | 794474|

# refined stats
# n contigs: 1,079,947
assm_ref_stats |> nrow()
summary(assm_ref_stats$length) |> enframe() |> knitr::kable()

# |name    |  value|
# |:-------|------:|
# |Min.    |   1000|
# |1st Qu. |   1218|
# |Median  |   1608|
# |Mean    |   2251|
# |3rd Qu. |   2515|
# |Max.    | 181713|

assm_ref_stats |>
  bind_rows(assm_comb_stats) |>
  inner_join(kapk_cdata) |>
  select(length, member_unit, class) |>
  group_by(member_unit, class) |>
  skimr::skim()



ggthemr(palette = "fresh")
assm_ref_stats |>
  bind_rows(assm_comb_stats) |>
  inner_join(kapk_cdata) |>
  ggplot(aes(class, length, fill = class)) +
  #geom_violin(alpha = 0.5) +
  geom_lv(size = 0.5, width.method = "height", color = "#404040", width = 0.5, alpha = 0.5) +
  scale_y_log10(labels = scales::comma) +
  facet_wrap(~member_unit) +
  ylab("Contig length (bp)") +
  xlab("") +
  theme_light() +
  theme(legend.position = "top")

# Taxonomic annotation ----------------------------------------------------

contig_tax_nr <- read_tsv(file = "data/taxonomy/assm.refined.tax-summary.NR.tsv.gz") |>
  tidyr::separate(col = tax_string, sep = ";", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), remove = FALSE)

contig_tax_gtdb <- read_tsv(file = "data/taxonomy/assm.refined.tax-summary.GTDB.tsv.gz") |>
  tidyr::separate(col = tax_string, sep = ";", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), remove = FALSE)


# Total contigs: 1,063,700
# Contigs annotated: 1,025,572
contig_tax_nr |>
  select(contig) |>
  bind_rows(contig_tax_gtdb |> select(contig)) |>
  distinct() |>
  nrow()


contig_tax_ov_sel <- contig_tax_nr |>
  select(contig, support, Kingdom) |>
  rename(support_nr = support,
         Kingdom_nr = Kingdom) |>
  full_join(contig_tax_gtdb |> select(contig, support, Kingdom)) |>
  mutate(class = case_when(is.na(Kingdom_nr) & is.na(Kingdom) ~ "unclassified",
                           Kingdom_nr == "unknown" & Kingdom == "unknown" ~ "unclassified",
                           Kingdom_nr == "unknown" & is.na(Kingdom) ~ "unclassified",
                           is.na(Kingdom_nr) & Kingdom == "unknown" ~ "unclassified",
                           Kingdom_nr != "unknown" & Kingdom == "unknown" ~ "NR",
                           Kingdom_nr == "unknown" & Kingdom != "unknown" ~ "GTDB",
                           is.na(Kingdom_nr) & !is.na(Kingdom) ~ "GTDB",
                           !is.na(Kingdom_nr) & is.na(Kingdom) ~ "NR",
                           Kingdom_nr == "Bacteria" & Kingdom == "Bacteria" ~ "GTDB",
                           Kingdom_nr == "Archaea" & Kingdom == "Archaea" ~ "GTDB",
                           Kingdom_nr == "Bacteria" & Kingdom == "Archaea" ~ "GTDB",
                           Kingdom_nr == "Archaea" & Kingdom == "Bacteria" ~ "GTDB",
                           support >= support_nr ~ "GTDB",
                           TRUE ~ "NR"))

contig_tax_nr_sel <- contig_tax_nr |>
  filter(contig %in% (contig_tax_ov_sel |> filter(class == "NR") |> .$contig)) |>
  mutate(Phylum = ifelse(Phylum == "unknown", paste0("uc_", Kingdom), Phylum))
contig_tax_gtdb_sel <- contig_tax_gtdb |>
  filter(contig %in% (contig_tax_ov_sel |> filter(class == "GTDB") |> .$contig))
contig_tax_uncl <- contig_tax_gtdb |>
  bind_rows(contig_tax_nr) |>
  filter(contig %in% (contig_tax_ov_sel |> filter(class == "unclassified") |> .$contig)) |>
  select(contig) |>
  distinct()
nrow(contig_tax_nr_sel) + nrow(contig_tax_gtdb_sel) + nrow(contig_tax_uncl)

dom_h <- c("unclassified", "root", "Eukaryota",  "Archaea", "Bacteria", "Viruses")
contig_tax_nr_sel |>
  mutate(db = "NCBI-nr") |>
  bind_rows(contig_tax_gtdb_sel |> mutate(db = "GTDB")) |>
  group_by(Kingdom, db) |>
  count() |>
  ungroup() |>
  mutate(Kingdom = fct_relevel(Kingdom, rev(dom_h)),
         db = fct_rev(db)) |>
  ggplot(aes(Kingdom, n)) +
  geom_col(color = "black", fill = "#333333") +
  facet_wrap(~db, scales = "free") +
  scale_y_continuous(labels = scales::comma) +
  ggpubr::rotate() +
  xlab("") +
  ylab("Number of contigs")


contig_tax_euk <- contig_tax_nr_sel |>
  filter(Kingdom == "Eukaryota")
contig_tax_noneuk <- contig_tax_nr_sel |>
  filter(Kingdom != "Eukaryota") |>
  bind_rows(contig_tax_gtdb_sel)

top_phylum <- contig_tax_euk |>
  group_by(label, Phylum) |>
  count() |>
  ungroup() |>
  group_by(label) |>
  mutate(prop = n/sum(n)) |>
  ungroup() |>
  select(Phylum, prop) |>
  group_by(Phylum) |>
  summarise(mean_prop = mean(prop)) |>
  ungroup() |>
  arrange(desc(mean_prop))

library(wesanderson)
pal <- wes_palette("Zissou1", 100, type = "continuous")

contig_tax_euk |>
  group_by(label, Phylum) |>
  count() |>
  #filter(Phylum %in% top_phylum$Phylum) |>
  #mutate(name = ifelse(name %in% top_phylum$name, name, "Other")) |>
  #droplevels() |>
  ungroup() |>
  inner_join(kapk_cdata) |>
  inner_join(site_order) |>
  mutate(Phylum = fct_relevel(Phylum, rev(top_phylum$Phylum)),
         figure_names = fct_reorder(figure_names, site_rnk)) |>
  select(Phylum, figure_names, member_unit, n) |>
  group_by(member_unit) |>
  complete(Phylum, figure_names, fill = list(n = NA)) |>
  ungroup() |>
  inner_join(kapk_cdata |> select(figure_names, member_unit)) |>
  ggplot(aes(figure_names, Phylum,  fill = n)) +
  geom_tile()+
  geom_tile(colour="gray30",size=0.3, show.legend = FALSE, na.rm = TRUE)+
  labs(x="",y="",title="")+
  scale_y_discrete(expand=c(0,0), position = "left")+
  scale_x_discrete(expand=c(0,0))+
  scale_fill_gradientn(colours = pal, labels = scales::comma_format(accuracy = 1), trans = "log10") + #, name="Number of restaurant", guide = guide_legend( keyheight = unit(3, units= "mm"), keywidth=unit(12, units = "mm"), label.position = "bottom", title.position = 'top', nrow=1) ) +
  theme_grey(base_size=10)+
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
  facet_grid(~member_unit, scales = "free_x", space = "free_x")





top_phylum <- contig_tax_gtdb_sel |>
  group_by(label, Phylum) |>
  count() |>
  ungroup() |>
  group_by(label) |>
  mutate(prop = n/sum(n)) |>
  ungroup() |>
  select(Phylum, prop) |>
  group_by(Phylum) |>
  summarise(mean_prop = mean(prop)) |>
  ungroup() |>
  arrange(desc(mean_prop)) |>
  head(50)



library(wesanderson)
pal <- wes_palette("Zissou1", 100, type = "continuous")

contig_tax_gtdb_sel |>
  filter(Phylum %in% top_phylum$Phylum) |>
  group_by(label, Phylum) |>
  count() |>
  ungroup() |>
  inner_join(kapk_cdata) |>
  inner_join(site_order) |>
  mutate(Phylum = fct_relevel(Phylum, rev(top_phylum$Phylum)),
         figure_names = fct_reorder(figure_names, site_rnk)) |>
  select(Phylum, figure_names, member_unit, n) |>
  group_by(member_unit) |>
  complete(Phylum, figure_names, fill = list(n = NA)) |>
  ungroup() |>
  inner_join(kapk_cdata |> select(figure_names, member_unit)) |>
  ggplot(aes(figure_names, Phylum,  fill = n)) +
  geom_tile()+
  geom_tile(colour="gray30",size=0.3, show.legend = FALSE, na.rm = TRUE)+
  labs(x="",y="",title="")+
  scale_y_discrete(expand=c(0,0), position = "left")+
  scale_x_discrete(expand=c(0,0))+
  scale_fill_gradientn(colours = pal, labels = scales::comma_format(accuracy = 1), trans = "log10") + #, name="Number of restaurant", guide = guide_legend( keyheight = unit(3, units= "mm"), keywidth=unit(12, units = "mm"), label.position = "bottom", title.position = 'top', nrow=1) ) +
  theme_grey(base_size=10)+
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
  facet_grid(~member_unit, scales = "free_x", space = "free_x")



# Damage ------------------------------------------------------------------

contig_dmg <- read_tsv("data/pydamage/pydamage-filt-summary.tsv.gz") |>
  rename(contig = reference)

contig_dmg |>
  select(label, contig, predicted_accuracy, qvalue, pvalue) |>
  inner_join(assm_ref_stats) |>
  left_join(bind_rows(contig_tax_gtdb_sel, contig_tax_nr_sel)) |>
  inner_join(kapk_cdata) |>
  mutate(Kingdom = ifelse(is.na(Kingdom), "Unclassified", Kingdom),
         Kingdom = fct_relevel(Kingdom, (dom_h))) |> arrange(desc(length)) |>
  ggplot(aes(Kingdom, length)) +
  geom_lv(size = 0.5, width.method = "height", color = "#404040", width = 0.5, alpha = 0.5) +
  scale_y_log10(labels = scales::comma_format(accuracy = 1))+
  theme_light() +
  ylab("Contig length (bp)") +
  xlab("")


contig_dmg |>
  select(label, contig, predicted_accuracy, qvalue, pvalue) |>
  inner_join(assm_ref_stats) |>
  left_join(bind_rows(contig_tax_gtdb_sel, contig_tax_nr_sel)) |>
  group_by(Kingdom) |> count()


assm_cov <- read_tsv("data/pydamage/assm-contig-cov-summary.tsv.gz") |>
  clean_names() |>
  rename(contig = number_id)

organelle_seqs <- read_tsv("data/organelles/assm-organelles-seqs.tsv.gz", col_names = c("contig", "seq", "X1")) |>
  select(-X1)

contig_organelle_old <- contig_organelle
contig_organelle <- read_tsv("data/organelles/organelles-summary-arctic.cov-60.id-90.tsv.gz") |>
  rename(contig = query) |>
  select(contig, theader, pident, qlen, tlen, qcov, bits, organelle, label) |>
  group_by(contig) |>
  arrange(desc(bits)) |>
  do(head(., n = 1)) |>
  ungroup()


organelle_accs <- read_tsv("data/organelles/organelles-acc.tsv")

contig_organelle <- contig_organelle |>
  left_join(contig_dmg |> select(contig) |> mutate(damaged = "damaged")) |>
  #inner_join(organelle_accs |> rename(theader = accession)) |>
  mutate(damaged = ifelse(is.na(damaged), "non_damaged", damaged))

organelle_seqs |> filter(!(contig %in% contig_organelle$contig))


contig_organelle |>
  inner_join(kapk_cdata) |>
  ggplot(aes(organelle, pident/100)) +
  geom_quasirandom(shape = 21, dodge.width=0.5, color = "black", fill = "#F1A747", alpha = 0.5, size = 1.5) +
  facet_wrap(~member_unit) +
  ggpubr::rotate() +
  theme_light() +
  scale_y_continuous(labels = scales::percent) +
  xlab("") +
  ylab("Identity")

contig_organelle |>
  inner_join(kapk_cdata) |>
  ggplot(aes(organelle, qlen)) +
  geom_quasirandom(shape = 21, dodge.width=0.5, color = "black", fill = "#F1A747", alpha = 0.5, size = 1.5) +
  facet_wrap(~member_unit) +
  ggpubr::rotate() +
  theme_light() +
  scale_y_continuous(labels = scales::comma) +
  xlab("") +
  ylab("Contig length (bp)")

colors <- c("#252525", "#C02942", "#4477A8", "#F1A747")
names(colors) <- c("Other", "Salix", "Betula", "Populus")

data <- contig_organelle |>
  filter(pident >= 95) |>
  mutate(class = case_when(grepl("salix", organism_name, ignore.case = TRUE) ~ "Salix",
                           grepl("betula", organism_name, ignore.case = TRUE) ~ "Betula",
                           grepl("populus", organism_name, ignore.case = TRUE) ~ "Populus",
                           TRUE ~ "Other")) |>
  inner_join(kapk_cdata) |>
  inner_join(assm_cov)
ggplot() +
  geom_point(data = data |> filter(class == "Other"), aes(pident/100, qlen, fill = class), shape = 21, color = "black", alpha = 0.3, size = 2) +
  geom_point(data = data |> filter(class == "Salix"), aes(pident/100, qlen, fill = class), shape = 21, color = "black", alpha = 0.8, size = 2) +
  geom_point(data = data |> filter(class == "Populus"), aes(pident/100, qlen, fill = class), shape = 21, color = "black", alpha = 0.8, size = 2) +
  geom_point(data = data |> filter(class == "Betula"), aes(pident/100, qlen, fill = class), shape = 21, color = "black", alpha = 0.8, size = 2) +
  facet_grid(organelle~member_unit) +
  theme_light() +
  scale_y_log10(labels = scales::comma) +
  scale_x_log10(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = colors) +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  xlab("Identity") +
  ylab("Contig length (bp)")


ggplot() +
  geom_point(data = data |> filter(class == "Other"), aes(avg_fold, qlen, fill = class), shape = 21, color = "black", alpha = 0.3, size = 2) +
  geom_point(data = data |> filter(class == "Salix"), aes(avg_fold, qlen, fill = class), shape = 21, color = "black", alpha = 0.8, size = 2) +
  geom_point(data = data |> filter(class == "Populus"), aes(avg_fold, qlen, fill = class), shape = 21, color = "black", alpha = 0.8, size = 2) +
  geom_point(data = data |> filter(class == "Betula"), aes(avg_fold, qlen, fill = class), shape = 21, color = "black", alpha = 0.8, size = 2) +
  facet_grid(organelle~member_unit) +
  theme_light() +
  scale_y_log10(labels = scales::comma) +
  scale_x_log10() +
  scale_fill_manual(values = colors) +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  xlab("Mean coverage") +
  ylab("Contig length (bp)")


contig_organelle_salix <- contig_organelle |>
  filter(grepl("salix", organism_name, ignore.case = TRUE)) |>
  arrange(desc(qlen))

contig_organelle_betula <- contig_organelle |>
  filter(grepl("betula", organism_name, ignore.case = TRUE)) |>
  arrange(desc(qlen))

contig_organelle_populus <- contig_organelle |>
  filter(grepl("populus", organism_name, ignore.case = TRUE)) |>
  arrange(desc(qlen))

contig_organelle |>
  inner_join(kapk_cdata) |>
  group_by(organelle, damaged) |>
  dplyr::count() |>
  ungroup() |>
  pivot_wider(names_from = "damaged", values_from = n, values_fill = 0) |>
  mutate(total = damaged + non_damaged) |>
  knitr::kable()

contig_organelle |> inner_join(assm_cov) |>
  select(-length) |>
  setNames(c("contig", "best_hit", "identity", "contig_length",
             "best_hit_length", "query_coverage", "bits",
             "organelle", "label", "is_damaged", "organism_name",
             "organism_groups", "contig_gc_content", "coverage_mean", "coverage_breadth")) |>
  mutate(coverage_breadth = coverage_breadth/100,
         identity = identity/100) |>
  inner_join(kapk_cdata) |>
  filter(grepl("betula", organism_name, ignore.case = TRUE)) |>
  write_tsv("results/kapk-assm-organelles-betula.tsv")


customFun  = function(DF) {
  label <- DF$label |> unique()
  DF |>
    select(contig, seq) |>
    write_tsv(paste0("~/Downloads/KapK_organelles/",label,"-organelles.tsv"), col_names = FALSE)
  return(DF)
}

organelle_seqs |>
  inner_join(contig_organelle |> select(contig, label)) |>
  group_by(label) |>
  do(customFun(.))

contig_organelle |>
  inner_join(kapk_cdata) |>
  inner_join(site_order) |>
  mutate(figure_names = fct_reorder(figure_names, site_rnk)) |>
  group_by(organelle,figure_names, member_unit) |>
  count() |>
  select(organelle, figure_names, member_unit, n) |>
  group_by(member_unit) |>
  complete(organelle, figure_names, fill = list(n = NA)) |>
  ungroup() |>
  inner_join(kapk_cdata |> select(figure_names, member_unit)) |>
  ggplot(aes(figure_names, organelle,  fill = n)) +
  geom_tile()+
  geom_tile(colour="gray30",size=0.3, show.legend = FALSE, na.rm = TRUE)+
  labs(x="",y="",title="")+
  scale_y_discrete(expand=c(0,0), position = "left")+
  scale_x_discrete(expand=c(0,0))+
  scale_fill_gradientn(colours = pal, labels = scales::comma_format(accuracy = 1), trans = "log10") + #, name="Number of restaurant", guide = guide_legend( keyheight = unit(3, units= "mm"), keywidth=unit(12, units = "mm"), label.position = "bottom", title.position = 'top', nrow=1) ) +
  theme_grey(base_size=10)+
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
  facet_grid(~member_unit, scales = "free_x", space = "free_x")



contig_organelle |>
  filter(pident >= 95) |>
  inner_join(contig_dmg) |>
  mutate(class = case_when(grepl("salix", organism_name, ignore.case = TRUE) ~ "Salix",
                           grepl("betula", organism_name, ignore.case = TRUE) ~ "Betula",
                           grepl("populus", organism_name, ignore.case = TRUE) ~ "Populus",
                           TRUE ~ "Other")) |>
  inner_join(kapk_cdata) |>
  inner_join(assm_cov) |>
  select(contig, member_unit, class, starts_with("Ctot")) |>
  pivot_longer(names_to = "CtoT", values_to = "prop", cols = -!starts_with("CtoT")) |>
  mutate(CtoT = as.numeric(gsub("CtoT_", "", CtoT)) + 1) |>
  ggplot(aes(CtoT, prop, color = class, fill = class)) +
  stat_summary(fun.data = mean_sd, geom = "ribbon", alpha = .3, size = 0) +
  stat_summary(fun = mean, geom = "line", size = 0.8) +
  facet_grid(class~member_unit) +
  xlab("5'-end position") +
  ylab("C-to-T frequency") +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors)




# For tom -----------------------------------------------------------------

"55b07932586644efda03b438b26d0313"
"df264c7a33522a9fac14690c78ef62b4"
"f6c6573d3db8c8e59a61b0b37504093c"


# "55b07932586644efda03b438b26d0313"
# ncontigs: 14,204
# n_contigs damaged:
# ncontigs > 2000: 2,640
X <- "55b07932586644efda03b438b26d0313"
get_stats_tom <- function(X){
  n_contigs <- assm_ref_stats |>
    filter(label == X)
  n_contigs_dmg <- assm_ref_stats |>
    filter(label == X) |>
    inner_join(contig_dmg)
  n_contigs_filt <- assm_ref_stats |>
    filter(label == X, length >= 2000)
  n_contigs_filt_dmg <- assm_ref_stats |>
    filter(label == X, length >= 2000) |>
    inner_join(contig_dmg)
  data <- assm_ref_stats |>
    filter(label == X) |>
    mutate(is_damaged = ifelse(contig %in% contig_dmg$contig, TRUE, FALSE)) |>
    left_join(bind_rows(contig_tax_gtdb_sel, contig_tax_nr_sel)) |>
    separate(contig, into="short_label", remove=FALSE, extra="drop", sep = "_")
  res <- tibble("class" = c("n_contigs", "n_contigs_dmg", "n_contigs_filt", "n_contigs_filt_dmg"),
                "n" = c(nrow(n_contigs), nrow(n_contigs_dmg), nrow(n_contigs_filt), nrow(n_contigs_filt_dmg))) |>
    mutate(label = X, short_label = data$short_label |> unique())

  return(list(res = res,data = data))
}

df1 <- get_stats_tom("55b07932586644efda03b438b26d0313")
df2 <- get_stats_tom("df264c7a33522a9fac14690c78ef62b4")
df3 <- get_stats_tom("f6c6573d3db8c8e59a61b0b37504093c")
df4 <- get_stats_tom("d12d52cc3a95f4d0714eba529c00ea57")
map_df(list(df1, df2, df3, df4), "res") |>
  knitr::kable()

map_df(list(df1, df2, df3, df4), "data") |> write_tsv("~/Downloads/kapk4Tom/kapk_cdata.tsv")



# Viral contigs -----------------------------------------------------------

vs2_contigs <- read_tsv("data/viral/vs2/final-viral-score.tsv") |>
  clean_names() |>
  separate(seqname, into = c("contig", "cat"), sep = "\\|\\|", extra = "drop", remove = FALSE)



vs2_contigs |>
  left_join(contig_dmg) |>
  mutate(damage = ifelse(is.na(predicted_accuracy), "non-damaged", "damaged")) |>
  ggplot(aes(length, fill = damage)) +
  geom_density(color = "black", alpha = 0.5) +
  scale_x_log10(labels=scales::comma_format()) +
  xlab("Contig length (bp)")

vs2_contigs |>
  inner_join(contig_dmg) |>
  arrange(desc(length)) |>
  select(contig,starts_with("Ctot")) |>
  pivot_longer(names_to = "CtoT", values_to = "prop", cols = -!starts_with("CtoT")) |>
  mutate(CtoT = as.numeric(gsub("CtoT_", "", CtoT)) + 1) |>
  ggplot(aes(CtoT, prop)) +
  stat_summary(fun.data = mean_sd, geom = "ribbon", alpha = .3, size = 0) +
  stat_summary(fun = mean, geom = "line", size = 0.8) +
  xlab("5'-end position") +
  ylab("C-to-T frequency")
