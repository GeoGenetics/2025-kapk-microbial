estimate_mode <- function(x) {
  if (length(x) == 2) {
    mean(x)
  } else if (length(x) == 1) {
    x
  } else if (length(x) == 0) {
    NA
  } else {
    d <- density(x)
    d$x[which.max(d$y)]
  }
}


# Fixed coords ------------------------------------------------------------

get_ratio <- function(x = x, y = y, display = 4 / 3) {
  ratio_display <- display
  ratio_values <- ((max(x) - min(x))) / ((max(y) - min(y)))
  ratio_values / ratio_display
}

get_ratio_rev <- function(x = x, y = y, display = 4 / 3) {
  ratio_display <- display
  ratio_values <- ((min(x) - max(x))) / ((max(y) - min(y)))
  ratio_values / ratio_display
}


get_ratio_log10 <- function(x = x, y = y, display = 4 / 3) {
  ratio_display <- display
  ratio_values <- (log10(max(x) - min(x))) / (log10(max(y) - min(y)))
  ratio_values / ratio_display
}

myround <- function(x) {
  trunc(x + 0.5)
}



# Get core data (PHYLOSEQ) ------------------------------------------------


get_core_data <- function(core, bdata, cdata, threshold = 0, p_threshold = 0) {
  cat(paste(core, "\n"))
  labels <- cdata %>%
    filter(site_abbreviation == core)
  data <- bdata %>% filter(label %in% labels$label)
  data_df <- data %>%
    filter(fraction_total_reads > p_threshold) %>%
    select(label, name, new_est_reads) %>%
    pivot_wider(names_from = label, values_from = new_est_reads, values_fill = 0) %>%
    as.data.frame() %>%
    column_to_rownames("name")
  cdata_df <- cdata %>%
    filter(site_abbreviation == core) %>%
    mutate(rname = label) %>%
    as.data.frame() %>%
    column_to_rownames("rname")
  tax_data <- data %>%
    select(name, tax_string) %>%
    distinct() %>%
    tidyr::separate(col = tax_string, sep = ";", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) %>%
    as.data.frame() %>%
    column_to_rownames("name")
  ps <- phyloseq(otu_table(data_df, taxa_are_rows = TRUE), sample_data(cdata_df), tax_table(as.matrix(tax_data)))

  tax_filt <- qpsmelt(ps) %>%
    dplyr::select(label, OTU, Abundance) %>%
    dplyr::group_by(label) %>%
    mutate(prop = Abundance / sum(Abundance)) %>%
    ungroup() %>%
    group_by(OTU) %>%
    dplyr::summarise(mean_prop = sum(prop) / nsamples(ps)) %>%
    ungroup() %>%
    filter(mean_prop >= threshold) %>%
    .$OTU

  ps_filt <- prune_taxa(tax_filt, ps)
  ps_filt_phylum <- tax_glom(ps_filt, taxrank = "Phylum")
  ps_prop <- transform_sample_counts(ps_filt, function(X) {
    X / sum(X)
  })
  ps_prop_phylum <- transform_sample_counts(ps_filt_phylum, function(X) {
    X / sum(X)
  })
  ps_norm <- MDPpreprocess(phyloObject = ps_filt, processDir = "/tmp/")
  list(ps = ps_filt, ps_norm = ps_norm, ps_prop = ps_prop, ps_prop_phylum = ps_prop_phylum, ps_filt_phylum = ps_filt_phylum)
}



get_core_data_st <- function(core, bdata, cdata, rdata, threshold = 0, p_threshold = 0) {
  cat(paste(core, "\n"))
  cores_st <- rdata[[core]] %>% .$site_abbreviation
  labels <- cdata %>%
    filter(site_abbreviation %in% cores_st)
  data <- bdata %>% filter(label %in% labels$label)
  data_df <- data %>%
    filter(fraction_total_reads > p_threshold) %>%
    select(label, name, new_est_reads) %>%
    pivot_wider(names_from = label, values_from = new_est_reads, values_fill = 0) %>%
    as.data.frame() %>%
    column_to_rownames("name")
  cdata_df <- cdata %>%
    filter(site_abbreviation %in% cores_st) %>%
    mutate(rname = label) %>%
    as.data.frame() %>%
    column_to_rownames("rname")
  tax_data <- data %>%
    select(name, tax_string) %>%
    distinct() %>%
    tidyr::separate(col = tax_string, sep = ";", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) %>%
    as.data.frame() %>%
    column_to_rownames("name")

  ps <- phyloseq(otu_table(data_df, taxa_are_rows = TRUE), sample_data(cdata_df), tax_table(as.matrix(tax_data)))

  tax_filt <- qpsmelt(ps) %>%
    dplyr::select(label, OTU, Abundance) %>%
    dplyr::group_by(label) %>%
    mutate(prop = Abundance / sum(Abundance)) %>%
    ungroup() %>%
    group_by(OTU) %>%
    dplyr::summarise(mean_prop = sum(prop) / nsamples(ps)) %>%
    ungroup() %>%
    filter(mean_prop >= threshold) %>%
    .$OTU

  ps_filt <- prune_taxa(tax_filt, ps)
  ps_filt_phylum <- tax_glom(ps_filt, taxrank = "Phylum")
  ps_prop <- transform_sample_counts(ps_filt, function(X) {
    X / sum(X)
  })
  ps_prop_phylum <- transform_sample_counts(ps_filt_phylum, function(X) {
    X / sum(X)
  })
  ps_norm <- MDPpreprocess(phyloObject = ps_filt, processDir = "/tmp/")
  list(ps = ps_filt, ps_norm = ps_norm, ps_prop = ps_prop, ps_prop_phylum = ps_prop_phylum, ps_filt_phylum = ps_filt_phylum)
}



# get prevalent taxa ------------------------------------------------------

get_prev <- function(ps, prev = 50 / 100) {
  rare_taxa <- microbiome::rare_members(ps, prevalence = prev)
  ps_norare <- microbiome::remove_taxa(rare_taxa, ps)
  ps_norare_phylum <- tax_glom(ps_norare, taxrank = "Phylum")
  ps_prop <- transform_sample_counts(ps_norare, function(X) {
    X / sum(X)
  })
  ps_prop_phylum <- transform_sample_counts(ps_norare_phylum, function(X) {
    X / sum(X)
  })
  ps_norm <- MDPpreprocess(phyloObject = ps_norare, processDir = "/tmp/")
  list(ps = ps_norare, ps_norm = ps_norm, ps_prop = ps_prop, ps_prop_phylum = ps_prop_phylum, ps_phylum = ps_norare_phylum)
}


smooth_sp <- function(X, data, groups) {
  gr <- groups[X, ]
  data <- data %>% filter(group == gr$group, core == gr$core)
  dat2 <- as.data.frame(spline(data$x, data$y))
  dat2$group <- gr$group
  dat2$core <- gr$core
  names(dat2) <- c("x", "y", "group", "core")
  dat2
}
library(data.table)

get_stream <- function(df) {
  data <- df %>%
    select(pos, Abundance, Phylum, core) %>%
    setNames(c("x", "y", "group", "core"))

  d_groups <- data %>%
    dplyr::select(group, core) %>%
    distinct()

  data <- map_dfr(1:nrow(d_groups), smooth_sp, data = data, groups = d_groups)

  setDT(data)
  datasdorder <- data[,
    list(sdy = sd(y)),
    by = list(group)
  ][
    ,
    ranksdy := rank(sdy, ties.method = "random")
  ][
    ranksdy %% 2 == 0,
    ranksdy := -ranksdy
  ]

  data <- merge(data,
    datasdorder,
    "group",
    all = T
  )
  setkey(data, x, ranksdy)
  data[, ymax := cumsum(y) - (sum(y) / 2), by = list(x, core)]
  data[, ymin := ymax - y, by = list(x, core)]
  data %>% as_tibble()
}

get_stream_path <- function(df) {
  data <- df %>%
    select(pos, Abundance, Var, core) %>%
    setNames(c("x", "y", "group", "core"))

  d_groups <- data %>%
    dplyr::select(group, core) %>%
    distinct()

  data <- map_dfr(1:nrow(d_groups), smooth_sp, data = data, groups = d_groups)

  setDT(data)
  datasdorder <- data[,
    list(sdy = sd(y)),
    by = list(group)
  ][
    ,
    ranksdy := rank(sdy, ties.method = "random")
  ][
    ranksdy %% 2 == 0,
    ranksdy := -ranksdy
  ]

  data <- merge(data,
    datasdorder,
    "group",
    all = T
  )
  setkey(data, x, ranksdy)
  data[, ymax := cumsum(y) - (sum(y) / 2), by = list(x, core)]
  data[, ymin := ymax - y, by = list(x, core)]
  data %>% as_tibble()
}

qpsmelt_zeroes <- function(X) {
  if (taxa_are_rows(X)) {
    count_table <- as(otu_table(X), "matrix") %>%
      as_tibble(rownames = "OTU") %>%
      gather(label, Abundance, -OTU)
  } else {
    count_table <- as(otu_table(X), "matrix") %>%
      as_tibble(rownames = "label") %>%
      gather(OTU, Abundance, -label)
  }
  sample_table <- as(sample_data(X), "matrix") %>%
    as_tibble()
  taxa_table <- as(tax_table(X), "matrix") %>%
    as_tibble(rownames = "OTU")

  count_table %>%
    left_join(sample_table) %>%
    left_join(taxa_table)
}

qpsmelt <- function(X) {
  if (taxa_are_rows(X)) {
    count_table <- as(otu_table(X), "matrix") %>%
      as_tibble(rownames = "OTU") %>%
      pivot_longer(cols = -OTU, names_to = "Sample", values_to = "Abundance")
  } else {
    count_table <- as(otu_table(X), "matrix") %>%
      as_tibble(rownames = "Sample") %>%
      pivot_longer(cols = -Sample, names_to = "OTU", values_to = "Abundance")
  }
  sample_table <- as(sample_data(X), "data.frame") %>%
    as_tibble(rownames = "Sample")

  taxa_table <- as(tax_table(X), "matrix") %>%
    as_tibble(rownames = "OTU")

  count_table %>%
    left_join(sample_table) %>%
    left_join(taxa_table) %>%
    filter(Abundance > 0)
}


# Calculate horizon plot --------------------------------------------------


estimate_horizons <- function(X, data, combs) {
  comb <- combs[X, ]
  origin <- 0
  dat1 <- data %>%
    filter(core == comb$core, OTU == comb$OTU) %>%
    mutate(mc = Abundance - median(Abundance)) %>%
    arrange(pos)

  dat2 <- as.data.frame(spline(dat1$pos, dat1$mc))
  dat2$OTU <- dat1$OTU
  names(dat2) <- c("pos", "Abundance", "OTU")

  df.az <- dat2 %>% select(pos, OTU, Abundance)
  colnames(df.az) <- c("date", "grouping", "y")
  df.az <- as.data.table(df.az)

  max.y <- max(abs(df.az$y - origin))
  nbands <- 3
  # horizonscale <- max(abs(df.az$y-origin))/nbands
  horizonscale <- mad(abs(df.az$y - origin))

  h1 <- horizonscale
  h2 <- horizonscale * 2
  h3 <- horizonscale * 3
  h1n <- -horizonscale
  h2n <- -horizonscale * 2
  h3n <- -horizonscale * 3

  df5 <- df.az[, ":="(ypos1 = ifelse(y > 0, min(y, h1), 0),
    ypos2 = ifelse(y > h1, min(y, h2) - h1, 0),
    ypos3 = ifelse(y > h2, min(y, h3) - h2, 0),
    yneg1 = -ifelse(y < 0, max(y, h1n), 0),
    yneg2 = -ifelse(y < h1n, max(y, h2n) - h1n, 0),
    yneg3 = -ifelse(y < h2n, max(y, h3n) - h2n, 0)),
  by = c("date", "grouping")
  ]
  df6 <- df5 %>%
    select(-y) %>%
    gather(type, value, 3:8)
  colnames(df6) <- c("date", "grouping", "band", "value")
  df6$vmin <- 0
  df6$v2 <- ifelse(abs(df6$value) < abs(df6$vmin), df6$vmin, df6$value)
  df6$core <- comb$core
  return(df6 %>% as_tibble())
}



estimate_horizons_st <- function(X, data, combs) {
  comb <- combs[X, ]
  origin <- 0
  dat1 <- data %>%
    filter(region_st == comb$region_st, OTU == comb$OTU) %>%
    mutate(mc = Abundance - median(Abundance)) %>%
    arrange(pos)

  dat2 <- as.data.frame(spline(dat1$pos, dat1$mc))
  dat2$OTU <- dat1$OTU
  names(dat2) <- c("pos", "Abundance", "OTU")

  df.az <- dat2 %>% select(pos, OTU, Abundance)
  colnames(df.az) <- c("date", "grouping", "y")
  df.az <- as.data.table(df.az)

  max.y <- max(abs(df.az$y - origin))
  nbands <- 3
  # horizonscale <- max(abs(df.az$y-origin))/nbands
  horizonscale <- mad(abs(df.az$y - origin))

  h1 <- horizonscale
  h2 <- horizonscale * 2
  h3 <- horizonscale * 3
  h1n <- -horizonscale
  h2n <- -horizonscale * 2
  h3n <- -horizonscale * 3

  df5 <- df.az[, ":="(ypos1 = ifelse(y > 0, min(y, h1), 0),
    ypos2 = ifelse(y > h1, min(y, h2) - h1, 0),
    ypos3 = ifelse(y > h2, min(y, h3) - h2, 0),
    yneg1 = -ifelse(y < 0, max(y, h1n), 0),
    yneg2 = -ifelse(y < h1n, max(y, h2n) - h1n, 0),
    yneg3 = -ifelse(y < h2n, max(y, h3n) - h2n, 0)),
  by = c("date", "grouping")
  ]
  df6 <- df5 %>%
    select(-y) %>%
    gather(type, value, 3:8)
  colnames(df6) <- c("date", "grouping", "band", "value")
  df6$vmin <- 0
  df6$v2 <- ifelse(abs(df6$value) < abs(df6$vmin), df6$vmin, df6$value)
  df6$region_st <- comb$region_st
  return(df6 %>% as_tibble())
}


motus_expected <- tibble(
  number_of_reads = c(5e6, 8e6, 15e6, 25e6, 35e6, 50e6, 100e6),
  median_motus_count = c(600, 900, 1900, 3300, 5500, 8800, 13000)
)

find_closest_nmotus <- function(x) {
  val <- which.min(abs(motus_expected$number_of_reads - x))
  motus_expected[val, ]$median_motus_count
}



# Create phyloseq objects -------------------------------------------------

get_philr <- function(physeq, glom = NULL, min_counts = 100000, ncounts = 3, nsites = 0.2, vcoeff = 3,
                      part.weights = "enorm.x.gm.counts", ilr.weights = "uniform", colors = NULL) {
  suppressMessages({
    if (!is.null(glom)) {
      cat(paste("Tax glomming to:", glom, "\n"))
      physeq <- speedyseq::tax_glom(physeq, taxrank = glom)
    }

    if (is.null(colors)) {
      colors <- c(
        "#88CCAA", "#44AAAA", "#771155", "#DDDD77", "#DD7788", "#777711",
        "#44AA77", "#774411", "#77CCCC", "#4477AA", "#77AADD", "#AAAA44",
        "#CC99BB", "#AA7744", "#117744", "#DDAA77", "#771122", "#AA4455",
        "#AA4488", "#114477", "#117777"
      )
    }

    cat(paste("Filtering samples with at least", min_counts, "counts...\n"))
    physeq_filt <- prune_samples(sample_sums(physeq) >= min_counts, physeq)

    cat(paste("Filtering taxa not seen at least", ncounts, "in", nsites, "of the samples...\n"))
    physeq_filt <- prune_taxa(taxa_sums(physeq_filt) > 0, physeq_filt)
    physeq_filt <- filter_taxa(physeq_filt, function(x) sum(x > ncounts) > (nsites * length(x)), TRUE)

    total <- median(sample_sums(physeq_filt))
    standf <- function(x, t = total) round(t * (x / sum(x)))
    physeq_filt_m <- transform_sample_counts(physeq_filt, standf)

    # Filter ASVs with a CV <= 3.0
    # osd2014_dada2_phyloseq_alpha_norare = filter_taxa(osd2014_dada2_phyloseq_alpha_filt, function(x) sd(x)/mean(x) > 3.0, TRUE)

    if (is.null(vcoeff)) {
      vcoeff <- speedyseq::psmelt(physeq_filt) %>%
        as_tibble() %>%
        group_by(OTU) %>%
        summarise(cv = sd(Abundance) / mean(Abundance)) %>%
        .$cv %>%
        median()
      cat(paste("WARNING: Using median CV of all OGUs\n"))
      cat(paste("Filtering taxa with a lower CV than", vcoeff, "...\n"))
    } else {
      cat(paste("Filtering taxa with a lower CV than", vcoeff, "...\n"))
    }

    physeq_filt_v <- filter_taxa(physeq_filt_m, function(x) sd(x) / mean(x) > vcoeff, TRUE)
    physeq_filt <- prune_taxa(taxa_names(physeq_filt) %in% taxa_names(physeq_filt_v), physeq_filt)

    cat(paste0("Kept ", ntaxa(physeq_filt), " taxa (Original:", ntaxa(physeq), ") summing a total of ", sum(taxa_sums(physeq_filt)), " counts (Original:", sum(taxa_sums(physeq)), ")\n"))
    physeq_filt <- transform_sample_counts(physeq_filt, function(x) x + 1)
    ogus.outgroup <- row.names(subset(tax_table(physeq_filt), tax_table(physeq_filt)[, "Domain"] == "d__Archaea"))
    phy_tree(physeq_filt) <- root(phy_tree(physeq_filt), ogus.outgroup, resolve.root = TRUE)
    tax_table(physeq_filt) <- tax_table(physeq_filt)[, -ncol(tax_table(physeq_filt))]

    cat(paste("Is tree rooted?", is.rooted(phy_tree(physeq_filt))), "\n")
    cat(paste("Is tree binary?", is.binary.tree(phy_tree(physeq_filt))), "\n")

    cat("Preparing data for PhilR", "\n")
    phy_tree(physeq_filt) <- makeNodeLabel(phy_tree(physeq_filt), method = "number", prefix = "n")
    name.balance(phy_tree(physeq_filt), tax_table(physeq_filt), "n1")
    otu.table <- t(otu_table(physeq_filt))
    tree <- phy_tree(physeq_filt)
    metadata <- sample_data(physeq_filt)
    tax <- tax_table(physeq_filt)

    cat("Performing PhilR transformation", "\n")
    ps_philr <- philr(otu.table, tree,
      part.weights = part.weights,
      ilr.weights = ilr.weights
    )

    cat("Calculating PCoA", "\n")
    philr_dist <- dist(ps_philr, method = "euclidean")
    philr_pcoa <- ordinate(physeq_filt, "PCoA", distance = philr_dist)

    pcoa_sc <- philr_pcoa$vectors %>%
      as_tibble(rownames = "label") %>%
      inner_join(sample_sums(physeq_filt) %>% enframe(name = "label", value = "counts")) %>%
      inner_join(speedyseq::psmelt(physeq_filt))
    #
    ratio_values <- get_ratio(x = pcoa_sc$Axis.1, y = rev(pcoa_sc$Axis.2))

    p <- pcoa_sc %>%
      select(Axis.1, Axis.2, estimated_age_ka_bp, geographic_location, counts) %>%
      distinct() %>%
      ggplot(aes(x = Axis.1, y = Axis.2, fill = estimated_age_ka_bp, shape = geographic_location)) +
      # ggConvexHull::geom_convexhull(data = pcoa_sc %>%
      #                                select(Axis.1, Axis.2, site, member_unit), aes(Axis.1, Axis.2), fill = "grey80", alpha = 0.2, size = 0.5) +
      geom_point(size = 3, color = "black") +
      theme_bw() +
      xlab(paste0("PC1: ", round(philr_pcoa$values$Relative_eig[1], 3) * 100, "%")) +
      ylab(paste0("PC2: ", round(philr_pcoa$values$Relative_eig[2], 3) * 100, "%")) +
      theme(
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "top"
      ) +
      coord_fixed(ratio = ratio_values) +
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      # scale_shape_manual(values = c(21, 22, 24, 23)) +
      guides(
        fill = guide_legend(override.aes = list(shape = 21)),
        shape = guide_legend(override.aes = list(fill = "black"))
      )
    cat("All done", "\n")
  })
  list(ps = physeq_filt, philr = ps_philr, philr_dist = philr_dist, philr_pcoa = philr_pcoa, plot = p)
}



get_st <- function(physeq, glom = NULL, min_counts = 10000, ncounts = 3, nsites = 0.2, vcoeff = 3, mean_prop_thresh = 1e-5) {
  suppressMessages({
    if (!is.null(glom)) {
      cat(paste("Tax glomming to:", glom, "\n"))
      physeq <- speedyseq::tax_glom(physeq, taxrank = glom)
    }

    cat(paste("Filtering samples with at least", min_counts, "counts...\n"))
    physeq_filt <- prune_samples(sample_sums(physeq) >= min_counts, physeq)

    cat(paste("Filtering taxa not seen at least", ncounts, "in", nsites, "of the samples...\n"))
    physeq_filt <- prune_taxa(taxa_sums(physeq_filt) > 0, physeq_filt)
    physeq_filt <- filter_taxa(physeq_filt, function(x) sum(x > ncounts) > (nsites * length(x)), TRUE)
    physeq_filt <- prune_samples(sample_sums(physeq_filt) >= min_counts, physeq_filt)

    # Scale to the median of the sample depth
    cat("Scaling to median sample depth", "\n")
    total <- median(sample_sums(physeq_filt))
    standf <- function(x, t = total) round(t * (x / sum(x)))
    physeq_filt_m <- transform_sample_counts(physeq_filt, standf)
    # Filter ASVs with a CV <= 3.0
    # osd2014_dada2_phyloseq_alpha_norare = filter_taxa(osd2014_dada2_phyloseq_alpha_filt, function(x) sd(x)/mean(x) > 3.0, TRUE)

    if (is.null(vcoeff)) {
      vcoeff <- speedyseq::psmelt(physeq_filt_m) %>%
        as_tibble() %>%
        group_by(OTU) %>%
        summarise(cv = sd(Abundance) / mean(Abundance)) %>%
        .$cv %>%
        median()
      cat(paste("WARNING: Using median CV of all OGUs\n"))
      cat(paste("Filtering taxa with a lower CV than", vcoeff, "...\n"))
    } else {
      cat(paste("Filtering taxa with a lower CV than", vcoeff, "...\n"))
    }

    physeq_filt_v <- filter_taxa(physeq_filt_m, function(x) sd(x) / mean(x) > vcoeff, TRUE)
    physeq_filt <- prune_taxa(taxa_names(physeq_filt) %in% taxa_names(physeq_filt_v), physeq_filt)

    # remove those with very low abundance across samples
    cat(paste("Filtering taxa with a mean proportion <=", mean_prop_thresh, "\n"))
    tax_filt <- speedyseq::psmelt(physeq_filt) %>%
      dplyr::select(Sample, OTU, Abundance) %>%
      group_by(Sample) %>%
      mutate(prop = Abundance / sum(Abundance)) %>%
      ungroup() %>%
      group_by(OTU) %>%
      summarise(mean_prop = sum(prop) / nsamples(physeq_filt)) %>%
      ungroup() %>%
      filter(mean_prop >= mean_prop_thresh) %>%
      .$OTU
    physeq_filt <- prune_taxa(tax_filt, physeq_filt)
    physeq_filt <- prune_taxa(taxa_sums(physeq_filt) > 0, physeq_filt)
    physeq_filt <- prune_samples(sample_sums(physeq_filt) >= 0, physeq_filt)
    cat(paste0("Kept ", ntaxa(physeq_filt), " taxa (Original:", ntaxa(physeq), ") summing a total of ", sum(taxa_sums(physeq_filt)), " counts (Original:", sum(taxa_sums(physeq)), ")\n"))
    return(physeq_filt)
  })
}


# Get the LV values
# From https://mgimond.github.io/ES218/Week08b.html
lsum <- function(x, l = 6, all = TRUE) {
  # Limit max letter summaries to 9
  if (l > 10) {
    print("Limit level summary to 9")
    return()
  }
  # letter summary labels
  let <- c("M", "H", "E", "D", "C", "B", "A", "Z", "Y", "X")
  # Remove missing values
  x <- na.omit(x)
  # Sort values
  x <- sort(x)
  # Find depths from each end
  n <- length(x)
  Lrnk <- vector()
  Mrnk <- vector()
  Rrnk <- vector()
  Lrnk[1] <- n
  Mrnk[1] <- n
  Rrnk[1] <- n
  i <- 1
  while ((i <= l) & (Lrnk[i] > 1)) {
    i <- i + 1
    Lrnk[i] <- floor(Lrnk[i - 1] + 1) / 2
    Mrnk[i] <- floor(Lrnk[i])
    Rrnk[i] <- floor(Lrnk[i] + 0.5)
  }
  # Get final set of letters
  val <- factor(let[1:length(Lrnk[-1])], levels = let[1:length(Lrnk[-1])])
  # Find the summary values
  LO <- (x[Mrnk[-1]] + x[Rrnk[-1]]) / 2
  HI <- (x[n - Mrnk[-1] + 1] + x[n - Rrnk[-1] + 1]) / 2
  MD <- (LO + HI) / 2
  SP <- HI - LO
  # Generate output
  if (all == TRUE) {
    out <- data.frame(
      letter = val, depth = Lrnk[-1], lower = LO,
      mid = MD, upper = HI, spread = SP
    )
  } else {
    out <- data.frame(letter = val, mid = MD)
  }
  return(out)
}

read_gtdb_tree <- function(X) {
  gtdb_tree <- read.newick(X)
  gtdb_tree$tip.label <- gsub("__", "---", gtdb_tree$tip.label)
  gtdb_tree$tip.label <- gsub("_", " ", gtdb_tree$tip.label)
  gtdb_tree$tip.label <- gsub("---", "__", gtdb_tree$tip.label)
  return(gtdb_tree)
}

# Function to calculate the geometric mean
gm_mean <- function(x, na_rm = TRUE) {
    exp(sum(log(x[x > 0]), na.rm = na_rm) / length(x))
}
