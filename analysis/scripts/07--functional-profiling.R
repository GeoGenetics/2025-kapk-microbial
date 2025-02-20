# Required libraries
library(tidyverse)
library(janitor)
library(showtext)

# Initialize settings
source("./libs/lib.R")
showtext_auto()

#' Get functional information from KEGG and dbCAN data
#' @param kegg_data Path to KEGG data
#' @param dbcan_data Path to dbCAN data
#' @param completion_method Method for module completion
#' @param completion_threshold Threshold for completion
#' @return list containing functional information
get_functional_info <- function(kegg_data, dbcan_data,
                                completion_method = "stepwise_module_completeness",
                                completion_threshold = 0.9) {
    # Load KEGG data
    kegg_all_data <- read_tsv(kegg_data)

    # Filter based on completion method
    kegg_all_data_filt <- if (completion_method == "stepwise_module_completeness") {
        kegg_all_data %>% filter(stepwise_module_completeness >= completion_threshold)
    } else if (completion_method == "pathwise_module_completeness") {
        kegg_all_data %>% filter(pathwise_module_completeness >= completion_threshold)
    } else {
        kegg_all_data %>% filter(module_completeness >= completion_threshold)
    }

    # Get Woodcroft modules
    woodcroft2018_modules <- kegg_all_data_filt %>%
        filter(module_class %in% c("Woodcroft2018"))

    # Get Borrel modules
    borrel2022_modules <- kegg_all_data_filt %>%
        filter(module_class %in% c("Borrell2022")) %>%
        mutate(module_class = ifelse(module_class == "Borrell2022", "Borrel2023", module_class))

    # Get signature and pathway modules
    signature_modules <- kegg_all_data_filt %>%
        filter(module_class %in% c("Signature modules"))
    pathway_modules <- kegg_all_data_filt %>%
        filter(module_class %in% c("Pathway modules"))

    # Load dbCAN data
    dbcan_cdata <- read_delim("./data/function/dbcan/CAZyDB.08062022.fam-activities.txt",
        delim = "\t  ",
        col_names = c("family", "description"),
        comment = "#"
    ) %>%
        mutate(class = case_when(
            grepl("^GH", family) ~ "Glycoside Hydrolases",
            grepl("^GT", family) ~ "GlycosylTransferases",
            grepl("^PL", family) ~ "Polysaccharide Lyases",
            grepl("^CE", family) ~ "Carbohydrate Esterases",
            grepl("^AA", family) ~ "Auxiliary Activities",
            grepl("^CBM", family) ~ "Carbohydrate Binding Modules"
        ))

    dbcan_all_data <- read_tsv(dbcan_data)

    # Define degradation groups
    cellulose_degradation <- c("GH5", "GH9", "3.2.1.4", "GH51", "GH6", "GH7", "GH48", "3.2.1.91")
    hemicellulose_xylan_degradation <- c(
        "GH5", "GH8", "GH10", "GH11", "GH43", "3.2.1.8", "GH3", "GH30",
        "GH39", "GH52", "GH54", "GH116", "GH120", "3.2.1.37", "GH67", "GH115", "3.2.1.139", "CE1", "CE2",
        "CE3", "CE4", "CE5", "CE6", "CE7", "CE12", "3.1.1.72"
    )
    xylose_degradation <- c(
        "3.2.1.37", "GH3", "GH30", "GH39", "GH10", "GH43",
        "xylose degradation (isomerase pathway)",
        "xylose degradation (oxidoreductase pathway)",
        "xylose degradation (xylonate hydratase pathway)",
        "xylose degradation (weimburg/dahms)"
    )

    # Process degradation data
    dbcan_cellulose <- dbcan_all_data %>% filter(group %in% cellulose_degradation)
    dbcan_hemicellulose <- dbcan_all_data %>% filter(group %in% hemicellulose_xylan_degradation)
    dbcan_xylose <- dbcan_all_data %>% filter(group %in% xylose_degradation)

    # Get Woodcroft xylose degradation data
    woodcroft2018_xylose <- woodcroft2018_modules %>%
        filter(module_name %in% c(
            "xylose degradation (isomerase pathway)",
            "xylose degradation (oxidoreductase pathway)",
            "xylose degradation (xylonate hydratase pathway)",
            "xylose degradation (weimburg/dahms)"
        ))

    # Combine polysaccharide degradation data
    polysac_deg <- dbcan_cellulose %>%
        mutate(process = "Cellulose degradation") %>%
        bind_rows(dbcan_hemicellulose %>% mutate(process = "Hemicellulose xylan degradation")) %>%
        bind_rows(dbcan_xylose %>% mutate(process = "Xylose degradation"))

    # Process fermentation modules
    fermentation_modules <- c(
        "ethanol fermentation",
        "acetatogenesis",
        "lactate fermentation",
        "propionate fermentation"
    )

    fermentation <- woodcroft2018_modules %>%
        filter(module_name %in% fermentation_modules)

    central_carbon_meta <- fermentation %>%
        mutate(process = module_name) %>%
        select(annotation = module, abundance = avg_coverage, label, process)

    # Process methane metabolism modules
    pathway_modules_methanogen <- pathway_modules %>%
        filter(module_name %in% c(
            "Methanogenesis, CO2 => methane",
            "Methanogenesis, acetate => methane",
            "Methanogenesis, methanol => methane",
            "Methanogenesis, methylamine/dimethylamine/trimethylamine => methane"
        ))

    woodcroft2018_modules_methanogen <- woodcroft2018_modules %>%
        filter(module_name %in% c("acetoclastic methanogenesis", "hydrogenotrophic methanogenesis"))

    pathway_modules_methanotroph <- pathway_modules %>%
        filter(module_name == "Methane oxidation, methanotroph, methane => formaldehyde")

    # Combine methane metabolism data
    methane_metabolism <- pathway_modules_methanogen %>%
        mutate(process = module_name) %>%
        select(annotation = module, label, process, abundance = avg_coverage, module_class) %>%
        bind_rows(woodcroft2018_modules_methanogen %>%
            mutate(process = module_name) %>%
            select(annotation = module, label, process, abundance = avg_coverage, module_class)) %>%
        bind_rows(pathway_modules_methanotroph %>%
            mutate(process = module_name) %>%
            select(annotation = module, label, process, abundance = avg_coverage, module_class)) %>%
        bind_rows(borrel2022_modules %>%
            mutate(
                process = module_name,
                module_class = as.character(module_class)
            ) %>%
            select(annotation = module, label, process, abundance = avg_coverage, module_class))

    # Return results
    list(
        polysac_deg = polysac_deg,
        central_carbon_meta = central_carbon_meta,
        methane_metabolism = methane_metabolism,
        signature_modules = signature_modules
    )
}

#' Process functional data
#' @param func_all All functional data
#' @param func_dmg Damaged functional data
#' @param kapk_cdata Sample metadata
#' @param label_nobloom No bloom labels
#' @return list containing processed data
process_functional_data <- function(func_all, func_dmg, kapk_cdata, label_nobloom) {
    # Process signature paths
    signature_paths <- func_all$signature_modules %>%
        mutate(type = "All") %>%
        bind_rows(func_dmg$signature_modules %>% mutate(type = "Damaged")) %>%
        filter(
            module_category == "Module set",
            module_subcategory == "Metabolic capacity"
        ) %>%
        select(label, abundance = avg_coverage, type, process = module_name) %>%
        filter(process == "Acetogen") %>%
        mutate(
            prank = 2,
            class = "Central carbon metabolism"
        )

    # Process methane metabolism
    methane_metabolism <- func_all$methane_metabolism %>%
        mutate(type = "All") %>%
        bind_rows(func_dmg$methane_metabolism %>% mutate(type = "Damaged")) %>%
        filter(module_class == "Borrel2023" |
            process == "Methane oxidation, methanotroph, methane => formaldehyde" |
            process == "hydrogenotrophic methanogenesis" |
            process == "Methanogenesis, methanol => methane" |
            process == "Methanogenesis, methylamine/dimethylamine/trimethylamine => methane") %>%
        select(label, abundance, type, process) %>%
        mutate(
            process = sub("(.)", "\\U\\1", process, perl = TRUE),
            process = case_when(
                process == "Methane oxidation, methanotroph, methane => formaldehyde" ~ "Methanotrophy",
                process == "Methanogenesis, methanol => methane" ~ "Methylotrophic methanogenesis (methanol)",
                process == "Methanogenesis, methylamine/dimethylamine/trimethylamine => methane" ~ "Methylotrophic methanogenesis (methylamine)",
                TRUE ~ process
            )
        ) %>%
        mutate(
            prank = 3,
            class = "Methane metabolism"
        )

    # Process central carbon metabolism
    central_carbon_meta <- func_all$central_carbon_meta %>%
        mutate(type = "All") %>%
        bind_rows(func_dmg$central_carbon_meta %>% mutate(type = "Damaged")) %>%
        select(label, abundance, type, process) %>%
        mutate(process = sub("(.)", "\\U\\1", process, perl = TRUE)) %>%
        mutate(
            prank = 2,
            class = "Central carbon metabolism"
        ) %>%
        filter(process != "Acetatogenesis")

    # Process polysaccharide degradation
    polysac_deg <- func_all$polysac_deg %>%
        mutate(type = "All") %>%
        bind_rows(func_dmg$polysac_deg %>% mutate(type = "Damaged")) %>%
        select(label, abundance = mean, type, process) %>%
        filter(!grepl("Xylose", process)) %>%
        mutate(
            prank = 1,
            class = "Polysaccharide degradation"
        )

    # Define data order
    data_all_order <- c(
        "Hemicellulose xylan degradation",
        "Cellulose degradation",
        "Methanotrophy",
        "Methylotrophic methanogenesis (methanol)",
        "Methylotrophic methanogenesis (methylamine)",
        "Acetoclastic methanogenesis",
        "Hydrogenotrophic methanogenesis",
        "Ethanol fermentation",
        "Propionate fermentation",
        "Lactate fermentation",
        "Acetogen"
    )

    # Process all data
    data_all <- polysac_deg %>%
        bind_rows(central_carbon_meta) %>%
        bind_rows(methane_metabolism) %>%
        bind_rows(signature_paths) %>%
        filter(type == "All") %>%
        droplevels() %>%
        complete(label, nesting(process, class), type, fill = list(abundance = NA)) %>%
        inner_join(kapk_cdata) %>%
        select(label = figure_names, member_unit, process, class, abundance, type, site_rnk) %>%
        ungroup() %>%
        mutate(
            label = fct_reorder(label, site_rnk),
            member_unit = fct_relevel(member_unit, c("B1", "B2", "B3")),
            class = fct_relevel(class, c(
                "Polysaccharide degradation",
                "Central carbon metabolism",
                "Methane metabolism"
            )),
            process = fct_relevel(process, data_all_order)
        )

    # Process damaged data
    data_dmg <- polysac_deg %>%
        bind_rows(central_carbon_meta) %>%
        bind_rows(methane_metabolism) %>%
        bind_rows(signature_paths) %>%
        filter(type == "Damaged") %>%
        droplevels() %>%
        inner_join(kapk_cdata) %>%
        mutate(
            figure_names = fct_reorder(figure_names, site_rnk),
            member_unit = fct_relevel(member_unit, c("B1", "B2", "B3")),
            class = fct_relevel(class, c(
                "Polysaccharide degradation",
                "Central carbon metabolism",
                "Methane metabolism"
            )),
            process = fct_relevel(process, data_all_order)
        ) %>%
        select(label = figure_names, member_unit, process, class, abundance, type)

    list(
        all = data_all,
        dmg = data_dmg,
        order = data_all_order
    )
}

#' Create functional heatmap
#' @param data_all All data
#' @param data_dmg Damaged data
#' @return ggplot object
plot_functional_heatmap <- function(data_all, data_dmg) {
    ggplot(data_all, aes(x = label, y = process, fill = abundance)) +
        geom_tile(color = "black", size = 0.2) +
        geom_tile(
            data = data_dmg, aes(x = label, y = process, color = type),
            color = "#323232", size = 1
        ) +
        facet_grid(class ~ member_unit, scales = "free", space = "free") +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "top",
            legend.title = element_blank(),
            text = element_text(size = 10),
            strip.background = element_blank()
        ) +
        scale_fill_gradientn(
            colors = c("#D7E7F0", "#F7E2B5", "#F1AF82", "#E67B70", "#C94A6B", "#9B2471"),
            n.breaks = 6,
            guide = "colorbar",
            na.value = "#ffffff",
            trans = "log2"
        ) +
        xlab("") +
        ylab("")
}

#' Plot KEGG completion density
#' @param kegg_data KEGG data
#' @param kapk_cdata Sample metadata
#' @return list of ggplot objects
plot_kegg_completion <- function(kegg_data, kapk_cdata) {
    kegg_data <- kegg_data %>%
        inner_join(kapk_cdata %>% select(label, figure_names)) %>%
        select(-label) %>%
        rename(label = figure_names) %>%
        mutate(module_class = ifelse(module_class == "Borrell2022", "Borrel2023", module_class))

    p1 <- kegg_data %>%
        filter(grepl("acetoclastic", module_name) | module == "M00567") %>%
        ggplot(aes(pathwise_module_completeness, fill = module_class)) +
        geom_density(alpha = 0.5)

    p2 <- kegg_data %>%
        filter(grepl("acetoclastic", module_name) | module == "M00567") %>%
        ggplot(aes(stepwise_module_completeness, fill = module_class)) +
        geom_density(alpha = 0.5)

    list(
        pathwise = p1,
        stepwise = p2
    )
}

#' Main execution function
#' @return list containing results and plots
main <- function() {
    # Load metadata
    kapk_cdata <- read_tsv("./data/cdata/KapK-cdata-manuscript-20221211.tsv",
        show_col_types = FALSE
    )
    label_nobloom <- read_tsv("./data/cdata/KapK-cdata-manuscript-20221211-labels-nobloom.tsv",
        show_col_types = FALSE
    )

    # Filter metadata
    kapk_cdata <- kapk_cdata %>%
        filter(figure_names %in% label_nobloom$label)

    # Get functional data
    func_all <- get_functional_info(
        kegg_data = "./data/function/kegg/all/kegg-modules-summary.tsv.gz",
        dbcan_data = "./data/function/dbcan/all/dbcan.group-abundances-agg.tsv.gz",
        completion_method = "stepwise_module_completeness",
        completion_threshold = 1.0
    )

    func_dmg <- get_functional_info(
        kegg_data = "./data/function/kegg/damaged/kegg-modules-summary.tsv.gz",
        dbcan_data = "./data/function/dbcan/damaged/dbcan.damaged.group-abundances-agg.tsv.gz",
        completion_method = "stepwise_module_completeness",
        completion_threshold = 0.8
    )

    # Process functional data
    processed_data <- process_functional_data(
        func_all,
        func_dmg,
        kapk_cdata,
        label_nobloom
    )

    # Create heatmap
    heatmap_plot <- plot_functional_heatmap(
        processed_data$all,
        processed_data$dmg
    )

    # Get KEGG completion plots
    kegg_data <- read_tsv("./data/function/kegg/all/kegg-modules-summary.tsv.gz",
        show_col_types = FALSE
    )
    kegg_plots <- plot_kegg_completion(kegg_data, kapk_cdata)

    # Create combined KEGG plot
    kegg_combined <- ggpubr::ggarrange(
        kegg_plots$pathwise,
        kegg_plots$stepwise,
        nrow = 1,
        common.legend = TRUE,
        legend = "bottom"
    )

    # Save data

    # Return results
    list(
        data = list(
            all = processed_data$all,
            dmg = processed_data$dmg,
            kegg = kegg_data
        ),
        plots = list(
            heatmap = heatmap_plot,
            kegg = kegg_plots,
            kegg_combined = kegg_combined
        )
    )
}

# Run main function and store results
results <- main()
