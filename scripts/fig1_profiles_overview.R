# Figure 1 | Protein-profile overview

# 0. Setup ----
source("scripts/_project_setup.R")
use_packages(c("data.table", "tidyverse", "readxl", "ComplexHeatmap", "circlize", "ggridges", "ggpubr", "ggplotify", "openxlsx", "showtext"))
source("utils/figure_style.R")
source("utils/feature_mapping.R")
plasmix_theme <- if (is.function(theme_plasmix)) theme_plasmix() else theme_plasmix
showtext_auto(); showtext_opts(dpi = 600)
label_style <- list(size = 12, face = "bold")

# 1. Inputs ----
paths <- c(
    metadata = "data/study_metadata.xlsx", feature_metadata = "data/feature_metadata.tsv.gz",
    profiles = "data/protein_profiles_long.tsv.gz", detection = "results/detection_status.tsv.gz",
    physchem = "data/physchem_matrix.tsv.gz"
)
input_check <- tibble(Input = names(paths), Path = unname(paths), Exists = file.exists(paths))
print(input_check)
stopifnot(all(input_check$Exists))

meta_batch <- read_xlsx(paths["metadata"], sheet = "batch")
feat_meta <- fread(paths["feature_metadata"])

detection_status <- fread(paths["detection"])
physchem_matrix <- fread(paths["physchem"])
long_df <- fread(paths["profiles"])
long_df <- filter_batch_analysis_features(long_df, feat_meta, strict_platforms = character())
detection_status <- detection_status %>% semi_join(distinct(long_df, Platform, Batch, UniqueID), by = c("Platform", "Batch", "UniqueID"))
stopifnot(nrow(anti_join(distinct(long_df, Platform, Batch, UniqueID), detection_status, by = c("Platform", "Batch", "UniqueID"))) == 0)

required_hpa_columns <- c(
    "Entry", "BloodConc_log10_pgml", "Abundance_Source", "HPA_Protein_Class", "HPA_Subcellular"
)
stopifnot(all(required_hpa_columns %in% names(physchem_matrix)))

hpa_annotation <- physchem_matrix %>%
    select(all_of(required_hpa_columns)) %>%
    rename(UniProtID = Entry) %>%
    distinct(UniProtID, .keep_all = TRUE)

hpa_conc <- hpa_annotation %>%
    filter(Abundance_Source != "Unknown", is.finite(BloodConc_log10_pgml)) %>%
    select(UniProtID, BloodConc_log10_pgml, Abundance_Source)

long_df_filter <- long_df %>%
    filter(
        Sample %in% c("M", "Y", "P", "X", "F", "N"),
        DataTier == "Baseline" | Batch %in% c("OLK_P1_B1", "OLK_P1_B2")
    )

# 2. Panel b: HPA-referenced abundance distributions ----
protocol_proteins <- long_df_filter %>%
    distinct(Platform, Batch, UniProtID) %>%
    filter(!is.na(UniProtID), UniProtID != "") %>%
    left_join(meta_batch %>% select(Batch, Protocol), by = "Batch") %>%
    mutate(
        Protocol = case_when(
            Protocol %in% c("AAgAtlas IgA", "AAgAtlas IgG") ~ "AAgAtlas IgA/IgG",
            TRUE ~ Protocol
        )
    ) %>%
    select(-Batch) %>%
    distinct()

abundance_plot_data <- protocol_proteins %>%
    left_join(hpa_conc, by = "UniProtID") %>%
    filter(is.finite(BloodConc_log10_pgml)) %>%
    distinct(Platform, Protocol, UniProtID, BloodConc_log10_pgml)

protocol_medians <- abundance_plot_data %>%
    group_by(Protocol, Platform) %>%
    summarize(median_val = median(BloodConc_log10_pgml, na.rm = TRUE), .groups = "drop") %>%
    arrange(median_val)

abundance_plot_data <- abundance_plot_data %>%
    mutate(Protocol = factor(Protocol, levels = protocol_medians$Protocol))

p_abundance <- ggplot(abundance_plot_data, aes(x = BloodConc_log10_pgml, y = Protocol, fill = Protocol)) +
    geom_density_ridges(scale = 1.5, quantile_lines = TRUE, quantiles = 2, alpha = 0.9, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = protocol_color) +
    scale_x_continuous(limits = c(0, 9), breaks = c(0, 3, 6, 9), expand = c(0, 0)) +
    scale_y_discrete(expand = expansion(mult = c(0, 0))) +
    labs(x = "Estimated conc. (log10 pg/mL)", y = NULL) +
    theme_plasmix() +
    theme(
        legend.position = "none", panel.grid.major = element_blank(),
        axis.text.x = element_text(hjust = c(0.1, 0.5, 0.5, 0.9)), axis.title.x = element_text(hjust = 1)
    )

p_abundance

# 3. Shared annotation functions for the UpSet panels ----
subcellular_colors <- c(
    "Secreted" = "#E41A1C",
    "Membrane" = "#377EB8",
    "Intracellular" = "#4DAF4A",
    "Secreted & Membrane" = "#984EA3",
    "Unknown" = "#E0E0E0"
)
comb_colors <- c(
    "#1B9E77FF", "#D95F02FF", "#7570B3FF", "#E7298AFF",
    "#E6AB02FF", "#66A61EFF", "#A6761DFF", "#666666FF"
)
location_levels <- c("Secreted", "Secreted & Membrane", "Membrane", "Intracellular", "Unknown")

hpa_subcell <- hpa_annotation %>%
    transmute(
        UniProtID, Subcellular = replace_na(HPA_Subcellular, "Unknown"),
        Subcellular_Factor = factor(replace_na(HPA_Subcellular, "Unknown"), levels = location_levels)
    )

get_subcell_counts <- function(protein_lists, accession_ids = FALSE) {
    ids <- unique(unlist(protein_lists, use.names = FALSE))
    id_map <- if (accession_ids) tibble(UniqueID = ids, UniProtID = ids) else feat_meta %>% select(UniqueID, UniProtID) %>% distinct()
    composition <- lapply(names(protein_lists), function(name) {
        mapped <- tibble(UniqueID = protein_lists[[name]]) %>%
            left_join(id_map, by = "UniqueID") %>%
            left_join(
                hpa_subcell %>% select(UniProtID, Subcellular, Subcellular_Factor) %>% distinct(),
                by = "UniProtID"
            ) %>%
            mutate(
                Subcellular = replace_na(Subcellular, "Unknown"),
                Subcellular_Factor = factor(
                    replace_na(as.character(Subcellular_Factor), "Unknown"),
                    levels = location_levels
                )
            ) %>%
            group_by(UniqueID) %>%
            arrange(Subcellular_Factor) %>%
            slice(1) %>%
            ungroup()

        counts <- table(mapped$Subcellular)
        result <- setNames(rep(0, length(subcellular_colors)), names(subcellular_colors))
        shared <- intersect(names(result), names(counts))
        result[shared] <- counts[shared]
        result
    })
    matrix_result <- do.call(rbind, composition)
    rownames(matrix_result) <- names(protein_lists)
    matrix_result[, names(subcellular_colors), drop = FALSE]
}

intersection_height <- function(x) x^(1 / 3)
transform_stacked_counts <- function(x) {
    totals <- rowSums(x)
    sweep(x, 1, ifelse(totals > 0, intersection_height(totals) / totals, 0), "*")
}
# format_axis_counts <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

# 4. Panel c: DIA protocol target overlap ----
dia_metadata <- meta_batch %>% filter(Platform == "DIA")
dia_protocol_proteins <- list()

# Figure 1 shows assay coverage, so retain every anchor-tier feature without filtering IsDetected.
for (protocol in unique(dia_metadata$Protocol)) {
    target_batches <- dia_metadata %>%
        filter(Protocol == protocol) %>%
        pull(Batch)

    dia_protocol_proteins[[protocol]] <- detection_status %>%
        filter(Batch %in% target_batches) %>%
        pull(UniqueID) %>%
        unique()
}

dia_matrix <- list_to_matrix(dia_protocol_proteins)
dia_combination <- make_comb_mat(dia_matrix)

dia_combination_names <- comb_name(dia_combination)
dia_combination_lists <- lapply(dia_combination_names, function(name) extract_comb(dia_combination, name))
names(dia_combination_lists) <- dia_combination_names

dia_combination_composition <- get_subcell_counts(dia_combination_lists)
dia_combination_composition_plot <- transform_stacked_counts(dia_combination_composition)
dia_set_composition <- get_subcell_counts(dia_protocol_proteins)[set_name(dia_combination), , drop = FALSE]
dia_set_sizes <- set_size(dia_combination)
dia_axis_breaks <- c(0, 100, 800)
dia_axis_limit <- intersection_height(800) * 1.10

upset_dia <- UpSet(
    dia_combination,
    set_order = order(set_size(dia_combination), decreasing = TRUE),
    comb_order = order(
        comb_degree(dia_combination),
        comb_size(dia_combination),
        decreasing = TRUE
    ),
    comb_col = comb_colors[comb_degree(dia_combination)],
    pt_size = unit(2, "mm"),
    lwd = 1,
    row_names_gp = gpar(fontsize = 8),
    top_annotation = HeatmapAnnotation(
        "Intersection\ncount" = anno_barplot(
            dia_combination_composition_plot, ylim = c(0, dia_axis_limit), extend = 0, bar_width = 0.6,
            gp = gpar(fill = subcellular_colors[colnames(dia_combination_composition_plot)], lwd = 0.2),
            height = unit(1.5, "cm"), add_numbers = FALSE,
            axis_param = list(side = "left", at = intersection_height(dia_axis_breaks),
                              labels = dia_axis_breaks, gp = gpar(fontsize = 7)),
            border = FALSE
        ),
        annotation_name_side = "left",
        annotation_name_rot = 90,
        annotation_name_offset = unit(0.95, "cm"),
        annotation_name_gp = gpar(fontsize = 8)
    ),
    right_annotation = rowAnnotation(
        "Protein composition" = anno_barplot(
            dia_set_composition,
            bar_width = 0.6,
            gp = gpar(
                fill = subcellular_colors[colnames(dia_set_composition)],
                lwd = 0.2
            ),
            border = FALSE,
            width = unit(2, "cm"),
            axis_param = list(labels_rot = 0, gp = gpar(fontsize = 7), at = c(0, 1000, 2000))
        ),
        "Total labels" = anno_barplot(
            dia_set_sizes,
            gp = gpar(fill = "transparent", col = "transparent"),
            bar_width = 0.6,
            add_numbers = TRUE,
            numbers_rot = 0,
            numbers_gp = gpar(fontsize = 6, col = "black"),
            numbers_offset = unit(-1.3, "cm"),
            border = FALSE,
            width = unit(1.4, "cm"),
            axis = FALSE
        ),
        show_annotation_name = FALSE,
        annotation_name_gp = gpar(fontsize = c(9, 0)),
        gap = unit(-0.63, "cm")
    )
)

# 5. Panel d: cross-platform protein coverage by UniProt accession ----
platform_labels <- c(DIA = "MS-DIA", SOM = "SomaScan", OLK = "Olink", NLS = "NULISA", AAG = "AAgAtlas")

platform_proteins <- list()
for (platform_code in names(platform_labels)) {
    platform_proteins[[platform_labels[[platform_code]]]] <- detection_status %>%
        filter(Platform == platform_code) %>%
        inner_join(feat_meta %>% distinct(Platform, UniqueID, UniProtID), by = c("Platform", "UniqueID")) %>%
        filter(!is.na(UniProtID), UniProtID != "") %>%
        pull(UniProtID) %>%
        unique()
}

platform_matrix <- list_to_matrix(platform_proteins)
platform_combination <- make_comb_mat(platform_matrix)

platform_combination_names <- comb_name(platform_combination)
platform_combination_lists <- lapply(
    platform_combination_names, function(name) extract_comb(platform_combination, name)
)
names(platform_combination_lists) <- platform_combination_names

platform_combination_composition <- get_subcell_counts(platform_combination_lists, accession_ids = TRUE)
platform_combination_composition_plot <- transform_stacked_counts(platform_combination_composition)
platform_set_composition <- get_subcell_counts(platform_proteins, accession_ids = TRUE)[set_name(platform_combination), , drop = FALSE]
platform_set_sizes <- set_size(platform_combination)
platform_axis_breaks <- c(0, 100, 1000, 8000)
platform_axis_limit <- intersection_height(8000) * 1.15

upset_platform <- UpSet(
    platform_combination,
    set_order = order(set_size(platform_combination), decreasing = TRUE),
    comb_order = order(
        comb_degree(platform_combination),
        comb_size(platform_combination),
        decreasing = TRUE
    ),
    comb_col = comb_colors[comb_degree(platform_combination)],
    pt_size = unit(2, "mm"),
    lwd = 1,
    row_names_gp = gpar(fontsize = 8),
    top_annotation = HeatmapAnnotation(
        "Intersection\ncount" = anno_barplot(
            platform_combination_composition_plot, ylim = c(0, platform_axis_limit), extend = 0, bar_width = 0.6,
            gp = gpar(fill = subcellular_colors[colnames(platform_combination_composition_plot)], lwd = 0.2),
            height = unit(1.5, "cm"), add_numbers = FALSE,
            axis_param = list(side = "left", at = intersection_height(platform_axis_breaks),
                              labels = platform_axis_breaks, gp = gpar(fontsize = 7)),
            border = FALSE
        ),
        annotation_name_side = "left",
        annotation_name_rot = 90,
        annotation_name_offset = unit(0.95, "cm"),
        annotation_name_gp = gpar(fontsize = 8)
    ),
    right_annotation = rowAnnotation(
        "Protein composition" = anno_barplot(
            platform_set_composition,
            bar_width = 0.6,
            gp = gpar(
                fill = subcellular_colors[colnames(platform_set_composition)],
                lwd = 0.2
            ),
            border = FALSE,
            width = unit(2, "cm"),
            axis_param = list(labels_rot = 0, gp = gpar(fontsize = 7))
        ),
        "Total labels" = anno_barplot(
            platform_set_sizes,
            gp = gpar(fill = "transparent", col = "transparent"),
            bar_width = 0.6,
            add_numbers = TRUE,
            numbers_rot = 0,
            numbers_gp = gpar(fontsize = 6, col = "black"),
            numbers_offset = unit(-1.25, "cm"),
            border = FALSE,
            width = unit(1.5, "cm"),
            axis = FALSE
        ),
        show_annotation_name = FALSE,
        annotation_name_gp = gpar(fontsize = c(9, 0)),
        gap = unit(-0.65, "cm")
    )
)

# 6. Assemble Figure 1 panels b–d ----
subcellular_legend <- Legend(
    at = names(subcellular_colors), title = "Subcellular location", legend_gp = gpar(fill = subcellular_colors),
    title_gp = gpar(fontsize = 7, fontface = "bold"), labels_gp = gpar(fontsize = 7),
    grid_height = unit(2, "mm"), grid_width = unit(2, "mm")
)

dia_grob <- grid.grabExpr({
    drawn <- draw(upset_dia, newpage = FALSE, padding = unit(c(0, 0, 1.5, 0), "mm"))
    ordered_columns <- column_order(drawn)
    decorate_annotation("Intersection\ncount", {
        values <- comb_size(dia_combination)[ordered_columns]
        grid.text(values, x = seq_along(values),
                  y = unit(intersection_height(values), "native") + unit(0.45, "mm"),
                  default.units = "native", just = "left", rot = 90, gp = gpar(fontsize = 6))
    })
})

platform_grob <- grid.grabExpr({
    drawn <- draw(upset_platform, newpage = FALSE, padding = unit(c(0, 0, 3, 0), "mm"))
    ordered_columns <- column_order(drawn)
    decorate_annotation("Intersection\ncount", {
        values <- comb_size(platform_combination)[ordered_columns]
        grid.text(values, x = seq_along(values),
                  y = unit(intersection_height(values), "native") + unit(0.45, "mm"),
                  default.units = "native", just = "left", rot = 90, gp = gpar(fontsize = 6))
    })
})

subcellular_legend_grob <- grid.grabExpr(
    draw(subcellular_legend, x = unit(0, "npc"), y = unit(0.5, "npc"), just = c("left", "center"))
)

dia_plot <- as.ggplot(dia_grob) + theme(plot.margin = margin(0, -5, 0, 5))
platform_plot <- as.ggplot(platform_grob) + theme(plot.margin = margin(0, -5, 0, 16.5))
legend_plot <- as.ggplot(subcellular_legend_grob)

dia_panel <- ggarrange(
    dia_plot, nrow = 1, labels = "c", font.label = label_style,
    label.x = 0, label.y = 1, hjust = -0.2, vjust = 1)

platform_panel <- ggarrange(
    platform_plot, nrow = 1, labels = "d", font.label = label_style,
    label.x = 0, label.y = 1, hjust = -0.2, vjust = 1)

upset_stack <- ggarrange(dia_panel, platform_panel, nrow = 2, align = "v")
right <- ggarrange(upset_stack, legend_plot, nrow = 1, widths = c(1, 0.163), align = "h")

abundance_panel <- ggarrange(
    p_abundance + theme(plot.margin = margin(l = 5, t = 5, b = 5, r = 0)), nrow = 1, labels = "b",
    font.label = label_style, label.x = 0, label.y = 1, hjust = -0.2, vjust = 1
)

figure1_b_to_d <- ggarrange(abundance_panel, right, nrow = 1, widths = c(0.5, 1.8), align = "h")
ggsave("figures/fig1_profiles_overview.pdf", figure1_b_to_d, width = 10, height = 3, bg = "white")
ggsave("figures/fig1_profiles_overview.png", figure1_b_to_d, width = 10, height = 3, dpi = 600, bg = "white")

# 7. Source data ----
# The HPA summary reports mapping percentages across 10 protocol-specific native-intensity deciles; D1 and D10 are the lowest and highest deciles.
hpa_status <- hpa_annotation %>%
    transmute(
        UniProtID, HPA_mapped = (Abundance_Source != "Unknown" & is.finite(BloodConc_log10_pgml)),
        BloodConc_log10_pgml
    ) %>%
    distinct(UniProtID, .keep_all = TRUE)

# Derive one native-intensity rank per target protein and analytical protocol.
# Ranking is first performed within each batch and then summarized across
# replicate batches of the same protocol.
intensity_rank_data <- long_df_filter %>%
    filter(is.finite(Value), !is.na(UniProtID), UniProtID != "") %>%
    left_join(meta_batch %>% select(Batch, Protocol), by = "Batch") %>%
    mutate(
        Protocol = case_when(
            Protocol %in% c("AAgAtlas IgA", "AAgAtlas IgG") ~
                "AAgAtlas IgA/IgG",
            TRUE ~ Protocol
        )
    ) %>%
    group_by(Protocol, Batch, UniProtID) %>%
    summarize(Median_native_intensity = median(Value, na.rm = TRUE), .groups = "drop") %>%
    group_by(Protocol, Batch) %>%
    mutate(Batch_intensity_rank = percent_rank(Median_native_intensity)) %>%
    ungroup() %>%
    group_by(Protocol, UniProtID) %>%
    summarize(Native_intensity_rank = median(Batch_intensity_rank, na.rm = TRUE), .groups = "drop") %>%
    group_by(Protocol) %>%
    mutate(Intensity_decile = ntile(Native_intensity_rank, 10)) %>%
    ungroup() %>%
    left_join(hpa_status %>% select(UniProtID, HPA_mapped), by = "UniProtID") %>%
    mutate(HPA_mapped = replace_na(HPA_mapped, FALSE))

decile_summary_long <- intensity_rank_data %>%
    group_by(Protocol, Intensity_decile) %>%
    summarize(
        Total_targets = n_distinct(UniProtID), HPA_mapped_targets = n_distinct(UniProtID[HPA_mapped]),
        HPA_mapped_percent = 100 * HPA_mapped_targets / Total_targets, .groups = "drop"
    ) %>%
    complete(Protocol, Intensity_decile = 1:10)

decile_summary_wide <- decile_summary_long %>%
    select(Protocol, Intensity_decile, HPA_mapped_percent) %>%
    pivot_wider(
        names_from = Intensity_decile,
        values_from = HPA_mapped_percent,
        names_glue = "HPA mapped D{Intensity_decile} (%)"
    )

protocol_summary <- protocol_proteins %>%
    select(Protocol, UniProtID) %>%
    distinct() %>%
    left_join(hpa_status, by = "UniProtID") %>%
    group_by(Protocol) %>%
    summarize(
        `Total target proteins` = n_distinct(UniProtID),
        `HPA mapped (%)` = 100 * n_distinct(UniProtID[HPA_mapped]) / n_distinct(UniProtID),
        `Min conc. (log10 pg/mL)` = min(BloodConc_log10_pgml[HPA_mapped], na.rm = TRUE),
        `Q1 conc. (log10 pg/mL)` = quantile(BloodConc_log10_pgml[HPA_mapped], 0.25, na.rm = TRUE),
        `Median conc. (log10 pg/mL)` = median(BloodConc_log10_pgml[HPA_mapped], na.rm = TRUE),
        `Q3 conc. (log10 pg/mL)` = quantile(BloodConc_log10_pgml[HPA_mapped], 0.75, na.rm = TRUE),
        `Max conc. (log10 pg/mL)` = max(BloodConc_log10_pgml[HPA_mapped], na.rm = TRUE), .groups = "drop"
    ) %>%
    left_join(decile_summary_wide, by = "Protocol") %>%
    rename(`Analytical protocol` = Protocol)

expected_decile_columns <- paste0("HPA mapped D", 1:10, " (%)")
stopifnot(
    all(expected_decile_columns %in% names(protocol_summary)),
    all(vapply(protocol_summary[expected_decile_columns], is.numeric, logical(1)))
)

protocol_order <- protocol_medians$Protocol
protocol_summary <- protocol_summary %>%
    mutate(`Analytical protocol` = factor(`Analytical protocol`, levels = protocol_order)) %>%
    arrange(`Analytical protocol`) %>%
    mutate(`Analytical protocol` = as.character(`Analytical protocol`))

dia_membership <- enframe(dia_protocol_proteins, name = "Protocol", value = "UniqueID") %>%
    unnest_longer(UniqueID)

platform_membership <- enframe(platform_proteins, name = "Platform", value = "UniProtID") %>%
    unnest_longer(UniProtID)

output_file <- "tables/SourceData_Figure1.xlsx"
write.xlsx(
    list(
        HPA_protocol_summary = protocol_summary,
        HPA_target_intensity_deciles = intensity_rank_data,
        HPA_protocol_values = abundance_plot_data,
        DIA_membership = dia_membership,
        Platform_membership = platform_membership
    ),
    output_file,
    overwrite = TRUE,
    keepNA = TRUE,
    na.string = "NA"
)

# Format HPA_protocol_summary:
# column 2: integer with thousands separators
# columns 3 onward: two decimal places
wb <- loadWorkbook(output_file)

summary_rows <- seq_len(nrow(protocol_summary)) + 1
integer_style <- createStyle(numFmt = "#,##0")
decimal_style <- createStyle(numFmt = "0.00")

addStyle(
    wb, sheet = "HPA_protocol_summary", style = integer_style,
    rows = summary_rows, cols = 2, gridExpand = TRUE, stack = TRUE
)

addStyle(
    wb, sheet = "HPA_protocol_summary", style = decimal_style,
    rows = summary_rows, cols = 3:ncol(protocol_summary), gridExpand = TRUE, stack = TRUE
)

# ST3: protein coverage intersections and subcellular localization ----
build_intersection_summary <- function(combination, set_counts, combination_counts, set_order) {
    combination_ids <- comb_name(combination)
    members <- matrix(as.integer(unlist(strsplit(combination_ids, "", fixed = TRUE))), nrow = length(combination_ids), byrow = TRUE,
                      dimnames = list(combination_ids, set_name(combination)))[, set_order, drop = FALSE]
    sizes <- as.integer(comb_size(combination))
    total_sizes <- as.integer(set_size(combination)[match(set_order, set_name(combination))])
    set_counts <- set_counts[set_order, names(subcellular_colors), drop = FALSE]
    combination_counts <- combination_counts[combination_ids, names(subcellular_colors), drop = FALSE]
    stopifnot(all(rowSums(set_counts) == total_sizes), all(rowSums(combination_counts) == sizes),
              all(as.vector(crossprod(members, sizes)) == total_sizes))
    headers <- c(set_order, "Intersection size", "Secreted", "Membrane", "Intracellular", "Secreted & membrane", "Unknown")
    totals <- data.frame(ifelse(diag(length(set_order)) == 1, "\u25CF", NA_character_), total_sizes, set_counts, check.names = FALSE)
    intersections <- data.frame(ifelse(members == 1, "\u25CF", NA_character_), sizes, combination_counts, check.names = FALSE)
    names(totals) <- names(intersections) <- headers
    intersections <- intersections[order(-sizes, combination_ids), , drop = FALSE]
    rownames(totals) <- rownames(intersections) <- NULL
    list(totals = totals, intersections = intersections)
}

st3_dia <- build_intersection_summary(dia_combination, dia_set_composition, dia_combination_composition, names(dia_protocol_proteins))
st3_platform <- build_intersection_summary(platform_combination, platform_set_composition, platform_combination_composition, names(platform_proteins))
st3_sheet <- "Feature_intersections"
if (st3_sheet %in% names(wb)) removeWorksheet(wb, st3_sheet)
addWorksheet(wb, st3_sheet)
st3_body <- createStyle(fontName = "Aptos Narrow", fontSize = 12, halign = "center", valign = "center")
st3_header <- createStyle(fontName = "Aptos Narrow", fontSize = 12, textDecoration = "bold", halign = "center", valign = "center")
st3_title <- createStyle(fontName = "Aptos Narrow", fontSize = 12, textDecoration = "bold", halign = "left", valign = "center")
st3_band <- createStyle(fontName = "Aptos Narrow", fontSize = 12, textDecoration = "bold", halign = "center", valign = "center", fgFill = "#E8E8E8")

write_intersection_section <- function(summary, title, start_row) {
    total_rows <- start_row + 2L + seq_len(nrow(summary$totals))
    intersection_label <- max(total_rows) + 1L
    intersection_rows <- intersection_label + seq_len(nrow(summary$intersections))
    last_row <- max(intersection_rows)
    writeData(wb, st3_sheet, title, startRow = start_row, colNames = FALSE)
    writeData(wb, st3_sheet, t(names(summary$totals)), startRow = start_row + 1L, colNames = FALSE)
    writeData(wb, st3_sheet, "(Totals)", startRow = start_row + 2L, colNames = FALSE)
    writeData(wb, st3_sheet, summary$totals, startRow = min(total_rows), colNames = FALSE, rowNames = FALSE)
    writeData(wb, st3_sheet, "(Intersections)", startRow = intersection_label, colNames = FALSE)
    writeData(wb, st3_sheet, summary$intersections, startRow = min(intersection_rows), colNames = FALSE, rowNames = FALSE)
    addStyle(wb, st3_sheet, st3_body, rows = start_row:last_row, cols = 1:11, gridExpand = TRUE)
    addStyle(wb, st3_sheet, st3_title, rows = start_row, cols = 1, stack = TRUE)
    addStyle(wb, st3_sheet, st3_header, rows = start_row + 1L, cols = 1:11, gridExpand = TRUE, stack = TRUE)
    addStyle(wb, st3_sheet, st3_band, rows = c(start_row + 2L, intersection_label), cols = 1:11, gridExpand = TRUE, stack = TRUE)
    addStyle(wb, st3_sheet, createStyle(numFmt = "#,##0"), rows = c(total_rows, intersection_rows), cols = 6:11, gridExpand = TRUE, stack = TRUE)
    setRowHeights(wb, st3_sheet, rows = start_row:last_row, heights = 16)
    last_row
}

st3_section_a_end <- write_intersection_section(st3_dia, "Section A: DIA-MS pre-analytical protocol intersections", 1L)
st3_section_b_end <- write_intersection_section(st3_platform, "Section B: Technology platform protein intersections (distinct UniProt accessions)", st3_section_a_end + 4L)
setColWidths(wb, st3_sheet, cols = 1:11, widths = c(15, 13.5, 12.16, 11, 13.83, 17, 9.83, 11.33, 12.5, 23, 10))
# These source-data worksheets contain no drawing parts to reference.
for (i in seq_along(wb$worksheets)) {
    wb$worksheets[[i]]$drawing <- character(0)
    wb$worksheets_rels[[i]] <- wb$worksheets_rels[[i]][!grepl("/(drawing|vmlDrawing)\"", wb$worksheets_rels[[i]])]
}
saveWorkbook(wb, output_file, overwrite = TRUE)

print(protocol_summary)
