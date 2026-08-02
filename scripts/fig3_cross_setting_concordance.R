# Figure 3 | Cross-setting concordance of titration-valid measurements

# 0. Setup ----
source("scripts/_project_setup.R")
use_packages(c("data.table", "tidyverse", "readxl", "ggridges", "ggrepel", "ggpubr", "patchwork", "openxlsx", "showtext", "clue"))
source("utils/figure_style.R")
plasmix_theme <- if (is.function(theme_plasmix)) theme_plasmix() else theme_plasmix
set.seed(20260718)
label_style <- list(size = 12, face = "bold")
to_bool <- function(x) {
    if (is.logical(x)) return(replace(x, is.na(x), FALSE))
    !is.na(x) & toupper(as.character(x)) %in% c("TRUE", "T", "1")
}
first_text <- function(x) {
    x <- x[!is.na(x) & x != ""]
    if (length(x)) x[[1]] else NA_character_
}
paths <- c(
    dea = "results/dea_df_multi.tsv.gz",
    trc = "results/trc_feature_level.tsv.gz",
    external_sex = "results/external_sex_effects_long.tsv.gz",
    feature_metadata = "data/feature_metadata.tsv.gz",
    metadata = "data/study_metadata.xlsx",
    ckb_correlation = upstream_path("cohort", "nc_ckb", "41467_2025_56935_MOESM4_ESM.xlsx"),
    ckb_traits = upstream_path("cohort", "nc_ckb", "41467_2025_56935_MOESM8_ESM.xlsx")
)
missing_inputs <- paths[!file.exists(paths)]
if (length(missing_inputs)) stop("Missing Figure 3 inputs:\n", paste(missing_inputs, collapse = "\n"))
selected_batches <- c("DIA_P1_B1", "DIA_P2_B1", "DIA_P3_B1", "DIA_P4_B1", "DIA_P5_B1", "DIA_P5_B2", "OLK_P2_B1", "OLK_P2_B2", "SOM_P1_B1", "SOM_P1_B2", "SOM_P2_B1", "SOM_P2_B2")
meta_batch <- read_xlsx(paths["metadata"], sheet = "batch") %>% filter(Batch %in% selected_batches) %>% select(Batch, Platform, Protocol) %>% distinct()
if (n_distinct(meta_batch$Batch) != length(selected_batches)) stop("study_metadata.xlsx does not contain all 12 high-throughput batches.")

# 1. Inputs and external-effect assembly ----
dea <- fread(paths["dea"]) %>% as_tibble()
external_sex <- fread(paths["external_sex"]) %>% as_tibble()
feature_metadata <- fread(paths["feature_metadata"]) %>% as_tibble() %>%
    filter(!to_bool(Is_Protein_Group), !to_bool(Is_Unknown), !is.na(UniProtID), UniProtID != "") %>%
    select(Platform, UniqueID, UniProtID) %>% distinct()
unique_feature_map <- feature_metadata %>% group_by(Platform, UniqueID) %>% filter(n_distinct(UniProtID) == 1) %>% ungroup()

plasmix_effect <- dea %>%
    filter(Pair == "M/F", Batch %in% c("OLK_P2_B1", "OLK_P2_B2"), Platform == "OLK", DataTier == "Calibrated", str_detect(ProcessLevel, regex("NPX", ignore_case = TRUE))) %>%
    mutate(UniProtID = tstrsplit(UniqueID, "_", fixed = TRUE)[[1]]) %>%
    group_by(UniProtID) %>% filter(n_distinct(Batch) == 2, all(logFC > 0, na.rm = TRUE) | all(logFC < 0, na.rm = TRUE)) %>%
    summarise(Effect = median(logFC, na.rm = TRUE), .groups = "drop")
tier_map <- external_sex %>% transmute(UniProtID = as.character(UniProt), Tier, Gene = as.character(Gene.Symbol)) %>%
    group_by(UniProtID) %>% summarise(Tier = first_text(Tier), Gene = first_text(Gene), .groups = "drop")
plasmix_effect <- plasmix_effect %>% left_join(tier_map, by = "UniProtID")
external_effect <- external_sex %>% filter(Cohort %in% c("UKBiobank", "Iceland", "Wellness", "BAMSE"), is.finite(Effect_Size)) %>%
    transmute(Cohort, UniProtID = as.character(UniProt), Effect = as.numeric(Effect_Size)) %>%
    group_by(Cohort, UniProtID) %>% summarise(Effect = median(Effect, na.rm = TRUE), .groups = "drop")

cohort_labels <- c(BAMSE = "BAMSE\n(Olink HT)", Wellness = "Wellness\n(Olink HT)", Iceland = "Iceland\n(SomaScan 5K)", UKBiobank = "UK Biobank\n(Olink 3072)", Plasmix = "Plasmix\n(Olink HT)")
cohort_levels <- unname(cohort_labels[c("BAMSE", "Wellness", "Iceland", "UKBiobank", "Plasmix")])
magnitude_data <- bind_rows(external_effect %>% mutate(Cohort = recode(Cohort, !!!cohort_labels)), plasmix_effect %>% transmute(Cohort = cohort_labels[["Plasmix"]], UniProtID, Effect)) %>%
    filter(is.finite(Effect)) %>% mutate(Cohort = factor(Cohort, levels = cohort_levels))
tier1_plasmix <- plasmix_effect %>% filter(Tier == "Tier1", is.finite(Effect)) %>% mutate(Cohort = factor(cohort_labels[["Plasmix"]], levels = cohort_levels))
tier1_labels <- tier1_plasmix %>% mutate(Label = coalesce(na_if(Gene, ""), UniProtID)) %>% arrange(desc(abs(Effect))) %>% distinct(Label, .keep_all = TRUE) %>% slice_head(n = 6)
cohort_colors <- c("Plasmix\n(Olink HT)" = "#56106E", "UK Biobank\n(Olink 3072)" = "#065EAD", "Iceland\n(SomaScan 5K)" = "#4DBBD5", "Wellness\n(Olink HT)" = "#D85170", "BAMSE\n(Olink HT)" = "#00A087")

# 2. Panel a: M/F effects in external cohorts and Plasmix ----
priority_genes <- c("SP3")
tier1_label_pool <- tier1_plasmix %>%
    mutate(Label = coalesce(na_if(Gene, ""), UniProtID), Priority = Label %in% priority_genes) %>%
    arrange(desc(Priority), desc(abs(Effect))) %>% distinct(Label, .keep_all = TRUE)
tier1_labels <- tier1_label_pool %>% slice_head(n = 7)

tail_cut <- 1
tail_scale <- .1
compress_tail_trans <- scales::trans_new(
    name = "compress_tail",
    transform = function(x) ifelse(x < -tail_cut, -tail_cut + tail_scale * (x + tail_cut),
                                   ifelse(x > tail_cut, tail_cut + tail_scale * (x - tail_cut), x)),
    inverse = function(x) ifelse(x < -tail_cut, -tail_cut + (x + tail_cut) / tail_scale,
                                 ifelse(x > tail_cut, tail_cut + (x - tail_cut) / tail_scale, x)),
    domain = c(-Inf, Inf)
)
tail_coord <- if ("coord_transform" %in% getNamespaceExports("ggplot2")) {
    coord_transform(x = compress_tail_trans, clip = "off")
} else {
    coord_trans(x = compress_tail_trans, clip = "off")
}
p_a <- ggplot(magnitude_data, aes(Effect, Cohort, fill = Cohort)) +
    stat_density_ridges(quantile_lines = TRUE, quantiles = 4, alpha = .75, scale = 1.5, linewidth = .2, color = "grey20", rel_min_height = .001) +
    geom_point(data = tier1_plasmix, aes(Effect, Cohort), inherit.aes = FALSE, shape = 124, color = "#D85170",
               size = 2, alpha = .65, position = position_nudge(y = -.1)) +
    geom_text_repel(data = tier1_labels, aes(Effect, Cohort, label = Label), inherit.aes = FALSE,
                             nudge_y = -.34, direction = "both", seed = 20260718, size = 2.25, color = "#B93858",
                             segment.color = "#D98AA0", segment.size = .25, min.segment.length = 0,
                             box.padding = .2, point.padding = .12, max.overlaps = Inf) +
    scale_fill_manual(values = cohort_colors) +
    plasmix_theme +
    scale_x_continuous(limits = c(-4, 8), breaks = c(-6, -4, -1, 0, 1, 4, 8), expand = expansion(mult = c(0, 0))) +
    tail_coord + scale_y_discrete(expand = expansion(mult = c(.06, .17))) +
    labs(x = expression("Male–female effect size (" * plain(log)[2] * " scale)"), y = NULL) +
    theme(legend.position = "none", panel.grid.major.y = element_blank())

# 3. Protein-batch and batch-pair results ----
trc <- fread(paths["trc"]) %>% as_tibble()
required_trc_fields <- c("Platform", "Batch", "UniqueID", "DataTier", "ProcessLevel", "Is_Detected", "ExpectedResponseValid", "M", "Y", "P", "X", "F")
missing_trc_fields <- setdiff(required_trc_fields, names(trc))
if (length(missing_trc_fields)) stop("trc_feature_level.tsv.gz is missing Figure 3 fields: ", paste(missing_trc_fields, collapse = ", "))
trc <- trc %>%
    mutate(across(c(M, Y, P, X, F), as.numeric), Is_Detected = to_bool(Is_Detected), ExpectedResponseValid = to_bool(ExpectedResponseValid)) %>%
    filter(Batch %in% selected_batches) %>%
    filter((Platform == "DIA" & DataTier == "Baseline" & ProcessLevel == "Intensity") |
           (Platform == "SOM" & DataTier == "Calibrated" & ProcessLevel == "Calibrate") |
           (Platform == "OLK" & DataTier == "Calibrated" & str_detect(ProcessLevel, regex("NPX", ignore_case = TRUE))))
if (!setequal(unique(trc$Batch), selected_batches)) stop("trc_feature_level.tsv.gz does not contain all 12 default batches")
if (anyDuplicated(trc[c("Batch", "UniqueID")])) stop("Duplicate Batch × UniqueID keys remain after default-process filtering")

protein_batch <- trc %>% inner_join(unique_feature_map, by = c("Platform", "UniqueID"), relationship = "many-to-one") %>%
    group_by(Batch, UniProtID) %>% filter(n_distinct(UniqueID) == 1) %>% ungroup() %>%
    left_join(meta_batch, by = c("Batch", "Platform"), relationship = "many-to-one") %>%
    mutate(Figure3_eligible = (Is_Detected | ExpectedResponseValid) & is.finite(M) & is.finite(F), MF_effect = M - F)

pair_type <- function(platform1, protocol1, platform2, protocol2) {
    if (platform1 != platform2) return("Cross-platform")
    if (protocol1 != protocol2) return("Cross-protocol")
    "Inter-batch"
}
platform_combination <- function(platform1, platform2) {
    if (platform1 == platform2) return(platform1)
    platform_order <- c("DIA", "OLK", "SOM")
    paste(platform_order[platform_order %in% c(platform1, platform2)], collapse = "-")
}

pair_tables <- map(combn(selected_batches, 2, simplify = FALSE), function(pair) {
    batch1 <- pair[[1]]; batch2 <- pair[[2]]
    meta1 <- filter(meta_batch, Batch == batch1); meta2 <- filter(meta_batch, Batch == batch2)
    left <- protein_batch %>% filter(Batch == batch1, Figure3_eligible) %>% select(UniProtID, M1 = M, Y1 = Y, P1 = P, X1 = X, F1 = F, Valid1 = ExpectedResponseValid)
    right <- protein_batch %>% filter(Batch == batch2, Figure3_eligible) %>% select(UniProtID, M2 = M, Y2 = Y, P2 = P, X2 = X, F2 = F, Valid2 = ExpectedResponseValid)
    inner_join(left, right, by = "UniProtID", relationship = "one-to-one") %>% mutate(
        pair_id = paste(batch1, batch2, sep = "__"), Batch1 = batch1, Batch2 = batch2,
        Pair_type = pair_type(meta1$Platform, meta1$Protocol, meta2$Platform, meta2$Protocol),
        Platform_combination = platform_combination(meta1$Platform, meta2$Platform),
        Effect1 = M1 - F1, Effect2 = M2 - F2,
        Validity_status = case_when(Valid1 & Valid2 ~ "Joint-valid", Valid1 | Valid2 ~ "One-valid", TRUE ~ "Neither-valid")
    )
})
names(pair_tables) <- map_chr(pair_tables, ~first(.x$pair_id))
pair_inventory <- map_dfr(pair_tables, ~summarise(.x, pair_id = first(pair_id), Batch1 = first(Batch1), Batch2 = first(Batch2), Pair_type = first(Pair_type), Platform_combination = first(Platform_combination), Common_evaluable = n(), Joint_valid = sum(Validity_status == "Joint-valid"), One_valid = sum(Validity_status == "One-valid"), Neither_valid = sum(Validity_status == "Neither-valid")))

# 4. Panel b: Concordance at the top ----
dea_cat <- dea %>% filter(Pair %in% c("M/F", "N/P"), Batch %in% selected_batches) %>%
    filter((Platform == "DIA" & DataTier == "Baseline" & ProcessLevel == "Intensity") |
           (Platform == "SOM" & DataTier == "Calibrated" & ProcessLevel == "Calibrate") |
           (Platform == "OLK" & DataTier == "Calibrated" & str_detect(ProcessLevel, regex("NPX", ignore_case = TRUE)))) %>%
    inner_join(unique_feature_map, by = c("Platform", "UniqueID"), relationship = "many-to-one") %>%
    group_by(Batch, Pair, UniProtID) %>% filter(n_distinct(UniqueID) == 1) %>% slice(1) %>% ungroup() %>%
    filter(is.finite(logFC))

cat_dea_tables <- unlist(map(combn(selected_batches, 2, simplify = FALSE), function(pair) {
    batch1 <- pair[[1]]; batch2 <- pair[[2]]
    meta1 <- filter(meta_batch, Batch == batch1); meta2 <- filter(meta_batch, Batch == batch2)
    map(c("M/F", "N/P"), function(contrast) {
        left <- dea_cat %>% filter(Batch == batch1, Pair == contrast) %>% select(UniProtID, Effect1 = logFC)
        right <- dea_cat %>% filter(Batch == batch2, Pair == contrast) %>% select(UniProtID, Effect2 = logFC)
        joined <- inner_join(left, right, by = "UniProtID", relationship = "one-to-one")
        if (!nrow(joined)) return(tibble())
        joined %>% mutate(
            pair_id = paste(batch1, batch2, sep = "__"), Contrast = contrast,
            Pair_type = pair_type(meta1$Platform, meta1$Protocol, meta2$Platform, meta2$Protocol),
            Color_Label = platform_combination(meta1$Platform, meta2$Platform)
        )
    })
}), recursive = FALSE)

calc_cat_pair <- function(data, step = 10, max_k = 10000) {
    n <- min(nrow(data), max_k)
    if (n < step) return(tibble())
    ranked1 <- arrange(data, desc(abs(Effect1)), UniProtID); ranked2 <- arrange(data, desc(abs(Effect2)), UniProtID)
    direction1 <- setNames(sign(data$Effect1), data$UniProtID); direction2 <- setNames(sign(data$Effect2), data$UniProtID)
    map_dfr(seq(step, n, by = step), function(k) {
        top1 <- head(ranked1$UniProtID, k); top2 <- head(ranked2$UniProtID, k); overlap <- intersect(top1, top2)
        tibble(TopN = k, Consistency = 100 * sum(direction1[overlap] == direction2[overlap]) / k)
    })
}
cat_pair_curves <- map_dfr(cat_dea_tables, function(data) {
    if (!nrow(data)) return(tibble())
    calc_cat_pair(data) %>% mutate(pair_id = first(data$pair_id), Pair_type = first(data$Pair_type), Color_Label = first(data$Color_Label), Contrast = first(data$Contrast))
})
# cat_curve_support <- cat_pair_curves %>%
#     group_by(Pair_type, Color_Label, Contrast, pair_id) %>% summarise(Curve_max = max(TopN), .groups = "drop") %>%
#     group_by(Pair_type, Color_Label, Contrast) %>% summarise(Common_max = min(Curve_max), Total_pairs = n_distinct(pair_id), .groups = "drop")
cat_summary <- cat_pair_curves %>%
    # left_join(cat_curve_support, by = c("Pair_type", "Color_Label", "Contrast")) %>% filter(TopN <= Common_max) %>%
    group_by(Pair_type, Color_Label, Contrast, TopN) %>%
    summarise(Consistency = mean(Consistency, na.rm = TRUE), Batch_pairs = n_distinct(pair_id), .groups = "drop")

cols_platform <- c(DIA = "#155289", SOM = "#B33E90", OLK = "#489FA7")
combo_colors <- c("DIA-OLK" = "#7f404a", "DIA-SOM" = "#5b4080", "OLK-SOM" = "#408073")
pair_type_labels <- c("Inter-batch" = "Within protocol", "Cross-protocol" = "Across protocols", "Cross-platform" = "Across platforms")

plot_cat_panel <- function(facet_name, show_y = FALSE) {
    sub_data <- filter(cat_summary, Pair_type == facet_name)
    if (!nrow(sub_data)) return(ggplot() + theme_void() + labs(title = facet_name))
    max_x <- if (facet_name == "Cross-platform") 2000 else 10000
    axis_breaks <- if (facet_name == "Cross-platform") c(10, 100, 2000) else c(10, 100, 1000, 10000)
    axis_hjust <- c(.1, rep(.5, length(axis_breaks) - 2), .9)
    k_vals <- seq(10, max_x, by = 10)
    direction_pool <- bind_rows(map(cat_dea_tables, ~filter(.x, Pair_type == facet_name, Contrast == "M/F") %>% select(Effect1, Effect2))) %>%
        pivot_longer(everything(), values_to = "Effect")
    p_pos <- mean(direction_pool$Effect > 0, na.rm = TRUE)
    p_match <- p_pos^2 + (1 - p_pos)^2
    random_ribbon <- tibble(TopN = k_vals) %>%
        mutate(p_rand = p_match * TopN / max_x, expected = 100 * p_rand,
               sd = 100 * sqrt(TopN * p_rand * (1 - p_rand)) / TopN,
               lower = pmax(0, expected - 1.96 * sd), upper = pmin(100, expected + 1.96 * sd))
    palette <- if (facet_name == "Cross-platform") combo_colors else cols_platform
    legend_pos <- if (facet_name == "Cross-platform") c(0, 1) else c(1, 0)
    legend_just <- if (facet_name == "Cross-platform") c(0, 1) else c(1, 0)
    ggplot(sub_data, aes(TopN, Consistency, color = Color_Label, linetype = Contrast)) +
        geom_ribbon(data = random_ribbon, aes(TopN, ymin = lower, ymax = upper), inherit.aes = FALSE, fill = "grey90", alpha = .6) +
        geom_line(data = random_ribbon, aes(TopN, expected), inherit.aes = FALSE, linetype = "dashed", color = "black", linewidth = .45) +
        geom_line(linewidth = .6) +
        scale_color_manual(values = palette, breaks = names(palette)) +
        scale_linetype_manual(values = c("M/F" = "solid", "N/P" = "dotted")) +
        scale_x_log10(limits = c(10, max_x), breaks = axis_breaks, expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
        labs(x = expression("Top " * italic(k) * " proteins"), y = if (show_y) "Concordance at top (%)" else NULL,
            title = unname(pair_type_labels[facet_name]), color = NULL, linetype = NULL) +
        plasmix_theme +
        theme(plot.title = element_text(size = 8.5, face = "bold", hjust = .5), panel.border = element_rect(color = "black", fill = NA, linewidth = .35),
              axis.text.x = element_text(hjust = axis_hjust),
              axis.line = element_blank(), panel.grid = element_blank(), legend.position = legend_pos, legend.justification = legend_just,
              legend.margin = margin(5, 5, 5, 5), legend.box.spacing = unit(0, "pt"), legend.spacing = unit(0, "pt"),
              legend.background = element_blank(), legend.key = element_blank())
}

p_b <- plot_cat_panel("Inter-batch", TRUE) | plot_cat_panel("Cross-protocol") | plot_cat_panel("Cross-platform")
p_b <- p_b & theme(plot.margin = margin(3, 2, 2.5, 5))
# p_b <- ggarrange(plot_cat_panel("Inter-batch", TRUE), plot_cat_panel("Cross-protocol"), plot_cat_panel("Cross-platform"), nrow = 1)

# 5. Panel c: Matched profile sMAPE ----
calculate_profile_smape <- function(M1, Y1, P1, X1, F1, M2, Y2, P2, X2, F2) {
    intermediates1 <- c(Y = Y1, P = P1, X = X1); intermediates2 <- c(Y = Y2, P = P2, X = X2)
    common_names <- names(intermediates1)[is.finite(intermediates1) & is.finite(intermediates2)]
    if (!is.finite(M1) || !is.finite(F1) || !is.finite(M2) || !is.finite(F2) || length(common_names) < 2) return(c(Common_intermediates = length(common_names), Profile_sMAPE = NA_real_))
    profile1 <- 2^(c(M = M1, intermediates1[common_names]) - F1); profile2 <- 2^(c(M = M2, intermediates2[common_names]) - F2)
    error <- abs(profile1 - profile2) / ((abs(profile1) + abs(profile2)) / 2)
    c(Common_intermediates = length(common_names), Profile_sMAPE = 100 * mean(error))
}
profile_tables <- map(pair_tables, function(data) {
    profile_metrics <- pmap_dfr(select(data, M1, Y1, P1, X1, F1, M2, Y2, P2, X2, F2), function(...) as_tibble_row(as.list(calculate_profile_smape(...))))
    bind_cols(data, profile_metrics) %>% mutate(Mean_abs_MF = (abs(Effect1) + abs(Effect2)) / 2) %>% filter(Common_intermediates >= 2, is.finite(Profile_sMAPE), is.finite(Mean_abs_MF))
})

standardized_mean_difference <- function(x, y) {
    if (length(x) < 2 || length(y) < 2) return(NA_real_)
    pooled <- sqrt((var(x) + var(y)) / 2)
    if (!is.finite(pooled) || pooled == 0) return(NA_real_)
    (mean(x) - mean(y)) / pooled
}
match_profile_error <- function(data) {
    joint <- filter(data, Validity_status == "Joint-valid") %>% arrange(UniProtID)
    neither <- filter(data, Validity_status == "Neither-valid") %>% arrange(UniProtID)
    if (!nrow(joint) || !nrow(neither)) return(NULL)
    x_joint <- log1p(joint$Mean_abs_MF); x_neither <- log1p(neither$Mean_abs_MF)
    if (nrow(joint) <= nrow(neither)) {
        assignment <- as.integer(solve_LSAP(abs(outer(x_joint, x_neither, "-"))))
        matched_joint <- joint; matched_neither <- neither[assignment, , drop = FALSE]
    } else {
        assignment <- as.integer(solve_LSAP(abs(outer(x_neither, x_joint, "-"))))
        matched_joint <- joint[assignment, , drop = FALSE]; matched_neither <- neither
    }
    matched_pairs <- tibble(pair_id = first(data$pair_id), match_id = seq_len(nrow(matched_joint)), Joint_UniProtID = matched_joint$UniProtID, Neither_UniProtID = matched_neither$UniProtID, Joint_Mean_abs_MF = matched_joint$Mean_abs_MF, Neither_Mean_abs_MF = matched_neither$Mean_abs_MF, Joint_sMAPE = matched_joint$Profile_sMAPE, Neither_sMAPE = matched_neither$Profile_sMAPE)
    list(
        matched = matched_pairs,
        balance = tibble(pair_id = first(data$pair_id), Pair_type = first(data$Pair_type), Platform_combination = first(data$Platform_combination), Joint_before = nrow(joint), Neither_before = nrow(neither), Matched_each = nrow(matched_joint), Absolute_SMD_before = abs(standardized_mean_difference(x_joint, x_neither)), Absolute_SMD_after = abs(standardized_mean_difference(log1p(matched_joint$Mean_abs_MF), log1p(matched_neither$Mean_abs_MF)))),
        result = tibble(pair_id = first(data$pair_id), Batch1 = first(data$Batch1), Batch2 = first(data$Batch2), Pair_type = first(data$Pair_type), Platform_combination = first(data$Platform_combination), Matched_each = nrow(matched_joint), Joint_median_sMAPE = median(matched_joint$Profile_sMAPE), Neither_median_sMAPE = median(matched_neither$Profile_sMAPE), Neither_minus_joint_sMAPE = median(matched_neither$Profile_sMAPE) - median(matched_joint$Profile_sMAPE))
    )
}
profile_matches <- compact(map(profile_tables, match_profile_error))
profile_matched_proteins <- map_dfr(profile_matches, "matched")
profile_matching_balance <- map_dfr(profile_matches, "balance")
profile_pair_results <- map_dfr(profile_matches, "result")
if (!nrow(profile_pair_results)) stop("No batch pair contains both joint-valid and neither-valid proteins for Panel c")
profile_summary <- profile_pair_results %>% group_by(Pair_type) %>% summarise(Eligible_pairs = n(), Pairs_favoring_joint = sum(Neither_minus_joint_sMAPE > 0), Median_difference = median(Neither_minus_joint_sMAPE), .groups = "drop")

pair_levels <- c("Inter-batch", "Cross-protocol", "Cross-platform")
pair_display_labels_c <- c("Inter-batch" = "Within\nprotocol", "Cross-protocol" = "Across\nprotocols", "Cross-platform" = "Across\nplatforms")
pair_display_levels_c <- rev(unname(pair_display_labels_c[pair_levels]))

setting_levels <- c("DIA", "SOM", "OLK", "DIA-OLK", "DIA-SOM", "OLK-SOM")
setting_colors <- c(cols_platform, combo_colors)
setting_shapes <- setNames(c(21, 24, 22, 21, 24, 22), setting_levels)

profile_pair_results <- profile_pair_results %>%
    mutate(Pair_display = factor(unname(pair_display_labels_c[Pair_type]), levels = pair_display_levels_c),
           Platform_combination = factor(Platform_combination, levels = setting_levels))

profile_summary <- profile_pair_results %>%
    group_by(Pair_type) %>%
    summarise(Eligible_pairs = n(), Positive_pairs = sum(Neither_minus_joint_sMAPE > 0),
              Median_difference = median(Neither_minus_joint_sMAPE), .groups = "drop") %>%
    mutate(Pair_display = factor(unname(pair_display_labels_c[Pair_type]), levels = pair_display_levels_c),
           Count_label = sprintf("%d/%d pairs", Positive_pairs, Eligible_pairs))

profile_medians <- profile_pair_results %>%
    group_by(Pair_display) %>%
    summarise(Neither_minus_joint_sMAPE = median(Neither_minus_joint_sMAPE), .groups = "drop")

count_header <- tibble(Pair_display = factor(tail(pair_display_levels_c, 1), levels = pair_display_levels_c), Header = "Positive difference")

p_c <- ggplot(profile_pair_results,
              aes(Neither_minus_joint_sMAPE, Pair_display, fill = Platform_combination, shape = Platform_combination)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = .4) +
    geom_jitter(height = .22, width = 0, alpha = 0.95, size = 2.5, stroke = .5, color = "white") +
    geom_segment(data = profile_medians, aes(x = 0, xend = Neither_minus_joint_sMAPE, y = Pair_display, yend = Pair_display),
                 inherit.aes = FALSE, color = "black", linewidth = .45, position = position_nudge(y = .27)) +
    geom_point(data = profile_medians, aes(Neither_minus_joint_sMAPE, Pair_display), inherit.aes = FALSE,
               shape = 18, size = 3.2, color = "black", position = position_nudge(y = .27)) +
    geom_text(data = profile_summary, aes(Inf, Pair_display, label = Count_label), inherit.aes = FALSE,
              hjust = 1.04, size = 2.8, position = position_nudge(y = .27)) +
    geom_text(data = count_header, aes(Inf, Pair_display, label = Header), inherit.aes = FALSE,
              hjust = 1.04, fontface = "bold", size = 2.8, position = position_nudge(y = .62)) +
    scale_fill_manual(values = setting_colors, breaks = setting_levels, name = NULL) +
    scale_shape_manual(values = setting_shapes, breaks = setting_levels, name = NULL) +
    scale_y_discrete(expand = expansion(add = c(.5, .75))) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(alpha = 1, size = 2.3, stroke = .5, color = "white", shape = unname(setting_shapes[setting_levels]))), shape = "none") +
    labs(x = "Profile sMAPE difference: neither pass − joint pass", y = NULL, color = NULL, shape = NULL) +
    plasmix_theme +
    theme(
        legend.position = "top", legend.direction = "horizontal", legend.background = element_blank(),
        legend.margin = margin(-3.5, 5, 5, -10), legend.box.spacing = unit(0, "pt"), legend.spacing.x = unit(.02, "cm"),
        legend.key.width = unit(.3, "cm"), legend.key.height = unit(.3, "cm"),
        legend.text = element_text(vjust = .5, lineheight = 1, margin = margin(r = .1)),
        axis.text.y = element_text(lineheight = .85),
        panel.grid.major.y = element_blank()
    )

# 6. Panels d-e: CKB analyses ----
id_columns <- c("uniprot_id", "olink_id", "somascan_id")
ckb_correlation_raw <- read_excel(paths["ckb_correlation"], sheet = "Supplementary_Data_1") %>% distinct(across(all_of(id_columns)), .keep_all = TRUE)
ckb_traits_raw <- read_excel(paths["ckb_traits"], sheet = "Supplementary_Data_5") %>% distinct(across(all_of(id_columns)), .keep_all = TRUE)
strict_pairs <- ckb_correlation_raw %>% select(all_of(id_columns)) %>% drop_na() %>% distinct() %>%
    group_by(uniprot_id) %>% filter(n_distinct(olink_id) == 1, n_distinct(somascan_id) == 1) %>% ungroup() %>%
    group_by(olink_id) %>% filter(n_distinct(uniprot_id) == 1) %>% ungroup() %>%
    group_by(somascan_id) %>% filter(n_distinct(uniprot_id) == 1) %>% ungroup() %>% distinct(uniprot_id, .keep_all = TRUE)
ckb_raw <- strict_pairs %>% inner_join(ckb_correlation_raw, by = id_columns, relationship = "one-to-one") %>%
    inner_join(ckb_traits_raw, by = id_columns, relationship = "one-to-one") %>% rename(UniProtID = uniprot_id)

ckb_protein_batch <- protein_batch %>% filter(Batch %in% c("OLK_P2_B1", "OLK_P2_B2", "SOM_P1_B1", "SOM_P1_B2"), Figure3_eligible)
collapse_ckb_platform <- function(data, prefix) {
    data %>% group_by(UniProtID) %>% summarise(Valid_count = sum(ExpectedResponseValid), Batch_count = n_distinct(Batch), MF_effect = median(MF_effect, na.rm = TRUE), .groups = "drop") %>% rename_with(~paste0(prefix, "_", .x), -UniProtID)
}
ckb_membership <- inner_join(collapse_ckb_platform(filter(ckb_protein_batch, Platform == "OLK"), "Olink"), collapse_ckb_platform(filter(ckb_protein_batch, Platform == "SOM"), "Soma"), by = "UniProtID", relationship = "one-to-one") %>%
    filter(Olink_Batch_count == 2, Soma_Batch_count == 2) %>%
    mutate(Joint_pass_count = Olink_Valid_count * Soma_Valid_count, Replicated_joint_pass = Joint_pass_count >= 2, Zero_joint_pass = Joint_pass_count == 0, MF_magnitude = (abs(Olink_MF_effect) + abs(Soma_MF_effect)) / 2)

ckb_base <- ckb_membership %>% inner_join(ckb_raw, by = "UniProtID", relationship = "one-to-one") %>% filter(is.finite(MF_magnitude), is.finite(rho_olink_soma_non_ANML)) %>% arrange(UniProtID)
ckb_joint <- filter(ckb_base, Replicated_joint_pass) %>% arrange(UniProtID)
ckb_zero_pool <- filter(ckb_base, Zero_joint_pass) %>% arrange(UniProtID)
if (!nrow(ckb_joint) || nrow(ckb_joint) > nrow(ckb_zero_pool)) stop("CKB has insufficient replicated-pass and zero-joint controls for 1:1 matching")
ckb_cost <- abs(outer(log1p(ckb_joint$MF_magnitude), log1p(ckb_zero_pool$MF_magnitude), "-"))
ckb_zero <- ckb_zero_pool[as.integer(solve_LSAP(ckb_cost)), , drop = FALSE]
ckb_measurement <- tibble(match_id = seq_len(nrow(ckb_joint)), joint_UniProtID = ckb_joint$UniProtID, zero_UniProtID = ckb_zero$UniProtID, joint_mf_magnitude = ckb_joint$MF_magnitude, zero_mf_magnitude = ckb_zero$MF_magnitude, match_log1p_distance = abs(log1p(ckb_joint$MF_magnitude) - log1p(ckb_zero$MF_magnitude)), rho_zero = ckb_zero$rho_olink_soma_non_ANML, rho_joint = ckb_joint$rho_olink_soma_non_ANML, delta_rho = rho_joint - rho_zero)
bootstrap_median_ci <- function(x, n_boot = 1000, seed = 20260718) {
    x <- x[is.finite(x)]; set.seed(seed)
    unname(quantile(replicate(n_boot, median(sample(x, length(x), replace = TRUE))), c(.025, .975)))
}
rho_ci <- bootstrap_median_ci(ckb_measurement$delta_rho)
ckb_summary <- tibble(n_replicated_joint_pass = nrow(ckb_joint), n_zero_joint_pass_original = nrow(ckb_zero_pool), n_matched_controls = nrow(ckb_zero), median_mf_joint = median(ckb_joint$MF_magnitude), median_mf_zero_before = median(ckb_zero_pool$MF_magnitude), median_mf_zero_after = median(ckb_zero$MF_magnitude), median_match_distance = median(ckb_measurement$match_log1p_distance), median_rho_zero = median(ckb_measurement$rho_zero), median_rho_joint = median(ckb_measurement$rho_joint), median_paired_delta_rho = median(ckb_measurement$delta_rho), bootstrap_ci_low = rho_ci[1], bootstrap_ci_high = rho_ci[2])

safe_cor <- function(x, y, method) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
    cor(x[ok], y[ok], method = method)
}
trait_names <- setdiff(sub("^olink_es_", "", grep("^olink_es_", names(ckb_raw), value = TRUE)), "is_female")
if (length(trait_names) != 18) stop("The CKB supplement does not contain exactly 18 non-sex traits: ", paste(trait_names, collapse = ", "))
trait_labels <- c(age = "Age", alcohol_regular_vs_occasion = "Regular alcohol use", ambient_temperature = "Ambient temperature", bmi = "BMI", cancer_diagnosis = "Cancer diagnosis", dbp = "Diastolic blood pressure", diabetes_diagnosis = "Diabetes diagnosis", heart_rate = "Heart rate", hours_since_last_ate = "Hours since last meal", kidney_dis_diagnosis = "Kidney disease diagnosis", married = "Married", physical_activity = "Physical activity", poor_health = "Poor health", random_glucose = "Random glucose", region_is_urban = "Urban region", sbp = "Systolic blood pressure", school = "School education", smoking_ever_regular = "Ever regular smoking")
ckb_traits <- map_dfr(sort(trait_names), function(trait) {
    joint_data <- filter(ckb_raw, UniProtID %in% ckb_joint$UniProtID); zero_data <- filter(ckb_raw, UniProtID %in% ckb_zero$UniProtID)
    olink_col <- paste0("olink_es_", trait); soma_col <- paste0("soma_non_ANML_es_", trait)
    tibble(trait = trait, trait_label = unname(trait_labels[trait]), n_joint = sum(is.finite(joint_data[[olink_col]]) & is.finite(joint_data[[soma_col]])), n_zero = sum(is.finite(zero_data[[olink_col]]) & is.finite(zero_data[[soma_col]])), pearson_zero = safe_cor(zero_data[[olink_col]], zero_data[[soma_col]], "pearson"), pearson_joint = safe_cor(joint_data[[olink_col]], joint_data[[soma_col]], "pearson"), spearman_zero = safe_cor(zero_data[[olink_col]], zero_data[[soma_col]], "spearman"), spearman_joint = safe_cor(joint_data[[olink_col]], joint_data[[soma_col]], "spearman"))
}) %>% mutate(pearson_delta = pearson_joint - pearson_zero, spearman_delta = spearman_joint - spearman_zero)

measurement_long <- ckb_measurement %>%
    select(match_id, rho_zero, rho_joint) %>%
    pivot_longer(c(rho_zero, rho_joint), names_to = "Group", values_to = "Correlation") %>%
    mutate(Group = factor(Group, c("rho_zero", "rho_joint"), c("Absent", "Supported")))

measurement_annotation <- sprintf("Median paired Δρ = %.2f\n95%% CI %.2f–%.2f",
                                  ckb_summary$median_paired_delta_rho, ckb_summary$bootstrap_ci_low, ckb_summary$bootstrap_ci_high)

p_d <- ggplot(measurement_long, aes(Group, Correlation, group = match_id)) +
    geom_line(color = "grey75", alpha = .4, linewidth = .3) +
    geom_boxplot(aes(group = Group), width = .42, outlier.shape = NA, fill = "white", color = "black", linewidth = .35) +
    geom_point(aes(color = Group), alpha = .8, size = 1.15) +
    annotate("text", x = 1.5, y = Inf, label = measurement_annotation, hjust = .5, vjust = 0, size = 2.8, lineheight = .95) +
    scale_color_manual(values = c("Absent" = "grey55", "Supported" = "#D85170"), guide = "none") +
    labs(x = NULL, y = "Participant-level Spearman ρ") +
    plasmix_theme +
    scale_y_continuous(expand = expansion(mult = c(.05, .08))) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          plot.margin = margin(17.5, 5, 5, 5))

trait_plot <- ckb_traits %>%
    mutate(Trait_label = factor(trait_label, levels = trait_label[order(pearson_delta)]))

p_e <- ggplot(trait_plot, aes(y = Trait_label)) +
    geom_segment(aes(x = pearson_zero, xend = pearson_joint, yend = Trait_label), color = "grey72", linewidth = .55) +
    geom_point(aes(x = pearson_zero, color = "Absent"), size = 2) +
    geom_point(aes(x = pearson_joint, color = "Supported"), size = 2) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("%d/%d traits higher; Median Δr = %.2f",
                             sum(ckb_traits$pearson_delta > 0, na.rm = TRUE),
                             nrow(ckb_traits), median(ckb_traits$pearson_delta, na.rm = TRUE)),
             hjust = 1.2, vjust = 0, size = 2.8) +
    scale_color_manual(values = c("Absent" = "grey55", "Supported" = "#D85170"), breaks = c("Absent", "Supported")) +
    guides(color = guide_legend(title.position = "left", title.hjust = 0, nrow = 1, byrow = TRUE,
                                override.aes = list(alpha = 1, size = 2))) +
    labs(x = "Olink–SomaScan phenotype-effect Pearson correlation", y = NULL, color = "Joint pass status") +
    plasmix_theme +
    scale_y_discrete(expand = expansion(add = c(.35, 0.9))) +
    theme(
        legend.position = "top", legend.justification = "left",
        legend.direction = "horizontal", legend.background = element_blank(),
        legend.margin = margin(-3.5, 0, 10, 0), legend.box.spacing = unit(0, "pt"),
        legend.spacing.x = unit(.02, "cm"), legend.key.width = unit(.3, "cm"), legend.key.height = unit(.3, "cm"),
        legend.title = element_text(margin = margin(r = 1)), legend.text = element_text(margin = margin(r = .1)),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(hjust = 1),
        plot.margin = margin(5, 5, 5, 5)
    )

# 7. Assemble and export ----
showtext_auto()
showtext_opts(dpi = 600)
row1 <- ggarrange(p_a, p_b, nrow = 1, widths = c(1, 2), labels = c("a", "b"), font.label = label_style, label.x = 0, label.y = 1, hjust = -.2, vjust = 1)
row2 <- ggarrange(p_c, p_d, p_e, nrow = 1, widths = c(1, 0.7, 1), labels = c("c", "d", "e"), font.label = label_style, label.x = 0, label.y = 1, hjust = -.2, vjust = 1)
figure3 <- ggarrange(row1, row2, ncol = 1, heights = c(1, 1.1))
ggsave("figures/fig3_cross_setting_concordance.pdf", figure3, width = 10, height = 6)
ggsave("figures/fig3_cross_setting_concordance.png", figure3, width = 10, height = 6, dpi = 600, bg = "white")

write.xlsx(list(
    Cohort_magnitude = magnitude_data, Plasmix_Tier1 = tier1_plasmix, Pair_inventory = pair_inventory,
    CAT_pair_curves = cat_pair_curves, CAT_summary = cat_summary,
    Profile_sMAPE_pair_results = profile_pair_results, Profile_sMAPE_matching_balance = profile_matching_balance, Profile_sMAPE_matched_proteins = profile_matched_proteins,
    CKB_measurement_pairs = ckb_measurement, CKB_measurement_summary = ckb_summary,
    CKB_trait_Pearson = ckb_traits %>% select(trait, trait_label, n_joint, n_zero, pearson_zero, pearson_joint, pearson_delta),
    CKB_trait_Spearman = ckb_traits %>% select(trait, trait_label, n_joint, n_zero, spearman_zero, spearman_joint, spearman_delta)
), "tables/SourceData_Figure3.xlsx", overwrite = TRUE, keepNA = TRUE, na.string = "NA")

figure3
message("Figure 3 completed: figures/fig3_cross_setting_concordance.pdf")
