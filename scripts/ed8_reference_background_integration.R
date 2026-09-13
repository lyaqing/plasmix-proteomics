# Extended Data Figure 8 | Reference choice and background subtraction

# 0. Setup ----
source("scripts/_project_setup.R")
use_packages(c("data.table", "tidyverse", "readxl", "patchwork", "ggpubr", "openxlsx", "showtext"))
source("utils/figure_style.R")
source("utils/benchmark_metrics.R")
source("utils/feature_mapping.R")
plasmix_theme <- if (is.function(theme_plasmix)) theme_plasmix() else theme_plasmix
showtext_auto(); showtext_opts(dpi = 600)
label_style <- list(size = 12, face = "bold")
paths <- c(metadata = "data/study_metadata.xlsx", feature_metadata = "data/feature_metadata.tsv.gz", profiles = "data/protein_profiles_long.tsv.gz",
           fig6_inputs = "results/fig6_extended_data_inputs.rds")
missing_inputs <- paths[!file.exists(paths)]
if (length(missing_inputs)) stop("Extended Data Figure 8 is missing input files:\n", paste(missing_inputs, collapse = "\n"))

# 1. Inputs and shared definitions ----
meta_batch <- read_xlsx(paths["metadata"], sheet = "batch")
feature_metadata <- fread(paths["feature_metadata"]) %>% as_tibble()
single_accession_df <- fread(paths["profiles"]) %>% as_tibble() %>% filter_batch_analysis_features(feature_metadata, strict_platforms = character())
long_df <- aggregate_som_profiles(single_accession_df) %>% filter_batch_analysis_features(analysis_feature_metadata(feature_metadata))
fig6_ed_inputs <- readRDS(paths["fig6_inputs"])
required_cache_objects <- c("srr_reference_feature_results", "srr_reference_pair_summary")
missing_cache_objects <- setdiff(required_cache_objects, names(fig6_ed_inputs))
if (length(missing_cache_objects)) stop("fig6_extended_data_inputs.rds lacks: ", paste(missing_cache_objects, collapse = ", "))

srr_feature_results <- fig6_ed_inputs$srr_reference_feature_results
srr_pair_summary <- fig6_ed_inputs$srr_reference_pair_summary
smape_cutoff <- if (!is.null(fig6_ed_inputs$thresholds$sMAPE)) fig6_ed_inputs$thresholds$sMAPE else 30
trc_cutoff <- if (!is.null(fig6_ed_inputs$thresholds$TRC_deviation)) fig6_ed_inputs$thresholds$TRC_deviation else 0.25
if (smape_cutoff <= 1) smape_cutoff <- 100 * smape_cutoff

detailed_type_levels <- c("Intra-DIA", "Intra-OLK", "Intra-SOM", "DIA-OLK", "DIA-SOM", "OLK-SOM")
batch_levels <- c("DIA_P1_B1", "DIA_P2_B1", "DIA_P3_B1", "DIA_P4_B1", "DIA_P5_B1", "DIA_P5_B2",
                  "OLK_P2_B1", "OLK_P2_B2", "SOM_P1_B1", "SOM_P1_B2", "SOM_P2_B1", "SOM_P2_B2")
selected_batches <- intersect(batch_levels, unique(c(srr_pair_summary$Batch1, srr_pair_summary$Batch2)))
affinity_batches <- intersect(selected_batches, c("OLK_P2_B1", "OLK_P2_B2", "SOM_P1_B1", "SOM_P1_B2", "SOM_P2_B1", "SOM_P2_B2"))

rbb_colors <- c("OLK_P2_B1"="#389191", "OLK_P2_B2"="#226b6b", "SOM_P1_B1"="#d162b5", "SOM_P1_B2"="#830062", "SOM_P2_B1"="#f0a22e", "SOM_P2_B2"="#c85203")
plot_batch_colors <- batch_color[selected_batches]
plot_batch_shapes <- setNames(c(16, 17, 15, 18, 3, 4, 8, 7, 0, 1, 2, 5)[seq_along(selected_batches)], selected_batches)

parse_reference <- function(method, method_reference = NULL) {
    parsed <- case_when(str_detect(method, "\\(P\\)$") ~ "P", str_detect(method, "\\(N\\)$") ~ "N", TRUE ~ NA_character_)
    if (!is.null(method_reference)) parsed <- coalesce(as.character(method_reference), parsed)
    parsed
}

metric_levels <- c("Rate_sMAPE", "Rate_Response", "Rate_Harmonized")
metric_labels <- c("Quantitative agreement", "Expected response", "Harmonized")
metric_colors <- c("Quantitative agreement" = "#d4ae93", "Expected response" = "#d3605e", "Harmonized" = "#531E1B")
reference_shapes <- c("P" = 23, "N" = 9)

# Separate text labels only; points and numerical summaries retain their original values.
spread_reference_labels <- function(y, gap = 0.045) {
    if (length(y) < 2L || any(!is.finite(y))) return(y)
    ord <- order(y); offset <- seq_along(y) * gap
    placed <- stats::isoreg(seq_along(y), y[ord] - offset)$yf + offset
    placed <- placed + max(0, 0.015 - min(placed)) - max(0, max(placed) - 0.985)
    out <- y; out[ord] <- placed; out
}

make_reference_plot <- function(pair_values, show_legend = TRUE) {
    pair_values <- pair_values %>%
        mutate(Detailed_Type = factor(Detailed_Type, levels = detailed_type_levels), Reference = factor(Reference, levels = c("P", "N")))
    summary_values <- pair_values %>%
        group_by(Detailed_Type, Reference) %>%
        summarize(across(all_of(metric_levels), ~median(.x, na.rm = TRUE)), N_Pairs = n(), .groups = "drop")
    pair_long <- pair_values %>% pivot_longer(all_of(metric_levels), names_to = "Metric", values_to = "Proportion") %>%
        mutate(Metric = factor(Metric, levels = metric_levels, labels = metric_labels))
    summary_long <- summary_values %>% pivot_longer(all_of(metric_levels), names_to = "Metric", values_to = "Proportion") %>%
        mutate(Metric = factor(Metric, levels = metric_levels, labels = metric_labels), X_Point = as.numeric(Reference) - 0.2) %>%
        group_by(Detailed_Type, Reference) %>% mutate(Label_Y = spread_reference_labels(Proportion)) %>% ungroup()
    connector_data <- summary_long %>% group_by(Detailed_Type, Reference, X_Point) %>%
        summarize(Y_Min = min(Proportion, na.rm = TRUE), Y_Max = max(Proportion, na.rm = TRUE), .groups = "drop")

    p <- ggplot(summary_long, aes(y = Proportion, color = Metric, fill = Metric, shape = Reference)) +
        geom_point(data = pair_long, aes(x = as.numeric(Reference) - 0.2, y = Proportion, color = Metric, fill = Metric),
                   alpha = 0.15, size = 0.5, position = position_jitter(width = 0.12)) +
        geom_segment(data = connector_data, aes(x = X_Point, xend = X_Point, y = Y_Min, yend = Y_Max), inherit.aes = FALSE,
                     color = "grey70", linewidth = 0.5) +
        geom_point(aes(x = X_Point), size = 2) +
        geom_text(aes(x = as.numeric(Reference), y = Label_Y, label = scales::percent(Proportion, accuracy = 0.1)),
                  hjust = 0, size = 2.5, fontface = "bold", show.legend = FALSE) +
        facet_wrap(~Detailed_Type, scales = "free_x", nrow = 1, drop = TRUE) +
        scale_x_continuous(breaks = c(1, 2), labels = c("P", "N"), limits = c(0.5, 2.5), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = c("0", "25", "50", "75", "100"),
                           expand = expansion(mult = c(0, 0.02))) +
        scale_fill_manual(values = metric_colors) + scale_color_manual(values = metric_colors) +
        scale_shape_manual(values = reference_shapes, labels = c("P" = "Matched (P)", "N" = "Unmatched (N)")) +
        labs(x = NULL, y = "Proportion of proteins (%)", color = "Metric", shape = "Matrix") +
        guides(shape = guide_legend(order = 1), color = guide_legend(order = 2), fill = "none") + plasmix_theme +
        theme(panel.grid.major.x = element_blank(), legend.position = if (show_legend) "bottom" else "none",
              legend.box = "horizontal", legend.box.just = "left", legend.direction = "horizontal",
              legend.spacing.x = unit(0.05, "cm"), legend.title = element_text(size = 7.5, face = "bold", vjust = 0.5),
              legend.text = element_text(size = 7.5, vjust = 0.5, margin = margin(0, 0, 0, 3)),
              legend.box.margin = margin(t = -10, b = -5, l = 0))
    list(plot = p, pair = pair_long, summary = summary_long, connectors = connector_data)
}

# 2. Panel a: standard SRR using P or N, retaining the former Figure 6 panel-a design ----
standard_pair_values <- srr_pair_summary %>% filter(Design == "Balanced", Method %in% c("SRR (P)", "SRR (N)"))
if ("MethodReference" %in% names(standard_pair_values)) {
    standard_pair_values <- standard_pair_values %>% mutate(Reference = parse_reference(Method, MethodReference))
} else {
    standard_pair_values <- standard_pair_values %>% mutate(Reference = parse_reference(Method))
}
standard_pair_values <- standard_pair_values %>% filter(Reference %in% c("P", "N")) %>%
    select(Detailed_Type, Macro_Type, Batch1, Batch2, Design, Method, Reference, N_Features,
           Rate_sMAPE, Rate_Response, Rate_Harmonized)
if (!nrow(standard_pair_values)) stop("No Balanced SRR(P/N) pair summaries were found in fig6_extended_data_inputs.rds.")
standard_plot <- make_reference_plot(standard_pair_values, show_legend = FALSE)
p_standard_srr <- standard_plot$plot +
    theme(legend.position = "bottom", legend.box.margin = margin(t = -5, b = 0, r = 0, l = 0),
          legend.text = element_text(size = 7.5, vjust = 0.5, margin = margin(0, 0, 0, 3))) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE), shape = guide_legend(nrow = 1, byrow = TRUE))

# 3. Panel b: three representative protein trajectories, retaining the former Figure 6 panel-b design ----
target_genes_map <- c("INSL3" = "P51460", "SPINT3" = "P49223", "PZP" = "P20742") #"LEP" = "P41159")
representative_batches <- c("DIA_P3_B1", "DIA_P4_B1", "OLK_P2_B1", "OLK_P2_B2", "SOM_P1_B2", "SOM_P2_B1")
representative_values_all <- long_df %>%
    filter(Batch %in% representative_batches, UniqueID %in% unname(target_genes_map), Sample %in% c("M", "Y", "P", "X", "F"),
           (Platform %in% c("OLK", "SOM") & ProcessLevel == "HybNorm") | (Platform == "DIA" & ProcessLevel == "Intensity")) %>%
    mutate(Normalized_Linear = 2^Value, Fraction = recode(Sample, M = 1, Y = 0.75, P = 0.50, X = 0.25, F = 0),
           GeneSymbol = names(target_genes_map)[match(UniqueID, unname(target_genes_map))],
           GeneSymbol = factor(GeneSymbol, levels = names(target_genes_map)), Batch = factor(Batch, levels = selected_batches))

missing_representatives <- setdiff(names(target_genes_map), unique(as.character(na.omit(representative_values_all$GeneSymbol))))
if (length(missing_representatives)) stop("Representative-protein input is missing: ", paste(missing_representatives, collapse = ", "))
message("Representative proteins found before fitting: ", paste(names(target_genes_map), collapse = ", "))

linear_fit_quality <- representative_values_all %>% filter(is.finite(Normalized_Linear), is.finite(Fraction)) %>%
    group_by(Platform, Batch, UniqueID, GeneSymbol) %>%
    summarize(R2_Lin = if (n() >= 3 && n_distinct(Fraction) >= 2) summary(lm(Normalized_Linear ~ Fraction))$r.squared else NA_real_,
              .groups = "drop")

representative_values <- representative_values_all %>%
    inner_join(linear_fit_quality %>% filter(is.finite(R2_Lin), R2_Lin > 0.6), by = c("Platform", "Batch", "UniqueID", "GeneSymbol")) %>%
    group_by(Platform, Batch, GeneSymbol, UniqueID) %>%
    mutate(P_Reference = if (any(Sample == "P" & is.finite(Normalized_Linear))) mean(Normalized_Linear[Sample == "P" & is.finite(Normalized_Linear)]) else NA_real_,
           SRR = Normalized_Linear / P_Reference) %>% ungroup() %>% filter(is.finite(SRR), SRR > 0)
if (!nrow(representative_values)) stop("No representative-protein observations passed the linear-fit filter.")
# Display two batches per platform; individual proteins still require R2 > 0.6.
missing_after_fit <- setdiff(names(target_genes_map), unique(as.character(na.omit(representative_values$GeneSymbol))))
if (length(missing_after_fit)) stop("Representative proteins removed by the R2 > 0.6 filter: ", paste(missing_after_fit, collapse = ", "))
message("Representative proteins retained after fitting: ", paste(names(target_genes_map), collapse = ", "))
displayed_batches <- unique(as.character(representative_values$Batch))
representative_values <- representative_values %>% mutate(Batch = factor(as.character(Batch), levels = displayed_batches))

fit_representative_curves <- function(df) {
    tryCatch({
        prediction_grid <- tibble(Fraction = seq(0.001, 1, length.out = 100))
        prediction_grid %>% mutate(Fitted_Raw = predict(lm(Normalized_Linear ~ Fraction, data = df), newdata = prediction_grid),
                                   Fitted_SRR = predict(lm(SRR ~ Fraction, data = df), newdata = prediction_grid))
    }, error = function(e) tibble(Fraction = numeric(), Fitted_Raw = numeric(), Fitted_SRR = numeric()))
}

representative_fits <- representative_values %>% group_by(Platform, Batch, GeneSymbol, UniqueID) %>% nest() %>%
    mutate(CurveData = map(data, fit_representative_curves)) %>% unnest(CurveData) %>% select(-data) %>% ungroup() %>% mutate(Batch = factor(as.character(Batch), levels = displayed_batches))

p_rep_raw <- ggplot() +
    geom_point(data = representative_values, aes(Fraction, Normalized_Linear, color = Batch, shape = Batch), size = 1.25) +
    geom_line(data = representative_fits, aes(Fraction, Fitted_Raw, color = Batch, group = Batch), linewidth = 0.5) +
    facet_wrap(~GeneSymbol, scales = "free_y", nrow = 3) +
    scale_color_manual(values = plot_batch_colors[displayed_batches]) + scale_shape_manual(values = plot_batch_shapes[displayed_batches]) +
    scale_x_reverse(breaks = c(1, 0.75, 0.5, 0.25, 0), labels = c("M", "Y", "P", "X", "F")) +
    scale_y_continuous(labels = scales::label_number(scale = 1e-3), limits = c(0, NA), expand = c(0, 0)) +
    labs(y = expression("Before SRR (Baseline intensity," ~ "\u00D7" ~ 10^3 * ")"), x = NULL, color = "Batch", shape = "Batch") +
    plasmix_theme +
    theme(legend.position = "bottom", legend.box.margin = margin(t = -5, b = 0, r = 0, l = 0),
          legend.text = element_text(size = 7.5, vjust = 0.5, margin = margin(0, 0, 0, 3)),
          strip.text = element_text(size = 7.5, face = "bold", margin = margin(0, 0, 0, 0)), strip.background = element_blank()) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE))

p_rep_srr <- ggplot() +
    geom_point(data = representative_values, aes(Fraction, SRR, color = Batch, shape = Batch), size = 1.25) +
    geom_line(data = representative_fits, aes(Fraction, Fitted_SRR, color = Batch, group = Batch), linewidth = 0.5) +
    facet_wrap(~GeneSymbol, scales = "free_y", nrow = 3) +
    scale_color_manual(values = plot_batch_colors[displayed_batches]) + scale_shape_manual(values = plot_batch_shapes[displayed_batches]) +
    scale_x_reverse(breaks = c(1, 0.75, 0.5, 0.25, 0), labels = c("M", "Y", "P", "X", "F")) +
    scale_y_continuous(n.breaks = 5, limits = c(0, NA), expand = c(0, 0)) +
    labs(y = "After SRR (Ratio to P)", x = NULL, color = "Batch", shape = "Batch") + plasmix_theme +
    theme(legend.position = "none", strip.text = element_text(size = 7.5, face = "bold", margin = margin(0, 0, 0, 0)),
          strip.background = element_blank())

p_representative <- ggarrange(p_rep_raw, p_rep_srr, nrow = 1, widths = c(1, 1), common.legend = TRUE, legend = "bottom")

# 4. Panel c: relative blank burden ----
# Within-batch background burden is evaluated per single-accession assay, without cross-platform matching.
raw_affinity <- single_accession_df %>%
    filter(Batch %in% affinity_batches, Platform %in% c("OLK", "SOM"), ProcessLevel == "Raw", Sample %in% c("P", "BLK")) %>%
    mutate(Linear_Value = 2^Value)
raw_blank <- raw_affinity %>% filter(Sample == "BLK") %>% group_by(Platform, Batch, UniqueID) %>%
    summarize(Median_BLK_Raw = median(Linear_Value, na.rm = TRUE), MAD_BLK_Raw = mad(Linear_Value, na.rm = TRUE), .groups = "drop")
raw_p_detection <- raw_affinity %>% filter(Sample == "P") %>% inner_join(raw_blank, by = c("Platform", "Batch", "UniqueID")) %>%
    group_by(Platform, Batch, UniqueID) %>%
    summarize(P_Replicates = sum(is.finite(Linear_Value)), P_Above_BLK = sum(is.finite(Linear_Value) & Linear_Value > Median_BLK_Raw + 3 * MAD_BLK_Raw),
              Median_P_Raw = median(Linear_Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(P_Detected_Raw = P_Replicates > 0 & P_Above_BLK / P_Replicates > 0.5)
rbb_values <- raw_p_detection %>% inner_join(raw_blank %>% select(Platform, Batch, UniqueID, Median_BLK_Raw), by = c("Platform", "Batch", "UniqueID")) %>%
    filter(P_Detected_Raw, is.finite(Median_P_Raw), Median_P_Raw > 0) %>% mutate(RBB = Median_BLK_Raw / Median_P_Raw) %>%
    filter(is.finite(RBB), RBB >= 0, RBB <= 1)
rbb_medians <- rbb_values %>% group_by(Platform, Batch) %>%
    summarize(Median_RBB = median(RBB, na.rm = TRUE), N_Features = n(), .groups = "drop") %>%
    mutate(Label = paste0(Batch, ": ", scales::percent(Median_RBB, accuracy = 0.1)))

rbb_labels <- rbb_medians %>% arrange(desc(Median_RBB)) %>% mutate(y_pos = seq(6 - 0.8 * (n() - 1), 6, by = 0.8))

p_rbb <- ggplot(rbb_values, aes(RBB, fill = Batch, color = Batch)) +
    geom_density(alpha = 0.2, linewidth = 0.4, show.legend = FALSE) +
    geom_segment(data = rbb_labels, aes(x = Median_RBB, xend = Median_RBB, y = 0, yend = y_pos, color = Batch),
                 linetype = "dotted", linewidth = 0.5, show.legend = FALSE) +
    geom_label(data = rbb_labels, aes(x = Median_RBB, y = y_pos, label = Label, color = Batch), fill = "white", size = 2.5,
               fontface = "bold", hjust = 0, vjust = 0, label.padding = unit(0.4, "lines"), label.size = 0.2, show.legend = FALSE) +
    scale_fill_manual(values = rbb_colors) + scale_color_manual(values = rbb_colors) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 7)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.2), labels = scales::percent) +
    labs(x = "Relative blank burden", y = "Density") + plasmix_theme +
    theme(legend.position = "none", panel.grid.major = element_blank(),
          axis.text.x = element_text(hjust = c(0.1, 0.5, 0.5, 0.5, 0.5, 0.9)))

# 5. Panel d: BLK-subtracted SRR using P or N; Intra-DIA is not evaluated ----
long_baseline <- long_df %>%
    filter(Batch %in% selected_batches, Sample %in% c("M", "Y", "P", "X", "F", "N"),
           (Platform %in% c("OLK", "SOM") & ProcessLevel == "HybNorm") | (Platform == "DIA" & ProcessLevel == "Intensity")) %>%
    mutate(Baseline_Linear = 2^Value) %>% select(Batch, Platform, UniqueID, Sample, ColName, Baseline_Linear)
raw_blank_all <- single_accession_df %>%
    filter(Batch %in% affinity_batches, Platform %in% c("OLK", "SOM"), Sample == "BLK", ProcessLevel == "Raw") %>%
    mutate(Raw_Linear = 2^Value) %>% group_by(Batch, Platform, UniqueID) %>%
    summarize(Median_BLK_Raw = median(Raw_Linear, na.rm = TRUE), .groups = "drop")
subtracted_affinity <- single_accession_df %>%
    filter(Batch %in% affinity_batches, Platform %in% c("OLK", "SOM"), Sample %in% c("M", "Y", "P", "X", "F", "N"),
           ProcessLevel %in% c("Raw", "HybNorm")) %>%
    mutate(Linear_Value = 2^Value) %>% select(Batch, Platform, UniqueID, Sample, ColName, ProcessLevel, Linear_Value) %>%
    pivot_wider(names_from = ProcessLevel, values_from = Linear_Value) %>% rename(Raw_Linear = Raw, Baseline_Linear = HybNorm) %>%
    inner_join(raw_blank_all, by = c("Batch", "Platform", "UniqueID")) %>%
    mutate(Well_Scale = if_else(is.finite(Raw_Linear) & Raw_Linear > 0, Baseline_Linear / Raw_Linear, NA_real_),
           Subtracted_Linear = pmax(Raw_Linear - Median_BLK_Raw, 1e-5) * Well_Scale) %>%
    select(Batch, Platform, UniqueID, Sample, ColName, Subtracted_Linear) %>%
    left_join(feature_metadata %>% distinct(Platform, UniqueID, UniProtID), by = c("Platform", "UniqueID"), relationship = "many-to-one") %>%
    mutate(UniqueID = if_else(Platform == "SOM", UniProtID, UniqueID)) %>%
    group_by(Batch, Platform, UniqueID, Sample, ColName) %>%
    summarize(Subtracted_Linear = exp(mean(log(Subtracted_Linear))), .groups = "drop")
subtracted_dia <- long_baseline %>% filter(Platform == "DIA") %>%
    transmute(Batch, Platform, UniqueID, Sample, ColName, Subtracted_Linear = Baseline_Linear)
subtracted_all <- bind_rows(subtracted_affinity, subtracted_dia)
reference_means <- subtracted_all %>% filter(Sample %in% c("P", "N")) %>% group_by(Batch, UniqueID, Sample) %>%
    summarize(Reference_Mean = if (any(is.finite(Subtracted_Linear))) mean(Subtracted_Linear[is.finite(Subtracted_Linear)]) else NA_real_,
              .groups = "drop") %>% pivot_wider(names_from = Sample, values_from = Reference_Mean, names_prefix = "Reference_")
ratio_replicates <- subtracted_all %>% left_join(reference_means, by = c("Batch", "UniqueID")) %>%
    mutate(P = Subtracted_Linear / Reference_P, N = Subtracted_Linear / Reference_N) %>%
    select(Batch, Platform, UniqueID, Sample, ColName, P, N) %>% pivot_longer(c(P, N), names_to = "Reference", values_to = "Value") %>%
    filter(is.finite(Value), Value > 0)

eligibility_keys <- srr_feature_results %>%
    filter(Design == "Balanced", Detailed_Type != "Intra-DIA", Method %in% c("SRR (P)", "SRR (N)"))
if ("MethodReference" %in% names(eligibility_keys)) {
    eligibility_keys <- eligibility_keys %>% mutate(Reference = parse_reference(Method, MethodReference))
} else {
    eligibility_keys <- eligibility_keys %>% mutate(Reference = parse_reference(Method))
}
eligibility_keys <- eligibility_keys %>% filter(Reference %in% c("P", "N")) %>%
    distinct(Detailed_Type, Macro_Type, Batch1, Batch2, Design, Reference, UniqueID)

calculate_smape <- function(df_pair, batch1, batch2) {
    samples <- c("M", "Y", "X", "F")
    df_pair %>% filter(Sample %in% samples) %>% group_by(UniqueID, Sample) %>%
        summarize(Batch1_Value = if (any(Batch == batch1)) mean(Value[Batch == batch1], na.rm = TRUE) else NA_real_,
                  Batch2_Value = if (any(Batch == batch2)) mean(Value[Batch == batch2], na.rm = TRUE) else NA_real_, .groups = "drop") %>%
        mutate(Denominator = (abs(Batch1_Value) + abs(Batch2_Value)) / 2,
               Sample_sMAPE = if_else(is.finite(Denominator) & Denominator > 0, abs(Batch1_Value - Batch2_Value) / Denominator, NA_real_)) %>%
        group_by(UniqueID) %>%
        summarize(sMAPE = if (sum(is.finite(Sample_sMAPE)) == length(samples)) 100 * mean(Sample_sMAPE) else NA_real_, .groups = "drop")
}

calculate_expected_response <- function(df_pair) {
    calc_feature_titration_metrics(df_pair %>% filter(Sample %in% c("M", "Y", "X", "F")), method = "Mean", is_log2 = FALSE,
                                   trc_deviation_cutoff = trc_cutoff, min_valid_replicates = 2) %>%
        transmute(UniqueID, TitrationMonoValid, TRC_N_finite, MeanTRCDev, TRCDevValid, ExpectedResponseValid)
}

pair_reference_grid <- eligibility_keys %>% distinct(Detailed_Type, Macro_Type, Batch1, Batch2, Design, Reference)
blank_feature_results <- pmap_dfr(pair_reference_grid, function(Detailed_Type, Macro_Type, Batch1, Batch2, Design, Reference) {
    eligible <- eligibility_keys %>%
        filter(.data$Batch1 == .env$Batch1, .data$Batch2 == .env$Batch2, .data$Reference == .env$Reference) %>% pull(UniqueID)
    df_pair <- ratio_replicates %>%
        filter(Batch %in% c(.env$Batch1, .env$Batch2), .data$Reference == .env$Reference, UniqueID %in% eligible) %>%
        group_by(Batch, UniqueID, Sample) %>% slice_head(n = 3) %>% ungroup()
    if (!nrow(df_pair)) return(tibble())
    calculate_smape(df_pair, Batch1, Batch2) %>% full_join(calculate_expected_response(df_pair), by = "UniqueID") %>%
        mutate(QuantAgreementValid = is.finite(sMAPE) & sMAPE <= smape_cutoff,
               OverallSuccess = QuantAgreementValid & coalesce(ExpectedResponseValid, FALSE), Detailed_Type = .env$Detailed_Type,
               Macro_Type = .env$Macro_Type, Batch1 = .env$Batch1, Batch2 = .env$Batch2, Design = .env$Design,
               Reference = .env$Reference, Calculation = "BLK-subtracted SRR")
})
blank_pair_values <- blank_feature_results %>% group_by(Detailed_Type, Macro_Type, Batch1, Batch2, Design, Reference) %>%
    summarize(N_Features = n(), Rate_sMAPE = mean(coalesce(QuantAgreementValid, FALSE)),
              Rate_Response = mean(coalesce(ExpectedResponseValid, FALSE)), Rate_Harmonized = mean(coalesce(OverallSuccess, FALSE)), .groups = "drop")
if (!nrow(blank_pair_values)) stop("No BLK-subtracted SRR results were generated.")
blank_plot <- make_reference_plot(blank_pair_values, show_legend = FALSE)
p_blank_srr <- blank_plot$plot + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = c("0", "25", "50", "75", "100"), expand = c(0, 0))

# 6. Assemble and export ----
# Panels a and b occupy the same row and therefore have exactly the same outer height.
top_row <- ggarrange(p_standard_srr, p_representative, nrow = 1, widths = c(1.65, 1), labels = c("a", "b"),
                             label.x = 0, label.y = 1, hjust = -0.2, vjust = 1.2, font.label = label_style)
bottom_row <- ggarrange(p_rbb, p_blank_srr, nrow = 1, widths = c(1, 1.35), labels = c("c", "d"),
                                label.x = 0, label.y = 1, hjust = -0.2, vjust = 1.2, font.label = label_style)
figure_ed8 <- ggarrange(top_row, bottom_row, ncol = 1, heights = c(1, 0.9))
figure_pdf <- file.path("figures", "ed8_reference_background_integration.pdf")
figure_png <- file.path("figures", "ed8_reference_background_integration.png")
ggsave(figure_pdf, figure_ed8, width = 10, height = 6.27)
ggsave(figure_png, figure_ed8, width = 10, height = 6.27, dpi = 600, bg = "white")

source_data <- list(
    ED8a_standard_pair_values = standard_pair_values %>% arrange(Detailed_Type, Batch1, Batch2, Reference),
    ED8a_standard_medians = standard_plot$summary %>% mutate(across(where(is.factor), as.character)),
    ED8b_representative_values = representative_values %>% mutate(across(where(is.factor), as.character)) %>%
        select(Platform, Batch, GeneSymbol, UniqueID, Sample, ColName, Fraction, Normalized_Linear, P_Reference, SRR, R2_Lin),
    ED8b_linear_fit_quality = linear_fit_quality %>% mutate(across(where(is.factor), as.character)),
    ED8b_fitted_curves = representative_fits %>% mutate(across(where(is.factor), as.character)),
    ED8c_relative_blank_burden = rbb_values %>% arrange(Platform, Batch, UniqueID),
    ED8c_RBB_medians = rbb_medians %>% arrange(Platform, Batch),
    ED8d_BLK_subtracted_features = blank_feature_results %>% arrange(Detailed_Type, Batch1, Batch2, Reference, UniqueID),
    ED8d_BLK_subtracted_pairs = blank_pair_values %>% arrange(Detailed_Type, Batch1, Batch2, Reference),
    ED8d_BLK_subtracted_medians = blank_plot$summary %>% mutate(across(where(is.factor), as.character))
)
source_data_file <- file.path("tables", "SourceData_EDFigure8.xlsx")
write.xlsx(source_data, source_data_file, overwrite = TRUE, keepNA = TRUE, na.string = "NA")
message("Extended Data Figure 8 exported to: ", figure_pdf)
message("Extended Data Figure 8 source data exported to: ", source_data_file)
