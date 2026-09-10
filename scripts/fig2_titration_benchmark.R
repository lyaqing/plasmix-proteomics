# Figure 2 | Titration response, signal-to-noise ratio and replicate CV

# 0. Setup ----
source("scripts/_project_setup.R")
use_packages(c("data.table", "tidyverse", "readxl", "ggpubr", "ggpp", "openxlsx", "showtext", "RColorBrewer"), "matrixStats")
source("utils/benchmark_metrics.R")
source("utils/imputation.R")
source("utils/figure_style.R")
source("utils/feature_mapping.R")
plasmix_theme <- if (is.function(theme_plasmix)) theme_plasmix() else theme_plasmix
label_style <- list(size = 12, face = "bold")
set.seed(2026)
trc_deviation_cutoff <- 0.25

# 1. Inputs ----
paths <- c(metadata = "data/study_metadata.xlsx", feature_metadata = "data/feature_metadata.tsv.gz",
           profiles = "data/protein_profiles_long.tsv.gz", detection = "results/analyte_detection_status.tsv.gz")
missing_inputs <- paths[!file.exists(paths)]
if (length(missing_inputs) > 0) stop("Missing input files: ", paste(missing_inputs, collapse = ", "), call. = FALSE)

meta_sample <- read_xlsx(paths["metadata"], sheet = "sample") %>%
    filter(Sample %in% c("M", "Y", "P", "X", "F", "N"))

lod_status <- fread(paths["detection"]) %>% mutate(Is_Detected = M | Y | P | X | F)
feature_metadata <- analysis_feature_metadata(fread(paths["feature_metadata"]))

long_df <- aggregate_som_profiles(fread(paths["profiles"]))
long_df <- filter_batch_analysis_features(long_df, feature_metadata, strict_platforms = character())
lod_status <- lod_status %>% semi_join(distinct(long_df, Platform, Batch, UniqueID), by = c("Platform", "Batch", "UniqueID"))
stopifnot(nrow(anti_join(distinct(long_df, Platform, Batch, UniqueID), lod_status, by = c("Platform", "Batch", "UniqueID"))) == 0)
long_df_filter <- long_df %>%
    filter(Sample %in% c("M", "Y", "P", "X", "F", "N"),
           DataTier %in% c("Baseline", "Calibrated", "Reshaped"))

target_samples <- c("M", "F", "P", "X", "Y")
task_combinations <- long_df_filter %>% distinct(Platform, Batch, ProcessLevel, DataTier)

filter_primary_process <- function(data) {
    data %>% filter(
        (Platform == "AAG" & DataTier == "Baseline" & ProcessLevel == "SNR") |
        (Platform == "DIA" & DataTier == "Baseline" & ProcessLevel == "Intensity") |
        (Platform == "SOM" & DataTier == "Calibrated" & ProcessLevel == "Calibrate") |
        # (Platform == "SOM" & DataTier == "Reshaped" & ProcessLevel %in% c("ANML-SMP", "MedNormExt")) |
        (Platform == "NLS" & DataTier == "Calibrated" & grepl("NPQ", ProcessLevel, ignore.case = TRUE)) |
        (Platform == "OLK" & DataTier == "Calibrated" & grepl("NPX", ProcessLevel, ignore.case = TRUE))
    )
}

# 2. Feature- and batch-level titration metrics ----
titration_iteration_results <- vector("list", nrow(task_combinations))
for (index in seq_len(nrow(task_combinations))) {
    task <- task_combinations[index, ]
    subset_data <- long_df_filter %>%
        filter(Batch == task$Batch, ProcessLevel == task$ProcessLevel, DataTier == task$DataTier)
    if (n_distinct(subset_data$UniqueID) < 3) next

    sample_metadata <- meta_sample %>%
        filter(Batch == task$Batch, Sample %in% target_samples)
    available_columns <- subset_data %>%
        group_by(ColName) %>%
        summarize(Available = any(is.finite(Value)), .groups = "drop") %>%
        filter(Available) %>%
        pull(ColName)
    if (nrow(sample_metadata) == 0 || length(available_columns) == 0) next

    replicate_plan <- make_replicate_plan(
        sample_metadata, available_columns, n_replicates = 3,
        required_samples = c("M", "F")
    )
    if (length(replicate_plan) == 0) {
        warning(sprintf("Could not generate a three-replicate plan: %s / %s / %s", task$Batch,
                        task$ProcessLevel, task$DataTier), call. = FALSE)
        next
    }
    iteration_results <- vector("list", length(replicate_plan))

    for (iteration in seq_along(replicate_plan)) {
        selected_metadata <- replicate_plan[[iteration]]$Metadata
        sampled_columns <- selected_metadata$ColName
        iteration_data <- subset_data %>% filter(ColName %in% sampled_columns)

        result <- tryCatch(
            calc_feature_titration_metrics(
                iteration_data, is_log2 = TRUE, method = "Mean",
                trc_deviation_cutoff = trc_deviation_cutoff,
                min_valid_replicates = 2
            ),
            error = function(error) {
                warning(sprintf("Titration metric calculation failed: %s / %s / %s / combination %s；%s", task$Batch,
                                task$ProcessLevel, task$DataTier, iteration, conditionMessage(error)), call. = FALSE)
                NULL
            }
        )
        if (is.null(result) || nrow(result) == 0) next

        iteration_results[[iteration]] <- result %>%
            mutate(Batch = task$Batch, Platform = task$Platform,
                   ProcessLevel = task$ProcessLevel, DataTier = task$DataTier,
                   Iteration = iteration,
                   Replicate_set = paste(replicate_plan[[iteration]]$TargetReplicates, collapse = ","),
                   Substitution_map = replicate_plan[[iteration]]$SubstitutionMap)
    }
    titration_iteration_results[[index]] <- bind_rows(iteration_results)
}

safe_mean <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) NA_real_ else mean(x)
}

trc_iteration_level <- bind_rows(titration_iteration_results) %>%
    left_join(lod_status %>% select(Batch, UniqueID, Is_Detected), by = c("Batch", "UniqueID")) %>%
    mutate(Is_Detected = replace_na(Is_Detected, FALSE))
if (nrow(trc_iteration_level) == 0) stop("No feature-level titration results were generated.", call. = FALSE)

titration_summary <- summarize_titration_subsamples(
    trc_iteration_level,
    group_cols = c("Batch", "Platform", "ProcessLevel", "DataTier"),
    majority_cutoff = 0.5,
    trc_deviation_cutoff = trc_deviation_cutoff
)
trc_feature_level <- titration_summary$feature %>%
    left_join(lod_status %>% select(Batch, UniqueID, Is_Detected), by = c("Batch", "UniqueID")) %>%
    mutate(Is_Detected = replace_na(Is_Detected, FALSE))
trc_relation_level <- titration_summary$relation %>%
    left_join(lod_status %>% select(Batch, UniqueID, Is_Detected), by = c("Batch", "UniqueID")) %>%
    mutate(Is_Detected = replace_na(Is_Detected, FALSE))

gradient_fit_iteration <- calc_gradient_fit(
    trc_iteration_level,
    group_cols = c("Batch", "Platform", "ProcessLevel", "DataTier", "Iteration")
)
gradient_fit_level <- gradient_fit_iteration %>%
    group_by(Batch, Platform, ProcessLevel, DataTier, Gradient) %>%
    summarize(Nominal_TRC = first(Nominal_TRC),
              All_R2_SD = if (sum(is.finite(All_R2)) > 1) sd(All_R2, na.rm = TRUE) else NA_real_,
              Detected_R2_SD = if (sum(is.finite(Detected_R2)) > 1) sd(Detected_R2, na.rm = TRUE) else NA_real_,
              All_R2 = safe_mean(All_R2), Detected_R2 = safe_mean(Detected_R2),
              All_N = round(safe_mean(All_N)), Detected_N = round(safe_mean(Detected_N)),
              Iterations = n_distinct(Iteration), .groups = "drop")

replicate_plan_audit <- trc_iteration_level %>%
    distinct(Platform, Batch, ProcessLevel, DataTier, Iteration, Replicate_set, Substitution_map)

fwrite(as.data.table(trc_iteration_level), "results/trc_iteration_level.tsv.gz", sep = "\t", na = "NA")
fwrite(as.data.table(trc_feature_level), "results/trc_feature_level.tsv.gz", sep = "\t", na = "NA")
fwrite(as.data.table(trc_relation_level), "results/trc_relation_level.tsv.gz", sep = "\t", na = "NA")
fwrite(as.data.table(replicate_plan_audit), "results/replicate_plan_audit.tsv.gz", sep = "\t", na = "NA")
fwrite(as.data.table(gradient_fit_level), "results/gradient_fit_level.tsv.gz", sep = "\t", na = "NA")

# 3. Panel a: Effect size and titration fidelity ----
default_trc <- trc_feature_level %>% filter_primary_process()
platform_levels <- c("DIA", "SOM", "OLK", "NLS", "AAG")
default_trc <- trc_feature_level %>%
    filter_primary_process() %>%
    mutate(MF_log2 = M - F, Absolute_MF = abs(MF_log2), Platform = factor(Platform, levels = platform_levels))

trc_scope <- bind_rows(
    default_trc %>% mutate(TargetSpace = "All"),
    default_trc %>% filter(Is_Detected %in% TRUE) %>% mutate(TargetSpace = "Detected")
) %>%
    mutate(TargetSpace = factor(TargetSpace, levels = c("All", "Detected")))

bin_width <- 0.2
titration_by_effect <- trc_scope %>%
    filter(is.finite(MF_log2), !is.na(ExpectedResponseValid)) %>%
    mutate(MF_bin = floor(MF_log2 / bin_width) * bin_width + bin_width / 2) %>%
    group_by(Platform, TargetSpace, MF_bin) %>%
    summarize(Titration_fidelity = 100 * mean(ExpectedResponseValid, na.rm = TRUE), Features = n(), .groups = "drop") %>%
    filter(Features >= 2, is.finite(MF_bin), is.finite(Titration_fidelity)) %>%
    mutate(Platform = factor(Platform, levels = platform_levels))

platform_pass_n <- default_trc %>%
    group_by(Platform) %>%
    summarize(N_pass = n_distinct(UniqueID[ExpectedResponseValid %in% TRUE]), .groups = "drop") %>%
    mutate(Platform = factor(Platform, levels = platform_levels),
          Label = paste0("plain('Pass')~italic(n) == '", scales::comma(N_pass), "'"))

integer_breaks <- function(n = 5) {
    function(x) { brks <- scales::pretty_breaks(n = n)(x); brks[brks == floor(brks)] }
}
nls_axis_anchor <- tibble(Platform = factor("NLS", levels = platform_levels), MF_bin = c(-1, 1), Titration_fidelity = c(0, 0))

p_u <- ggplot(titration_by_effect, aes(x = MF_bin, y = Titration_fidelity, color = Platform)) +
    geom_blank(data = nls_axis_anchor, aes(x = MF_bin, y = Titration_fidelity), inherit.aes = FALSE) +
    geom_point(data = filter(titration_by_effect, TargetSpace == "All"), aes(size = Features, fill = Platform), shape = 21, color = "grey70", stroke = 0.5, alpha = 0.8) +
    geom_point(data = filter(titration_by_effect, TargetSpace == "Detected"), aes(size = Features, fill = Platform), shape = 21, color = "black", stroke = 0.5, alpha = 0.2) +
    geom_smooth(data = filter(titration_by_effect, TargetSpace == "All"), aes(linetype = TargetSpace, alpha = TargetSpace), method = "loess", formula = y ~ x, span = 0.8, method.args = list(degree = 2), se = FALSE, linewidth = 0.5) +
    geom_smooth(data = filter(titration_by_effect, TargetSpace == "Detected"), aes(linetype = TargetSpace, alpha = TargetSpace), method = "loess", formula = y ~ x, span = 0.8, method.args = list(degree = 2), se = FALSE, linewidth = 0.5) +
    geom_text_npc(data = platform_pass_n, aes(npcx = 0.5, npcy = 0.96, label = Label), inherit.aes = FALSE, parse = TRUE, family = "serif", size = 3) +
    facet_wrap(~Platform, nrow = 1, scales = "free_x") +
    scale_x_continuous(breaks = integer_breaks(n = 4)) +
    scale_y_continuous(limits = c(-5, 118), breaks = seq(0, 100, 25), expand = c(0, 0)) +
    scale_color_manual(values = platform_color, limits = platform_levels, name = "Platform") +
    scale_fill_manual(values = platform_color, limits = platform_levels, guide = "none") +
    scale_size_continuous(range = c(1.5, 6), breaks = c(10, 50, 100, 1000), name = "Analytes per bin") +
    scale_linetype_manual(values = c("All" = "dotted", "Detected" = "solid"), name = "Analyte scope") +
    scale_alpha_manual(values = c("All" = 0.65, "Detected" = 1), guide = "none") +
    labs(x = expression(log[2]*"(M/F)"), y = "Expected-response analytes (%)") +
    plasmix_theme +
    theme(panel.grid.major.x = element_blank(), legend.position = "bottom",
          legend.box.margin = margin(-8, 0, 0, 0), panel.spacing.x = unit(0.25, "lines")) +
    guides(color = guide_legend(order = 1, override.aes = list(linetype = "solid", alpha = 1, linewidth = 0.8)),
           size = guide_legend(order = 2),
           linetype = guide_legend(order = 3, override.aes = list(color = "black", alpha = 1, linewidth = 0.6)))

# 4. Panel b: Cumulative fidelity by effect size ----
batch_ranked <- trc_scope %>%
    group_by(Batch, Platform, TargetSpace) %>%
    arrange(desc(Absolute_MF), .by_group = TRUE) %>%
    mutate(Rank = row_number(), Cumulative_fidelity = 100 * cummean(replace_na(ExpectedResponseValid, FALSE))) %>%
    ungroup()

cumulative_batch <- batch_ranked %>%
    group_by(Platform, TargetSpace) %>%
    mutate(Platform_max_rank = max(Rank, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(Batch, Platform, TargetSpace) %>%
    complete(Rank = seq_len(max(Platform_max_rank, na.rm = TRUE))) %>%
    fill(Cumulative_fidelity, .direction = "down") %>%
    ungroup()

cumulative_fidelity <- cumulative_batch %>%
    group_by(Platform, TargetSpace, Rank) %>%
    summarize(Mean = mean(Cumulative_fidelity, na.rm = TRUE),
              Q1 = quantile(Cumulative_fidelity, 0.25, na.rm = TRUE),
              Q3 = quantile(Cumulative_fidelity, 0.75, na.rm = TRUE),
              Batches = n_distinct(Batch), .groups = "drop")
topk_breaks <- 10^(0:floor(log10(max(cumulative_fidelity$Rank, na.rm = TRUE))))

p_topk <- ggplot(cumulative_fidelity,
                 aes(x = Rank, y = Mean, color = Platform, linetype = TargetSpace,
                     alpha = TargetSpace, group = interaction(Platform, TargetSpace))) +
    geom_ribbon(data = filter(cumulative_fidelity, TargetSpace == "Detected"),
                aes(x = Rank, ymin = Q1, ymax = Q3, fill = Platform, group = Platform),
                inherit.aes = FALSE, alpha = 0.12, color = NA, show.legend = FALSE) +
    geom_line(linewidth = 0.55) +
    scale_color_manual(values = platform_color, limits = platform_levels, name = "Platform") +
    scale_fill_manual(values = platform_color, limits = platform_levels) +
    scale_linetype_manual(values = c("All" = "dotted", "Detected" = "solid"), name = "Analytes") +
    scale_alpha_manual(values = c("All" = 0.65, "Detected" = 1), guide = "none") +
    plasmix_theme +
    scale_x_log10(limits = c(1, NA), breaks = topk_breaks, labels = scales::comma, expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0,0)) +
    labs(x = expression(Top~italic(k)~analytes), y = "Cumulative expected-response rate (%)") +
    theme(panel.grid.major.x = element_blank(), legend.position = "bottom",
          axis.title.y = element_text(hjust = 1.5), axis.text.x = element_text(hjust = c(rep(0.5, length(topk_breaks) - 1L), 0.95))) +
    guides(color = guide_legend(order = 1, override.aes = list(linetype = "solid", alpha = 1, linewidth = 0.7)),
          linetype = guide_legend(order = 3, override.aes = list(color = "black", alpha = 1, linewidth = 0.6)))

# 5. Panel e: Batch-level gradient fit ----
global_plot_data <- gradient_fit_level %>%
    filter_primary_process() %>%
    filter(is.finite(Detected_R2)) %>%
    group_by(Batch) %>% mutate(Average_R2 = mean(Detected_R2, na.rm = TRUE)) %>% ungroup() %>%
    mutate(Gradient = factor(Gradient, levels = c("Y", "P", "X")))

p_r2 <- ggplot(global_plot_data, aes(x = Gradient, y = reorder(Batch, Average_R2), fill = Detected_R2)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.2f", Detected_R2), color = Detected_R2 > 0.60), size = 2.8, angle = 0) +
    scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    scale_fill_gradientn(colours = brewer.pal(9, "YlGnBu"), limits = c(0, 1), oob = scales::squish, name = expression(bolditalic(R)^bold(2))) +
    guides(fill = guide_colorbar(title.position = "left", title.vjust = 1, barwidth = unit(3, "cm"), barheight = unit(0.25, "cm"))) +
    labs(x = "Mixture gradient", y = NULL) +
    plasmix_theme +
    theme(axis.line.x = element_blank(), axis.line.y = element_blank(),
          axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
          legend.position = "bottom", legend.title = element_text(face = "bold", margin = margin(r = 5)),
          legend.box.spacing = unit(1, "pt"), legend.margin = margin(t = 3, l = -10, b = 0),
          legend.text = element_text(margin = margin(t = 2, r = 0, b = 0, l = 0)),
          panel.grid.major = element_blank(), panel.border = element_blank(),
          plot.background = element_blank(), panel.background = element_blank())

# 6. Panel f: Expected-response feature yield ----
calc_titration_yield <- function(data) {
    data %>%
        group_by(Platform, Batch) %>%
        summarize(Detected = n_distinct(UniqueID[Is_Detected %in% TRUE]),
                  Valid_detected = n_distinct(UniqueID[Is_Detected %in% TRUE & ExpectedResponseValid %in% TRUE]),
                  Valid_lod_excluded = n_distinct(UniqueID[Is_Detected %in% FALSE & ExpectedResponseValid %in% TRUE]), .groups = "drop") %>%
        mutate(Denominator = Detected + Valid_lod_excluded, Valid_total = Valid_detected + Valid_lod_excluded,
               Detected_valid_pct = 100 * Valid_detected / Denominator,
               Lod_excluded_valid_pct = 100 * Valid_lod_excluded / Denominator,
               Fidelity = 100 * Valid_total / Denominator)
}

yield_primary <- trc_feature_level %>%
    filter_primary_process() %>%
    calc_titration_yield() %>%
    arrange(Fidelity) %>%
    mutate(Batch = factor(Batch, levels = Batch), Primary_label = if_else(Platform == "SOM", paste0("italic(n)[C] == '", scales::comma(Valid_total), "'"), paste0("italic(n) == '", scales::comma(Valid_total), "'")))

yield_reshaped <- trc_feature_level %>%
    filter(Platform == "SOM", DataTier == "Reshaped",
           (Batch %in% c("SOM_P1_B1", "SOM_P1_B2") & ProcessLevel == "ANML-SMP") |
           (Batch %in% c("SOM_P2_B1", "SOM_P2_B2") & ProcessLevel == "MedNormExt")) %>%
    calc_titration_yield() %>%
    left_join(yield_primary %>% select(Batch, Primary_fidelity = Fidelity), by = "Batch") %>%
    mutate(Batch = factor(Batch, levels = levels(yield_primary$Batch)),
           Reshaped_label = paste0("italic(n)[R] == '", scales::comma(Valid_total), "'"))

yield_segments <- yield_primary %>%
    select(Platform, Batch, Detected_valid_pct, Lod_excluded_valid_pct) %>%
    pivot_longer(c(Detected_valid_pct, Lod_excluded_valid_pct), names_to = "Detection_status", values_to = "Percentage") %>%
    mutate(Detection_status = factor(Detection_status, levels = c("Detected_valid_pct", "Lod_excluded_valid_pct"), labels = c("Detected", "LoD-excluded")))

p_pass_num <- ggplot() +
    geom_col(data = yield_segments, aes(x = Percentage, y = Batch, fill = Detection_status), width = 0.70,
             color = "white", linewidth = 0.1, position = position_stack(reverse = TRUE)) +
    geom_segment(data = yield_reshaped, aes(x = Primary_fidelity, xend = Fidelity, y = Batch, yend = Batch, color = Platform), linewidth = 0.4) +
    geom_point(data = yield_reshaped, aes(x = Fidelity, y = Batch, color = Platform, shape = "Reshaped"), size = 2.2, fill = "white", stroke = 0.65) +
    geom_text(data = yield_reshaped, aes(x = Fidelity, y = Batch, label = Reshaped_label, color = Platform),
              hjust = -0.18, size = 2.8, parse = TRUE) +
    geom_tile(data = yield_primary, aes(x = 90, y = Batch, color = Platform), width = 34, height = 0.75, fill = "grey98", linewidth = 0.45) +
    geom_text(data = yield_primary, aes(x = 90, y = Batch, label = Primary_label, color = Platform), size = 2.8, parse = TRUE) +
    scale_color_manual(values = platform_color, guide = "none") +
    scale_fill_manual(values = c("Detected" = "#397DB7", "LoD-excluded" = "#A9D5EA"), name = NULL) +
    scale_shape_manual(values = c(Reshaped = 23), name = NULL) +
    scale_x_continuous(breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(xlim = c(0, 119), clip = "off") +
    labs(x = "Detection-adjusted response fraction (%)", y = NULL) +
    plasmix_theme +
    theme(panel.grid.major.y = element_blank(), legend.position = "bottom",
          legend.margin = margin(-5, 10, 0, -20),
          axis.title.x = element_text(hjust = 1),
          plot.margin = margin(5, 5, 5, 5)) +
    guides(fill = guide_legend(order = 1), shape = guide_legend(order = 2, override.aes = list(fill = "white")))

# 7. Panel c: PCA-based signal-to-noise ratio ----
snr_results <- list()
for (index in seq_len(nrow(task_combinations))) {
    task <- task_combinations[index, ]
    subset_data <- long_df_filter %>%
        filter(Batch == task$Batch, ProcessLevel == task$ProcessLevel, DataTier == task$DataTier)
    if (nrow(subset_data) == 0) next

    expression_matrix <- subset_data %>%
        select(UniqueID, ColName, Value) %>%
        pivot_wider(names_from = ColName, values_from = Value) %>%
        column_to_rownames("UniqueID") %>% as.matrix()
    sample_metadata <- meta_sample %>%
        filter(Batch == task$Batch, Sample %in% target_samples)
    available_columns <- colnames(expression_matrix)[colSums(is.finite(expression_matrix)) > 0]
    replicate_plan <- make_replicate_plan(sample_metadata, available_columns, n_replicates = 3)
    if (length(replicate_plan) == 0) {
        warning(sprintf("SNRCould not generate a three-replicate plan: %s / %s / %s", task$Batch,
                        task$ProcessLevel, task$DataTier), call. = FALSE)
        next
    }
    valid_columns <- unique(unlist(lapply(replicate_plan, function(plan) plan$Metadata$ColName)))
    available_metadata <- bind_rows(lapply(replicate_plan, function(plan) plan$Metadata)) %>%
        distinct(Sample, ColName, .keep_all = TRUE)
    iterations <- length(replicate_plan)

    detected_features <- lod_status %>%
        filter(Batch == task$Batch, Is_Detected %in% TRUE) %>% pull(UniqueID)
    feature_sets <- list("All features" = rownames(expression_matrix),
                         "Detected features" = intersect(detected_features, rownames(expression_matrix)))

    for (subset_name in names(feature_sets)) {
        common_features <- feature_sets[[subset_name]]
        if (length(common_features) < 3) next
        task_matrix <- expression_matrix[common_features, valid_columns, drop = FALSE]
        task_matrix <- as.matrix(impute_lod_noise(task_matrix, available_metadata, seed = 42))

        for (scale_parameter in c(TRUE, FALSE)) {
            iteration_values <- numeric()
            retained_features <- integer()
            for (iteration in seq_len(iterations)) {
                selected_metadata <- replicate_plan[[iteration]]$Metadata
                sampled_columns <- selected_metadata$ColName
                sampled_groups <- selected_metadata$Sample
                if (length(unique(sampled_groups)) < 2) next

                iteration_matrix <- task_matrix[, sampled_columns, drop = FALSE]
                variable <- matrixStats::rowVars(iteration_matrix, na.rm = TRUE) > 1e-9
                iteration_matrix <- iteration_matrix[variable, , drop = FALSE]
                if (nrow(iteration_matrix) < 3) next

                value <- tryCatch(calc_pca_snr(iteration_matrix, sampled_groups, scale. = scale_parameter),
                                  error = function(error) NA_real_)
                if (is.finite(value)) {
                    iteration_values <- c(iteration_values, value)
                    retained_features <- c(retained_features, nrow(iteration_matrix))
                }
            }
            if (length(iteration_values) == 0) next

            snr_results[[length(snr_results) + 1]] <- tibble(
                Platform = task$Platform, Batch = task$Batch, ProcessLevel = task$ProcessLevel,
                DataTier = task$DataTier, Subset = subset_name, Scale_Param = scale_parameter,
                SNR_mean = mean(iteration_values),
                SNR_SD = if (length(iteration_values) > 1) sd(iteration_values) else NA_real_,
                Assays = round(mean(retained_features)), Iterations = length(iteration_values)
            )
        }
    }
}

snr_all_tiers <- bind_rows(snr_results) %>%
    filter((Platform == "DIA" & Scale_Param) | (Platform %in% c("SOM", "NLS", "AAG", "OLK") & !Scale_Param))

snr_main <- snr_all_tiers %>%
    filter(Subset == "Detected features") %>%
    group_by(Platform, Batch, DataTier) %>%
    slice_max(SNR_mean, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(DataTier = factor(DataTier, levels = c("Baseline", "Calibrated", "Reshaped")))

snr_order <- snr_main %>%
    group_by(Batch) %>%
    summarize(Max_SNR = max(SNR_mean, na.rm = TRUE), .groups = "drop") %>%
    arrange(Max_SNR) %>%
    pull(Batch)

snr_main <- snr_main %>%
    mutate(Batch = factor(Batch, levels = snr_order),
           Lower = SNR_mean - SNR_SD,
           Upper = SNR_mean + SNR_SD)

snr_error <- snr_main %>%
    group_by(Batch) %>%
    slice_max(SNR_mean, n = 1, with_ties = FALSE) %>%
    ungroup()

p_snr <- ggplot(snr_main, aes(x = SNR_mean, y = Batch)) +
    geom_col(aes(fill = Platform, alpha = DataTier), position = "identity", width = 0.7, color = "white", linewidth = 0.1, orientation = "y") +
    geom_errorbar(aes(xmin = SNR_mean, xmax = SNR_mean), width = 0.7,
                  color = "black", linewidth = 0.3, orientation = "y") +
    geom_errorbar(data = snr_error, aes(xmin = Lower, xmax = Upper), width = 0.25,
                  linewidth = 0.3, orientation = "y", na.rm = TRUE) +
    scale_fill_manual(values = platform_color, guide = "none") +
    scale_alpha_manual(values = c(Baseline = 1, Calibrated = 0.7, Reshaped = 0.3), breaks = c("Baseline", "Calibrated", "Reshaped"), name = "Data tier") +
    labs(x = "Signal-to-noise ratio", y = NULL) +
    plasmix_theme +
    scale_x_continuous(limits = c(0, 26), breaks = seq(0, 25, 5), expand = expansion(mult = c(0, 0))) +
    theme(panel.grid.major.y = element_blank(), legend.position = c(1, 0), legend.justification = c(1, 0), legend.box = "horizontal",
          legend.key.size = unit(0.3, "cm"), legend.title = element_text(margin = margin(b = 3)),
          legend.margin = margin(5, 0, 5, 5))

# 8. Panel d: Technical-replicate CV ----
cv_results <- list()
for (index in seq_len(nrow(task_combinations))) {
    task <- task_combinations[index, ]
    subset_data <- long_df_filter %>%
        filter(Batch == task$Batch, ProcessLevel == task$ProcessLevel, DataTier == task$DataTier)
    if (nrow(subset_data) == 0) next

    expression_matrix <- subset_data %>%
        select(UniqueID, ColName, Value) %>%
        pivot_wider(names_from = ColName, values_from = Value) %>%
        column_to_rownames("UniqueID") %>% as.matrix()
    expression_linear <- 2^expression_matrix
    sample_metadata <- meta_sample %>%
        filter(Batch == task$Batch, Sample %in% target_samples)
    available_columns <- colnames(expression_linear)[colSums(is.finite(expression_linear)) > 0]
    replicate_plan <- make_replicate_plan(sample_metadata, available_columns, n_replicates = 3)
    if (length(replicate_plan) == 0) {
        warning(sprintf("CVCould not generate a three-replicate plan: %s / %s / %s", task$Batch,
                        task$ProcessLevel, task$DataTier), call. = FALSE)
        next
    }
    valid_columns <- unique(unlist(lapply(replicate_plan, function(plan) plan$Metadata$ColName)))
    iterations <- length(replicate_plan)

    detected_features <- lod_status %>%
        filter(Batch == task$Batch, Is_Detected %in% TRUE) %>% pull(UniqueID)
    feature_sets <- list("All features" = rownames(expression_linear),
                         "Detected features" = intersect(detected_features, rownames(expression_linear)))

    for (subset_name in names(feature_sets)) {
        common_features <- feature_sets[[subset_name]]
        if (length(common_features) < 3) next
        task_matrix <- expression_linear[common_features, valid_columns, drop = FALSE]
        iteration_cv <- matrix(NA_real_, nrow = nrow(task_matrix), ncol = iterations)

        for (iteration in seq_len(iterations)) {
            selected_metadata <- replicate_plan[[iteration]]$Metadata
            sampled_groups <- split(selected_metadata$ColName, selected_metadata$Sample)
            sampled_groups <- sampled_groups[lengths(sampled_groups) >= 2]
            if (length(sampled_groups) == 0) next

            group_cv <- vapply(sampled_groups, function(columns) {
                values <- task_matrix[, columns, drop = FALSE]
                valid_n <- rowSums(is.finite(values))
                means <- matrixStats::rowMeans2(values, na.rm = TRUE)
                cvs <- matrixStats::rowSds(values, na.rm = TRUE) / means
                cvs[valid_n < 3 | !is.finite(cvs)] <- NA_real_
                cvs
            }, numeric(nrow(task_matrix)))
            if (is.null(dim(group_cv))) group_cv <- matrix(group_cv, ncol = 1)
            iteration_cv[, iteration] <- rowMeans(group_cv, na.rm = TRUE)
        }

        mean_cv <- rowMeans(iteration_cv, na.rm = TRUE)
        mean_cv[!is.finite(mean_cv)] <- NA_real_
        cv_results[[length(cv_results) + 1]] <- tibble(
            Platform = task$Platform, Batch = task$Batch, ProcessLevel = task$ProcessLevel,
            DataTier = task$DataTier, Subset = subset_name, UniqueID = common_features, CV = mean_cv
        )
    }
}

cv_all_tiers <- bind_rows(cv_results)
fwrite(as.data.table(cv_all_tiers), "results/cv_feature_level.tsv.gz", sep = "\t", na = "NA")

cv_main <- cv_all_tiers %>% filter(Subset == "Detected features") %>% filter_primary_process()
cv_sort <- cv_main %>% group_by(Batch, Platform) %>% summarize(Default_CV = median(CV, na.rm = TRUE), .groups = "drop")
batch_order <- cv_sort %>% arrange(desc(Default_CV)) %>% pull(Batch)

cv_main <- cv_main %>%
    left_join(cv_sort, by = c("Batch", "Platform")) %>%
    mutate(Batch = factor(Batch, levels = batch_order))

cv_unfiltered <- cv_all_tiers %>%
    filter(Subset == "All features", Platform %in% c("AAG", "OLK"),
           !Batch %in% c("OLK_P1_B1", "OLK_P1_B2")) %>%
    filter_primary_process() %>%
    group_by(Batch, Platform) %>%
    summarize(Marker_CV = median(CV, na.rm = TRUE), .groups = "drop") %>%
    mutate(Processing = "Unfiltered")

cv_reshaped <- cv_all_tiers %>%
    filter(Subset == "Detected features", Platform == "SOM", DataTier == "Reshaped") %>%
    group_by(Batch, Platform) %>%
    summarize(Marker_CV = median(CV, na.rm = TRUE), .groups = "drop") %>%
    mutate(Processing = "Reshaped")

cv_markers <- bind_rows(cv_unfiltered, cv_reshaped) %>%
    left_join(cv_sort, by = c("Batch", "Platform")) %>%
    filter(is.finite(Marker_CV), is.finite(Default_CV)) %>%
    mutate(Batch = factor(Batch, levels = rev(batch_order)),
           Start_CV = if_else(Processing == "Unfiltered", Marker_CV, Default_CV),
           End_CV = if_else(Processing == "Unfiltered", Default_CV, Marker_CV))

p_cv <- ggplot(cv_main, aes(x = CV * 100, y = Batch)) +
    geom_boxplot(aes(color = Platform), width = 0.65, linewidth = 0.4, outlier.size = 0.1, outlier.alpha = 0.1, fill = "white", orientation = "y") +
    geom_segment(data = cv_markers,
                 aes(x = Start_CV * 100, xend = End_CV * 100, yend = Batch, color = Platform),
                 arrow = arrow(length = unit(0.1, "cm"), type = "closed"), linewidth = 0.3) +
    geom_point(data = cv_markers, aes(x = Marker_CV * 100, shape = Processing, fill = Platform),
               size = 2, color = "black", stroke = 0.65) +
    scale_color_manual(values = platform_color, guide = "none") +
    scale_fill_manual(values = platform_color, guide = "none") +
    scale_shape_manual(values = c(Unfiltered = 21, Reshaped = 23), name = "Processing") +
    labs(x = "Coefficient of variation (%)", y = NULL) +
    plasmix_theme +
    coord_cartesian(xlim = c(0, 40)) +
    scale_x_continuous(breaks = seq(0, 40, 10), expand = expansion(mult = c(0, 0))) +
    theme(panel.grid.major.y = element_blank(), legend.position = c(1, 1), legend.justification = c(1, 1), legend.key.size = unit(0.3, "cm"),
          legend.title = element_text(margin = margin(b = 3)), legend.margin = margin(5, 0, 5, 5)) +
    guides(shape = guide_legend(override.aes = list(fill = "grey50")))

# 9. Assemble and export Figure 2 ----
showtext_auto()
showtext_opts(dpi = 600)
row1 <- ggarrange(p_u, p_topk, nrow = 1, widths = c(2.25, 1), labels = c("a", "b"),
                  font.label = label_style, label.x = 0, label.y = 1, hjust = -0.2, vjust = 1,
                  common.legend = TRUE, legend = "bottom")
row2 <- ggarrange(p_snr, p_cv, p_r2, p_pass_num, nrow = 1, widths = c(1, 1, 1, 1.1), labels = c("c", "d", "e", "f"),
                  font.label = label_style, label.x = 0, label.y = 1, hjust = -0.2, vjust = 1)
figure2 <- ggarrange(row1, row2, ncol = 1, heights = c(1, 1.5))
ggsave("figures/fig2_titration_benchmark.pdf", figure2, width = 10, height = 6.6)
ggsave("figures/fig2_titration_benchmark.png", figure2, width = 10, height = 6.6, dpi = 600, bg = "white")

# 10. Batch-tier performance source data ----
fmt_num <- function(x, digits = 2) ifelse(is.finite(x), formatC(x, format = "f", digits = digits), NA_character_)
fmt_snr <- function(mean_value, sd_value) {
    ifelse(is.finite(mean_value) & is.finite(sd_value),
           paste0(fmt_num(mean_value), " ± ", fmt_num(sd_value)), fmt_num(mean_value))
}
fmt_median_iqr <- function(median_value, q1, q3, digits = 3) {
    ifelse(is.finite(median_value) & is.finite(q1) & is.finite(q3),
           paste0(fmt_num(median_value, digits), " [", fmt_num(q1, digits), "–", fmt_num(q3, digits), "]"), NA_character_)
}
fmt_pct_frac <- function(numerator, denominator) {
    ifelse(is.finite(denominator) & denominator > 0,
           paste0(fmt_num(100 * numerator / denominator, 1), "% (", format(numerator, big.mark = ",", scientific = FALSE, trim = TRUE), "/", format(denominator, big.mark = ",", scientific = FALSE, trim = TRUE), ")"), NA_character_)
}
safe_median <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) NA_real_ else median(x)
}
safe_quantile <- function(x, probability) {
    x <- x[is.finite(x)]
    if (length(x) == 0) NA_real_ else unname(quantile(x, probability))
}
ensure_columns <- function(data, columns, value = NA) {
    for (column in columns) if (!column %in% colnames(data)) data[[column]] <- value
    data
}

format_quantification_output <- function(platform, batch, data_tier, process_level) {
    case_when(
        platform == "OLK" & data_tier == "Baseline" ~ "ExtNPX",
        platform == "NLS" & data_tier == "Baseline" ~ "ICNorm",
        batch %in% c("SOM_P1_B1", "SOM_P1_B2") & data_tier == "Reshaped" ~ "ANML",
        batch %in% c("SOM_P2_B1", "SOM_P2_B2") & data_tier == "Baseline" ~ "ReadoutNorm",
        batch %in% c("SOM_P2_B1", "SOM_P2_B2") & data_tier == "Calibrated" ~ "PlateNorm",
        batch %in% c("SOM_P2_B1", "SOM_P2_B2") & data_tier == "Reshaped" ~ "SampleNorm",
        TRUE ~ process_level
    )
}

st4_keys <- bind_rows(
    trc_feature_level %>% distinct(Platform, Batch, DataTier, ProcessLevel),
    snr_all_tiers %>% distinct(Platform, Batch, DataTier, ProcessLevel),
    cv_all_tiers %>% distinct(Platform, Batch, DataTier, ProcessLevel),
    gradient_fit_level %>% distinct(Platform, Batch, DataTier, ProcessLevel)
) %>%
    distinct() %>%
    filter(!is.na(Platform), !is.na(Batch), !is.na(DataTier), !is.na(ProcessLevel))

st4_base <- st4_keys %>%
    mutate(`Data tier` = DataTier,
           `Quantification output` = format_quantification_output(Platform, Batch, DataTier, ProcessLevel))

feature_counts <- trc_feature_level %>%
    group_by(Platform, Batch, DataTier, ProcessLevel) %>%
    summarize(`Analytes (All evaluable)` = format(n_distinct(UniqueID), big.mark = ",", scientific = FALSE, trim = TRUE),
              `Analytes (Detected)` = format(n_distinct(UniqueID[Is_Detected %in% TRUE]), big.mark = ",", scientific = FALSE, trim = TRUE), .groups = "drop")

snr_wide <- snr_all_tiers %>%
    mutate(SNR_stat = fmt_snr(SNR_mean, SNR_SD),
           Subset = recode(Subset, `All features` = "All evaluable", `Detected features` = "Detected")) %>%
    select(Platform, Batch, DataTier, ProcessLevel, Subset, SNR_stat) %>%
    pivot_wider(id_cols = c(Platform, Batch, DataTier, ProcessLevel), names_from = Subset,
                values_from = SNR_stat, names_glue = "PCA-based SNR ({Subset})") %>%
    ensure_columns(c("PCA-based SNR (All evaluable)", "PCA-based SNR (Detected)"), NA_character_)

cv_wide <- cv_all_tiers %>%
    mutate(Subset = recode(Subset, `All features` = "All evaluable", `Detected features` = "Detected")) %>%
    group_by(Platform, Batch, DataTier, ProcessLevel, Subset) %>%
    summarize(Median_CV = safe_median(CV), Q1_CV = safe_quantile(CV, 0.25),
              Q3_CV = safe_quantile(CV, 0.75), .groups = "drop") %>%
    mutate(CV_stat = fmt_median_iqr(100 * Median_CV, 100 * Q1_CV, 100 * Q3_CV, digits = 1)) %>%
    select(Platform, Batch, DataTier, ProcessLevel, Subset, CV_stat) %>%
    pivot_wider(id_cols = c(Platform, Batch, DataTier, ProcessLevel), names_from = Subset,
                values_from = CV_stat, names_glue = "Technical CV, % ({Subset})") %>%
    ensure_columns(c("Technical CV, % (All evaluable)", "Technical CV, % (Detected)"), NA_character_)

gradient_fit_wide <- gradient_fit_level %>%
    mutate(All_R2 = round(All_R2, 3), Detected_R2 = round(Detected_R2, 3)) %>%
    select(Platform, Batch, DataTier, ProcessLevel, Gradient, All_R2, Detected_R2) %>%
    pivot_wider(id_cols = c(Platform, Batch, DataTier, ProcessLevel), names_from = Gradient,
                values_from = c(All_R2, Detected_R2), names_sep = "_") %>%
    ensure_columns(c("All_R2_Y", "All_R2_P", "All_R2_X", "Detected_R2_Y", "Detected_R2_P", "Detected_R2_X"), NA_real_) %>%
    rename(`Gradient-fit R² Y (All evaluable)` = All_R2_Y,
           `Gradient-fit R² P (All evaluable)` = All_R2_P,
           `Gradient-fit R² X (All evaluable)` = All_R2_X,
           `Gradient-fit R² Y (Detected)` = Detected_R2_Y,
           `Gradient-fit R² P (Detected)` = Detected_R2_P,
           `Gradient-fit R² X (Detected)` = Detected_R2_X)

feature_metrics <- trc_feature_level %>%
    mutate(Is_Detected = replace_na(Is_Detected, FALSE), TitrationMonoValid = replace_na(TitrationMonoValid, FALSE),
           ExpectedResponseValid = replace_na(ExpectedResponseValid, FALSE),
           Mean_absolute_TRC_deviation = MeanTRCDev) %>%
    group_by(Platform, Batch, DataTier, ProcessLevel) %>%
    summarize(N_all = n(), N_detected = sum(Is_Detected), N_lod_excluded = sum(!Is_Detected),
              Mono_all = sum(TitrationMonoValid), Mono_detected = sum(Is_Detected & TitrationMonoValid),
              TRC_median_all = safe_median(Mean_absolute_TRC_deviation),
              TRC_q1_all = safe_quantile(Mean_absolute_TRC_deviation, 0.25),
              TRC_q3_all = safe_quantile(Mean_absolute_TRC_deviation, 0.75),
              TRC_median_detected = safe_median(Mean_absolute_TRC_deviation[Is_Detected]),
              TRC_q1_detected = safe_quantile(Mean_absolute_TRC_deviation[Is_Detected], 0.25),
              TRC_q3_detected = safe_quantile(Mean_absolute_TRC_deviation[Is_Detected], 0.75),
              Recovery_all = sum(ExpectedResponseValid),
              Recovery_detected = sum(Is_Detected & ExpectedResponseValid),
              Recovery_lod_excluded = sum(!Is_Detected & ExpectedResponseValid), .groups = "drop") %>%
    mutate(`Titration monotonicity (All evaluable)` = fmt_pct_frac(Mono_all, N_all),
           `Titration monotonicity (Detected)` = fmt_pct_frac(Mono_detected, N_detected),
           `Mean absolute TRC deviation (All evaluable)` = fmt_median_iqr(TRC_median_all, TRC_q1_all, TRC_q3_all),
           `Mean absolute TRC deviation (Detected)` = fmt_median_iqr(TRC_median_detected, TRC_q1_detected, TRC_q3_detected),
           `Expected response (All evaluable)` = fmt_pct_frac(Recovery_all, N_all),
           `Expected response (Detected)` = fmt_pct_frac(Recovery_detected, N_detected),
           `Expected response (LoD-excluded)` = fmt_pct_frac(Recovery_lod_excluded, N_lod_excluded)) %>%
    select(Platform, Batch, DataTier, ProcessLevel, `Titration monotonicity (All evaluable)`, `Titration monotonicity (Detected)`,
           `Mean absolute TRC deviation (All evaluable)`, `Mean absolute TRC deviation (Detected)`,
           `Expected response (All evaluable)`,
           `Expected response (Detected)`,
           `Expected response (LoD-excluded)`)

join_keys <- c("Platform", "Batch", "DataTier", "ProcessLevel")
st4 <- st4_base %>%
    left_join(feature_counts, by = join_keys) %>%
    left_join(snr_wide, by = join_keys) %>%
    left_join(cv_wide, by = join_keys) %>%
    left_join(gradient_fit_wide, by = join_keys) %>%
    left_join(feature_metrics, by = join_keys) %>%
    mutate(Platform = factor(Platform, levels = platform_levels),
           DataTier = factor(DataTier, levels = c("Baseline", "Calibrated", "Reshaped"))) %>%
    arrange(Platform, Batch, DataTier, `Quantification output`) %>%
    select(Batch, `Data tier`, `Quantification output`,
           `Analytes (All evaluable)`, `Analytes (Detected)`,
           `PCA-based SNR (All evaluable)`, `PCA-based SNR (Detected)`,
           `Technical CV, % (All evaluable)`, `Technical CV, % (Detected)`,
           `Gradient-fit R² Y (All evaluable)`, `Gradient-fit R² P (All evaluable)`,
           `Gradient-fit R² X (All evaluable)`, `Gradient-fit R² Y (Detected)`,
           `Gradient-fit R² P (Detected)`, `Gradient-fit R² X (Detected)`,
           `Titration monotonicity (All evaluable)`, `Titration monotonicity (Detected)`,
           `Mean absolute TRC deviation (All evaluable)`, `Mean absolute TRC deviation (Detected)`,
           `Expected response (All evaluable)`,
           `Expected response (Detected)`,
           `Expected response (LoD-excluded)`)

st4_definitions <- tribble(
    ~Column, ~Definition,
    "Batch", "Unique batch identifier. The platform is encoded in the batch name.",
    "Data tier", "Data processing tier evaluated: Baseline, Calibrated or Reshaped.",
    "Quantification output", "Platform-specific quantification or processing output.",
    "Analytes (All evaluable)", "Number of analytes included in the gradient evaluation.",
    "Analytes (Detected)", "Number of analytes meeting the platform-specific detection criterion.",
    "PCA-based SNR (All evaluable)", "PCA-based signal-to-noise ratio across all evaluable analytes, shown as mean ± SD across three-replicate combinations.",
    "PCA-based SNR (Detected)", "PCA-based signal-to-noise ratio across detected analytes, shown as mean ± SD across three-replicate combinations.",
    "Technical CV, % (All evaluable)", "Analyte-level technical replicate coefficient of variation across all evaluable analytes, shown as median [Q1–Q3] in percent.",
    "Technical CV, % (Detected)", "Analyte-level technical replicate coefficient of variation across detected analytes, shown as median [Q1–Q3] in percent.",
    "Gradient-fit R² Y (All evaluable)", "Cross-analyte coefficient of determination for recovery of the expected Y-mixture response among all evaluable analytes.",
    "Gradient-fit R² P (All evaluable)", "Cross-analyte coefficient of determination for recovery of the expected P-mixture response among all evaluable analytes.",
    "Gradient-fit R² X (All evaluable)", "Cross-analyte coefficient of determination for recovery of the expected X-mixture response among all evaluable analytes.",
    "Gradient-fit R² Y (Detected)", "Cross-analyte coefficient of determination for recovery of the expected Y-mixture response among detected analytes.",
    "Gradient-fit R² P (Detected)", "Cross-analyte coefficient of determination for recovery of the expected P-mixture response among detected analytes.",
    "Gradient-fit R² X (Detected)", "Cross-analyte coefficient of determination for recovery of the expected X-mixture response among detected analytes.",
    "Titration monotonicity (All evaluable)", "Percentage and fraction of all evaluable analytes for which every evaluable adjacent titration relation received greater than 50% support across three-replicate combinations.",
    "Titration monotonicity (Detected)", "Percentage and fraction of detected analytes for which every evaluable adjacent titration relation received greater than 50% support across three-replicate combinations.",
    "Mean absolute TRC deviation (All evaluable)", "Analyte-level mean absolute deviation from the nominal Y, P and X TRCs, shown as median [Q1–Q3] among all evaluable analytes.",
    "Mean absolute TRC deviation (Detected)", "Analyte-level mean absolute deviation from the nominal Y, P and X TRCs, shown as median [Q1–Q3] among detected analytes.",
    "Expected response (All evaluable)", "Percentage and fraction of all evaluable analytes satisfying titration monotonicity and having a mean absolute TRC deviation <0.25 across at least two intermediate gradients.",
    "Expected response (Detected)", "Percentage and fraction of detected analytes satisfying titration monotonicity and having a mean absolute TRC deviation <0.25 across at least two intermediate gradients.",
    "Expected response (LoD-excluded)", "Percentage and fraction of non-detected evaluable analytes satisfying titration monotonicity and having a mean absolute TRC deviation <0.25 across at least two intermediate gradients."
)

# 11. Source data ----
fig2a_source <- titration_by_effect %>%
    mutate(`Analyte set` = recode(as.character(TargetSpace), All = "All evaluable", Detected = "Detected")) %>%
    transmute(Platform = as.character(Platform), `Analyte set`, `log2(M/F) bin` = MF_bin,
              `Analytes per bin` = Features, `Expected response (%)` = Titration_fidelity)

fig2a_counts_source <- platform_pass_n %>%
    transmute(Platform = as.character(Platform), `Expected-response analytes, n` = N_pass)

fig2b_source <- cumulative_fidelity %>%
    mutate(`Analyte set` = recode(as.character(TargetSpace), All = "All evaluable", Detected = "Detected")) %>%
    transmute(Platform = as.character(Platform), `Analyte set`, `Analyte rank` = Rank,
              `Mean expected-response (%)` = Mean, `Q1 (%)` = Q1, `Q3 (%)` = Q3, Batches)

fig2c_source <- snr_main %>%
    mutate(`Data tier` = as.character(DataTier),
           `Quantification output` = format_quantification_output(Platform, as.character(Batch), as.character(DataTier), ProcessLevel)) %>%
    transmute(Platform, Batch = as.character(Batch), `Data tier`, `Quantification output`,
              `PCA-based SNR` = SNR_mean, `SNR SD` = SNR_SD,
              `Analytes retained` = Assays, Iterations)

fig2d_source <- cv_main %>%
    mutate(`Data tier` = as.character(DataTier),
           `Quantification output` = format_quantification_output(Platform, as.character(Batch), as.character(DataTier), ProcessLevel)) %>%
    transmute(Platform, Batch = as.character(Batch), `Data tier`, `Quantification output`,
              UniqueID, `Technical CV (%)` = 100 * CV)

fig2d_markers_source <- cv_markers %>%
    transmute(Platform, Batch = as.character(Batch), Processing,
              `Marker CV (%)` = 100 * Marker_CV, `Primary CV (%)` = 100 * Default_CV,
              `Arrow start (%)` = 100 * Start_CV, `Arrow end (%)` = 100 * End_CV)

fig2e_source <- global_plot_data %>%
    mutate(`Data tier` = as.character(DataTier),
           `Quantification output` = format_quantification_output(Platform, Batch, as.character(DataTier), ProcessLevel)) %>%
    transmute(Platform, Batch, `Data tier`, `Quantification output`,
              Gradient = as.character(Gradient), `Gradient-fit R²` = Detected_R2)

fig2f_source <- yield_primary %>%
    transmute(Platform, Batch = as.character(Batch), Detected,
              `Expected-response detected analytes` = Valid_detected,
              `Expected-response LoD-excluded analytes` = Valid_lod_excluded,
              Denominator, `Expected-response analytes, total` = Valid_total,
              `Expected-response detected analytes (%)` = Detected_valid_pct,
              `Expected-response LoD-excluded analytes (%)` = Lod_excluded_valid_pct,
              `Detection-adjusted recovery yield (%)` = Fidelity)

fig2f_segments_source <- yield_segments %>%
    transmute(Platform, Batch = as.character(Batch),
              `Detection status` = as.character(Detection_status), Percentage)

fig2f_reshaped_source <- yield_reshaped %>%
    transmute(Platform, Batch = as.character(Batch), Detected,
              `Expected-response detected analytes` = Valid_detected,
              `Expected-response LoD-excluded analytes` = Valid_lod_excluded,
              Denominator, `Expected-response analytes, total` = Valid_total,
              `Detection-adjusted recovery yield (%)` = Fidelity,
              `Primary-tier recovery yield (%)` = Primary_fidelity)

write.xlsx(
    list(Fig2a_effect_response = fig2a_source,
         Fig2a_recovered_counts = fig2a_counts_source,
         Fig2b_cumulative_recovery = fig2b_source,
         Fig2c_SNR = fig2c_source,
         Fig2d_CV = fig2d_source,
         Fig2d_CV_markers = fig2d_markers_source,
         Fig2e_gradient_fit_R2 = fig2e_source,
         Fig2f_recovery_yield = fig2f_source,
         Fig2f_recovery_segments = fig2f_segments_source,
         Fig2f_reshaped = fig2f_reshaped_source,
         Batch_tier_performance = st4,
         Batch_tier_definitions = st4_definitions),
    "tables/SourceData_Figure2.xlsx", overwrite = TRUE, keepNA = TRUE, na.string = "NA"
)

invisible(figure2)
