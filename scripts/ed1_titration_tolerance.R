# Extended Data Figure 1 | TRC deviation and technical precision

# 0. Setup ----
source("scripts/_project_setup.R")
use_packages(c("data.table", "tidyverse", "ggpubr", "openxlsx", "showtext"), "scales")
source("utils/figure_style.R")
plasmix_theme <- if (is.function(theme_plasmix)) theme_plasmix() else theme_plasmix
label_style <- list(size = 12, face = "bold")

# 1. Inputs ----
trc_path <- "results/trc_feature_level.tsv.gz"
trc_iteration_path <- "results/trc_iteration_level.tsv.gz"
cv_path <- "results/cv_feature_level.tsv.gz"
missing_inputs <- c(trc_path, trc_iteration_path, cv_path)[!file.exists(c(trc_path, trc_iteration_path, cv_path))]
if (length(missing_inputs) > 0) stop("Missing input files: ", paste(missing_inputs, collapse = ", "), call. = FALSE)

trc_feature_level <- fread(trc_path)
trc_iteration_level <- fread(trc_iteration_path)
cv_feature_level <- fread(cv_path)

filter_primary_process <- function(data) {
    data %>% filter(
        (Platform == "AAG" & DataTier == "Baseline" & ProcessLevel == "SNR") |
            (Platform == "DIA" & DataTier == "Baseline" & ProcessLevel == "Intensity") |
            (Platform == "SOM" & DataTier == "Calibrated" & ProcessLevel == "Calibrate") |
            (Platform == "NLS" & DataTier == "Calibrated" & grepl("NPQ", ProcessLevel, ignore.case = TRUE)) |
            (Platform == "OLK" & DataTier == "Calibrated" & grepl("NPX", ProcessLevel, ignore.case = TRUE))
    )
}

platform_levels <- c("DIA", "SOM", "OLK", "NLS", "AAG")

default_trc <- trc_feature_level %>%
    filter_primary_process() %>%
    mutate(Platform = factor(Platform, levels = platform_levels),
          Is_Detected = replace_na(Is_Detected, FALSE),
          TitrationMonoValid = replace_na(TitrationMonoValid, FALSE),
          ExpectedResponseValid = replace_na(ExpectedResponseValid, FALSE))

detected_trc <- default_trc %>% filter(Is_Detected %in% TRUE)

default_trc_iteration <- trc_iteration_level %>%
    filter_primary_process() %>%
    mutate(Platform = factor(Platform, levels = platform_levels),
          Is_Detected = replace_na(Is_Detected, FALSE),
          TitrationMonoValid = replace_na(TitrationMonoValid, FALSE),
          ExpectedResponseValid = replace_na(ExpectedResponseValid, FALSE))

detected_trc_iteration <- default_trc_iteration %>% filter(Is_Detected %in% TRUE)

# 2. Panel a: TRC-deviation tolerance ----
calc_tolerance_summary <- function(data) {
    long_data <- data %>%
        select(Platform, Batch, UniqueID, Iteration, TRC_Y, TRC_P, TRC_X) %>%
        pivot_longer(c(TRC_Y, TRC_P, TRC_X), names_to = "Gradient", values_to = "TRC") %>%
        mutate(Gradient = sub("TRC_", "", Gradient),
              Expected = recode(Gradient, Y = 0.25, P = 0.50, X = 0.75),
              Absolute_deviation = abs(TRC - Expected))
    map_dfr(c(0.05, 0.10, 0.15, 0.20, 0.25), function(cutoff) {
        long_data %>%
            group_by(Platform, Batch, Gradient, Iteration) %>%
            summarize(Features = n_distinct(UniqueID), Evaluable = sum(is.finite(Absolute_deviation)),
                      Within_cutoff = sum(Absolute_deviation < cutoff, na.rm = TRUE),
                      Percentage = if_else(Evaluable > 0, 100 * Within_cutoff / Evaluable, NA_real_), .groups = "drop") %>%
            group_by(Platform, Batch, Gradient) %>%
            summarize(Cutoff = cutoff, Features = max(Features), Evaluable = mean(Evaluable),
                      Within_cutoff = mean(Within_cutoff),
                      Percentage_SD = if (sum(is.finite(Percentage)) > 1) sd(Percentage, na.rm = TRUE) else NA_real_,
                      Percentage = mean(Percentage, na.rm = TRUE), Iterations = n_distinct(Iteration), .groups = "drop")
    })
}

tolerance_data <- calc_tolerance_summary(detected_trc_iteration)

batch_label_order <- tolerance_data %>%
    transmute(Platform, Batch, Batch_label = sprintf("%s (n=%s)", Batch, scales::comma(Features))) %>%
    distinct() %>%
    arrange(factor(Platform, levels = platform_levels), Batch) %>%
    pull(Batch_label) %>% unique() %>% rev()

heatmap_data <- tolerance_data %>%
    mutate(Platform = factor(Platform, levels = platform_levels),
          Gradient = factor(Gradient, levels = c("Y", "P", "X")),
          Batch_label = factor(sprintf("%s (n=%s)", Batch, scales::comma(Features)), levels = batch_label_order),
          Cutoff_label = factor(sprintf("%.2f", Cutoff), levels = sprintf("%.2f", c(0.05, 0.10, 0.15, 0.20, 0.25))),
          Label = if_else(is.finite(Percentage), sprintf("%.1f", Percentage), ""),
          Text_color = if_else(is.finite(Percentage) & Percentage >= 50, "white", "black"))

p_a <- ggplot(heatmap_data, aes(x = Cutoff_label, y = Batch_label, fill = Percentage)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = Label, color = Text_color), size = 2.5) +
    scale_color_identity() +
    scale_fill_gradientn(colours = colorRampPalette(c("#F7FBFF", "#6BAED6", "#08306B"))(100),
                        limits = c(0, 75), breaks = seq(0, 75, 25), oob = scales::squish,
                        na.value = "grey90", name = "Features within cutoff (%)") +
    facet_wrap(~Gradient, ncol = 3) +
    labs(x = "TRC absolute-deviation cutoff", y = "Batch (n = detected features)") +
    plasmix_theme +
    theme(panel.grid.major = element_blank(),
          legend.position = "bottom",
          legend.title.position = "left",
          legend.title = element_text(face = "bold", margin = margin(r = 5)),
          legend.box.spacing = unit(1, "pt"),
          legend.margin = margin(t = 3, l = -10, b = 0),
          legend.text = element_text(margin = margin(t = 2, r = 0, b = 0, l = 0))) +
    guides(fill = guide_colorbar(title.position = "left", title.vjust = 0.75, barwidth = unit(3, "cm"), barheight = unit(0.25, "cm")))

# 3. Panel b: TRC-deviation distributions ----
deviation_data <- detected_trc %>%
    filter(TitrationMonoValid %in% TRUE) %>%
    select(Platform, Batch, UniqueID, TRC_abs_err_Y, TRC_abs_err_P, TRC_abs_err_X) %>%
    pivot_longer(starts_with("TRC_abs_err_"), names_to = "Gradient", values_to = "Absolute_deviation") %>%
    mutate(Platform = factor(Platform, levels = platform_levels),
          Gradient = factor(sub("TRC_abs_err_", "", Gradient), levels = c("Y", "P", "X"))) %>%
    filter(is.finite(Absolute_deviation))

threshold_grid <- seq(0, 0.75, by = 0.0025)

deviation_survival <- deviation_data %>%
    group_by(Platform, Gradient) %>%
    group_modify(~{
        values <- .x$Absolute_deviation
        tibble(Absolute_deviation = threshold_grid,
               Exceeding_fraction = 100 * vapply(threshold_grid, function(cutoff) mean(values >= cutoff), numeric(1)))
    }) %>%
    ungroup()

p_b <- ggplot(deviation_survival, aes(x = Absolute_deviation, y = Exceeding_fraction, color = Platform)) +
    geom_line(linewidth = 0.55) +
    facet_wrap(~Gradient, ncol = 1) +
    scale_color_manual(values = platform_color, limits = platform_levels, name = "Platform") +
    scale_x_continuous(limits = c(0, 0.75), breaks = seq(0, 0.75, 0.15), labels = function(x) formatC(x, format = "f", digits = 2), expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0))) +
    labs(x = "Absolute TRC deviation threshold", y = "Features exceeding threshold (%)") +
    plasmix_theme +
    theme(plot.margin = margin(5, 8, 5, 5),
          legend.position = "bottom", legend.box = "horizontal",
          legend.key.width = unit(0.4, "cm"),
          legend.box.margin = margin(-8, 20, 0, -10)) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(linewidth = 0.7)))

# cv_primary <- cv_feature_level %>%
#     filter(Subset == "Detected features") %>%
#     filter_primary_process() %>%
#     select(Platform, Batch, ProcessLevel, DataTier, UniqueID, CV) %>%
#     filter(is.finite(CV), CV >= 0)

# cv_outcome_data <- detected_trc %>%
#     select(Platform, Batch, ProcessLevel, DataTier, UniqueID, TRC_N_finite, MonoRelationN, MonoRelationEvaluableN, TitrationMonoValid, TRCDevValid) %>%
#     inner_join(cv_primary, by = c("Platform", "Batch", "ProcessLevel", "DataTier", "UniqueID")) %>%
#     mutate(ExpectedResponseEvaluable = TRC_N_finite >= 2 & MonoRelationN >= 3 & MonoRelationEvaluableN == MonoRelationN,
#            CV_percent = 100 * CV, CV_interval = cut(CV_percent, breaks = c(-Inf, 10, 20, 30, Inf), labels = c("≤10", "10–20", "20–30", ">30"), right = TRUE)) %>%
#     filter(ExpectedResponseEvaluable, !is.na(CV_interval)) %>%
#     mutate(Outcome = case_when(
#               !TitrationMonoValid ~ "Monotonicity failure",
#               TitrationMonoValid & !TRCDevValid ~ "TRC-deviation failure",
#               TitrationMonoValid & TRCDevValid ~ "Expected response",
#               TRUE ~ NA_character_
#           ),
#           Platform = factor(Platform, levels = platform_levels),
#           CV_interval = factor(CV_interval, levels = c("≤10", "10–20", "20–30", ">30")),
#           Outcome = factor(Outcome, levels = c("Monotonicity failure", "TRC-deviation failure", "Expected response"))) %>%
#     filter(!is.na(Outcome))

# cv_composition <- cv_outcome_data %>%
#     count(Platform, CV_interval, Outcome, name = "Feature_batch_observations") %>%
#     complete(Platform, CV_interval, Outcome, fill = list(Feature_batch_observations = 0)) %>%
#     group_by(Platform, CV_interval) %>%
#     mutate(Total = sum(Feature_batch_observations),
#            Percentage = if_else(Total > 0, 100 * Feature_batch_observations / Total, NA_real_)) %>%
#     ungroup()

# p_c <- ggplot(cv_composition, aes(x = CV_interval, y = Percentage, fill = Outcome)) +
#     geom_col(width = 0.75, color = "white", linewidth = 0.2) +
#     facet_wrap(~Platform, nrow = 1) +
#     scale_fill_manual(values = c("Monotonicity failure" = "#BDBDBD", "TRC-deviation failure" = "#F4B183", "Expected response" = "#397DB7"), name = "Outcome") +
#     scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0))) +
#     labs(x = "Technical CV (%)", y = "Feature–batch observations (%)") +
#     plasmix_theme +
#     theme(panel.grid.major.x = element_blank(), legend.position = "bottom", legend.box.margin = margin(-8, 0, 0, 0)) +
#     guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# 4. Panel c: Precision and response validity ----
cv_bin_width <- 2.5
cv_plot_min <- 0
cv_plot_max <- 25
y_plot_min <- 0
y_plot_max <- 53
min_bin_fraction <- 0.025
trc_cutoffs <- c(0.05, 0.10, 0.15, 0.20, 0.25)

metric_levels <- c("Titration monotonicity", sprintf("Mean TRC deviation < %.2f", trc_cutoffs))
metric_labels <- c("Monotonicity", sprintf("TRC < %.2f", trc_cutoffs))
metric_colors <- setNames(c("#0ABAB5", "#FDD49E", "#FDBB84", "#FC8D59", "#E34A33", "#B30000"), metric_levels)

cv_primary <- cv_feature_level %>%
    filter(Subset == "Detected features") %>% filter_primary_process() %>%
    select(Platform, Batch, ProcessLevel, DataTier, UniqueID, CV) %>%
    filter(is.finite(CV), CV >= 0)

cv_outcome_data <- detected_trc %>%
    select(Platform, Batch, ProcessLevel, DataTier, UniqueID, TRC_N_finite, MonoRelationN, MonoRelationEvaluableN,
           TitrationMonoValid, MeanTRCDev) %>%
    inner_join(cv_primary, by = c("Platform", "Batch", "ProcessLevel", "DataTier", "UniqueID")) %>%
    mutate(ExpectedResponseEvaluable = TRC_N_finite >= 2 & MonoRelationN >= 3 & MonoRelationEvaluableN == MonoRelationN,
           TitrationMonoValid = replace_na(TitrationMonoValid, FALSE), CV_percent = 100 * CV,
           CV_bin_lower = floor(CV_percent / cv_bin_width) * cv_bin_width, CV_bin_upper = CV_bin_lower + cv_bin_width,
           CV_bin = CV_bin_lower + cv_bin_width / 2, Platform = factor(Platform, levels = platform_levels)) %>%
    filter(ExpectedResponseEvaluable, is.finite(MeanTRCDev), is.finite(CV_percent), CV_percent >= 0) %>%
    group_by(Platform, Batch) %>% mutate(Batch_total_features = n()) %>% ungroup()

cv_batch_composition_all <- bind_rows(
    cv_outcome_data %>%
        group_by(Platform, Batch, CV_bin_lower, CV_bin_upper, CV_bin) %>%
        summarize(Features = n(), Batch_total_features = first(Batch_total_features),
                  Percentage = 100 * mean(TitrationMonoValid), .groups = "drop") %>%
        mutate(Metric = "Titration monotonicity"),
    map_dfr(trc_cutoffs, function(cutoff) {
        cv_outcome_data %>%
            group_by(Platform, Batch, CV_bin_lower, CV_bin_upper, CV_bin) %>%
            summarize(Features = n(), Batch_total_features = first(Batch_total_features),
                      Percentage = 100 * mean(MeanTRCDev < cutoff), .groups = "drop") %>%
            mutate(Metric = sprintf("Mean TRC deviation < %.2f", cutoff))
    })
) %>%
    mutate(Bin_fraction = Features / Batch_total_features, Included_by_fraction = Bin_fraction >= min_bin_fraction)

platform_batch_n <- cv_outcome_data %>% distinct(Platform, Batch) %>% count(Platform, name = "Total_batches")

cv_platform_composition_all <- cv_batch_composition_all %>%
    filter(Included_by_fraction) %>%
    group_by(Platform, CV_bin_lower, CV_bin_upper, CV_bin, Metric) %>%
    summarize(Percentage = mean(Percentage), Percentage_SD = if (n() > 1) sd(Percentage) else NA_real_,
              Batch_support = n_distinct(Batch), Mean_bin_fraction = mean(Bin_fraction),
              Min_bin_fraction = min(Bin_fraction), Max_bin_fraction = max(Bin_fraction), .groups = "drop") %>%
    left_join(platform_batch_n, by = "Platform") %>%
    mutate(Minimum_batches = pmax(1L, ceiling(Total_batches / 2)),
           Included_by_batch_support = Batch_support >= Minimum_batches,
           Metric = factor(Metric, levels = metric_levels))

cv_composition <- cv_platform_composition_all %>%
    filter(Included_by_batch_support) %>%
    group_by(Platform, Metric) %>% arrange(CV_bin, .by_group = TRUE) %>%
    mutate(Segment = cumsum(c(TRUE, diff(CV_bin) > cv_bin_width * 1.5))) %>% ungroup()

cv_plot_data <- cv_composition %>%
    filter(CV_bin_lower >= cv_plot_min, CV_bin_upper <= cv_plot_max)

if (nrow(cv_plot_data) == 0) stop("Panel c has no data within the configured CV display range.", call. = FALSE)
if (max(cv_plot_data$Percentage, na.rm = TRUE) > y_plot_max)
    warning(sprintf("Panel c contains points above the %.1f%% y-axis limit; statistics are unchanged, but these points are clipped.", y_plot_max), call. = FALSE)

plasmix_theme_no_coord <- plasmix_theme[!vapply(plasmix_theme, function(x) inherits(x, "Coord"), logical(1))]

p_c <- ggplot(cv_plot_data, aes(CV_bin, Percentage, color = Metric, group = interaction(Metric, Segment))) +
    geom_line(linewidth = 0.55, alpha = 0.9) +
    geom_point(size = 1.25, alpha = 0.9) +
    facet_wrap(~Platform, nrow = 1) +
    scale_color_manual(values = metric_colors, breaks = metric_levels, labels = metric_labels, name = NULL) +
    scale_x_continuous(breaks = seq(cv_plot_min, cv_plot_max, 5), expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50), expand = expansion(mult = c(0, 0))) +
    labs(x = "Technical CV (%)", y = "Features meeting criterion (%)") +
    plasmix_theme_no_coord +
    coord_cartesian(xlim = c(cv_plot_min, cv_plot_max), ylim = c(y_plot_min, y_plot_max), clip = "off") +
    theme(panel.grid.major.x = element_blank(), panel.spacing.x = unit(0.6, "lines"),
          legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal",
          legend.box.margin = margin(-6, 0, 0, 0)) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(linewidth = 0.7, size = 1.6, alpha = 1)))

p_c

# 5. Assemble and export ----
showtext_auto()
showtext_opts(dpi = 600)

row_top <- ggarrange(p_a, p_b, nrow = 1, widths = c(2.2, 1), labels = c("a", "b"), font.label = label_style,
                     label.x = 0, label.y = 1, hjust = -0.2, vjust = 1)

ed1_figure <- ggarrange(row_top, p_c, ncol = 1, heights = c(2, 1), labels = c("", "c"), font.label = label_style,
                        label.x = 0, label.y = 1, hjust = -0.2, vjust = 1)

ggsave("figures/ed1_trc_error_and_precision.pdf", ed1_figure, width = 10, height = 8, bg = "white")
ggsave("figures/ed1_trc_error_and_precision.png", ed1_figure, width = 10, height = 8, dpi = 600, bg = "white")

panel_c_parameters <- tibble(
    Parameter = c("CV bin width (%)", "Displayed CV range (%)", "Displayed y-axis range (%)",
                  "Minimum batch-bin fraction", "Minimum batch support", "TRC cutoffs"),
    Value = c(cv_bin_width, sprintf("%.1f to %.1f", cv_plot_min, cv_plot_max),
              sprintf("%.1f to %.1f", y_plot_min, y_plot_max), min_bin_fraction,
              "ceiling(total batches / 2), minimum 1", paste(sprintf("%.2f", trc_cutoffs), collapse = ", "))
)

write.xlsx(
    list(Gradient_cutoff_summary = tolerance_data, TRC_deviation_feature_level = deviation_data,
         TRC_deviation_survival = deviation_survival, CV_TRC_feature_level = cv_outcome_data,
         CV_TRC_batch_bins = cv_batch_composition_all, CV_TRC_platform_bins = cv_platform_composition_all,
         CV_TRC_plot_data = cv_plot_data, Panel_c_parameters = panel_c_parameters),
    "tables/SourceData_EDFigure1.xlsx", overwrite = TRUE, keepNA = TRUE, na.string = "NA"
)

invisible(ed1_figure)
