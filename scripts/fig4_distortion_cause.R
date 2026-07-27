# Figure 4 | Analytical distortion and physicochemical drivers

# 0. Setup ----
source("scripts/_project_setup.R")
use_packages(c("data.table", "tidyverse", "readxl", "ggh4x", "ggpubr", "patchwork", "ggtext", "openxlsx", "showtext", "grid", "xgboost", "lightgbm"), c("pdp", "randomForest"))
source("utils/figure_style.R")
plasmix_theme <- if (is.function(theme_plasmix)) theme_plasmix() else theme_plasmix
showtext_auto(); showtext_opts(dpi = 600)
label_style <- list(size = 12, face = "bold")
paths <- c(metadata = "data/study_metadata.xlsx", profiles = "data/protein_profiles_long.tsv.gz", detection = "results/detection_status.tsv.gz",
           physchem = "data/physchem_matrix.tsv.gz", physchem_dictionary = "data/physchem_dictionary.tsv")
missing_inputs <- paths[!file.exists(paths)]
if (length(missing_inputs)) stop("Figure 4 is missing the following release inputs:\n", paste(missing_inputs, collapse = "\n"))
meta_batch <- read_xlsx(paths["metadata"], sheet = "batch")
long_df <- fread(paths["profiles"]) %>% as_tibble()
lod_status <- fread(paths["detection"]) %>% as_tibble()
physchem_matrix <- fread(paths["physchem"]) %>% as_tibble()
physchem_dict <- fread(paths["physchem_dictionary"]) %>% as_tibble()
meta_batch_ht <- meta_batch %>% filter(Platform %in% c("SOM", "OLK", "DIA"), !Batch %in% c("OLK_P1_B1", "OLK_P1_B2"))

# Physicochemical annotations ----
cat_colors <- c("Structure" = "#E64B35", "Surface" = "#4DBBD5", "Charge" = "#00A087", "Disorder" = "#F39B7F", "Secretory" = "#8491B4", "Abundance" = "#91D1C2")
physchem_dict <- physchem_dict %>% mutate(Category = factor(Category, levels = names(cat_colors)))
final_physchem <- physchem_dict %>% filter(Retained == "Yes") %>% pull(Feature)
# Lookup vectors for property labels and categories.
name_map <- setNames(physchem_dict$Property, physchem_dict$Feature)
cat_map  <- setNames(physchem_dict$Category, physchem_dict$Feature)

# 1. Empirical response-envelope parameters ----
# HPA abundance is taken directly from the release annotation matrix.
hpa_sel <- physchem_matrix %>%
    transmute(UniProtID = Entry, BloodConc_log10_pgml, Abundance_Source) %>%
    filter(is.finite(BloodConc_log10_pgml), Abundance_Source != "Unknown") %>%
    distinct()

# Estimate replicate CV from normalized measurements.
ic_intensity <- long_df %>%
  filter(Batch %in% unique(meta_batch_ht$Batch), Sample %in% c("M", "Y", "X", "P", "F", "N")) %>%
  filter(
    (Platform %in% c("SOM", "OLK") & ProcessLevel == "HybNorm") |
    (Platform == "DIA" & ProcessLevel == "Intensity")
  ) %>%
  group_by(Batch, Sample, UniqueID, UniProtID, Platform) %>%
  summarize(
    Intensity = mean(2^Value, na.rm = TRUE),
    Protein_CV = sd(2^Value, na.rm = TRUE) / mean(2^Value, na.rm = TRUE),
    .groups = 'drop'
  )
ic_intensity_conc <- left_join(ic_intensity, hpa_sel, by = "UniProtID")

# Summarize each intensity interval by the median replicate CV of its measurements.
probs <- c(0, 0.001, seq(0.01, 0.99, by = 0.01), 0.999, 1)
calc_interval_stats <- function(df) {
  plat <- unique(df$Platform)[1]
  bat <- unique(df$Batch)[1]
  # Intensity quantile boundaries pooled across samples.
  q_bounds <- quantile(df$Intensity, probs = probs, na.rm = TRUE)
  res <- map_df(1:(length(probs)-1), function(i) {
    low_p <- probs[i]
    high_p <- probs[i+1]
    low_val <- q_bounds[i]
    high_val <- q_bounds[i+1]
    # Measurements within the current interval.
    if (i == 1) {
      sub_df <- df %>% filter(Intensity >= low_val & Intensity <= high_val)
    } else {
      sub_df <- df %>% filter(Intensity > low_val & Intensity <= high_val)
    }
    n_prot <- nrow(sub_df)
    if (n_prot > 0) {
      mean_int <- mean(sub_df$Intensity, na.rm = TRUE)
      # Median replicate CV limits the influence of failed replicates.
      cv <- median(sub_df$Protein_CV, na.rm = TRUE)
      n_hpa_valid <- sum(!is.na(sub_df$BloodConc_log10_pgml))
      hpa_cov <- n_hpa_valid / n_prot
      mean_hpa <- mean(sub_df$BloodConc_log10_pgml, na.rm = TRUE)
    } else {
      mean_int <- NA
      cv <- NA
      hpa_cov <- NA
      mean_hpa <- NA
    }
    interval_name <- paste0(low_p*100, "%-", high_p*100, "%")
    tibble(
      Platform = plat,
      Batch = bat,
      Interval = interval_name,
      Interval_Order = i,
      Interval_Lower_Intensity = low_val,
      Interval_Upper_Intensity = high_val,
      N_Proteins = n_prot,
      Mean_Intensity = mean_int,
      CV = cv,
      Mean_HPA_Conc_Log10 = mean_hpa,
      HPA_Coverage = hpa_cov
    )
  })
  return(res)
}

# Build the interval-level summary table.
final_long_df <- ic_intensity_conc %>%
  filter(!is.na(Intensity)) %>%
  group_by(Platform, Batch) %>%
  group_split() %>%
  map_dfr(calc_interval_stats) %>%
  pivot_longer(
    cols = c(Interval_Lower_Intensity, Interval_Upper_Intensity, N_Proteins, Mean_Intensity, CV, Mean_HPA_Conc_Log10, HPA_Coverage),
    names_to = "Metric", values_to = "Value")

# Estimate empirical sigmoidal response-envelope parameters.
estimate_response_envelope_parameters <- function(sub_df) {
  plat_data <- sub_df %>%
    group_by(Interval_Order, Metric) %>%
    summarize(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Metric, values_from = Value) %>%
    arrange(Interval_Order)

  min_order <- min(plat_data$Interval_Order)
  max_order <- max(plat_data$Interval_Order)

  # Empirical lower and upper response endpoints.
  A_val <- plat_data %>% filter(Interval_Order == min_order) %>% pull(Mean_Intensity)
  D_val <- plat_data %>% filter(Interval_Order == max_order) %>% pull(Mean_Intensity)

  # Inner-tail anchors used only for background burden.
  A_bg <- plat_data %>% filter(Interval_Order == min_order + 1) %>% pull(Mean_Intensity)
  D_bg <- plat_data %>% filter(Interval_Order == max_order - 1) %>% pull(Mean_Intensity)
  background_burden_pct <- if(is.finite(A_bg) && is.finite(D_bg) && D_bg > A_bg) 100 * A_bg / (D_bg - A_bg) else NA_real_

  # HPA abundance represented by the two intervals surrounding the 50% intensity boundary.
  central_hpa <- plat_data %>%
    filter(Interval_Order %in% c(51, 52), is.finite(Mean_HPA_Conc_Log10)) %>%
    mutate(HPA_n = N_Proteins * HPA_Coverage)

  abundance_mid_log <- weighted.mean(central_hpa$Mean_HPA_Conc_Log10, central_hpa$HPA_n, na.rm = TRUE)

  # Response slope based on the HPA abundance span covering 1%–99% of the A–D response range.
  clean_hpa <- plat_data %>% filter(is.finite(Mean_HPA_Conc_Log10), is.finite(Mean_Intensity)) %>% arrange(Mean_Intensity)
  S_1 <- A_val + 0.01 * (D_val - A_val)
  S_99 <- A_val + 0.99 * (D_val - A_val)
  hpa_1 <- approx(clean_hpa$Mean_Intensity, clean_hpa$Mean_HPA_Conc_Log10, xout = S_1, rule = 2, ties = mean)$y
  hpa_99 <- approx(clean_hpa$Mean_Intensity, clean_hpa$Mean_HPA_Conc_Log10, xout = S_99, rule = 2, ties = mean)$y
  delta_hpa <- hpa_99 - hpa_1
  response_slope <- if(is.finite(delta_hpa) && delta_hpa > 0) 2 * log10(99) / delta_hpa else NA_real_

  # Lower 20%, middle 60% and upper 20% CV components.
  cv_bottom <- plat_data %>% filter(Interval_Order <= 21) %>% summarize(med_cv = median(CV, na.rm = TRUE)) %>% pull(med_cv)
  cv_tech <- plat_data %>% filter(Interval_Order > 21, Interval_Order <= 81) %>%
    summarize(med_cv = median(CV, na.rm = TRUE)) %>% pull(med_cv)
  cv_top <- plat_data %>% filter(Interval_Order > 81) %>% summarize(med_cv = median(CV, na.rm = TRUE)) %>% pull(med_cv)

  # Noise-adjusted response bounds.
  bg_sd <- A_val * cv_bottom
  I_low <- A_val + 3 * bg_sd
  D_sd <- D_val * sqrt(cv_tech^2 + cv_top^2)
  I_high <- D_val - 3 * D_sd

  abundance_mid <- 10^abundance_mid_log
  abundance_low <- if(is.finite(response_slope) && I_low > A_val && I_low < D_val) {
    abundance_mid / ((D_val - A_val) / (I_low - A_val) - 1)^(1 / response_slope)
  } else NA_real_
  abundance_high <- if(is.finite(response_slope) && I_high > A_val && I_high < D_val) {
    abundance_mid / ((D_val - A_val) / (I_high - A_val) - 1)^(1 / response_slope)
  } else NA_real_

  tibble(
    A_log = log10(A_val), D_log = log10(D_val), Abundance_mid_log = abundance_mid_log,
    Response_slope = response_slope, Background_burden_pct = background_burden_pct,
    bg_sd = bg_sd, cv_tech = cv_tech, cv_bot = cv_bottom, cv_top = cv_top,
    Lower_bound_log = log10(I_low), Upper_bound_log = log10(I_high),
    Abundance_min = log10(abundance_low), Abundance_max = log10(abundance_high),
    HPA_lower_anchor_log = hpa_1, HPA_upper_anchor_log = hpa_99
  )
}

# Estimate response envelopes at platform and batch levels.
envelope_by_platform <- final_long_df %>% group_by(Platform) %>% group_modify(~ estimate_response_envelope_parameters(.x)) %>% ungroup()
envelope_by_batch <- final_long_df %>% group_by(Platform, Batch) %>% group_modify(~ estimate_response_envelope_parameters(.x)) %>% ungroup()
print(envelope_by_platform)
print(envelope_by_batch)

# Panel a: Empirical envelope of reference abundance and measured response ----
format_10_exp <- function(x) { parse(text = paste0("10^", x)) }
set.seed(42)
abundance_values <- 10^seq(-10, 10, length.out = 600)
generate_response_envelope_curve_data <- function(abundance, params) {
  A <- 10^params$A_log
  D <- 10^params$D_log
  abundance_mid <- 10^params$Abundance_mid_log
  response_slope <- params$Response_slope
  bg_sd <- params$bg_sd
  cv_tech <- params$cv_tech
  cv_top <- params$cv_top

  signal <- A + (D - A) / (1 + (abundance_mid / abundance)^response_slope)
  signal_ratio <- pmax(0, pmin(1, (signal - A) / (D - A)))
  cv_total <- sqrt((bg_sd / signal)^2 + cv_tech^2 + (signal_ratio^6) * cv_top^2)
  band_width_log <- log10(1 + 3 * cv_total)
  ideal_log <- log10(signal)

  tibble(
    log_abundance = log10(abundance),
    ideal = ideal_log,
    noisy = ideal_log + rnorm(length(abundance), 0, band_width_log * 0.3),
    upper = ideal_log + band_width_log,
    lower = ideal_log - band_width_log
  )
}

build_response_envelope_curve_dataset <- function(envelope_df, group_col = "Batch") {
  pmap_dfr(envelope_df, function(...) {
    row_data <- list(...)
    params <- list(A_log = row_data$A_log, D_log = row_data$D_log, Abundance_mid_log = row_data$Abundance_mid_log,
                   Response_slope = row_data$Response_slope, bg_sd = row_data$bg_sd,
                   cv_tech = row_data$cv_tech, cv_top = row_data$cv_top)
    if(is.na(params$Response_slope)) return(NULL)
    curve_data <- generate_response_envelope_curve_data(abundance_values, params)
    curve_data$Platform <- row_data$Platform
    curve_data[[group_col]] <- row_data[[group_col]]
    curve_data
  })
}
response_curve_platform <- build_response_envelope_curve_dataset(envelope_by_platform, group_col = "Platform")
response_curve_batch <- build_response_envelope_curve_dataset(envelope_by_batch, group_col = "Batch")

# Plot one platform response envelope.
plot_response_envelope <- function(curve_df, envelope_row, prefix, title, color, y_label = FALSE) {
    A_log <- envelope_row$A_log
    D_log <- envelope_row$D_log
    abundance_min <- envelope_row$Abundance_min
    abundance_max <- envelope_row$Abundance_max
    lower_bound_y <- envelope_row$Lower_bound_log
    upper_bound_y <- envelope_row$Upper_bound_log
    abundance_span <- round(abundance_max - abundance_min, 1)
    intensity_span <- round(D_log - A_log, 1)
    arr_x <- 6
    arr_y_horiz <- -0.9

    p <- ggplot(curve_df, aes(x = log_abundance)) +
        geom_hline(yintercept = A_log, color = "gray30", linetype = "dotted", linewidth = 0.5, alpha = 0.5) +
        geom_hline(yintercept = D_log, color = "gray30", linetype = "dotted", linewidth = 0.5, alpha = 0.5) +
        geom_ribbon(aes(ymin = lower, ymax = upper), fill = color, alpha = 0.3) +
        geom_line(aes(y = noisy), color = color, linewidth = 0.3, alpha = 0.8) +
        geom_line(aes(y = ideal), color = color, linewidth = 0.8) +
        annotate("segment", x = abundance_min, xend = abundance_min, y = -1.1, yend = lower_bound_y, color = color,
                 linetype = "dotted", linewidth = 0.5, alpha = 0.5) +
        annotate("segment", x = abundance_max, xend = abundance_max, y = -1.1, yend = upper_bound_y, color = color,
                 linetype = "dotted", linewidth = 0.5, alpha = 0.5) +
        annotate("point", x = abundance_min, y = lower_bound_y, size = 3.2, color = color, alpha = 0.3) +
        annotate("point", x = abundance_min, y = lower_bound_y, size = 1.8, color = color) +
        annotate("text", x = abundance_min + ifelse(prefix == "som", 0.8, 0.1), y = lower_bound_y + 0.85,
                 label = "LB", color = color, size = 2.8, fontface = "bold", hjust = 1) +
        annotate("point", x = abundance_max, y = upper_bound_y, size = 3.2, color = color, alpha = 0.3) +
        annotate("point", x = abundance_max, y = upper_bound_y, size = 1.8, color = color) +
        annotate("text", x = abundance_max - ifelse(prefix == "olk", 0, ifelse(prefix == "dia", 0.7, 0.5)),
                 y = upper_bound_y + ifelse(prefix == "olk", 1.2, ifelse(prefix == "dia", -0.2, 0.6)),
                 label = "UB", color = color, size = 2.8, fontface = "bold", hjust = 1) +
        annotate("richtext", x = 5.4, y = A_log + 0.5, label = "<i><b>A</b></i>", color = color, size = 2.8, label.padding = unit(0.1, "lines")) +
        annotate("richtext", x = 5.4, y = D_log - 0.7, label = "<i><b>D</b></i>", color = color, size = 2.8, label.padding = unit(0.1, "lines")) +
        annotate("segment", x = arr_x, xend = arr_x, y = A_log, yend = D_log, arrow = arrow(ends = "both", length = unit(0.045, "inches"), type = "closed"), color = "gray10", linewidth = 0.3) +
        annotate("text", x = arr_x - 0.2, y = A_log + (D_log - A_log) / 2,
                 label = paste0("Intensity span: ~", intensity_span, " logs"), hjust = 1, color = "gray10", size = 2.8) +
        annotate("segment", x = abundance_min, xend = abundance_max, y = arr_y_horiz, yend = arr_y_horiz,
                 arrow = arrow(ends = "both", length = unit(0.045, "inches"), type = "closed"), color = "gray10", linewidth = 0.3) +
        annotate("text", x = abundance_min, y = arr_y_horiz + 0.6,
                 label = paste0("Abundance span: ", abundance_span, " logs"), hjust = 0, color = "gray10", size = 2.8) +
        scale_x_continuous(limits = c(-3, 6), breaks = seq(-4, 8, 2), labels = format_10_exp, expand = c(0,0)) +
        scale_y_continuous(limits = c(-1.3, 10), breaks = seq(0, 10, 2), labels = format_10_exp, expand = c(0,0)) +
        coord_cartesian(clip = "off") +
        labs(title = title, x = "HPA abundance (pg/mL)", y = "Observed intensity") +
        plasmix_theme +
        theme(panel.grid.major = element_blank(), axis.text.x = element_text(hjust = c(0.8, 0.8, 0.8, 0.8, 0.9)), plot.title = element_text(margin = margin(b = -2.5)))
    if(!y_label) p <- p + theme(axis.title.y = element_blank())
    p
}

plot_envelope_dia <- plot_response_envelope(response_curve_platform %>% filter(Platform == "DIA"),
                                             envelope_by_platform %>% filter(Platform == "DIA"),
                                             "dia", "DIA", platform_color["DIA"], y_label = TRUE)
plot_envelope_olk <- plot_response_envelope(response_curve_platform %>% filter(Platform == "OLK"),
                                             envelope_by_platform %>% filter(Platform == "OLK"),
                                             "olk", "OLK", platform_color["OLK"])
plot_envelope_som <- plot_response_envelope(response_curve_platform %>% filter(Platform == "SOM"),
                                             envelope_by_platform %>% filter(Platform == "SOM"),
                                             "som", "SOM", platform_color["SOM"])
print(plot_envelope_dia + plot_envelope_olk + plot_envelope_som + plot_layout(ncol = 3))

envelope_metrics_platform <- envelope_by_platform %>%
    mutate(`Intensity span` = D_log - A_log, `Abundance span` = Abundance_max - Abundance_min,
           `Response slope` = Response_slope, `Background burden` = Background_burden_pct) %>%
    select(Platform, `Intensity span`, `Abundance span`, `Response slope`, `Background burden`)

# Panel b: Batch-level response-envelope metrics ----
envelope_metrics_batch <- envelope_by_batch %>%
    mutate(`Intensity span` = D_log - A_log, `Abundance span` = Abundance_max - Abundance_min,
           `Response slope` = Response_slope, `Background burden` = Background_burden_pct) %>%
    select(Platform, Batch, `Intensity span`, `Abundance span`, `Response slope`, `Background burden`)

# Intensity span uses A and D, whereas background burden uses the inner-tail anchors A_bg and D_bg.
envelope_metrics_long <- envelope_metrics_batch %>%
    pivot_longer(cols = c(`Abundance span`, `Response slope`, `Background burden`), names_to = "Metric", values_to = "Value") %>%
    mutate(Batch = factor(Batch, levels = unique(Batch[order(Platform)])))

# Plot batch-level metrics.
plot_macro_bar <- function(data, metric_name, use_10_exp = FALSE, show_x_axis = FALSE, use_percent = FALSE) {
    p <- ggplot(data %>% filter(Metric == metric_name), aes(x = Batch, y = Value, fill = Platform)) +
        geom_bar(stat = "identity", width = 0.6) +
        scale_fill_manual(values = platform_color) +
        labs(x = NULL, y = metric_name) +
        plasmix_theme +
        theme(
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
            axis.line.x = element_blank(), axis.line.y = element_blank(),
            panel.grid.major.x = element_blank(),
            legend.position = "none",
            plot.margin = margin(t = 0, r = 5, b = 5, l = 5)
        )
    if (show_x_axis) {
        p <- p + theme(
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
            axis.ticks.x = element_line(color = "black")
        )
    } else {
        p <- p + theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        )
    }
    if (use_10_exp) {
        p <- p + scale_y_continuous(labels = format_10_exp, expand = expansion(mult = c(0, 0.15)), breaks = function(x) {b <- pretty(x, n = 4); b[b %% 1 == 0]})
    } else if (use_percent) {
        p <- p + scale_y_continuous(labels = scales::label_number(accuracy = 0.1, suffix = "%"), expand = expansion(mult = c(0, 0.15)),
                                    trans = "sqrt", breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))
    } else {
        p <- p + scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
    }
    return(p)
}

p_abundance_span <- plot_macro_bar(envelope_metrics_long, "Abundance span", use_10_exp = TRUE)
p_response_slope <- plot_macro_bar(envelope_metrics_long, "Response slope")
p_baseline_ratio <- plot_macro_bar(envelope_metrics_long, "Background burden", show_x_axis = TRUE, use_percent = TRUE) +
    theme(plot.margin = margin(t = 0, r = 5, b = 0, l = 5))
p_envelope_metrics <- (p_abundance_span / p_response_slope / p_baseline_ratio) + plot_layout(heights = c(1, 1, 1))
print(p_envelope_metrics)

# Simulate within-platform distortion and its cross-platform projection.
# Panel c: Stochastic vise plot ----
generate_vise_data <- function(envelope_df, group_col = "Platform") {
    set.seed(123)
    n_sims <- 1000
    matrix_load_seq <- c(0, 10^seq(-3, 3, length.out = 100))
    pmap_dfr(envelope_df, function(...) {
        row_data <- list(...)
        A <- 10^row_data$A_log; D <- 10^row_data$D_log
        abundance_mid <- 10^row_data$Abundance_mid_log
        response_slope <- row_data$Response_slope
        # Apply the input ratios from the same relative position within each platform envelope.
        abundance_1 <- abundance_mid / 10
        abundance_up <- abundance_1 * 2
        abundance_down <- abundance_1 / 2
        signal_1 <- A + (D - A) / (1 + (abundance_mid / abundance_1)^response_slope)
        signal_up <- A + (D - A) / (1 + (abundance_mid / abundance_up)^response_slope)
        signal_down <- A + (D - A) / (1 + (abundance_mid / abundance_down)^response_slope)
        map_dfr(matrix_load_seq, function(matrix_load) {
            # Convert the common latent matrix load into platform-specific signal interference.
            mean_interference <- matrix_load * (row_data$Background_burden_pct / 100) * (D - A)
            interference_1 <- abs(rnorm(n_sims, mean_interference, mean_interference * 0.30))
            interference_2 <- abs(rnorm(n_sims, mean_interference, mean_interference * 0.30))
            obs_up <- log2((signal_up + interference_2) / (signal_1 + interference_1))
            obs_down <- log2((signal_down + interference_2) / (signal_1 + interference_1))
            plot_x <- if(matrix_load == 0) -4 else log10(matrix_load)
            tibble(
                Platform = row_data$Platform, Group_Var = row_data[[group_col]],
                MatrixLoad = matrix_load, PlotX = plot_x,
                Up_Mean = mean(obs_up), Up_Lower = quantile(obs_up, 0.05), Up_Upper = quantile(obs_up, 0.95),
                Down_Mean = mean(obs_down), Down_Lower = quantile(obs_down, 0.05), Down_Upper = quantile(obs_down, 0.95)
            )
        })
    })
}
plot_vise <- function(vise_df, is_batch = FALSE) {
    vise_df$Platform <- factor(vise_df$Platform, levels = c("OLK", "SOM", "DIA"))
    native_df <- vise_df %>% filter(MatrixLoad == 0)
    p <- ggplot(vise_df, aes(x = PlotX)) +
        geom_hline(yintercept = c(-1, 1), color = "grey50", linetype = "dotted", linewidth = 0.5) +
        geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3) +
        geom_ribbon(aes(ymin = Up_Lower, ymax = Up_Upper, fill = Platform), alpha = 0.15) +
        geom_line(aes(y = Up_Mean, color = Platform), linewidth = 0.8) +
        geom_ribbon(aes(ymin = Down_Lower, ymax = Down_Upper, fill = Platform), alpha = 0.15) +
        geom_line(aes(y = Down_Mean, color = Platform), linewidth = 0.8) +
        geom_point(data = native_df, aes(y = Up_Mean, color = Platform), size = 1.7) +
        geom_point(data = native_df, aes(y = Down_Mean, color = Platform), size = 1.7) +
        scale_color_manual(values = platform_color) +
        scale_fill_manual(values = platform_color) +
        scale_x_continuous(
            limits = c(-4.5, 3),
            breaks = c(-4, -3, -2, -1, 0, 1, 2, 3),
            labels = expression(0, 10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3),
            expand = c(0, 0)
        ) +
        scale_y_continuous(limits = c(-2.2, 2), expand = c(0, 0)) +
        labs(x = "Added matrix load", y = expression(Observed~log[2]~ratio)) +
        annotate("text", x = -3.85, y = 0.8, label = "Input ratio = 2", hjust = 0, vjust = 0.5, size = 2.8) +
        annotate("text", x = -3.85, y = -0.8, label = "Input ratio = 0.5", hjust = 0, vjust = 0.5, size = 2.8) +
        plasmix_theme +
        theme(panel.grid.major = element_blank())
    if(is_batch) {
        p <- p + facet_wrap(~ Group_Var, ncol = 4) + theme(legend.position = "bottom")
    } else {
        p <- p + theme(legend.position = "right", legend.key.size = unit(0.4, "lines"), legend.margin = margin(5, 0, 5, -5))
    }
    return(p)
}
vise_df_plat <- generate_vise_data(envelope_by_platform, "Platform")
p_vise <- plot_vise(vise_df_plat)

# Panel d: Cross-platform distortion grid ----
build_grid_plot <- function(param_x, param_y, name_x = "SOM", name_y = "OLK") {
    true_lfc_levels <- c(-3, -2, -1, 0, 1, 2, 3)
    matrix_load_levels <- c(0, 1, 3, 10, 30)
    # Apply the same latent matrix load at matched relative positions within the two platform envelopes.
    grid_df <- expand_grid(TrueLFC = true_lfc_levels, MatrixLoad = matrix_load_levels) %>%
        mutate(
            PlatX_Obs = map2_dbl(TrueLFC, MatrixLoad, ~ {
                A <- 10^param_x$A_log; D <- 10^param_x$D_log
                abundance_mid <- 10^param_x$Abundance_mid_log; response_slope <- param_x$Response_slope
                abundance_1 <- abundance_mid / 10; abundance_2 <- abundance_1 * 2^.x
                signal_1 <- A + (D - A) / (1 + (abundance_mid / abundance_1)^response_slope)
                signal_2 <- A + (D - A) / (1 + (abundance_mid / abundance_2)^response_slope)
                interference <- .y * (param_x$Background_burden_pct / 100) * (D - A)
                log2((signal_2 + interference) / (signal_1 + interference))
            }),
            PlatY_Obs = map2_dbl(TrueLFC, MatrixLoad, ~ {
                A <- 10^param_y$A_log; D <- 10^param_y$D_log
                abundance_mid <- 10^param_y$Abundance_mid_log; response_slope <- param_y$Response_slope
                abundance_1 <- abundance_mid / 10; abundance_2 <- abundance_1 * 2^.x
                signal_1 <- A + (D - A) / (1 + (abundance_mid / abundance_1)^response_slope)
                signal_2 <- A + (D - A) / (1 + (abundance_mid / abundance_2)^response_slope)
                interference <- .y * (param_y$Background_burden_pct / 100) * (D - A)
                log2((signal_2 + interference) / (signal_1 + interference))
            }),
            TrueLFC_Group = factor(TrueLFC, levels = true_lfc_levels),
            MatrixLoad_Group = factor(MatrixLoad, levels = matrix_load_levels)
        )
    p_grid <- ggplot(grid_df, aes(x = PlatX_Obs, y = PlatY_Obs)) +
        geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
        geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
        geom_abline(slope = 1, intercept = 0, color = "grey40", linetype = "dotted", linewidth = 0.5) +
        geom_line(aes(group = TrueLFC), color = "grey40", linewidth = 0.5) +
        geom_line(aes(group = MatrixLoad), color = "grey60", linetype = "dotted", linewidth = 0.5) +
        geom_point(
            data = grid_df %>% arrange(desc(MatrixLoad)),
            aes(fill = TrueLFC_Group, size = MatrixLoad_Group),
            shape = 21, color = "black", stroke = 0.5
        ) +
        # scale_fill_brewer(palette = "RdBu", direction = -1, name = expression(bold(Input~log[2]*ratio))) +
        scale_fill_brewer(palette = "RdBu", direction = -1, name = "Input ratio", breaks = c("-3", "-2", "-1", "0", "1", "2", "3"), labels = c("1/8", "1/4", "1/2", "1", "2", "4", "8")) +
        scale_size_manual(values = c("0" = 1, "1" = 1.2, "3" = 1.6, "10" = 3, "30" = 5), breaks = c("0", "1", "3", "10", "30"), name = "Matrix load") +
        labs(x = bquote("Simulated"~log[2]~"ratio ("*.(name_x)*")"), y = bquote("Simulated"~log[2]~"ratio ("*.(name_y)*")")) +
        plasmix_theme +
        coord_cartesian(xlim = c(-3.5, 3), ylim = c(-6.5, 6), expand = c(0, 0)) +
        theme(
            panel.grid.major = element_blank(),
            legend.key.size = unit(0.4, "lines"),
            legend.position = "right",
            legend.margin = margin(5, 0, 5, -5),
            legend.box.margin = margin(b = -10)
        )
    return(p_grid)
}
param_som_plat <- envelope_by_platform %>% filter(Platform == "SOM") %>% slice(1)
param_olk_plat <- envelope_by_platform %>% filter(Platform == "OLK") %>% slice(1)
p_grid <- build_grid_plot(param_som_plat, param_olk_plat, "SOM", "OLK")

# Part 2: BLK/N perturbation in observed data ----
# Unified processing with well-level mapping and the 1.1-fold gate.
process_batch_compression <- function(batch_name, lod_df, k = 1.1) {
    # Features passing the release detection rule in this batch.
    valid_probes <- lod_df %>% filter(Batch == batch_name, IsAboveLoD == TRUE) %>% pull(UniqueID)
    # M, F and N signals on the HybNorm or Intensity scale.
    df_std <- long_df %>%
        filter(
            Batch == batch_name,
            UniqueID %in% valid_probes,
            Sample %in% c("M", "F", "N"),
            (Platform %in% c("SOM", "OLK") & ProcessLevel == "HybNorm") | (Platform %in% c("DIA") & ProcessLevel == "Intensity")
        ) %>%
        group_by(UniqueID, UniProtID, Sample) %>%
        summarize(Val_Lin = median(2^Value, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = Sample, values_from = Val_Lin, names_prefix = "Std_")
    # Add absent sample columns when needed.
    for (col in c("Std_M", "Std_F", "Std_N")) {
        if (!col %in% colnames(df_std)) df_std[[col]] <- NA_real_
    }
    # Require both M and F measurements.
    df_std <- df_std %>% filter(!is.na(Std_M) & !is.na(Std_F))
    # BLK subtraction and physical gating for affinity platforms.
    plat <- unique(meta_batch_ht$Platform[meta_batch_ht$Batch == batch_name])
    if(plat %in% c("SOM", "OLK")) {
        # Batch-level median raw BLK.
        blk_raw <- long_df %>%
            filter(Batch == batch_name, UniqueID %in% valid_probes, Sample == "BLK", ProcessLevel == "Raw") %>%
            group_by(UniqueID, UniProtID) %>%
            summarize(Med_BLK_Raw = median(2^Value, na.rm = TRUE), .groups = "drop")
        # Map raw and HybNorm measurements by well before gating.
        df_sub <- long_df %>%
            filter(Batch == batch_name, UniqueID %in% valid_probes, Sample %in% c("M", "F"), ProcessLevel %in% c("Raw", "HybNorm")) %>%
            mutate(Lin_Val = 2^Value) %>%
            select(UniqueID, UniProtID, Sample, ColName, ProcessLevel, Lin_Val) %>%
            pivot_wider(names_from = ProcessLevel, values_from = Lin_Val) %>%
            left_join(blk_raw) %>%
            mutate(
                Med_BLK_Raw = coalesce(Med_BLK_Raw, 0),
                # Require the raw well signal to exceed 1.1 times the BLK median.
                Valid_BLK_Well = (Raw > k * Med_BLK_Raw),
                Scale = ifelse(Raw > 0, HybNorm / Raw, 1),
                Sub_Raw = pmax(Raw - Med_BLK_Raw, 1e-5),
                Sub_Hyb = Sub_Raw * Scale
            ) %>%
            group_by(UniqueID, UniProtID, Sample) %>%
            summarize(
                Sub_Lin = median(Sub_Hyb, na.rm = TRUE),
                # Require all M and F replicates to pass the gate.
                Valid_BLK_Sample = all(!is.na(Valid_BLK_Well) & Valid_BLK_Well),
                .groups = "drop"
            ) %>%
            pivot_wider(names_from = Sample, values_from = c(Sub_Lin, Valid_BLK_Sample))
    } else {
        # DIA has no BLK measurement.
        df_sub <- df_std %>%
            select(UniqueID, UniProtID) %>%
            mutate(Sub_Lin_M = NA_real_, Sub_Lin_F = NA_real_, Valid_BLK_Sample_M = FALSE, Valid_BLK_Sample_F = FALSE)
    }
    # Merge signals and calculate perturbation ratios.
    df_calc <- df_std %>%
        left_join(df_sub) %>%
        mutate(
            Batch = batch_name,
            # Native M/F ratio.
            Ratio_Raw = log2(Std_M) - log2(Std_F),
            # BLK-subtracted M/F ratio.
            Valid_BLK = Valid_BLK_Sample_M & Valid_BLK_Sample_F,
            Ratio_BLK = if_else(Valid_BLK & !is.na(Sub_Lin_M) & !is.na(Sub_Lin_F), log2(Sub_Lin_M) - log2(Sub_Lin_F), NA_real_),
            # N-subtracted M/F ratio.
            Valid_N = (Std_M > k * Std_N) & (Std_F > k * Std_N),
            Ratio_N = if_else(Valid_N & !is.na(Std_N), log2(Std_M - Std_N) - log2(Std_F - Std_N), NA_real_),
            # Relative expansion after N subtraction.
            Ratio_Raw_OK = abs(Ratio_Raw) >= 0.01, # Exclude near-zero native ratios for which expansion is undefined.
            Y_RelExpansion = if_else(
                Valid_N & !is.na(Ratio_N) & Ratio_Raw_OK,
                Ratio_N / Ratio_Raw,
                NA_real_
            )
        )
    return(df_calc)
}

# Process all batches into one table.
df_unified_all <- map_dfr(unique(meta_batch_ht$Batch), function(b) {
    process_batch_compression(b, lod_df = lod_status, k = 1.1)
})

# QC1: Near-zero native M/F ratios.
quantile(abs(df_unified_all$Ratio_Raw), c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.99, 1), na.rm=TRUE)
#         0%          1%          5%         10%         25%         50%         75%         99%        100%
# 0.000000000 0.002295175 0.010959420 0.022154614 0.057365031 0.131979837 0.279302657 2.281298547 9.907053583

# QC2: Sign reversals after perturbation.
df_unified_all %>%
  filter(Valid_N, !is.na(Ratio_Raw), !is.na(Ratio_N)) %>%
  summarize(n = n(), n_flip = sum(sign(Ratio_Raw) != sign(Ratio_N)), frac_flip = mean(sign(Ratio_Raw) != sign(Ratio_N)))
#       n n_flip frac_flip
# 1 10277      0         0

# QC3: Distribution of the raw model outcome.
quantile(df_unified_all$Y_RelExpansion, c(0, 0.01, 0.25, 0.5, 0.75, 0.99, 1), na.rm=TRUE)
#       0%        1%       25%       50%       75%       99%      100%
# 1.002090  1.141996  2.193674  3.332472  4.911410  8.942242 10.414996

# Batch- and platform-level expansion for batches containing P, N, M and F.
expansion_input <- df_unified_all %>%
  mutate(Platform = sub("_.*", "", Batch)) %>%
  filter(
    !is.na(Ratio_Raw), Valid_N, !is.na(Ratio_N),
    Platform == "DIA" | (Valid_BLK & !is.na(Ratio_BLK))
  )
batch_expansion <- expansion_input %>%
    group_by(Batch) %>%
    summarize(
        Features_Count = n(),
        Raw_d25 = quantile(Ratio_Raw, 0.025, na.rm = TRUE),
        Raw_d975 = quantile(Ratio_Raw, 0.975, na.rm = TRUE),
        BLK_d25 = quantile(Ratio_BLK, 0.025, na.rm = TRUE),
        BLK_d975 = quantile(Ratio_BLK, 0.975, na.rm = TRUE),
        N_d25 = quantile(Ratio_N, 0.025, na.rm = TRUE),
        N_d975 = quantile(Ratio_N, 0.975, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        Delta95_Raw = Raw_d975 - Raw_d25,
        Delta95_BLK = BLK_d975 - BLK_d25,
        Delta95_N = N_d975 - N_d25,
        Expansion_BLK = Delta95_BLK / Delta95_Raw,
        Expansion_N = Delta95_N / Delta95_Raw
    ) %>%
    arrange(desc(Expansion_N)) %>%
    select(Batch, Features_Count, Delta95_Raw, Delta95_BLK, Delta95_N, Expansion_BLK, Expansion_N)

platform_expansion <- expansion_input %>%
    group_by(Platform) %>%
    summarize(
        Features_Count = n(),
        Raw_d25 = quantile(Ratio_Raw, 0.025, na.rm = TRUE),
        Raw_d975 = quantile(Ratio_Raw, 0.975, na.rm = TRUE),
        BLK_d25 = quantile(Ratio_BLK, 0.025, na.rm = TRUE),
        BLK_d975 = quantile(Ratio_BLK, 0.975, na.rm = TRUE),
        N_d25 = quantile(Ratio_N, 0.025, na.rm = TRUE),
        N_d975 = quantile(Ratio_N, 0.975, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        Delta95_Raw = Raw_d975 - Raw_d25,
        Delta95_BLK = BLK_d975 - BLK_d25,
        Delta95_N = N_d975 - N_d25,
        Expansion_BLK = Delta95_BLK / Delta95_Raw,
        Expansion_N = Delta95_N / Delta95_Raw
    ) %>%
    arrange(desc(Expansion_N)) %>%
    select(Platform, Features_Count, Delta95_Raw, Delta95_BLK, Delta95_N, Expansion_BLK, Expansion_N)

# Panel e: Bowknot plots ----
# Expansion fold = perturbed 95% range / native 95% range.
plot_bowtie_from_unified_merged <- function(df_model, target_platform, bs_type, pt_color) {
    bs_col <- sym(paste0("Ratio_", bs_type))
    # Select the requested platform and perturbation.
    res_df <- df_model %>%
        mutate(Platform_tmp = sub("_.*", "", Batch)) %>%
        filter(Platform_tmp == target_platform) %>%
        filter(!is.na(Ratio_Raw), !is.na(!!bs_col)) %>%
        select(Ratio_Orig = Ratio_Raw, Ratio_BS = !!bs_col)
    if(nrow(res_df) == 0) return(ggplot() + theme_void() + ggtitle(paste(target_platform, "- No Data")))
    # Distribution summaries.
    st <- list(
        x_q025 = quantile(res_df$Ratio_Orig, 0.025, na.rm = TRUE),
        x_q975 = quantile(res_df$Ratio_Orig, 0.975, na.rm = TRUE),
        y_q025 = quantile(res_df$Ratio_BS, 0.025, na.rm = TRUE),
        y_q975 = quantile(res_df$Ratio_BS, 0.975, na.rm = TRUE),
        color = pt_color)
    st$d95_orig <- st$x_q975 - st$x_q025
    st$d95_bs   <- st$y_q975 - st$y_q025
    st$expansion <- st$d95_bs / st$d95_orig
    l1 <- bquote(Delta[95]*"Native" == .(sprintf("%.2f", st$d95_orig)))
    l2 <- bquote(Delta[95]*"Perturbed" == .(sprintf("%.2f", st$d95_bs)))
    l3 <- bquote(Expansion == .(sprintf("%.2f", st$expansion))*"×")
    plot_title <- sprintf("%s (%s-perturbed)", target_platform, bs_type)
    axis_min <- -6.04; axis_max <- 6.04
    # Bow-tie envelope based on the estimated expansion.
    ribbon_df <- data.frame(x_ribbon = seq(axis_min, axis_max, length.out = 100))
    # Identity line.
    ribbon_df$y1 <- ribbon_df$x_ribbon
    # Expansion boundary.
    ribbon_df$y2 <- ribbon_df$x_ribbon * st$expansion
    # Use pmin/pmax because ordering reverses for negative x.
    ribbon_df$ymin <- pmin(ribbon_df$y1, ribbon_df$y2)
    ribbon_df$ymax <- pmax(ribbon_df$y1, ribbon_df$y2)
    # Plot.
    p <- ggplot(res_df, aes(x = Ratio_Orig, y = Ratio_BS)) +
        geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
        geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
        geom_ribbon(data = ribbon_df, aes(x = x_ribbon, ymin = ymin, ymax = ymax), fill = st$color, alpha = 0.05, inherit.aes = FALSE) +
        geom_abline(slope = 1, intercept = 0, color = "grey50", linetype = "dotted", linewidth = 0.5) +
        geom_abline(slope = st$expansion, intercept = 0, color = st$color, linetype = "dotted", linewidth = 0.5) +
        annotate("segment", x = st$x_q025, xend = st$x_q975, y = -Inf, yend = -Inf, color = "grey50", linewidth = 2, alpha = 0.5) +
        annotate("segment", x = -Inf, xend = -Inf, y = st$y_q025, yend = st$y_q975, color = st$color, linewidth = 2, alpha = 0.5) +
        geom_point(alpha = 0.3, size = 0.5, color = st$color) +
        annotate("text", x = axis_min + 0.5, y = axis_max - 0.5, label = deparse(l1), hjust = 0, vjust = 1, parse = TRUE, size = 2.5) +
        annotate("text", x = axis_min + 0.5, y = axis_max - 1.5, label = deparse(l2), hjust = 0, vjust = 1, parse = TRUE, size = 2.5) +
        annotate("text", x = axis_min + 0.5, y = axis_max - 2.5, label = deparse(l3), hjust = 0, vjust = 1, parse = TRUE, size = 2.5) +
        labs(x = NULL, y = NULL, title = plot_title) +
        plasmix_theme +
        coord_cartesian(
            xlim = c(axis_min, axis_max), ylim = c(axis_min, axis_max),
            expand = FALSE,
            clip = "on"
        ) +
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35*2),
              panel.grid.major = element_blank(),
              axis.line.x = element_blank(), axis.line.y = element_blank(),
              plot.title = element_text(size = 8, face = "bold", hjust = 0.5))
    return(p)
}
x_expr <- expression(Native~log[2]*(M/F))
y_expr <- expression(Perturbed~log[2]*(M/F))

p_som_blk <- plot_bowtie_from_unified_merged(df_unified_all %>% filter(Valid_N == TRUE, Valid_BLK == TRUE), target_platform = "SOM", bs_type = "BLK", pt_color = "#B33E90")
p_som_n   <- plot_bowtie_from_unified_merged(df_unified_all %>% filter(Valid_N == TRUE, Valid_BLK == TRUE), target_platform = "SOM", bs_type = "N", pt_color = "#B33E90")
p_olk_blk   <- plot_bowtie_from_unified_merged(df_unified_all %>% filter(Valid_N == TRUE, Valid_BLK == TRUE), target_platform = "OLK", bs_type = "BLK", pt_color = "#489FA7")
p_olk_n   <- plot_bowtie_from_unified_merged(df_unified_all %>% filter(Valid_N == TRUE, Valid_BLK == TRUE), target_platform = "OLK", bs_type = "N", pt_color = "#489FA7")
p_dia_n   <- plot_bowtie_from_unified_merged(df_unified_all, target_platform = "DIA", bs_type = "N", pt_color = "#155289")

# Align the five square scatter panels.
p_row <- (p_som_blk + labs(y = y_expr)) + p_som_n + p_olk_blk + p_olk_n + p_dia_n + plot_layout(nrow = 1)
pb_scatter <- ggarrange(p_row,
                        textGrob(x_expr, vjust = -0.26, gp = gpar(fontsize = 9, color = "black")),
                        nrow = 2, heights = c(25, 1)) + theme(plot.margin = margin(-5, 0, 0, -5))

# Physicochemical modeling of perturbation susceptibility ----
platforms_to_run <- c("SOM", "OLK", "DIA")
# Join physicochemical annotations.
df_model <- df_unified_all %>% left_join(physchem_matrix, by = c("UniProtID" = "Entry")) %>% filter(Valid_N == TRUE) %>%
    mutate(Y_Log2Expansion = if_else(Y_RelExpansion > 0, log2(Y_RelExpansion), NA_real_))
print(dim(df_model))

# Platform-wide and platform-shared datasets.
df_full <- df_model %>% mutate(Platform = sub("_.*", "", Batch)) %>% filter(Platform %in% platforms_to_run, is.finite(Y_Log2Expansion))

# Shared UniProt targets.
common_prots <- df_full %>%
    group_by(Platform) %>%
    summarize(UP = list(unique(UniProtID))) %>%
    pull(UP) %>% Reduce(intersect, .)
length(common_prots)
df_intersect <- df_full %>% filter(UniProtID %in% common_prots)
datasets <- list(Full = df_full, Intersect = df_intersect)

# NA impact check before ML
na_by_feature <- imap_dfr(datasets, function(df_current, ds_name) {
  map_dfr(platforms_to_run, function(plat) {
    df_sub <- df_current %>% filter(Platform == plat)
    # Match the missingness denominator to rows eligible for modeling.
    df_sub <- df_sub[!is.na(df_sub[["Y_RelExpansion"]]), , drop = FALSE]
    n <- nrow(df_sub)
    map_dfr(final_physchem, function(f) {
      v <- df_sub[[f]]
      tibble(
        Dataset = ds_name, Platform = plat, Feature = f,
        N_model_rows = n,
        N_missing = sum(is.na(v)),
        Missing_Frac = sum(is.na(v)) / n,
        # Fraction assigned a single imputed value in RF.
        Imputed_To_Single_Value_Frac = sum(is.na(v)) / n
      )
    })
  })
}) %>% arrange(desc(Missing_Frac))
write.csv(na_by_feature, "results/fig4_na_impact_before_ML.csv", row.names = FALSE)

# Informative-missingness check ----
# Test whether property missingness itself is associated with raw Y_RelExpansion within each platform.
informative_na_check <- function(df, platform, y = "Y_RelExpansion", feats = final_physchem) {
    d <- df %>% filter(Platform == platform, !is.na(.data[[y]]))
    if (nrow(d) < 50) return(NULL)
    yv <- d[[y]]
    out <- map_dfr(feats, function(f) {
        if (!f %in% colnames(d)) return(NULL)
        v <- d[[f]]
        n_na <- sum(is.na(v))
        # No test is needed for complete properties.
        if (n_na == 0) {
            return(tibble(Feature = f, Missing_Frac = 0,
                          Y_mean_missing = NA_real_, Y_mean_present = mean(yv, na.rm = TRUE),
                          Cliffs_delta = NA_real_, p_raw = NA_real_, Test = "none"))
        }
        is_na <- is.na(v)
        # Skip tests with fewer than ten rows in either group.
        if (sum(is_na) < 10 || sum(!is_na) < 10) {
            return(tibble(Feature = f, Missing_Frac = mean(is_na),
                          Y_mean_missing = mean(yv[is_na], na.rm = TRUE),
                          Y_mean_present = mean(yv[!is_na], na.rm = TRUE),
                          Cliffs_delta = NA_real_, p_raw = NA_real_, Test = "too_few"))
        }
        y_miss <- yv[is_na]; y_pres <- yv[!is_na]
        # Wilcoxon test for the skewed raw outcome.
        p <- tryCatch(wilcox.test(y_miss, y_pres)$p.value, error = function(e) NA_real_)
        # Cliff's delta: <0.15 negligible, 0.15-0.33 small, 0.33-0.47 moderate, >0.47 large.
        delta <- tryCatch({
            m <- outer(y_miss, y_pres, FUN = function(a, b) sign(a - b))
            mean(m)
        }, error = function(e) NA_real_)
        tibble(
            Feature = f,
            Missing_Frac   = mean(is_na),
            Y_mean_missing = mean(y_miss, na.rm = TRUE),
            Y_mean_present = mean(y_pres, na.rm = TRUE),
            Cliffs_delta   = delta,
            p_raw          = p,
            Test           = "wilcoxon"
        )
    })
    out %>%
        # BH adjustment within platform.
        mutate(p_adj = p.adjust(p_raw, method = "BH")) %>%
        arrange(p_adj, desc(abs(Cliffs_delta)))
}

# Run separately by platform using the raw outcome.
na_info <- map_dfr(c("SOM", "OLK", "DIA"), function(p) {
    res <- informative_na_check(df_full, platform = p)
    if (!is.null(res)) res$Platform <- p
    res
})

# Report properties with missing values.
na_info_flagged <- na_info %>%
    filter(Missing_Frac > 0, Test == "wilcoxon") %>%
    mutate(
        Informative = case_when(
            is.na(p_adj)                             ~ "untested",
            p_adj < 0.05 & abs(Cliffs_delta) >= 0.15 ~ "YES (effect non-trivial)",
            p_adj < 0.05 & abs(Cliffs_delta) <  0.15 ~ "sig but tiny effect",
            TRUE                                     ~ "no"
        )
    ) %>%
    select(Platform, Feature, Missing_Frac, Cliffs_delta, p_adj, Informative, Y_mean_missing, Y_mean_present) %>%
    arrange(Platform, desc(Missing_Frac))

print(as.data.frame(na_info_flagged), digits = 3)

# Model training ----
ml_physchem <- setdiff(final_physchem, c("Log10_Abundance", "Helix_Fraction", "Beta_Fraction"))
model_physchem_names <- make.names(ml_physchem, unique = TRUE)
physchem_name_map <- c(setNames(ml_physchem, model_physchem_names), Random_Noise = "Random_Noise")
physchem_model_map <- setNames(model_physchem_names, ml_physchem)
rf_mtry <- min(5, length(ml_physchem)); rf_nodesize <- 5
k_folds <- 5; n_repeats <- 50

xgb_params <- list(objective = "reg:squarederror", eta = 0.05, max_depth = 3, min_child_weight = 3,
                   colsample_bytree = 0.8, subsample = 0.8)
lgb_params <- list(objective = "regression", metric = "rmse", learning_rate = 0.05, max_depth = 3,
                   num_leaves = 7, min_data_in_leaf = 5, feature_fraction = 0.8, bagging_fraction = 0.8,
                   bagging_freq = 1, use_missing = TRUE, zero_as_missing = FALSE, verbosity = -1)

calc_rf_impute_values <- function(x_train) {
    apply(x_train, 2, function(v) {
        v <- v[is.finite(v)]
        if(!length(v)) return(0)
        if(all(v %in% c(0, 1))) return(as.numeric(mean(v) >= 0.5))
        median(v)
    })
}

apply_rf_imputation <- function(x, impute_values) {
    x_imp <- x
    for(j in seq_along(impute_values)) x_imp[!is.finite(x_imp[, j]), j] <- impute_values[j]
    x_imp
}

make_predictor_matrix <- function(df, seed) {
    x_phys <- as.data.frame(lapply(df[, ml_physchem, drop = FALSE], function(x) {
        if(is.factor(x)) x <- as.character(x)
        suppressWarnings(as.numeric(x))
    }))
    names(x_phys) <- model_physchem_names
    set.seed(seed)
    x <- cbind(as.matrix(x_phys), Random_Noise = rnorm(nrow(df))); storage.mode(x) <- "numeric"
    x
}

# Balanced, reproducible row-wise folds shared by all three algorithms.
make_row_folds <- function(n, k = 5, seed = 1) {
    if(n < k) stop("The number of feature-batch observations is smaller than k_folds.", call. = FALSE)
    set.seed(seed)
    sample(rep(seq_len(k), length.out = n))
}

# Primary estimand: prediction across observed feature–batch contexts.
calc_r2 <- function(observed, predicted) {
    valid <- is.finite(observed) & is.finite(predicted)
    observed <- observed[valid]; predicted <- predicted[valid]
    denominator <- sum((observed - mean(observed))^2)
    if(!is.finite(denominator) || denominator <= 0) return(NA_real_)
    100 * (1 - sum((observed - predicted)^2) / denominator)
}

normalize_importance <- function(x) {
    x <- pmax(as.numeric(x), 0); maximum <- max(x, na.rm = TRUE)
    if(is.finite(maximum) && maximum > 0) x / maximum * 100 else rep(0, length(x))
}

format_importance <- function(raw_scores, algorithm, dataset, platform) {
    predictor_names <- c(model_physchem_names, "Random_Noise")
    scores <- setNames(rep(0, length(predictor_names)), predictor_names)
    shared <- intersect(names(raw_scores), predictor_names)
    scores[shared] <- raw_scores[shared]
    tibble(Dataset = dataset, Platform = platform, Algorithm = algorithm,
           Feature = unname(physchem_name_map[names(scores)]), Relative_Imp = normalize_importance(scores), Raw_score = as.numeric(scores))
}

model_input_rows <- imap_dfr(datasets, function(df_current, ds_name) {
    df_current %>% filter(Platform %in% platforms_to_run, !is.na(UniProtID), is.finite(Y_Log2Expansion)) %>%
        mutate(Dataset = ds_name, .before = 1) %>%
        select(Dataset, Platform, Batch, UniqueID, UniProtID, Y_RelExpansion, Y_Log2Expansion, all_of(ml_physchem))
})
fwrite(as.data.table(model_input_rows), "results/fig4_rowwise_model_input.tsv.gz", sep = "\t", na = "NA")

table_model_list <- list(); table_vip_list <- list(); global_models <- list(); cv_fold_list <- list(); cv_metric_list <- list()
total_steps <- length(datasets) * length(platforms_to_run) * n_repeats
pb <- txtProgressBar(min = 0, max = total_steps, style = 3); step_current <- 0

for(ds_name in names(datasets)) {
    global_models[[ds_name]] <- list()
    for(plat in platforms_to_run) {
        df_train <- datasets[[ds_name]] %>%
            filter(Platform == plat, !is.na(UniProtID), is.finite(Y_Log2Expansion)) %>% arrange(Batch, UniqueID, UniProtID)
        if(nrow(df_train) < k_folds) next
        model_seed <- 100000 + match(ds_name, names(datasets)) * 10000 + match(plat, platforms_to_run) * 1000
        if(any(df_train$Y_RelExpansion <= 0)) stop("Y_RelExpansion must be positive.")
        x_matrix <- make_predictor_matrix(df_train, model_seed); y_vector <- df_train$Y_Log2Expansion
        n_complete_feature_obs <- sum(complete.cases(df_train[, ml_physchem, drop = FALSE]))
        n_complete_feature_proteins <- n_distinct(df_train$UniProtID[complete.cases(df_train[, ml_physchem, drop = FALSE])])

        for(r in seq_len(n_repeats)) {
            fold_ids <- make_row_folds(nrow(df_train), k = k_folds, seed = model_seed + r)
            cv_fold_list[[length(cv_fold_list) + 1]] <- tibble(
                Dataset = ds_name, Platform = plat, Repeat = r, Row = seq_len(nrow(df_train)), Fold = fold_ids,
                Batch = df_train$Batch, UniqueID = df_train$UniqueID, UniProtID = df_train$UniProtID)
            pred_rf <- pred_xgb <- pred_lgb <- rep(NA_real_, length(y_vector))
            for(fold in seq_len(k_folds)) {
                idx_test <- which(fold_ids == fold); idx_train <- which(fold_ids != fold)
                # Estimate RF imputation values from the training fold only.
                rf_impute_values <- calc_rf_impute_values(x_matrix[idx_train, , drop = FALSE])
                x_train_rf <- apply_rf_imputation(x_matrix[idx_train, , drop = FALSE], rf_impute_values)
                x_test_rf <- apply_rf_imputation(x_matrix[idx_test, , drop = FALSE], rf_impute_values)
                fold_seed <- model_seed + r * 100 + fold
                set.seed(fold_seed)
                rf_model <- randomForest::randomForest(x = x_train_rf, y = y_vector[idx_train], ntree = 500,
                                                       mtry = rf_mtry, nodesize = rf_nodesize)
                pred_rf[idx_test] <- predict(rf_model, x_test_rf)

                xgb_params_fold <- xgb_params; xgb_params_fold$seed <- fold_seed + 1000000
                dtrain_xgb <- xgb.DMatrix(x_matrix[idx_train, , drop = FALSE], label = y_vector[idx_train], missing = NA)
                xgb_model <- xgb.train(params = xgb_params_fold, data = dtrain_xgb, nrounds = 100, verbose = 0)
                pred_xgb[idx_test] <- predict(xgb_model, xgb.DMatrix(x_matrix[idx_test, , drop = FALSE], missing = NA))

                lgb_params_fold <- lgb_params
                lgb_params_fold$seed <- fold_seed + 2000000
                lgb_params_fold$bagging_seed <- fold_seed + 2000001
                lgb_params_fold$feature_fraction_seed <- fold_seed + 2000002
                dtrain_lgb <- lgb.Dataset(x_matrix[idx_train, , drop = FALSE], label = y_vector[idx_train])
                lgb_model <- lgb.train(params = lgb_params_fold, data = dtrain_lgb, nrounds = 100, verbose = -1)
                pred_lgb[idx_test] <- predict(lgb_model, x_matrix[idx_test, , drop = FALSE])
            }
            repeat_metrics <- tibble(Algorithm = c("RandomForest", "XGBoost", "LightGBM"),
                                     R_Squared = c(calc_r2(y_vector, pred_rf), calc_r2(y_vector, pred_xgb), calc_r2(y_vector, pred_lgb))) %>%
                mutate(Dataset = ds_name, Platform = plat, Repeat = r, .before = 1)
            cv_metric_list[[length(cv_metric_list) + 1]] <- repeat_metrics
            step_current <- step_current + 1; setTxtProgressBar(pb, step_current)
        }

        platform_metrics <- bind_rows(cv_metric_list) %>% filter(Dataset == ds_name, Platform == plat)
        table_model_list[[length(table_model_list) + 1]] <- platform_metrics %>% group_by(Algorithm) %>%
            summarize(R_Squared_Mean = mean(R_Squared, na.rm = TRUE), R_Squared_SD = sd(R_Squared, na.rm = TRUE), .groups = "drop") %>%
            mutate(Platform = plat, Dataset = ds_name, CV_Unit = "Feature-batch row", Estimand = "Observed feature-batch context",
                   Target_Features = n_distinct(df_train$UniqueID),
                   Target_Proteins = n_distinct(df_train$UniProtID),
                   Total_Observations = nrow(df_train),
                   Complete_Feature_Proteins = n_complete_feature_proteins,
                   Complete_Feature_Observations = n_complete_feature_obs, .before = 1)

        # Full-data models are used only for importance and partial dependence; predictive performance remains cross-validated.
        rf_impute_values_full <- calc_rf_impute_values(x_matrix)
        x_matrix_rf_full <- apply_rf_imputation(x_matrix, rf_impute_values_full)
        set.seed(model_seed + 3000000)
        rf_final <- randomForest::randomForest(x = x_matrix_rf_full, y = y_vector, ntree = 500,
                                               mtry = rf_mtry, nodesize = rf_nodesize, importance = TRUE)
        xgb_params_final <- xgb_params; xgb_params_final$seed <- model_seed + 4000000
        xgb_final <- xgb.train(params = xgb_params_final,
            data = xgb.DMatrix(x_matrix, label = y_vector, missing = NA), nrounds = 100, verbose = 0)
        lgb_params_final <- lgb_params
        lgb_params_final$seed <- model_seed + 5000000
        lgb_params_final$bagging_seed <- model_seed + 5000001
        lgb_params_final$feature_fraction_seed <- model_seed + 5000002
        lgb_final <- lgb.train(params = lgb_params_final,
            data = lgb.Dataset(x_matrix, label = y_vector), nrounds = 100, verbose = -1)

        global_models[[ds_name]][[plat]] <- list(RandomForest = rf_final, XGBoost = xgb_final, LightGBM = lgb_final,
            X_train_raw = x_matrix, X_train_rf = x_matrix_rf_full, RF_Impute_Values = rf_impute_values_full,
            Y_train = y_vector, Outcome = "log2(Y_RelExpansion)",
            Protein_ID = df_train$UniProtID, Batch = df_train$Batch, CV_Unit = "Feature-batch row")

        rf_importance <- randomForest::importance(rf_final)[, "%IncMSE"]
        xgb_importance <- xgb.importance(model = xgb_final); xgb_importance <- setNames(xgb_importance$Gain, xgb_importance$Feature)
        lgb_importance <- lgb.importance(lgb_final, percentage = TRUE); lgb_importance <- setNames(lgb_importance$Gain, lgb_importance$Feature)
        table_vip_list[[length(table_vip_list) + 1]] <- bind_rows(
            format_importance(rf_importance, "RandomForest", ds_name, plat),
            format_importance(xgb_importance, "XGBoost", ds_name, plat),
            format_importance(lgb_importance, "LightGBM", ds_name, plat))
    }
}
close(pb)

cv_fold_assignments <- bind_rows(cv_fold_list); cv_repeat_metrics <- bind_rows(cv_metric_list)
fwrite(as.data.table(cv_fold_assignments), "results/fig4_rowwise_cv_folds.tsv.gz", sep = "\t", na = "NA")
fwrite(as.data.table(cv_repeat_metrics), "results/fig4_rowwise_cv_repeats.tsv.gz", sep = "\t", na = "NA")
model_support <- bind_rows(table_model_list) %>%
    select(Dataset, Platform, Algorithm, R_Squared_Mean, R_Squared_SD)
print(as.data.frame(model_support), digits = 3)
fwrite(as.data.table(model_support), "results/fig4_rowwise_model_support.tsv", sep = "\t", na = "NA")

# Panel f: Cross-validated model performance ----
df_panel_c <- bind_rows(table_model_list) %>%
    mutate(
        Platform = factor(Platform, levels = c("DIA", "OLK", "SOM")),
        Dataset = factor(Dataset, levels = c("Full", "Intersect")),
        Algorithm = factor(Algorithm, levels = c("RandomForest", "XGBoost", "LightGBM"))
    )
pc_bar <- ggplot(df_panel_c, aes(y = Platform, fill = Platform, alpha = Algorithm, group = Algorithm)) +
    geom_col(aes(x = R_Squared_Mean), position = position_dodge(width = 0.9), width = 0.85) +
    geom_errorbar(aes(xmin = R_Squared_Mean - R_Squared_SD, xmax = R_Squared_Mean + R_Squared_SD),
                  position = position_dodge(width = 0.9), width = 0.25, color = "black", alpha = 0.8, linewidth = 0.3) +
    facet_wrap(~ Dataset, ncol = 1, strip.position = "right", labeller = labeller(Dataset = c(Full = "Platform-wide", Intersect = "Platform-shared"))) +
    scale_fill_manual(values = platform_color, guide = "none") +
    scale_alpha_manual(values = c("RandomForest" = 1, "XGBoost" = 0.6, "LightGBM" = 0.3)) +
    labs(x = bquote("5-fold cross-validation model " ~ R^2 ~ "(%)"), y = NULL) +
    plasmix_theme +
    coord_cartesian(clip = "off") +
    guides(alpha = guide_legend(override.aes = list(label = ""))) +
    scale_x_continuous(limits = c(-2, 68), expand = c(0, 0)) +
    theme(panel.grid.major = element_blank(), legend.position = "bottom",
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(-7, 0, 0, 30))

# Panel g: Feature VIP lollipop chart (Top 5 union, dodged) ----
# Consensus relative importance across algorithms.
df_vip_consensus <- bind_rows(table_vip_list) %>%
    filter(Feature != "Random_Noise", Platform %in% c("SOM", "OLK", "DIA")) %>%
    group_by(Dataset, Platform, Feature) %>%
    summarize(Consensus_Score = mean(Relative_Imp, na.rm = TRUE), .groups = "drop")

# Rank properties in the platform-wide or platform-shared dataset.
df_wide_g <- df_vip_consensus %>% filter(Dataset == "Full") %>% group_by(Platform) %>% arrange(desc(Consensus_Score)) %>% mutate(Rank = row_number()) %>% ungroup()
df_shared_g <- df_vip_consensus %>% filter(Dataset == "Intersect") %>% group_by(Platform) %>% arrange(desc(Consensus_Score)) %>% mutate(Rank = row_number()) %>% ungroup()

# Define displayed properties using the platform-wide Top 5 only.
features_to_display <- df_wide_g %>%
    filter(Rank <= 5) %>%
    distinct(Platform, Feature) %>%
    pull(Feature) %>%
    unique()

# Order displayed properties by category and platform-wide importance.
feature_weights <- df_wide_g %>%
    filter(Feature %in% features_to_display) %>%
    group_by(Feature) %>%
    summarize(
        N_Platform_Top5 = sum(Rank <= 5),
        Mean_Wide_Score = mean(Consensus_Score, na.rm = TRUE),
        Global_Max_Score = max(Consensus_Score, na.rm = TRUE),
        .groups = "drop"
    )

# Final top-to-bottom property order.
unified_order <- tibble(Feature = features_to_display) %>%
    left_join(feature_weights, by = "Feature") %>%
    mutate(
        Category = factor(cat_map[Feature], levels = rev(names(cat_colors))),
        Property = unname(name_map[Feature])
    ) %>%
    arrange(Category, desc(N_Platform_Top5), desc(Mean_Wide_Score), desc(Global_Max_Score), Property) %>%
    pull(Property)

# Combine scopes and offset their plotting positions.
df_combined_g <- bind_rows(
    df_wide_g %>% mutate(Dataset_Type = "Platform-wide"),
    df_shared_g %>% mutate(Dataset_Type = "Platform-shared")
) %>%
    filter(Feature %in% features_to_display) %>%
    mutate(
        Property = factor(name_map[Feature], levels = rev(unified_order)),
        Category = factor(cat_map[Feature], levels = names(cat_colors)),
        Platform = factor(Platform, levels = c("SOM", "OLK", "DIA")),
        Dataset_Type = factor(Dataset_Type, levels = c("Platform-wide", "Platform-shared")),
        Y_Base = as.numeric(Property),
        Y_Plot = ifelse(Dataset_Type == "Platform-wide", Y_Base + 0.2, Y_Base - 0.2),
        Is_Top5 = ifelse(Rank <= 5, "Yes", "No")
    )

# Draw the offset lollipop chart.
pd_lollipop <- ggplot(df_combined_g) +
    geom_segment(aes(x = 0, xend = Consensus_Score - 5, y = Y_Plot, yend = Y_Plot, color = Category, linetype = Dataset_Type, alpha = Is_Top5), linewidth = 0.5) +
    geom_point(aes(x = Consensus_Score, y = Y_Plot, color = Category, alpha = Is_Top5, shape = Dataset_Type), size = 3.8, stroke = 0.8) +
    geom_text(data = filter(df_combined_g, Dataset_Type == "Platform-wide"), aes(x = Consensus_Score, y = Y_Plot, label = ifelse(Rank <= 5, Rank, "")), color = "white", size = 2.5, fontface = "plain") +
    geom_text(data = filter(df_combined_g, Dataset_Type == "Platform-shared"), aes(x = Consensus_Score, y = Y_Plot, label = ifelse(Rank <= 5, Rank, ""), color = Category), size = 2.5, fontface = "plain") +
    facet_wrap(~ Platform, nrow = 1) +
    scale_color_manual(values = cat_colors) +
    scale_linetype_manual(name = "Target scope", breaks = c("Platform-wide", "Platform-shared"), values = c("Platform-wide" = "solid", "Platform-shared" = "22")) +
    scale_shape_manual(name = "Target scope", breaks = c("Platform-wide", "Platform-shared"), values = c("Platform-wide" = 19, "Platform-shared" = 21)) +
    scale_alpha_manual(values = c("Yes" = 1, "No" = 0.3), guide = "none") +
    scale_y_continuous(breaks = 1:length(levels(df_combined_g$Property)), labels = levels(df_combined_g$Property), expand = expansion(add = 0.5)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) +
    labs(x = "Consensus relative importance (%)", y = NULL) +
    plasmix_theme +
    guides(
        color = guide_legend(title = "Category", override.aes = list(size = 0, linewidth = 2, nrow = 1, byrow = TRUE, keyheight = unit(0, "pt")), order = 1),
        linetype = guide_legend(title = "Target scope", override.aes = list(color = "black", linewidth = 0.8, shape = c(19, 21), fill = c("black", "white"), size = 2, alpha = 1), order = 2),
        shape = guide_legend(title = "Target scope", override.aes = list(color = "black", linewidth = 0.5, shape = c(19, 21), fill = c("black", "white"), size = 1, alpha = 1), order = 2)
    ) +
    theme(plot.background = element_rect(fill = "transparent", color = NA),
          panel.grid.major = element_blank(),
          panel.spacing.x = unit(10, "pt"),
          legend.key.width = unit(10, "pt"),
          legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
          legend.margin = margin(t = 2.5, b = 0, r = 0, l = 0),
          legend.title = element_text(size = 7.5, face = "bold", vjust = 0.5),
          legend.text = element_text(size = 7.5, vjust = 0.5, margin = margin(0, 0, 0, 3)),
          legend.box.margin = margin(t = -7.5, b = 0, r = 0, l = -20),
          axis.text.x = element_text(hjust = c(0.45, 0.5, 0.5, 0.5, 0.6)))

# Save Figure 4 analysis objects ----
fig4_analysis_results <- list(
    datasets = datasets,
    model_input = model_input_rows,
    models = global_models,
    model_tables = table_model_list,
    model_support = model_support,
    cv_repeat_metrics = cv_repeat_metrics,
    cv_fold_assignments = cv_fold_assignments,
    vip_tables = table_vip_list,
    target_properties = levels(df_combined_g$Property),
    prop_categories = df_combined_g %>% distinct(Property, Category) %>% arrange(Property)
)
saveRDS(fig4_analysis_results, "results/fig4_analysis_results.rds")

# Panel h: Two-dimensional partial dependence ----
pd1_som <- pdp::partial(global_models[["Intersect"]][["SOM"]]$RandomForest, pred.var = c("Delta_pI_7.4", "pLDDT_Fraction_Low"),
                  train = as.data.frame(global_models[["Intersect"]][["SOM"]]$X_train_rf), grid.resolution = 30, chull = TRUE)
pd2_olk <- pdp::partial(global_models[["Intersect"]][["OLK"]]$RandomForest, pred.var = c("Glyco_Density", "Sequence_Instability"),
                  train = as.data.frame(global_models[["Intersect"]][["OLK"]]$X_train_rf), grid.resolution = 30, chull = TRUE)
pd3_dia <- pdp::partial(global_models[["Intersect"]][["DIA"]]$RandomForest, pred.var = c("Glyco_Density", "Net_Charge_7.4"),
                  train = as.data.frame(global_models[["Intersect"]][["DIA"]]$X_train_rf), grid.resolution = 30, chull = TRUE)

pd1_som <- pd1_som %>% mutate(yhat = 2^yhat)
pd2_olk <- pd2_olk %>% mutate(yhat = 2^yhat)
pd3_dia <- pd3_dia %>% mutate(yhat = 2^yhat)

z_limits <- range(c(pd1_som$yhat, pd2_olk$yhat, pd3_dia$yhat), na.rm = TRUE)
z_breaks <- pretty(z_limits, n = 4)

draw_2d_pdp_shared <- function(pd_data, var_labs, platform_lab, z_limits, z_breaks) {
    ggplot(pd_data, aes_string(x = names(pd_data)[1], y = names(pd_data)[2], z = "yhat")) +
        geom_tile(aes(fill = yhat)) +
        geom_contour(color = "white", alpha = 0.4) +
        scale_fill_viridis_c(option = "magma", direction = 1, name = "Predicted\nexpansion", limits = z_limits, breaks = z_breaks,
                            guide = guide_colorbar(barheight = unit(3, "cm"), barwidth = unit(0.25, "cm"), title.position = "top")) +
        labs(title = platform_lab, x = var_labs[1], y = var_labs[2]) +
        plasmix_theme +
        theme(panel.grid.major = element_blank(),
              axis.line.x = element_blank(), axis.line.y = element_blank(),
              axis.text = element_blank(), axis.ticks = element_blank(),
              plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5, margin = margin(b = 2)),
              legend.title = element_text(angle = 90, hjust = 0, vjust = 1),
              legend.margin = margin(0, 5, 0, -5),
              axis.title.x = element_text(margin = margin(t = 2)),
              axis.title.y = element_text(margin = margin(r = 2), vjust = -1))
}

pe_pdp1_som <- draw_2d_pdp_shared(pd1_som, c(name_map["Delta_pI_7.4"], str_wrap(name_map["pLDDT_Fraction_Low"], width = 10)), "SOM", z_limits, z_breaks)
pe_pdp2_olk <- draw_2d_pdp_shared(pd2_olk, c(name_map["Glyco_Density"], str_wrap(name_map["Sequence_Instability"], width = 10)), "OLK", z_limits, z_breaks)
pe_pdp3_dia <- draw_2d_pdp_shared(pd3_dia, c(name_map["Glyco_Density"], str_wrap(name_map["Net_Charge_7.4"], width = 10)), "DIA", z_limits, z_breaks)
pe_pdp_combined <- ggarrange(pe_pdp1_som, pe_pdp2_olk, pe_pdp3_dia, ncol = 1, common.legend = TRUE, legend = "right")

# Assemble and export Figure 4 ----
row1_left_up <- ggarrange(plot_envelope_dia, plot_envelope_olk, plot_envelope_som, nrow = 1, widths = c(1.07, 1, 1),
                        labels = "a", font.label = label_style, label.x = 0, label.y = 1, hjust = -0.2, vjust = 1.2)
row1_left_bottom <- ggarrange(p_vise, p_grid, nrow = 1, widths = c(1, 1.1),
                            labels = c("c", "d"), font.label = label_style, label.x = 0, label.y = 1, hjust = -0.2, vjust = 1.2)
row1_left <- ggarrange(row1_left_up, row1_left_bottom, ncol = 1, heights = c(1.05, 1))
row1 <- ggarrange(row1_left, p_envelope_metrics, nrow = 1, widths = c(2, 1),
                labels = c("", "b"), font.label = label_style, label.x = 0, label.y = 1, hjust = -0.2, vjust = 1.2)
row2 <- ggarrange(pb_scatter, labels = "e", font.label = label_style,
                label.x = 0, label.y = 1, hjust = -0.2, vjust = 1.2)
row3 <- ggarrange(pc_bar, pd_lollipop, pe_pdp_combined, ncol = 3, widths = c(1, 2, 0.9),
                labels = c("f", "g", "h"), font.label = label_style, label.x = 0, label.y = 1, hjust = -0.2, vjust = 1.2)
fig_main <- ggarrange(row1, row2, row3, ncol = 1, heights = c(2, 1, 1.4))
ggsave("figures/fig4_distortion_cause.pdf", fig_main, width = 10, height = 10)
ggsave("figures/fig4_distortion_cause.png", fig_main, width = 10, height = 10, dpi = 600, bg = "white")

# Format source data ----
platform_order <- c("DIA", "OLK", "SOM")

format_response_section <- function(df, level, identifier) {
  df %>%
    transmute(
      Level = level,
      Identifier = as.character(.data[[identifier]]),
      `Feature observations` = as.integer(Features_Count),
      `Lower plateau, A (log10 intensity)` = A_log,
      `Response slope, B` = Response_slope,
      `Abundance midpoint, C (log10 pg/mL)` = Abundance_mid_log,
      `Upper plateau, D (log10 intensity)` = D_log,
      `Intensity span` = D_log - A_log,
      `Abundance span` = Abundance_max - Abundance_min,
      `Background burden (%)` = Background_burden_pct,
      `Mid-range CV` = cv_tech,
      `Lower-range CV` = cv_bot,
      `Upper-range CV` = cv_top,
      `Lower response bound (log10 intensity)` = Lower_bound_log,
      `Upper response bound (log10 intensity)` = Upper_bound_log,
      `Lower abundance bound (log10 pg/mL)` = Abundance_min,
      `Upper abundance bound (log10 pg/mL)` = Abundance_max,
      `Δ95 Native` = Delta95_Raw,
      `Δ95 BLK` = Delta95_BLK,
      `Δ95 N` = Delta95_N,
      `Expansion factor, BLK` = Expansion_BLK,
      `Expansion factor, N` = Expansion_N
    )
}

# Platform- and batch-level response envelopes and ratio expansion.
section_a <- envelope_by_platform %>%
  left_join(envelope_metrics_platform, by = "Platform") %>%
  left_join(platform_expansion, by = "Platform") %>%
  mutate(Platform = factor(Platform, levels = platform_order)) %>%
  arrange(Platform) %>%
  format_response_section(level = "Platform", identifier = "Platform")

section_b <- envelope_by_batch %>%
  left_join(envelope_metrics_batch, by = c("Platform", "Batch")) %>%
  left_join(batch_expansion, by = "Batch") %>%
  mutate(Platform = factor(Platform, levels = platform_order)) %>%
  arrange(Platform, Batch) %>%
  format_response_section(level = "Batch", identifier = "Batch")

st_response_envelope <- bind_rows(section_a, section_b)

# ST: Model R² and property-importance matrix.
model_r2_wide <- df_panel_c %>%
  mutate(`R squared` = sprintf("%.2f ± %.2f%%", R_Squared_Mean, R_Squared_SD)) %>%
  select(Platform, Dataset, Algorithm, `R squared`, Target_Features, Total_Observations) %>%
  pivot_wider(names_from = Algorithm, values_from = `R squared`, names_glue = "{Algorithm} R²")

# Rank properties within each algorithm.
df_vip_ensemble <- bind_rows(table_vip_list) %>%
  filter(Feature != "Random_Noise") %>%
  group_by(Dataset, Platform, Algorithm) %>%
  arrange(desc(Relative_Imp), .by_group = TRUE) %>%
  mutate(Rank = row_number()) %>%
  ungroup()

# Calculate consensus scores and ranks.
df_consensus <- df_vip_ensemble %>%
  group_by(Dataset, Platform, Feature) %>%
  summarize(`Consensus score` = mean(Relative_Imp, na.rm = TRUE), .groups = "drop") %>%
  group_by(Dataset, Platform) %>%
  arrange(desc(`Consensus score`), .by_group = TRUE) %>%
  mutate(`Consensus rank` = row_number()) %>%
  ungroup()

# Reshape raw importance scores and ranks.
df_metrics <- df_vip_ensemble %>%
  select(Dataset, Platform, Feature, Algorithm, Raw_score, Rank) %>%
  pivot_wider(names_from = Algorithm, values_from = c(Raw_score, Rank), names_glue = "{Algorithm}_{.value}")

# Assemble the final model-summary table.
vip_matrix <- df_consensus %>%
  left_join(df_metrics, by = c("Dataset", "Platform", "Feature")) %>%
  mutate(
    Category = factor(cat_map[Feature], levels = names(cat_colors)),
    Property = name_map[Feature],
    Platform = factor(Platform, levels = c("SOM", "OLK", "DIA")),
    Dataset = factor(Dataset, levels = c("Full", "Intersect"))
  )

st_model_vip <- vip_matrix %>%
  left_join(model_r2_wide, by = c("Platform", "Dataset")) %>%
  select(
    Platform, Dataset,
    `Target features` = Target_Features,
    `Feature–batch observations` = Total_Observations,
    `RF R² (Mean ± SD)` = `RandomForest R²`,
    `XGB R² (Mean ± SD)` = `XGBoost R²`,
    `LGBM R² (Mean ± SD)` = `LightGBM R²`,
    Category, Property,
    `RF %IncMSE` = RandomForest_Raw_score,
    `XGB gain` = XGBoost_Raw_score,
    `LGBM gain` = LightGBM_Raw_score,
    `RF rank` = RandomForest_Rank,
    `XGB rank` = XGBoost_Rank,
    `LGBM rank` = LightGBM_Rank,
    `Consensus score`, `Consensus rank`
  ) %>%
  arrange(Platform, Dataset, `Consensus rank`) %>%
  mutate(across(where(is.factor), as.character))

# Record the Excel rows occupied by each Platform–Dataset block.
model_merge_rows <- st_model_vip %>%
  mutate(.excel_row = row_number() + 1L) %>%
  group_by(Platform, Dataset) %>%
  summarize(Row_start = min(.excel_row), Row_end = max(.excel_row), .groups = "drop")

# Repeated values are blanked before the corresponding Excel cells are merged.
st_model_vip_display <- st_model_vip
duplicate_model_group <- duplicated(st_model_vip_display[c("Platform", "Dataset")])
st_model_vip_display[duplicate_model_group, 1:7] <- NA

# Write and format the source-data workbook ----
wb <- createWorkbook()
addWorksheet(wb, "Response_envelope", gridLines = FALSE)
addWorksheet(wb, "Model_summary", gridLines = FALSE)

header_style <- createStyle(
  fontName = "Aptos Narrow", fontSize = 12, textDecoration = "bold",
  valign = "center", wrapText = TRUE
)
body_style <- createStyle(
  fontName = "Aptos Narrow", fontSize = 12,
  halign = "left", valign = "center", wrapText = TRUE
)
level_style <- createStyle(
  fontName = "Aptos Narrow", fontSize = 12, textDecoration = "bold",
  halign = "left", valign = "center", wrapText = TRUE
)
integer_style <- createStyle(numFmt = "#,##0", halign = "left", valign = "center")
decimal_2_style <- createStyle(numFmt = "0.00", halign = "left", valign = "center")
decimal_3_style <- createStyle(numFmt = "0.000", halign = "left", valign = "center")
decimal_4_style <- createStyle(numFmt = "0.0000", halign = "left", valign = "center")
burden_style <- createStyle(numFmt = '0.0000"%"', halign = "left", valign = "center")

# Response_envelope.
writeData(
  wb, "Response_envelope", st_response_envelope,
  headerStyle = header_style, withFilter = FALSE, keepNA = FALSE
)
addStyle(
  wb, "Response_envelope", body_style,
  rows = 2:(nrow(st_response_envelope) + 1), cols = 1:ncol(st_response_envelope),
  gridExpand = TRUE, stack = TRUE
)
addStyle(
  wb, "Response_envelope", integer_style,
  rows = 2:(nrow(st_response_envelope) + 1), cols = 3,
  gridExpand = TRUE, stack = TRUE
)
addStyle(
  wb, "Response_envelope", decimal_3_style,
  rows = 2:(nrow(st_response_envelope) + 1), cols = c(4:9, 11:22),
  gridExpand = TRUE, stack = TRUE
)
addStyle(
  wb, "Response_envelope", burden_style,
  rows = 2:(nrow(st_response_envelope) + 1), cols = 10,
  gridExpand = TRUE, stack = TRUE
)

platform_excel_rows <- 2:(nrow(section_a) + 1)
batch_excel_rows <- (nrow(section_a) + 2):(nrow(st_response_envelope) + 1)
mergeCells(wb, "Response_envelope", cols = 1, rows = platform_excel_rows)
mergeCells(wb, "Response_envelope", cols = 1, rows = batch_excel_rows)
addStyle(wb, "Response_envelope", level_style, rows = c(platform_excel_rows[1], batch_excel_rows[1]), cols = 1, stack = TRUE)

setColWidths(
  wb, "Response_envelope", cols = 1:22,
  widths = c(9, 13, 18, 20, 15, 22, 20, 14, 15, 19, 14, 14, 14, 22, 22, 23, 23, 13, 11, 11, 17, 16)
)
setRowHeights(wb, "Response_envelope", rows = 1, heights = 52)
setRowHeights(wb, "Response_envelope", rows = 2:(nrow(st_response_envelope) + 1), heights = 18)
freezePane(wb, "Response_envelope", firstRow = TRUE)

# Model_summary.
writeData(
  wb, "Model_summary", st_model_vip_display,
  headerStyle = header_style, withFilter = FALSE, keepNA = FALSE
)
addStyle(
  wb, "Model_summary", body_style,
  rows = 2:(nrow(st_model_vip_display) + 1), cols = 1:ncol(st_model_vip_display),
  gridExpand = TRUE, stack = TRUE
)
addStyle(
  wb, "Model_summary", integer_style,
  rows = 2:(nrow(st_model_vip_display) + 1), cols = c(3, 4, 13:15, 17),
  gridExpand = TRUE, stack = TRUE
)
addStyle(
  wb, "Model_summary", decimal_2_style,
  rows = 2:(nrow(st_model_vip_display) + 1), cols = c(10, 16),
  gridExpand = TRUE, stack = TRUE
)
addStyle(
  wb, "Model_summary", decimal_4_style,
  rows = 2:(nrow(st_model_vip_display) + 1), cols = 11:12,
  gridExpand = TRUE, stack = TRUE
)

for(i in seq_len(nrow(model_merge_rows))) {
  merge_rows <- model_merge_rows$Row_start[i]:model_merge_rows$Row_end[i]
  for(j in 1:7) mergeCells(wb, "Model_summary", cols = j, rows = merge_rows)
}

setColWidths(
  wb, "Model_summary", cols = 1:17,
  widths = c(9, 10, 15, 21, 17, 17, 18, 12, 22, 12, 10, 11, 9, 10, 11, 16, 15)
)
setRowHeights(wb, "Model_summary", rows = 1, heights = 52)
setRowHeights(wb, "Model_summary", rows = 2:(nrow(st_model_vip_display) + 1), heights = 18)
freezePane(wb, "Model_summary", firstRow = TRUE)

saveWorkbook(wb, "tables/SourceData_Figure4.xlsx", overwrite = TRUE)
