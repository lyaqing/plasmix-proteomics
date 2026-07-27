# Extended Data Figure 6 | Physicochemical drivers of distortion

# 0. Setup ----
source("scripts/_project_setup.R")
use_packages(c("data.table", "tidyverse", "ggh4x", "ggpubr", "patchwork", "ggtext", "openxlsx", "showtext", "grid"))
source("utils/figure_style.R")
plasmix_theme <- if (is.function(theme_plasmix)) theme_plasmix() else theme_plasmix
showtext_auto(); showtext_opts(dpi = 600)
label_style <- list(size = 12, face = "bold")
paths <- c(
    physchem = "data/physchem_matrix.tsv.gz",
    physchem_dictionary = "data/physchem_dictionary.tsv",
    analysis_context = "results/fig4_analysis_results.rds"
)
missing_inputs <- paths[!file.exists(paths)]
if(length(missing_inputs)) stop("Missing Extended Data Figure 6 inputs:\n", paste(missing_inputs, collapse = "\n"))

# 1. Physicochemical features and model objects ----
physchem_matrix <- fread(paths["physchem"]) %>% as_tibble()
physchem_dict <- fread(paths["physchem_dictionary"]) %>% as_tibble()
ed6_context <- readRDS(paths["analysis_context"])

cat_colors <- c("Structure" = "#E64B35", "Surface" = "#4DBBD5", "Charge" = "#00A087", "Disorder" = "#F39B7F", "Secretory" = "#8491B4", "Abundance" = "#91D1C2")
category_order <- c("Structure", "Surface", "Charge", "Disorder", "Secretory", "Abundance")
physchem_dict <- physchem_dict %>% mutate(Category = factor(Category, levels = category_order))
final_physchem <- physchem_dict %>% filter(Retained == "Yes") %>% pull(Feature)

name_map <- setNames(physchem_dict$Property, physchem_dict$Feature)
cat_map <- setNames(physchem_dict$Category, physchem_dict$Feature)
name_map_chr <- setNames(as.character(name_map), names(name_map))
cat_map_chr <- setNames(as.character(cat_map), names(cat_map))

map_name <- function(x, map) {
    y <- unname(map[as.character(x)])
    ifelse(is.na(y), as.character(x), y)
}

datasets <- ed6_context$datasets
global_models <- ed6_context$models
target_properties <- as.character(ed6_context$target_properties)
prop_categories <- ed6_context$prop_categories %>%
    mutate(Property = factor(as.character(Property), levels = target_properties)) %>%
    arrange(Property)
platforms_to_run <- c("SOM", "OLK", "DIA")

# 2. Feature correlations ----
features_for_cor <- physchem_matrix %>%
    select(where(is.numeric)) %>%
    select(-any_of(c("CompBias_Count", "Transmembrane_Count", "BloodConc_log10_pgml")))

colnames(features_for_cor) <- sapply(colnames(features_for_cor), function(x) map_name(x, name_map_chr))

cor_matrix <- cor(features_for_cor, method = "spearman", use = "pairwise.complete.obs")

n_feat <- ncol(features_for_cor)
cor_pmat <- matrix(NA, n_feat, n_feat)
diag(cor_pmat) <- 0
for (i in 1:(n_feat - 1)) {
    for (j in (i + 1):n_feat) {
        tmp <- cor.test(features_for_cor[[i]], features_for_cor[[j]], method = "spearman", exact = FALSE)
        cor_pmat[i, j] <- cor_pmat[j, i] <- tmp$p.value
    }
}
colnames(cor_pmat) <- rownames(cor_pmat) <- colnames(cor_matrix)

dist_mat <- as.dist(1 - abs(cor_matrix))
hc <- hclust(dist_mat, method = "ward.D2")
ordered_features <- colnames(cor_matrix)[hc$order]

cor_matrix <- cor_matrix[ordered_features, ordered_features]
cor_pmat <- cor_pmat[ordered_features, ordered_features]

cor_matrix[upper.tri(cor_matrix)] <- NA
cor_pmat[upper.tri(cor_pmat)] <- NA

df_cor <- as.data.frame(as.table(cor_matrix)) %>%
    rename(Var1 = Var1, Var2 = Var2, cor_value = Freq) %>%
    mutate(
        Var1 = factor(Var1, levels = rev(ordered_features)),
        Var2 = factor(Var2, levels = ordered_features)
    ) %>%
    filter(!is.na(cor_value))

df_pmat <- as.data.frame(as.table(cor_pmat)) %>%
    rename(Var1 = Var1, Var2 = Var2, pval = Freq) %>%
    mutate(
        Var1 = factor(Var1, levels = rev(ordered_features)),
        Var2 = factor(Var2, levels = ordered_features)
    ) %>%
    filter(!is.na(pval))

plot_cor_df <- left_join(df_cor, df_pmat, by = c("Var1", "Var2"))
p_corr <- ggplot(plot_cor_df, aes(x = Var2, y = Var1)) +
    geom_tile(aes(fill = cor_value), color = "white", linewidth = 0.5) +
    geom_point(data = filter(plot_cor_df, pval > 0.05), shape = 4, size = 1.5, color = "grey50") +
    scale_fill_gradientn(colors = colorRampPalette(c("#053061", "#FFFFFF", "#67001F"))(200), limits = c(-1, 1), name = "Spearman\ncorrelation") +
    labs(x = NULL, y = NULL) +
    plasmix_theme +
    theme(
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = c(1, 1), legend.justification = c(1, 1),
        legend.direction = "horizontal",
        legend.title = element_text(vjust = 1),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.25, "cm")
    )
print(p_corr)

# 3. Wilcoxon heatmap ----
scan_thresholds <- seq(0.05, 0.50, by = 0.05)
wilcox_list <- list()
for (ds_name in names(datasets)) {
    df_current <- datasets[[ds_name]]
    for (plat in platforms_to_run) {
        df_plat <- df_current %>% filter(Platform == plat)
        y <- df_plat$Y_RelExpansion
        for (th in scan_thresholds) {
            y_low_cutoff <- quantile(y, th, na.rm = TRUE)
            y_high_cutoff <- quantile(y, 1 - th, na.rm = TRUE)
            idx_low_compression <- which(y <= y_low_cutoff)
            idx_high_compression <- which(y >= y_high_cutoff)
            for (feat in final_physchem) {
                val <- df_plat[[feat]]
                group_feat_uncompressed <- val[idx_low_compression]
                group_feat_highly_compressed <- val[idx_high_compression]
                group_feat_uncompressed <- group_feat_uncompressed[!is.na(group_feat_uncompressed)]
                group_feat_highly_compressed <- group_feat_highly_compressed[!is.na(group_feat_highly_compressed)]
                if(length(group_feat_highly_compressed) >= 10 && length(group_feat_uncompressed) >= 10) {
                    p_val <- wilcox.test(group_feat_highly_compressed, group_feat_uncompressed, exact = FALSE)$p.value
                    wilcox_list[[length(wilcox_list)+1]] <- data.frame(
                        Dataset = ds_name, Platform = plat, Feature = feat,
                        Threshold = paste0(th * 100, "%"), P_Value = p_val)
                }
            }
        }
    }
}
wilcox_wide_df <- bind_rows(wilcox_list) %>%
    mutate(P_Value = as.numeric(formatC(P_Value, format = "e", digits = 2))) %>%
    pivot_wider(names_from = Threshold, values_from = P_Value)

df_panel_e <- bind_rows(wilcox_list) %>%
    filter(Dataset == "Full") %>%
    mutate(
        Property = name_map[Feature],
        NegLog10_P = -log10(P_Value),
        Color_Scale = ifelse(NegLog10_P > 10, 10, NegLog10_P),
        Significance = case_when(
            P_Value < 0.001 ~ "***",
            P_Value < 0.01 ~ "**",
            P_Value < 0.05 ~ "*",
            TRUE ~ ""
        ),
        star_color = ifelse(Color_Scale > 3, "white", "black")
    )

df_panel_e$Threshold_Num <- as.numeric(gsub("%", "", df_panel_e$Threshold))
df_panel_e <- df_panel_e %>% arrange(Threshold_Num)
df_panel_e$Threshold <- factor(df_panel_e$Threshold, levels = unique(df_panel_e$Threshold))

df_panel_e <- df_panel_e %>%
    mutate(Category = factor(cat_map[Feature], levels = names(cat_colors)))

y_axis_order <- df_panel_e %>%
    distinct(Category, Property) %>%
    arrange(Category, Property) %>%
    pull(Property)

df_panel_e$Property <- factor(df_panel_e$Property, levels = rev(y_axis_order))
df_panel_e$Platform <- factor(df_panel_e$Platform, levels = c("SOM", "OLK", "DIA"))

df_cat <- df_panel_e %>% distinct(Property, Category) %>% mutate(Ghost_Panel = " ")

p_cat <- ggplot(df_cat, aes(x = "1", y = Property, fill = Category)) +
    geom_tile(color = "white", linewidth = 1, width = 0.6) +
    facet_wrap(~ Ghost_Panel) +
    scale_fill_manual(values = cat_colors, name = "Category") +
    plasmix_theme +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major = element_blank(),
        legend.key.size = unit(0.6, "lines"),
        legend.text = element_text(margin = margin(l = 3)),
        strip.background = element_blank(),
        strip.text = element_text(color = "transparent", size = 8.5),
        legend.margin = margin(0, 0, 0, -5),
        plot.margin = margin(5, 0, 5, 5)
    ) +
    guides(fill = guide_legend(order = 2, ncol = 1, keywidth = unit(0.6, "lines"), keyheight = unit(0.6, "lines"), override.aes = list(color = "white", linewidth = 0.25)))

p_main <- ggplot(df_panel_e, aes(x = Threshold, y = Property, fill = Color_Scale)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = Significance, color = star_color), size = 2.5, vjust = 0.7, angle = 90, show.legend = FALSE) +
    facet_wrap(~ Platform, ncol = 3) +
    scale_fill_gradientn(
        colors = c("white", "#FFD1D1", "#FF6B6B", "#CC0000", "#800000"),
        values = scales::rescale(c(0, 1.3, 3, 5, 10)),
        name = bquote(bold(-log[10](italic(P)))),
        limits = c(0, 10), oob = scales::squish
    ) +
    scale_color_manual(values = c("black" = "black", "white" = "white")) +
    labs(x = "Compression sampling threshold (top vs bottom)", y = NULL) +
    guides(fill = guide_colorbar(barheight = unit(2, "cm"), barwidth = unit(0.25, "cm"))) +
    plasmix_theme +
    theme(
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_text(margin = margin(5, 0, 0, 0)),
        legend.margin = margin(0, 0, 0, -5),
        plot.margin = margin(5, 5, 5, 0)
    )

sp_wilcox <- p_cat + p_main + plot_layout(widths = c(0.03, 1), guides = "collect") & theme(plot.margin = margin(2, 0, 2, 2))

# 4. One-dimensional ALE and PDP ----
## 1D-PDP Calculation ----
grid_n <- 30

pdp_grid <- function(x, n = 30) {
    x <- x[is.finite(x)]
    u <- sort(unique(x))
    if (length(u) <= 2) return(u)
    unique(as.numeric(quantile(x, seq(0.02, 0.98, length.out = n), na.rm = TRUE, type = 8)))
}

calc_pdp1d <- function(model, x_train, feature, n = 30) {
    x_train <- as.data.frame(x_train)
    g <- pdp_grid(x_train[[feature]], n)
    map_dfr(g, function(v) {
        x_new <- x_train
        x_new[[feature]] <- v
        tibble(
            Feature = feature,
            x = v,
            yhat = 2^mean(predict(model, newdata = x_new), na.rm = TRUE)
        )
    })
}

pdp_1d_all <- pmap_dfr(
    expand.grid(Dataset = c("Full", "Intersect"), Platform = c("SOM", "OLK", "DIA"), stringsAsFactors = FALSE),
    function(Dataset, Platform) {
        ds <- as.character(Dataset)
        plat <- as.character(Platform)
        x_train <- global_models[[ds]][[plat]]$X_train_rf
        model   <- global_models[[ds]][[plat]]$RandomForest
        feats   <- setdiff(colnames(x_train), "Random_Noise")
        map_dfr(feats, function(feat) {
            calc_pdp1d(model, x_train, feat, grid_n) %>%
                mutate(Dataset = ds, Platform = plat, Property = map_name(Feature, name_map_chr),
                       Category = map_name(Feature, cat_map_chr), Scale = "Original expansion")
        })
    }
)
write.csv(pdp_1d_all, "results/ed6_pdp_1d_full_intersect.csv", row.names = FALSE)

## 1D-ALE Calculation ----
K_ale <- 20

ale_feature_tbl <- pdp_1d_all %>%
    distinct(Platform, Feature) %>%
    filter(Platform %in% platforms_to_run) %>%
    crossing(Dataset = names(datasets)) %>%
    select(Dataset, Platform, Feature)

calc_ale1d <- function(model, x_train, feature, K = 20) {
    x_train <- as.data.frame(x_train)
    feature <- as.character(feature)
    x <- x_train[[feature]]
    ok <- is.finite(x)
    xv <- x[ok]
    u <- sort(unique(xv))
    if (length(u) < 2) {
        return(tibble(Feature = feature, x_mid = NA_real_, ale = NA_real_, n_bin = NA_integer_))
    }
    if (length(u) == 2) {
        x0 <- x_train; x1 <- x_train
        x0[[feature]] <- u[1]; x1[[feature]] <- u[2]
        d <- mean(predict(model, newdata = x1) - predict(model, newdata = x0), na.rm = TRUE)
        return(tibble(
            Feature = feature,
            x_mid = u,
            ale = c(-d / 2, d / 2),
            n_bin = c(sum(xv == u[1]), sum(xv == u[2]))
        ))
    }
    qs <- unique(as.numeric(quantile(xv, probs = seq(0, 1, length.out = K + 1), na.rm = TRUE, type = 8)))
    if (length(qs) < 3) return(tibble(Feature = feature, x_mid = NA_real_, ale = NA_real_, n_bin = NA_integer_))
    K_eff <- length(qs) - 1
    bin <- findInterval(x, qs, rightmost.closed = TRUE, all.inside = TRUE)
    d_mean <- numeric(K_eff)
    n_bin <- integer(K_eff)
    for (k in seq_len(K_eff)) {
        idx <- which(ok & bin == k)
        n_bin[k] <- length(idx)
        if (length(idx) == 0) next
        xl <- x_train[idx, , drop = FALSE]; xh <- x_train[idx, , drop = FALSE]
        xl[[feature]] <- qs[k]; xh[[feature]] <- qs[k + 1]
        d_mean[k] <- mean(predict(model, newdata = xh) - predict(model, newdata = xl), na.rm = TRUE)
    }
    ale_raw <- cumsum(d_mean)
    ale_center <- weighted.mean(ale_raw, w = pmax(n_bin, 1), na.rm = TRUE)
    tibble(
        Feature = feature,
        x_mid = (qs[-1] + qs[-length(qs)]) / 2,
        ale = ale_raw - ale_center,
        n_bin = n_bin
    )
}

ale_1d_all <- pmap_dfr(ale_feature_tbl, function(Dataset, Platform, Feature) {
    ds <- as.character(Dataset)
    plat <- as.character(Platform)
    feat <- as.character(Feature)
    calc_ale1d(model = global_models[[ds]][[plat]]$RandomForest, x_train = global_models[[ds]][[plat]]$X_train_rf,
               feature = feat, K = K_ale) %>%
        mutate(Dataset = ds, Platform = plat, Property = map_name(Feature, name_map_chr),
               Category = map_name(Feature, cat_map_chr), Scale = "log2 model scale")
})
write.csv(ale_1d_all, "results/ed6_ale_1d_full_intersect.csv", row.names = FALSE)

## Curve classification and summarization ----
safe_cor <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) return(NA_real_)
    suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

classify_curve <- function(x, y, min_effect = 0.15, cor_cut = 0.45) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]
    if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) return("Insufficient")
    amp <- diff(range(y, na.rm = TRUE))
    if (!is.finite(amp) || amp < min_effect) return("Flat/weak")
    if (length(unique(x)) <= 2) {
        y_low  <- mean(y[x == min(x)], na.rm = TRUE)
        y_high <- mean(y[x == max(x)], na.rm = TRUE)
        if (!is.finite(y_low) || !is.finite(y_high)) return("Insufficient")
        if (y_high > y_low) return("Increasing")
        if (y_high < y_low) return("Decreasing")
        return("Flat/weak")
    }
    rho <- suppressWarnings(cor(x, y, method = "spearman"))
    q <- quantile(x, probs = c(0.2, 0.4, 0.6, 0.8), na.rm = TRUE)
    y_low <- mean(y[x <= q[1]], na.rm = TRUE)
    y_mid <- mean(y[x >= q[2] & x <= q[3]], na.rm = TRUE)
    y_high <- mean(y[x >= q[4]], na.rm = TRUE)
    y_edge <- mean(c(y_low, y_high), na.rm = TRUE)
    if (is.finite(y_mid - y_edge) && (y_mid - y_edge) > min_effect) return("Middle high")
    if (is.finite(y_edge - y_mid) && (y_edge - y_mid) > min_effect) return("Middle low")
    if (is.finite(rho) && rho > cor_cut) return("Increasing")
    if (is.finite(rho) && rho < -cor_cut) return("Decreasing")
    return("Non-monotonic")
}

summarize_curve <- function(df, x_name, y_name, prefix) {
    df %>%
        transmute(
            Dataset, Platform, Feature, Property, Category,
            x = .data[[x_name]],
            y = .data[[y_name]]
        ) %>%
        filter(is.finite(x), is.finite(y)) %>%
        arrange(Dataset, Platform, Feature, x) %>%
        group_by(Dataset, Platform, Feature, Property, Category) %>%
        summarize(
            "{prefix}_min"       := min(y, na.rm = TRUE),
            "{prefix}_max"       := max(y, na.rm = TRUE),
            "{prefix}_range"     := max(y, na.rm = TRUE) - min(y, na.rm = TRUE),
            "{prefix}_start"     := first(y),
            "{prefix}_end"       := last(y),
            Spearman_rho         = safe_cor(x, y),
            Direction            = classify_curve(x, y),
            N_points             = n(),
            .groups = "drop"
        ) %>%
        mutate(across(where(is.numeric), ~ round(.x, 4)))
}

compare_full_intersect <- function(df, x_name, y_name, summary_tbl) {
    raw <- df %>%
        transmute(Dataset, Platform, Feature, x = .data[[x_name]], y = .data[[y_name]]) %>%
        filter(Dataset %in% c("Full", "Intersect"), is.finite(x), is.finite(y))
    shape_tbl <- raw %>%
        group_by(Platform, Feature) %>%
        group_modify(function(d, key) {
            f <- d %>% filter(Dataset == "Full") %>% arrange(x)
            i <- d %>% filter(Dataset == "Intersect") %>% arrange(x)
            if (nrow(f) < 3 || nrow(i) < 3) return(tibble(FI_Shape_rho = NA_real_))
            lo <- max(min(f$x), min(i$x))
            hi <- min(max(f$x), max(i$x))
            if (!is.finite(lo) || !is.finite(hi) || lo >= hi) return(tibble(FI_Shape_rho = NA_real_))
            g <- seq(lo, hi, length.out = 50)
            yf <- approx(f$x, f$y, xout = g, rule = 1)$y
            yi <- approx(i$x, i$y, xout = g, rule = 1)$y
            tibble(FI_Shape_rho = round(safe_cor(yf, yi), 4))
        }) %>%
        ungroup()
    dir_tbl <- summary_tbl %>%
        select(Platform, Feature, Dataset, Direction) %>%
        pivot_wider(names_from = Dataset, values_from = Direction, names_prefix = "Direction_")
    shape_tbl %>%
        left_join(dir_tbl, by = c("Platform", "Feature")) %>%
        mutate(
            Set_Consistency = case_when(
                is.na(FI_Shape_rho) ~ "Not comparable",
                Direction_Full == Direction_Intersect & FI_Shape_rho >= 0.6 ~ "Consistent",
                Direction_Full == Direction_Intersect ~ "Direction match",
                FI_Shape_rho >= 0.6 ~ "Shape match",
                TRUE ~ "Divergent"
            )
        )
}

pdp_summary <- summarize_curve(df = pdp_1d_all, x_name = "x", y_name = "yhat", prefix = "PDP")
pdp_consistency <- compare_full_intersect(df = pdp_1d_all, x_name = "x", y_name = "yhat", summary_tbl = pdp_summary)
pdp_summary <- pdp_summary %>% left_join(pdp_consistency, by = c("Platform", "Feature")) %>% mutate(Scale = "Original expansion")

ale_summary <- summarize_curve(df = ale_1d_all, x_name = "x_mid", y_name = "ale", prefix = "ALE")
ale_consistency <- compare_full_intersect(df = ale_1d_all, x_name = "x_mid", y_name = "ale", summary_tbl = ale_summary)
ale_summary <- ale_summary %>% left_join(ale_consistency, by = c("Platform", "Feature")) %>% mutate(Scale = "log2 model scale")

## Visualization: 1D-PDP & 1D-ALE (Top Features Patchwork) ----
df_pdp_plot <- pdp_1d_all %>%
    filter(Property %in% target_properties) %>%
    mutate(Property = factor(Property, levels = target_properties), Platform = factor(Platform, levels = platforms_to_run),
           Target_Scope = factor(Dataset, levels = c("Full", "Intersect"), labels = c("Platform-wide", "Platform-shared")))

df_ale_plot <- ale_1d_all %>%
    filter(Property %in% target_properties) %>%
    rename(x = x_mid) %>%
    mutate(Property = factor(Property, levels = target_properties), Platform = factor(Platform, levels = platforms_to_run),
           Target_Scope = factor(Dataset, levels = c("Full", "Intersect"), labels = c("Platform-wide", "Platform-shared")))

trim_x_tails <- function(df, x_col = "x", multiplier = 2.0) {
    df %>%
        group_by(Property) %>%
        mutate(
            q1 = quantile(.data[[x_col]], 0.25, na.rm = TRUE),
            q3 = quantile(.data[[x_col]], 0.75, na.rm = TRUE),
            iqr = q3 - q1,
            lower_bound = q1 - multiplier * iqr,
            upper_bound = q3 + multiplier * iqr
        ) %>%
        filter(.data[[x_col]] >= lower_bound & .data[[x_col]] <= upper_bound) %>%
        select(-q1, -q3, -iqr, -lower_bound, -upper_bound) %>%
        ungroup()
}

df_pdp_plot <- df_pdp_plot %>% trim_x_tails(x_col = "x", multiplier = 2.0)
df_ale_plot <- df_ale_plot %>% trim_x_tails(x_col = "x", multiplier = 2.0)

strip_bg_colors <- cat_colors[as.character(prop_categories$Category)]

plot_1d_curve <- function(df_data, y_col, y_label, show_legend = FALSE, show_strip_x = TRUE) {
    p <- ggplot(df_data, aes(x = .data[["x"]], y = .data[[y_col]])) +
        geom_line(aes(color = Platform, linetype = Target_Scope), linewidth = 0.5) +
        geom_point(aes(color = Platform, shape = Target_Scope), size = 1.2, alpha = 0.8) +
        facet_grid2(Platform ~ Property, scales = "free",
                          labeller = labeller(Property = function(x) str_wrap(x, width = 10)),
                          strip = strip_themed(background_x = elem_list_rect(fill = strip_bg_colors, color = "white"))) +
        scale_color_manual(values = platform_color) +
        scale_linetype_manual(values = c("Platform-wide" = "solid", "Platform-shared" = 22)) +
        scale_shape_manual(values = c("Platform-wide" = 16, "Platform-shared" = 1)) +
        labs(x = "Feature value", y = y_label) +
        plasmix_theme +
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35), panel.grid.major = element_blank(),
              axis.line = element_blank(), strip.text.x = element_text(size = 8),
              axis.title.x = element_text(margin = margin(t = -5)), axis.text.x = element_text(angle = 45, hjust = 1))
    if(show_strip_x) {
        p <- p + theme(strip.text.x = element_text(color = "white", face = "bold"))
    } else {
        p <- p + theme(strip.text.x = element_blank(), strip.background.x = element_blank())
    }
    if(!show_legend) {
        p <- p + theme(legend.position = "none")
    } else {
        p <- p + theme(legend.position = "bottom", legend.box = "horizontal", legend.margin = margin(t = -5, b = 0)) +
            guides(color = guide_legend(title = "Platform", nrow = 1, order = 1),
                   linetype = guide_legend(title = "Target scope", nrow = 1, order = 2, keywidth = unit(1, "lines"),
                                           override.aes = list(color = "black", linewidth = 0.5)),
                   shape = guide_legend(title = "Target scope", nrow = 1, order = 2))
    }
    p
}

p_pdp_all <- plot_1d_curve(df_pdp_plot, y_col = "yhat", y_label = "Partial dependence", show_legend = FALSE, show_strip_x = TRUE)
p_ale_all <- plot_1d_curve(df_ale_plot, y_col = "ale", y_label = "Accumulated local effects", show_legend = TRUE, show_strip_x = FALSE)

# 5. Assemble extended-data panels ----
top_row <- ggarrange(p_corr, sp_wilcox, ncol = 2, widths = c(1, 1.5),
                    labels = c("a", "b"), font.label = label_style, label.x = 0, label.y = 1, hjust = -0.2, vjust = 1.2)
fig_supp_combined <- ggarrange(top_row, p_pdp_all, p_ale_all, nrow = 3, heights = c(1.18, 1.05, 1),
                              labels = c("", "c", "d"), font.label = label_style, label.x = 0, label.y = 1, hjust = -0.2, vjust = 1.2)
ggsave("figures/ed6_distortion_features.pdf", fig_supp_combined, width = 10, height = 11)
ggsave("figures/ed6_distortion_features.png", fig_supp_combined, width = 10, height = 11, dpi = 600, bg = "white")

# 6. Source data ----
# ST: Wilcoxon Significance Matrix (P-values)
format_pvalue <- function(x) {
  x <- as.numeric(x)
  out <- character(length(x))
  for (i in seq_along(x)) {
    if (is.na(x[i])) {
      out[i] <- NA_character_
    } else if (x[i] < 0.001) {
      out[i] <- formatC(x[i], format = "e", digits = 2) # e.g., 1.23e-04
    } else {
      out[i] <- formatC(x[i], format = "f", digits = 3) # e.g., 0.045
    }
  }
  return(out)
}

st_wilcox_structured <- wilcox_wide_df %>%
  inner_join(physchem_dict %>% select(Feature, Property, Category), by = "Feature") %>%
  select(-Feature) %>%
  select(Platform, Dataset, Category, Property, everything()) %>%
  rename_with(~ ifelse(grepl("^0\\.", .), paste0("Threshold ", as.numeric(.) * 100, "%"), .),
              -c(Category, Property, Dataset, Platform)) %>%
  rename_with(~ ifelse(grepl("^\\d+%$", .), paste0("Threshold ", .), .),
              -c(Category, Property, Dataset, Platform)) %>%
  mutate(across(starts_with("Threshold"), format_pvalue))

ed6_source_data <- list(
    "Property_Correlation" = plot_cor_df,
    "Wilcoxon_P" = st_wilcox_structured,
    "PDP_Summary" = pdp_summary,
    "ALE_Summary" = ale_summary
)
write.xlsx(ed6_source_data, "tables/SourceData_EDFigure6.xlsx", overwrite = TRUE)

message("Extended Data Figure 6 completed: figures/ed6_distortion_features.pdf；source data: tables/SourceData_ExtendedDataFigure6.xlsx")
