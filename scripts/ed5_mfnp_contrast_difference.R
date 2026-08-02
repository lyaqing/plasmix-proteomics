# Extended Data Figure 5 | M/F and N/P contrast differences

# 0. Setup ----
source("scripts/_project_setup.R")
use_packages(c("data.table", "tidyverse", "readxl", "ggpubr", "ggridges", "smplot2", "openxlsx", "grid", "metap"), c("clusterProfiler", "org.Hs.eg.db"))
source("utils/figure_style.R")
plasmix_theme <- if (is.function(theme_plasmix)) theme_plasmix() else theme_plasmix
label_style <- list(size = 12, face = "bold")

# 1. Inputs ----
required_inputs <- c("results/dea_df_multi.tsv.gz", "data/feature_metadata.tsv.gz")
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) stop("Missing Extended Data Figure 5 inputs: ", paste(missing_inputs, collapse = ", "))
dea_df_multi <- fread("results/dea_df_multi.tsv.gz")
feat_meta <- fread("data/feature_metadata.tsv.gz") %>% filter(!Is_Protein_Group, !Is_Unknown)
uniprot <- read_tsv(upstream_path("references", "uniprotkb_AND_model_organism_9606_2026_04_01.tsv.gz"), show_col_types = FALSE) %>%
    transmute(Protein.ID = Entry, Gene.Symbol = `Gene Names (primary)`) %>% distinct(Protein.ID, .keep_all = TRUE)

# 2. Default M/F and N/P analysis tiers ----
dea_filter <- dea_df_multi %>%
    filter(
        Pair %in% c("M/F", "N/P"),
        (Platform %in% c("DIA", "AAG") & DataTier == "Baseline") |
            (Platform %in% c("OLK", "SOM", "NLS") & DataTier == "Calibrated")
    )

# 3. Panel a: M/F and N/P effect distributions ----
plat_order <- c("SOM", "NLS", "AAG", "DIA", "OLK")

ridge_plot_data <- dea_filter %>%
    mutate(Pair = factor(Pair, levels = c("M/F", "N/P")), Platform = factor(Platform, levels = plat_order))

ridge_order <- ridge_plot_data %>%
    group_by(Pair, Platform) %>%
    summarise(d025 = quantile(logFC, .025, na.rm = TRUE), d975 = quantile(logFC, .975, na.rm = TRUE), .groups = "drop") %>%
    mutate(delta95 = d975 - d025, Y_pos = as.numeric(Platform), label = paste0("Delta[95] == '", sprintf("%.2f", delta95), "'"))

p_ridge <- ggplot(ridge_plot_data, aes(logFC, Platform, fill = Platform, color = Platform)) +
    geom_density_ridges(alpha = .05, scale = 1.2, rel_min_height = .01, linewidth = .5) +
    geom_segment(data = ridge_order, aes(x = d025, xend = d025, y = Y_pos, yend = Y_pos + .6, color = Platform), inherit.aes = FALSE, linetype = "dotted", alpha = .8, linewidth = .5) +
    geom_segment(data = ridge_order, aes(x = d975, xend = d975, y = Y_pos, yend = Y_pos + .6, color = Platform), inherit.aes = FALSE, linetype = "dotted", alpha = .8, linewidth = .5) +
    geom_segment(data = ridge_order, aes(x = d025, xend = d975, y = Y_pos + .6, yend = Y_pos + .6, color = Platform), inherit.aes = FALSE, linewidth = .5, arrow = arrow(ends = "both", length = unit(.05, "inches"), type = "closed")) +
    geom_text(data = ridge_order, aes(x = (d025 + d975) / 2, y = Y_pos + .8, label = label, color = Platform), inherit.aes = FALSE, parse = TRUE, size = 3, fontface = "bold", family = "serif") +
    facet_grid(Pair ~ .) +
    scale_fill_manual(values = platform_color, drop = FALSE) +
    scale_color_manual(values = platform_color, drop = FALSE) +
    scale_x_continuous(limits = c(-4, 4), expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(drop = FALSE, expand = expansion(add = c(.2, 1.15))) +
    labs(x = expression(Observed~log[2]~fold~change), y = NULL) + plasmix_theme +
    theme(legend.position = "none", panel.grid.major.y = element_blank())

# 4. Panel b: Precision of differential features ----
area_colors <- c(
    "Non-significant" = "#EDEDED",
    "Precision-rejected" = "#DC0000",
    "Precision-verified" = "#00A087"
)

prepare_mixed_bins <- function(df) {
    df_nls <- df %>%
        filter(grepl("NLS", Batch)) %>%
        mutate(
            Ratio = 2^abs(logFC),
            Ratio_Bin = cut(Ratio, breaks = c(seq(1.0, 2.01, by = 0.1), Inf), include.lowest = TRUE)
        )

    df_others <- df %>%
        filter(!grepl("NLS", Batch)) %>%
        mutate(
            Ratio = 2^abs(logFC),
            Ratio_Bin = cut(Ratio, breaks = c(seq(1.0, 2.01, by = 0.025), Inf), include.lowest = TRUE)
        )

    bind_rows(df_nls, df_others) %>%
        group_by(Platform, Pair, Ratio_Bin, Classification) %>%
        summarize(Count = n(), .groups = "drop") %>%
        group_by(Platform, Pair) %>%
        complete(Ratio_Bin, Classification, fill = list(Count = 0)) %>%
        group_by(Platform, Pair, Ratio_Bin) %>%
        mutate(Pct = Count / sum(Count)) %>%
        ungroup() %>%
        mutate(
            Bin_Mid = as.numeric(gsub("\\(|\\[|\\,.*", "", Ratio_Bin)),
            Bin_Mid = ifelse(is.na(Bin_Mid), 2.01, Bin_Mid)
        )
}

global_retention <- dea_filter %>%
    filter(adj.P.Val < 0.05, abs(logFC) >= log2(1.2)) %>%
    group_by(Platform, Pair) %>%
    summarize(
        total_sig = n(),
        verified_n = sum(Classification %in% c("Small-magnitude", "Verified-DEP")),
        pct_val = verified_n / total_sig * 100,
        .groups = "drop"
    ) %>%
    mutate(
        label_text = paste0(sprintf("%.1f", pct_val), "% verified"),
        x_pos = 1.98,
        y_pos = 5,
        Platform = factor(Platform, levels = c("OLK", "SOM", "DIA", "AAG", "NLS"))
    )

pct_area_data <- dea_filter %>%
    mutate(Classification = case_when(
        Classification %in% c("Small-magnitude", "Verified-DEP") ~ "Precision-verified",
        TRUE ~ Classification
    )) %>%
    prepare_mixed_bins() %>%
    mutate(
        Platform = factor(Platform, levels = c("OLK", "SOM", "DIA", "AAG", "NLS")),
        Pair = factor(Pair, levels = c("M/F", "N/P")),
        Classification = factor(
            Classification,
            levels = c("Non-significant", "Precision-rejected", "Precision-verified")
        )
    )

p_pct_area <- ggplot(pct_area_data, aes(x = Bin_Mid, y = Pct * 100, fill = Classification)) +
    geom_area(alpha = 1, size = 0.05, color = "white") +
    facet_grid(Pair ~ Platform) +
    scale_fill_manual("Classification", values = area_colors) +
    geom_vline(xintercept = 1.2, linetype = "dotted", linewidth = 0.5, color = "black", alpha = 0.8) +
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100), expand = c(0, 0)) +
    scale_x_continuous(breaks = c(1.0, 1.2, 1.5, 2.0), limits = c(1.0, 2.0), expand = c(0, 0)) +
    labs(x = "Absolute fold change", y = "Classification (%)") +
    plasmix_theme +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major = element_blank(),
        strip.text = element_text(face = "bold", size = 8.5, color = "black", margin = margin(3, 3, 3, 3)),
        strip.background = element_rect(fill = "grey95", color = NA),
        panel.spacing.x = unit(0.55, "lines"),
        axis.text.x = element_text(hjust = c(0.1, 0.3, 0.3, 0.9)),
        axis.text.y = element_text(vjust = c(0.1, 0.5, 0.5, 0.5, 0.9)),
        legend.position = "bottom",
        legend.margin = margin(-7, 0, 0, 0),
        legend.key.spacing.x = unit(0.7, "lines"),
        legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 8, vjust = 0.5)
    ) +
    geom_text(
        data = global_retention,
        aes(x = x_pos, y = y_pos, label = label_text),
        inherit.aes = FALSE, color = "black", size = 3, hjust = 1, vjust = 0
    )

# 5. Panel c: Cross-platform consensus proteins ----
selected_pairs <- c(
    "M/F", "Y/F", "M/X", "M/P", "P/F", "Y/X", "X/F", "Y/P", "M/Y", "P/X",
    "N/M", "N/Y", "N/P", "N/X", "N/F"
)
fc_th <- 1.2
p_th <- 0.05

consensus_base <- dea_df_multi %>%
    filter(
        Pair %in% selected_pairs,
        Platform != "AAG",
        (Platform == "DIA" & DataTier == "Baseline") |
            (Platform %in% c("OLK", "SOM", "NLS") & DataTier == "Calibrated")
    )

df_verified <- consensus_base %>% filter(Classification == "Verified-DEP")

dea_merge <- df_verified %>%
    mutate(Direction_temp = ifelse(logFC > 0, "Up", "Down")) %>%
    group_by(Pair, UniqueID, Direction_temp) %>%
    mutate(
        n_plat_dir = n_distinct(Platform),
        n_batch_dir = n()
    ) %>%
    filter(n_plat_dir >= 2, n_batch_dir >= 4) %>%
    ungroup() %>%
    group_by(Pair, UniqueID) %>%
    filter(n_batch_dir == max(n_batch_dir)) %>%
    summarize(
        combined_logFC = median(logFC, na.rm = TRUE),
        combined_p = tryCatch({
            p_vals <- P.Value[!is.na(P.Value)]
            p_vals[p_vals == 0] <- 1e-300
            if (length(p_vals) > 0) sumlog(p_vals)$p else NA_real_
        }, error = function(e) NA_real_),
        n_batches = n(),
        .groups = "drop"
    ) %>%
    mutate(
        adj_p = p.adjust(combined_p, method = "fdr"),
        Label = case_when(
            combined_logFC >= log2(fc_th) & adj_p < p_th ~ "Up",
            combined_logFC <= -log2(fc_th) & adj_p < p_th ~ "Down",
            TRUE ~ "NotSig"
        )
    ) %>%
    ungroup()

fwrite(dea_merge, "results/ed5_dea_consensus.tsv.gz", sep = "\t", na = "NA")

avglogfc <- dea_merge %>% rename(avglogFC = combined_logFC, Freq = n_batches)

pair_order <- avglogfc %>%
    filter(Label != "NotSig") %>%
    count(Pair, name = "total_n") %>%
    arrange(desc(total_n)) %>%
    pull(Pair)
final_pair_levels <- c(pair_order, setdiff(selected_pairs, pair_order))

dd1 <- avglogfc %>%
    filter(Label != "NotSig", Pair %in% selected_pairs) %>%
    distinct() %>%
    merge(uniprot %>% select(Protein.ID, Gene.Symbol), by.x = "UniqueID", by.y = "Protein.ID", all.x = TRUE) %>%
    mutate(
        Gene.Symbol = ifelse(is.na(Gene.Symbol), UniqueID, Gene.Symbol),
        Pair = factor(Pair, levels = final_pair_levels)
    )

bar <- dd1 %>%
    group_by(Pair) %>%
    summarize(
        Q1 = quantile(avglogFC, 0.25, na.rm = TRUE),
        Q3 = quantile(avglogFC, 0.75, na.rm = TRUE),
        IQR_val = Q3 - Q1,
        upper_fence = Q3 + 3 * IQR_val,
        lower_fence = Q1 - 3 * IQR_val,
        max = {
            up_vals <- avglogFC[avglogFC > 0]
            normal_up <- up_vals[up_vals <= upper_fence]
            if (length(normal_up) > 0) max(normal_up) * 1.01 else if (length(up_vals) > 0) max(up_vals) * 1.01 else 0.1
        },
        min = {
            down_vals <- avglogFC[avglogFC < 0]
            normal_down <- down_vals[down_vals >= lower_fence]
            if (length(normal_down) > 0) min(normal_down) * 1.01 else if (length(down_vals) > 0) min(down_vals) * 1.01 else -0.1
        },
        .groups = "drop"
    ) %>%
    mutate(Pair = factor(Pair, levels = final_pair_levels))

p_consensus_base <- ggplot(dd1, aes(x = Pair, y = avglogFC)) +
    geom_col(data = bar, aes(x = Pair, y = min), fill = "#dcdcdc", alpha = 0.6, width = 0.8) +
    geom_col(data = bar, aes(x = Pair, y = max), fill = "#dcdcdc", alpha = 0.6, width = 0.8) +
    geom_jitter(aes(size = Freq, alpha = Freq, color = Label), width = 0.3) +
    scale_size_continuous(range = c(0.1, 3), breaks = c(3, 5, 7, 10)) +
    scale_alpha_continuous(trans = "reverse", range = c(0.9, 0.1), breaks = c(3, 5, 7, 10)) +
    scale_color_manual(values = c(Down = "#0c6399", Up = "#d21613")) +
    labs(y = expression(log[2]~FC), x = NULL, color = "Significant") +
    guides(color = guide_legend(override.aes = list(size = 3)))

count_data <- dd1 %>% count(Pair, Label, name = "count")
text_data <- bar %>%
    left_join(count_data, by = "Pair") %>%
    filter(!is.na(Label)) %>%
    mutate(
        y_pos = ifelse(Label == "Up", max, min),
        Pair = factor(Pair, levels = final_pair_levels)
    )

p_dea <- p_consensus_base +
    geom_tile(data = dd1 %>% distinct(Pair), aes(x = Pair, y = 0, fill = Pair), color = "black", height = 0.5, alpha = 1, show.legend = FALSE) +
    scale_fill_manual(values = sm_palette(15)) +
    geom_text(data = dd1 %>% distinct(Pair), aes(x = Pair, y = 0, label = Pair), color = "white", fontface = "bold", size = 2.5) +
    geom_text(data = text_data, aes(x = Pair, y = y_pos, label = paste0("n = ", count), color = Label, vjust = ifelse(Label == "Up", -0.5, 1.5)), size = 2, fontface = "bold") +
    plasmix_theme +
    theme(
        panel.grid.major = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.box = "vertical",
        legend.title = element_text(size = 7, face = "bold", vjust = 0.5),
        legend.text = element_text(size = 7, vjust = 0.5),
        plot.margin = margin(5, 5, 7, 5)
    ) +
    scale_y_continuous(limits = c(-5.2, 4.5), expand = c(0, 0), breaks = c(-4, -2, 0, 2, 4))

# 6. Panel d: GO enrichment of M/F and N/P consensus proteins ----
plot_input <- dea_merge %>%
    filter(Pair %in% c("M/F", "N/P"), Label != "NotSig") %>%
    transmute(Pair, UniProtID = sub("_.*", "", UniqueID)) %>%
    distinct()

# Define a shared measured-protein universe from the feature metadata.
universe_uniprot <- feat_meta %>%
    transmute(UniProtID = as.character(UniProtID)) %>%
    filter(!is.na(UniProtID), UniProtID != "") %>%
    pull(UniProtID) %>% unique()

all_uniprot_for_mapping <- unique(c(plot_input$UniProtID, universe_uniprot))
id_map <- clusterProfiler::bitr(all_uniprot_for_mapping, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db)

final_data <- plot_input %>%
    left_join(id_map, by = c("UniProtID" = "UNIPROT")) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(Pair, ENTREZID, .keep_all = TRUE)

universe_entrez <- id_map %>%
    filter(UNIPROT %in% universe_uniprot) %>%
    pull(ENTREZID) %>% unique()

foreground_outside_universe <- final_data %>% filter(!ENTREZID %in% universe_entrez)
if (nrow(foreground_outside_universe) > 0) {
    warning(nrow(foreground_outside_universe), " foreground Pair–Entrez records were not present in the measured-protein universe.")
}

message("GO foreground: M/F = ", n_distinct(final_data$ENTREZID[final_data$Pair == "M/F"]),
        "; N/P = ", n_distinct(final_data$ENTREZID[final_data$Pair == "N/P"]),
        "; shared measured-protein universe = ", length(universe_entrez))
# GO foreground: M/F = 100; N/P = 560; shared measured-protein universe = 12548

run_go <- function(ont_use) {
    comp_res <- clusterProfiler::compareCluster(
        ENTREZID ~ Pair, data = final_data, fun = "enrichGO", OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "ENTREZID", ont = ont_use, universe = universe_entrez, pvalueCutoff = 0.05, pAdjustMethod = "BH"
    )
    if (is.null(comp_res) || !nrow(as.data.frame(comp_res))) return(tibble())
    comp_res_sim <- clusterProfiler::simplify(comp_res, cutoff = 0.65, by = "p.adjust", select_fun = min)
    as_tibble(as.data.frame(comp_res_sim)) %>% mutate(Ontology = ont_use)
}

# Complete simplified BP, CC and MF results are retained for Source Data.
enrich_summary_all <- map_dfr(c("BP", "CC", "MF"), run_go)

enrich_summary <- enrich_summary_all %>%
    select(Ontology, Cluster, ID, Description, GeneRatio, BgRatio, p.adjust, Count) %>%
    mutate(
        Description = paste0(toupper(substr(Description, 1, 1)), substr(Description, 2, nchar(Description))),
        GR_val = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)),
        BR_val = as.numeric(sub("/.*", "", BgRatio)) / as.numeric(sub(".*/", "", BgRatio)),
        FoldEnrichment = GR_val / BR_val
    )

# Display up to seven BP terms per contrast after semantic-similarity simplification.
bp_n_show <- 7L
enrich_summary_plot <- enrich_summary %>%
    filter(Ontology == "BP") %>%
    group_by(Cluster) %>%
    arrange(p.adjust, desc(Count), desc(FoldEnrichment), .by_group = TRUE) %>%
    slice_head(n = bp_n_show) %>%
    ungroup() %>%
    mutate(Description_key = paste(Description, Cluster, sep = "|||"))

# Factor levels are reversed by adjusted P value so the most significant term appears at the top.
description_levels <- enrich_summary_plot %>%
    arrange(Cluster, desc(p.adjust), Count, FoldEnrichment) %>%
    pull(Description_key) %>% unique()

enrich_summary_plot <- enrich_summary_plot %>%
    mutate(Description_key = factor(Description_key, levels = description_levels))

ontology_labs <- c(BP = "Biological process", CC = "Cellular component", MF = "Molecular function")

p_go <- ggplot(enrich_summary_plot, aes(x = FoldEnrichment, y = reorder(Description, FoldEnrichment))) +
    geom_segment(aes(x = 0, xend = FoldEnrichment, y = Description, yend = Description), color = "gray90", linewidth = 0.8) +
    geom_point(aes(size = Count, fill = -log10(p.adjust)), alpha = 0.8, color = "black", shape = 21) +
    scale_fill_viridis_c(option = "rocket", name = expression(bold(-log[10](P[adj]))), direction = -1, end = 0.8, breaks = c(5, 15, 25)) +
    facet_grid(Ontology ~ Cluster, scales = "free_y", space = "free", labeller = labeller(Ontology = as_labeller(ontology_labs))) +
    scale_size_continuous(range = c(1, 5), name = "Protein\ncount", breaks = c(5, 10, 20, 40, 60, 80)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 90)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, NA)) +
    labs(x = "Fold enrichment", y = NULL) +
    plasmix_theme +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.background = element_rect(fill = "gray20"),
        strip.text = element_text(color = "white", face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(0.5, "lines"),
        legend.title = element_text(size = 7, face = "bold", vjust = 0.5),
        legend.text = element_text(size = 7, vjust = 0.5)
    )

# 7. Assemble and export ----
right_column <- ggarrange(p_pct_area, p_dea, ncol = 1, heights = c(1, 1), labels = c("b", "c"),
                          font.label = label_style, label.x = 0, label.y = 1, hjust = -0.2, vjust = 1.2)
top_block <- ggarrange(p_ridge, right_column, nrow = 1, widths = c(0.3, 0.7), labels = c("a", ""),
                      font.label = label_style, label.x = 0, label.y = 1, hjust = -0.2, vjust = 1.2)
final_figure <- ggarrange(top_block, p_go, ncol = 1, heights = c(1, 0.4), labels = c("", "d"),
                          font.label = label_style, label.x = 0, label.y = 1, hjust = -0.2, vjust = 1.2)
ggsave("figures/ed5_mfnp_contrast_difference.pdf", final_figure, width = 10, height = 8)
ggsave("figures/ed5_mfnp_contrast_difference.png", final_figure, width = 10, height = 8, dpi = 600, bg = "white")

# 8. Source data ----
format_pvalue <- function(x) {
    x <- as.numeric(x)
    ifelse(
        is.na(x),
        NA_character_,
        ifelse(x < 0.001, formatC(x, format = "e", digits = 2), formatC(x, format = "f", digits = 3))
    )
}

st_dep_batch_stat <- dea_df_multi %>%
    filter(Pair %in% c("M/F", "N/P")) %>%
    count(Pair, Batch, DataTier, ProcessLevel, Classification, Platform, name = "Count") %>%
    mutate(
        Pair = factor(Pair, levels = c("M/F", "N/P")),
        DataTier = factor(DataTier, levels = c("Baseline", "Calibrated", "Reshaped")),
        ProcessLevel = case_when(
            Platform == "OLK" & DataTier == "Baseline" ~ "ExtNPX",
            Platform == "NLS" & DataTier == "Baseline" ~ "ICNorm",
            Batch %in% c("SOM_P1_B1", "SOM_P1_B2") & DataTier == "Reshaped" ~ "ANML",
            Batch %in% c("SOM_P2_B1", "SOM_P2_B2") & DataTier == "Baseline" ~ "ReadoutNorm",
            Batch %in% c("SOM_P2_B1", "SOM_P2_B2") & DataTier == "Calibrated" ~ "PlateNorm",
            Batch %in% c("SOM_P2_B1", "SOM_P2_B2") & DataTier == "Reshaped" ~ "SampleNorm",
            TRUE ~ ProcessLevel
        )
    ) %>%
    arrange(Pair, Batch, DataTier) %>%
    rename(`Quantification metric` = ProcessLevel, `Data tier` = DataTier, Contrast = Pair) %>%
    pivot_wider(names_from = Classification, values_from = Count)

st_dep_consensus <- dea_merge %>%
    filter(Pair %in% c("M/F", "N/P")) %>%
    transmute(
        Contrast = Pair,
        Identifier = UniqueID,
        `Log2 FC` = round(combined_logFC, 3),
        `P-value` = format_pvalue(combined_p),
        `Adjusted P-value` = format_pvalue(adj_p),
        Direction = Label,
        Observations = n_batches
    )

st_enrich_summary <- enrich_summary %>%
    transmute(
        Ontology,
        Contrast = Cluster,
        `GO term ID` = ID,
        Description,
        `Gene ratio` = GeneRatio,
        `Background ratio` = BgRatio,
        `Adjusted P-value` = format_pvalue(p.adjust),
        Count,
        `Fold enrichment` = round(FoldEnrichment, 2)
    )

# Use one worksheet per panel-level dataset or supporting analysis table.
source_data <- list(
    "FC_distributions" = dea_filter,
    "FC_range_summary" = bind_rows(
        dea_filter %>% filter(Pair == "M/F") %>% group_by(Pair, Platform) %>%
            summarize(P2.5 = quantile(logFC, 0.025, na.rm = TRUE),
                      P97.5 = quantile(logFC, 0.975, na.rm = TRUE), .groups = "drop"),
        dea_filter %>% filter(Pair == "N/P") %>% group_by(Pair, Platform) %>%
            summarize(P2.5 = quantile(logFC, 0.025, na.rm = TRUE),
                      P97.5 = quantile(logFC, 0.975, na.rm = TRUE), .groups = "drop")
    ),
    "Precision_filtering" = pct_area_data,
    "Consensus_all_contrasts" = dea_merge,
    "DEP_classification" = st_dep_batch_stat,
    "Consensus_MF_NP" = st_dep_consensus,
    "GO_enrichment" = st_enrich_summary
)
write.xlsx(
    source_data,
    "tables/SourceData_EDFigure5.xlsx",
    overwrite = TRUE,
    keepNA = TRUE,
    na.string = "NA"
)
message("Extended Data Figure 5 completed: figures/ed5_mfnp_contrast_difference.pdf；source data: tables/SourceData_EDFigure5.xlsx")
