# Extended Data Figure 9 | Abundance by integration outcome

# 0. Setup ----
source("scripts/_project_setup.R")
use_packages(c("data.table", "tidyverse", "ggpubr", "openxlsx", "showtext"))
source("utils/figure_style.R")
plasmix_theme <- if (is.function(theme_plasmix)) theme_plasmix() else theme_plasmix
showtext_auto(); showtext_opts(dpi = 600)
label_style <- list(size = 12, face = "bold")
paths <- c(fig6_inputs = "results/fig6_extended_data_inputs.rds", profiles = "data/protein_profiles_long.tsv.gz",
           feature_metadata = "data/feature_metadata.tsv.gz", physchem = "data/physchem_matrix.tsv.gz")
check_inputs(paths)

# 1. Reconstruct final outcome groups and any annotations absent from the frozen Figure 6 object ----
fig6_ed_inputs <- readRDS(paths["fig6_inputs"])
if (is.null(fig6_ed_inputs$consensus_voting)) stop("fig6_extended_data_inputs.rds lacks consensus_voting.")
consensus_voting <- fig6_ed_inputs$consensus_voting

if (!("Native_Low_Side_Rank" %in% names(consensus_voting))) {
    if (is.null(fig6_ed_inputs$method_feature_results) || !file.exists(paths["profiles"])) {
        stop("Native_Low_Side_Rank is absent and cannot be reconstructed without method_feature_results and protein_profiles_long.tsv.gz.")
    }
    balanced_feature_pairs <- fig6_ed_inputs$method_feature_results %>% filter(Design == "Balanced") %>%
        distinct(Detailed_Type, Batch1, Batch2, UniqueID)
    selected_batches <- union(balanced_feature_pairs$Batch1, balanced_feature_pairs$Batch2)
    native_rank_by_batch <- fread(paths["profiles"]) %>% as_tibble() %>%
        filter(Batch %in% selected_batches,
               (Platform == "DIA" & DataTier == "Baseline" & ProcessLevel == "Intensity") |
               (Platform == "OLK" & DataTier == "Calibrated" & str_detect(ProcessLevel, regex("NPX", ignore_case = TRUE))) |
               (Platform == "SOM" & DataTier == "Calibrated" & ProcessLevel == "Calibrate"), is.finite(Value)) %>%
        group_by(Batch, UniqueID) %>% summarize(Median_Val = median(Value, na.rm = TRUE), .groups = "drop") %>%
        group_by(Batch) %>% mutate(Batch_Rank_pct = percent_rank(Median_Val) * 100) %>% ungroup() %>%
        select(Batch, UniqueID, Batch_Rank_pct)
    native_low_side_rank <- balanced_feature_pairs %>%
        left_join(native_rank_by_batch %>% rename(Batch1 = Batch, Rank_B1 = Batch_Rank_pct), by = c("Batch1", "UniqueID")) %>%
        left_join(native_rank_by_batch %>% rename(Batch2 = Batch, Rank_B2 = Batch_Rank_pct), by = c("Batch2", "UniqueID")) %>%
        mutate(Pair_Low_Rank = pmin(Rank_B1, Rank_B2, na.rm = TRUE), Pair_Low_Rank = if_else(is.infinite(Pair_Low_Rank), NA_real_, Pair_Low_Rank)) %>%
        group_by(Detailed_Type, UniqueID) %>%
        summarize(Native_Low_Side_Rank = median(Pair_Low_Rank, na.rm = TRUE), .groups = "drop") %>%
        mutate(Native_Low_Side_Rank = if_else(is.nan(Native_Low_Side_Rank), NA_real_, Native_Low_Side_Rank))
    consensus_voting <- consensus_voting %>% left_join(native_low_side_rank, by = c("Detailed_Type", "UniqueID"))
}

if (!("BloodConc_log10_pgml" %in% names(consensus_voting)) && !is.null(fig6_ed_inputs$consensus_voting_class)) {
    cached_hpa <- fig6_ed_inputs$consensus_voting_class %>%
        select(any_of(c("Detailed_Type", "UniqueID", "UniProtID", "BloodConc_log10_pgml"))) %>%
        distinct(Detailed_Type, UniqueID, .keep_all = TRUE)
    join_keys <- intersect(c("Detailed_Type", "UniqueID", "UniProtID"), names(consensus_voting))
    consensus_voting <- consensus_voting %>% left_join(cached_hpa, by = join_keys)
}

if (!("BloodConc_log10_pgml" %in% names(consensus_voting)) && file.exists(paths["physchem"])) {
    physchem_matrix <- fread(paths["physchem"]) %>% as_tibble()
    if (!("UniProtID" %in% names(consensus_voting)) && file.exists(paths["feature_metadata"])) {
        feature_metadata <- fread(paths["feature_metadata"]) %>% as_tibble()
        entry_column <- intersect(c("UniProtID", "UniProt_ID", "UniProt", "Uniprot", "Entry", "ProteinID"), names(feature_metadata))[1]
        if (length(entry_column) && !is.na(entry_column)) {
            id_map <- feature_metadata %>% transmute(UniqueID, UniProtID = as.character(.data[[entry_column]])) %>% distinct(UniqueID, .keep_all = TRUE)
            consensus_voting <- consensus_voting %>% left_join(id_map, by = "UniqueID")
        }
    }
    if (all(c("Entry", "BloodConc_log10_pgml") %in% names(physchem_matrix)) && "UniProtID" %in% names(consensus_voting)) {
        hpa_map <- physchem_matrix %>% transmute(UniProtID = as.character(Entry), BloodConc_log10_pgml = as.numeric(BloodConc_log10_pgml)) %>%
            distinct(UniProtID, .keep_all = TRUE)
        consensus_voting <- consensus_voting %>% left_join(hpa_map, by = "UniProtID")
    }
}

required_columns <- c("Detailed_Type", "UniqueID", "Success_Rate", "Native_Low_Side_Rank", "BloodConc_log10_pgml")
missing_columns <- setdiff(required_columns, names(consensus_voting))
if (length(missing_columns)) stop("consensus_voting lacks required columns after reconstruction: ", paste(missing_columns, collapse = ", "))

scenario_levels <- c("Intra-DIA", "Intra-OLK", "Intra-SOM", "DIA-OLK", "DIA-SOM", "OLK-SOM")
consensus_voting_class <- consensus_voting %>% filter(Detailed_Type %in% scenario_levels) %>% group_by(Detailed_Type) %>%
    mutate(Threshold_Q75 = quantile(Success_Rate, 0.75, na.rm = TRUE),
           Integration_Status = case_when(Success_Rate == 0 ~ "Failed", Success_Rate >= Threshold_Q75 ~ "Success", TRUE ~ "Intermediate")) %>%
    ungroup() %>% filter(Integration_Status %in% c("Success", "Failed")) %>%
    mutate(Detailed_Type = factor(Detailed_Type, levels = scenario_levels),
           Integration_Status = factor(Integration_Status, levels = c("Success", "Failed")))

# 2. Preserve the original distribution-plot form ----
plot_abundance_distribution <- function(data, metric_column, x_label) {
    plot_data <- data %>% filter(is.finite(.data[[metric_column]]))
    median_data <- plot_data %>% group_by(Detailed_Type, Integration_Status) %>%
        summarize(Median = median(.data[[metric_column]], na.rm = TRUE), N = n(), .groups = "drop")
    status_colors <- c("Success" = "#3171b8", "Failed" = "#DF6B6A")
    plot_object <- ggplot(plot_data, aes(x = .data[[metric_column]], fill = Integration_Status, color = Integration_Status)) +
        geom_histogram(aes(y = after_stat(density)), position = "identity", alpha = 0.3, bins = 30, color = NA) +
        geom_density(fill = NA, linewidth = 0.5) + geom_rug(alpha = 0.5, sides = "b", length = unit(0.04, "npc")) +
        geom_vline(data = median_data, aes(xintercept = Median, color = Integration_Status), linetype = "dashed", linewidth = 0.5) +
        facet_wrap(~Detailed_Type, nrow = 1) + scale_fill_manual(values = status_colors) + scale_color_manual(values = status_colors) +
        labs(x = x_label, y = "Density", fill = "Integration", color = "Integration") + plasmix_theme +
        theme(panel.grid.major.x = element_blank(), strip.background = element_blank(), legend.key.size = unit(0.6, "lines"))
    list(plot = plot_object, values = plot_data, medians = median_data)
}

rank_result <- plot_abundance_distribution(consensus_voting_class, "Native_Low_Side_Rank", "Pair-low-side native measurement rank (%)")
concentration_result <- plot_abundance_distribution(consensus_voting_class, "BloodConc_log10_pgml",
                                                     expression("HPA concentration (" * log[10] * " pg/mL)"))

# 3. Assemble and export ----
figure_ed9 <- ggarrange(rank_result$plot, concentration_result$plot, ncol = 1, nrow = 2, labels = c("a", "b"),
                                font.label = label_style, common.legend = TRUE, legend = "right")
ggsave("figures/ed9_abundance_integration_outcomes.pdf", figure_ed9, width = 10, height = 4)
ggsave("figures/ed9_abundance_integration_outcomes.png", figure_ed9, width = 10, height = 4, dpi = 600, bg = "white")

source_data <- list(
    ED9_classification = consensus_voting_class %>% mutate(across(where(is.factor), as.character)) %>%
        arrange(Detailed_Type, Integration_Status, desc(Success_Rate), UniqueID),
    ED9a_native_rank = rank_result$values %>% mutate(across(where(is.factor), as.character)) %>%
        select(Detailed_Type, UniqueID, any_of("UniProtID"), Success_Rate, Threshold_Q75, Integration_Status, Native_Low_Side_Rank),
    ED9a_medians = rank_result$medians %>% mutate(across(where(is.factor), as.character)),
    ED9b_HPA_concentration = concentration_result$values %>% mutate(across(where(is.factor), as.character)) %>%
        select(Detailed_Type, UniqueID, any_of("UniProtID"), Success_Rate, Threshold_Q75, Integration_Status, BloodConc_log10_pgml),
    ED9b_medians = concentration_result$medians %>% mutate(across(where(is.factor), as.character))
)
write.xlsx(source_data, "tables/SourceData_EDFigure9.xlsx", overwrite = TRUE, keepNA = TRUE, na.string = "NA")
message("Extended Data Figure 9 and its source data were exported.")
