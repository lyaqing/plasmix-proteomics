# Extended Data Figure 3 | SomaScan reshaping effects

# 0. Setup ----
source("scripts/_project_setup.R")
use_packages(c("data.table", "tidyverse", "readxl", "ggpubr", "openxlsx", "showtext", "limma"), "SomaDataIO")
source("utils/benchmark_metrics.R")
source("utils/imputation.R")
source("utils/differential_analysis.R")
source("utils/figure_style.R")
source("utils/feature_mapping.R")
plasmix_theme <- if (is.function(theme_plasmix)) theme_plasmix() else theme_plasmix
showtext_auto(); showtext_opts(dpi = 600)
label_style <- list(size = 12, face = "bold")

# 1. Inputs ----
target_batches <- c("SOM_P1_B1", "SOM_P1_B2", "SOM_P2_B1", "SOM_P2_B2")
batch_colors <- c("SOM_P1_B1" = "#d162b5", "SOM_P1_B2" = "#830062", "SOM_P2_B1" = "#f0a22e", "SOM_P2_B2" = "#c85203")

paths <- c(metadata = "data/study_metadata.xlsx", feature_metadata = "data/feature_metadata.tsv.gz",
           profiles = "data/protein_profiles_long.tsv.gz", detection = "results/detection_status.tsv.gz",
           replicate_selection = "results/dea_replicate_selection.tsv.gz",
           annotation_11k = upstream_path("references", "SomaScan_11K_Annotated_Content.xlsx"),
           annotation_9k = upstream_path("references", "Customer Facing - Protein Prep SOMAmer Annotations SEPT2025.xlsx"))
missing_inputs <- paths[!file.exists(paths)]
if (length(missing_inputs) > 0) stop("Missing input files: ", paste(missing_inputs, collapse = ", "), call. = FALSE)

meta_sample <- read_xlsx(paths[["metadata"]], sheet = "sample") %>% filter(Sample %in% c("M", "Y", "P", "X", "F", "N"))
feature_metadata <- fread(paths[["feature_metadata"]]) %>% filter(Platform == "SOM", !Is_Protein_Group, !Is_Unknown) %>%
    select(AssayID, UniqueID) %>% distinct()
long_som <- fread(paths[["profiles"]]) %>% filter(Platform == "SOM", Batch %in% target_batches,
    Sample %in% c("M", "Y", "P", "X", "F", "N")) %>%
    semi_join(feature_metadata %>% select(UniqueID), by = "UniqueID")
detection_som <- fread(paths[["detection"]]) %>% filter(Platform == "SOM", Batch %in% target_batches) %>%
    group_by(Batch, UniqueID) %>% summarize(across(c(M, Y, P, X, F, N), ~ any(.x %in% TRUE)), .groups = "drop")
replicate_selection <- fread(paths[["replicate_selection"]]) %>% filter(Batch %in% target_batches)

select_som_stage <- function(data) {
    data %>% filter(Platform == "SOM", Batch %in% target_batches) %>%
        mutate(Stage = case_when(DataTier == "Calibrated" & ProcessLevel == "Calibrate" ~ "Calibrated",
                                 DataTier == "Reshaped" & str_detect(Batch, "^SOM_P1_") & ProcessLevel == "ANML-SMP" ~ "Reshaped",
                                 DataTier == "Reshaped" & str_detect(Batch, "^SOM_P2_") & ProcessLevel == "MedNormExt" ~ "Reshaped",
                                 TRUE ~ NA_character_)) %>% filter(!is.na(Stage))
}

# Reconstruct only the assay-level inputs needed by this Extended Data figure.
# Shared protein-level TRC and DEA outputs are deliberately left untouched.
stage_tasks <- long_som %>% select_som_stage() %>% distinct(Batch, Platform, ProcessLevel, DataTier, Stage)
trc_iterations <- vector("list", nrow(stage_tasks))
for (i in seq_len(nrow(stage_tasks))) {
    task <- stage_tasks[i, ]
    subset_data <- long_som %>% filter(Batch == task$Batch, ProcessLevel == task$ProcessLevel, DataTier == task$DataTier)
    sample_metadata <- meta_sample %>% filter(Batch == task$Batch, Sample %in% c("M", "Y", "P", "X", "F"))
    available_columns <- subset_data %>% filter(is.finite(Value)) %>% distinct(ColName) %>% pull(ColName)
    replicate_plan <- make_replicate_plan(sample_metadata, available_columns, n_replicates = 3, required_samples = c("M", "F"))
    iter <- map(seq_along(replicate_plan), function(j) {
        result <- calc_feature_titration_metrics(subset_data %>% filter(ColName %in% replicate_plan[[j]]$Metadata$ColName),
                                                  is_log2 = TRUE, method = "Mean", trc_deviation_cutoff = 0.25,
                                                  min_valid_replicates = 2)
        result %>% mutate(Batch = task$Batch, Platform = task$Platform, ProcessLevel = task$ProcessLevel,
                          DataTier = task$DataTier, Iteration = j)
    })
    trc_iterations[[i]] <- bind_rows(iter)
}
trc_feature_level <- summarize_titration_subsamples(bind_rows(trc_iterations),
    group_cols = c("Batch", "Platform", "ProcessLevel", "DataTier"), majority_cutoff = 0.5,
    trc_deviation_cutoff = 0.25)$feature

dea_results <- list()
for (i in seq_len(nrow(stage_tasks))) {
    task <- stage_tasks[i, ]
    meta_task <- replicate_selection %>% filter(Batch == task$Batch)
    expr <- long_som %>% filter(Batch == task$Batch, ProcessLevel == task$ProcessLevel, DataTier == task$DataTier,
                                ColName %in% meta_task$ColName) %>%
        select(UniqueID, ColName, Value) %>% pivot_wider(names_from = ColName, values_from = Value) %>% column_to_rownames("UniqueID") %>% as.matrix()
    meta_task <- meta_sample %>% semi_join(meta_task, by = c("Batch", "Sample", "ColName")) %>% filter(ColName %in% colnames(expr))
    expr <- expr[, meta_task$ColName, drop = FALSE]
    expr <- expr[rowSums(!is.na(expr)) > 0, , drop = FALSE]
    set.seed(999L + i)
    expr <- as.matrix(impute_lod_noise(expr, meta_task, seed = 999L + i))
    for (contrast in c("M/F", "N/P")) {
        groups <- strsplit(contrast, "/", fixed = TRUE)[[1]]
        cols <- meta_task %>% filter(Sample %in% groups) %>% pull(ColName)
        valid_ids <- detection_som %>% filter(Batch == task$Batch, .data[[groups[1]]] | .data[[groups[2]]]) %>% pull(UniqueID)
        expr_contrast <- expr[rownames(expr) %in% valid_ids, cols, drop = FALSE]
        if (!nrow(expr_contrast)) next
        result <- dea_limma_flexible(expr_contrast, meta_task %>% filter(ColName %in% cols), contrast_pair = contrast)
        dea_results[[length(dea_results) + 1]] <- result %>% mutate(Batch = task$Batch, Platform = task$Platform,
            DataTier = task$DataTier, ProcessLevel = task$ProcessLevel)
    }
}
dea_df_multi <- bind_rows(dea_results)

# 2. Panel a: Expected-response transitions ----
normalize_seqid <- function(x) paste0("seq.", gsub("-", ".", sub("^seq\\.", "", as.character(x), ignore.case = TRUE)))
anno_11k <- SomaDataIO::read_annotations(paths[["annotation_11k"]]) %>% filter(Organism == "Human") %>%
    transmute(Annotation_panel = "11K", AssayID = normalize_seqid(SeqId),
              Dilution_Bin = case_when(Dilution == "20%" ~ "High conc", Dilution == "0.5%" ~ "Med conc", Dilution == "0.005%" ~ "Low conc", TRUE ~ NA_character_))
anno_9k <- read_excel(paths[["annotation_9k"]], range = "B7:I10464") %>% filter(Organism == "Human") %>%
    transmute(Annotation_panel = "9K", AssayID = normalize_seqid(SeqID), Dilution = as.numeric(Dilution),
              Dilution_Bin = case_when(Dilution == 0.2 ~ "High conc", Dilution == 0.005 ~ "Med conc", Dilution == 0.00005 ~ "Low conc", TRUE ~ NA_character_))
anno_combined <- bind_rows(anno_11k, anno_9k) %>% filter(!is.na(Dilution_Bin)) %>% inner_join(feature_metadata, by = "AssayID") %>%
    select(Annotation_panel, UniqueID, Dilution_Bin) %>% distinct()
trc_som <- trc_feature_level %>% select_som_stage() %>%
    mutate(Evaluable = TRC_N_finite >= 2 & MonoRelationN >= 3 & MonoRelationEvaluableN == MonoRelationN & !is.na(ExpectedResponseValid)) %>%
    select(Batch, UniqueID, Stage, Evaluable, ExpectedResponseValid) %>% distinct() %>%
    pivot_wider(names_from = Stage, values_from = c(Evaluable, ExpectedResponseValid), names_sep = "_") %>%
    filter(Evaluable_Calibrated %in% TRUE, Evaluable_Reshaped %in% TRUE) %>%
    mutate(Response_Transition = case_when(ExpectedResponseValid_Calibrated & ExpectedResponseValid_Reshaped ~ "Conserved",
                                           !ExpectedResponseValid_Calibrated & ExpectedResponseValid_Reshaped ~ "Recovered",
                                           ExpectedResponseValid_Calibrated & !ExpectedResponseValid_Reshaped ~ "Disrupted",
                                           !ExpectedResponseValid_Calibrated & !ExpectedResponseValid_Reshaped ~ "Unresolved"))

annotation_choice <- trc_som %>% distinct(Batch, UniqueID) %>% inner_join(anno_combined, by = "UniqueID") %>%
    count(Batch, Annotation_panel, name = "Matched_measurements") %>% arrange(Batch, desc(Matched_measurements), Annotation_panel) %>%
    group_by(Batch) %>% slice(1) %>% ungroup()
batch_annotation <- annotation_choice %>% select(Batch, Annotation_panel) %>% inner_join(anno_combined, by = "Annotation_panel")
df_transition_feature <- trc_som %>% inner_join(batch_annotation, by = c("Batch", "UniqueID")) %>%
    mutate(Batch = factor(Batch, levels = target_batches),
           Dilution_Bin = factor(Dilution_Bin, levels = c("High conc", "Med conc", "Low conc")),
           Response_Transition = factor(Response_Transition, levels = c("Conserved", "Recovered", "Disrupted", "Unresolved")))

df_transition_summary <- df_transition_feature %>% filter(Response_Transition != "Unresolved") %>%
    count(Batch, Dilution_Bin, Response_Transition, name = "Count") %>% group_by(Batch, Dilution_Bin) %>%
    complete(Response_Transition = factor(c("Conserved", "Recovered", "Disrupted"), levels = levels(df_transition_feature$Response_Transition)), fill = list(Count = 0)) %>%
    mutate(Total = sum(Count), Percentage = 100 * Count / Total, Label = ifelse(Count > 0, as.character(Count), "")) %>% ungroup()

transition_colors <- c("Conserved" = "#00A087", "Recovered" = "#4DBBD5", "Disrupted" = "#DC0000")
p_transition <- ggplot(df_transition_summary, aes(x = Batch, y = Count, fill = Response_Transition)) +
    geom_col(position = "fill", width = 0.7, color = "black", linewidth = 0.1) +
    geom_text(aes(label = Label), position = position_fill(vjust = 0.5), size = 2.5, fontface = "bold", color = "white") +
    facet_grid(Dilution_Bin ~ .) +
    scale_fill_manual(values = transition_colors, breaks = c("Conserved", "Recovered", "Disrupted")) +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    labs(x = NULL, y = "SOMAmer measurements (%)", fill = NULL) +
    plasmix_theme +
    theme(panel.grid.major.x = element_blank(), legend.position = "bottom", legend.margin = margin(t = -8),
          legend.title = element_text(size = 7.5, face = "bold", vjust = 0.5), legend.text = element_text(size = 7.5, vjust = 0.5),
          legend.box.margin = margin(5, 0, 5, 0),
          axis.text.x = element_text(angle = 30, hjust = 1), axis.text.y = element_text(vjust = c(0.1, 0.5, 0.5, 0.5, 0.9))) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# 3. Panel b: Top-k direction consistency ----
dea_som <- dea_df_multi %>% select_som_stage() %>% filter(Pair %in% c("M/F", "N/P"))

generate_cat_plot <- function(target_pair) {
    df_pair <- dea_som %>% filter(Pair == target_pair, is.finite(logFC))
    calc_batch_cat <- function(batch_name) {
        calibrated <- df_pair %>% filter(Batch == batch_name, Stage == "Calibrated") %>% select(UniqueID, logFC) %>% distinct(UniqueID, .keep_all = TRUE)
        reshaped <- df_pair %>% filter(Batch == batch_name, Stage == "Reshaped") %>% select(UniqueID, logFC) %>% distinct(UniqueID, .keep_all = TRUE)
        common <- intersect(calibrated$UniqueID, reshaped$UniqueID); max_k <- min(10000, length(common)); step <- 10
        if (max_k < step) return(NULL)
        calibrated <- calibrated %>% filter(UniqueID %in% common) %>% arrange(desc(abs(logFC)), UniqueID)
        reshaped <- reshaped %>% filter(UniqueID %in% common) %>% arrange(desc(abs(logFC)), UniqueID)
        direction_calibrated <- setNames(sign(calibrated$logFC), calibrated$UniqueID)
        direction_reshaped <- setNames(sign(reshaped$logFC), reshaped$UniqueID)
        map_dfr(seq(step, max_k, by = step), function(k) {
            overlap <- intersect(calibrated$UniqueID[1:k], reshaped$UniqueID[1:k])
            n_concordant <- if (length(overlap) > 0) sum(direction_calibrated[overlap] == direction_reshaped[overlap]) else 0
            tibble(TopN = k, Consistency = 100 * n_concordant / k)
        }) %>% mutate(Batch = batch_name, Pair = target_pair)
    }
    cat_data <- bind_rows(lapply(target_batches, calc_batch_cat))
    axis_breaks <- c(10, 100, 1000, 10000); axis_breaks <- axis_breaks[axis_breaks <= max(cat_data$TopN, na.rm = TRUE)]
    cat_plot <- ggplot(cat_data, aes(x = TopN, y = Consistency, color = Batch)) +
        geom_line(linewidth = 0.5, alpha = 0.8) +
        scale_color_manual(values = batch_colors) +
        scale_x_log10(breaks = axis_breaks, expand = c(0, 0)) +
        scale_y_continuous(limits = c(70, 100), breaks = seq(70, 100, 10), expand = c(0, 0)) +
        labs(x = expression(Top~italic(k)~SOMAmer~measurements), y = paste("Concordance of", target_pair, "(%)")) +
        plasmix_theme +
        theme(legend.margin = margin(-2, 5, 5, 5), legend.title = element_text(size = 7.5, face = "bold", vjust = 0.5),
              legend.text = element_text(size = 7.5, vjust = 0.5), legend.key.size = unit(0.6, "lines"),
              axis.text.x = element_text(hjust = 0.5), plot.title = element_text(face = "bold", size = 11, hjust = 0.5)) +
        guides(color = guide_legend(override.aes = list(linewidth = 1)))
    list(p = cat_plot, dat = cat_data)
}
cat_MF <- generate_cat_plot("M/F"); cat_NP <- generate_cat_plot("N/P")

# 4. Panel c: Adjusted P-value changes ----
df_pval_all <- dea_som %>% select(UniqueID, Batch, Pair, Stage, adj.P.Val) %>%
    pivot_wider(names_from = Stage, values_from = adj.P.Val) %>% filter(!is.na(Calibrated), !is.na(Reshaped)) %>%
    mutate(P_Pri = Calibrated, P_Out = Reshaped,
           Status = case_when(P_Pri < 0.05 & P_Out < 0.05 ~ "Consistent", P_Pri < 0.05 & P_Out >= 0.05 ~ "Squished",
                              P_Pri >= 0.05 & P_Out < 0.05 ~ "Inflated", TRUE ~ "Non-Sig"),
           Facet_Label = factor(paste0(Pair, " (", Batch, ")"), levels = c(paste0("M/F (", target_batches, ")"), paste0("N/P (", target_batches, ")"))))

anno_pval <- df_pval_all %>% group_by(Facet_Label) %>%
    summarize(total_n = n(), n_consist = sum(Status == "Consistent"), n_inflat = sum(Status == "Inflated"),
              n_squish = sum(Status == "Squished"), n_nosig = sum(Status == "Non-Sig"), .groups = "drop") %>%
    mutate(anno = "Pval", lab_consist = sprintf("Consistent\n%d (%.1f%%)", n_consist, n_consist / total_n * 100),
           lab_inflat = sprintf("Inflated\n%d (%.1f%%)", n_inflat, n_inflat / total_n * 100),
           lab_squish = sprintf("Squished\n%d (%.1f%%)", n_squish, n_squish / total_n * 100),
           lab_nosig = sprintf("Non-significant\n%d (%.1f%%)", n_nosig, n_nosig / total_n * 100))

pval_facet <- ggplot(df_pval_all, aes(x = -log10(P_Pri), y = -log10(P_Out))) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black", linewidth = 0.5) +
    geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "black", linewidth = 0.5) +
    geom_point(aes(color = Status), alpha = 0.66, size = 0.33, stroke = 0.05) +
    geom_text(data = anno_pval, aes(x = 14.5, y = 14.5, label = lab_consist), hjust = 1, vjust = 1, color = "#017765", fontface = "bold", size = 2.5, lineheight = 0.9) +
    geom_text(data = anno_pval, aes(x = 0.5, y = 14.5, label = lab_inflat), hjust = 0, vjust = 1, color = "#b90000", fontface = "bold", size = 2.5, lineheight = 0.9) +
    geom_text(data = anno_pval, aes(x = 14.5, y = 0.5, label = lab_squish), hjust = 1, vjust = 0, color = "#0072B2", fontface = "bold", size = 2.5, lineheight = 0.9) +
    geom_text(data = anno_pval, aes(x = 0.5, y = 0.5, label = lab_nosig), hjust = 0, vjust = 0, color = "grey20", fontface = "plain", size = 2.5, lineheight = 0.9) +
    facet_wrap(~Facet_Label, nrow = 2) +
    scale_color_manual(values = c("Consistent" = "#00A087", "Squished" = "#4DBBD5", "Inflated" = "#DC0000", "Non-Sig" = "grey85")) +
    scale_x_continuous(limits = c(0, 15), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 15), expand = c(0, 0)) +
    labs(x = expression("-log"[10]~italic(P)~"(Calibrated)"), y = expression("-log"[10]~italic(P)~"(Reshaped)")) +
    plasmix_theme +
    theme(legend.position = "none", panel.grid.major = element_blank(), axis.text.x = element_text(hjust = c(0.1, 0.5, 0.5, 0.9)))

# 5. Assemble and export ----
right_panel <- ggarrange(ggarrange(cat_MF$p, cat_NP$p, common.legend = TRUE, legend = "right"), pval_facet, nrow = 2, heights = c(1, 2),
                         labels = c("b", "c"), font.label = label_style, label.x = 0, label.y = 1, hjust = -0.2, vjust = 1.2)
fig_supp_final <- ggarrange(p_transition, right_panel, nrow = 1, widths = c(1, 2), labels = "a", font.label = label_style,
                            label.x = 0, label.y = 1, hjust = -0.2, vjust = 1.2)

ggsave("figures/ed3_som_reshaping.pdf", fig_supp_final, width = 10, height = 6)
ggsave("figures/ed3_som_reshaping.png", fig_supp_final, width = 10, height = 6, dpi = 600, bg = "white")

output_list <- list("a_feature" = df_transition_feature, "a_summary" = df_transition_summary,
                    "b" = bind_rows(cat_MF$dat, cat_NP$dat), "c_feature" = df_pval_all, "c_summary" = anno_pval,
                    "annotation_choice" = annotation_choice)
write.xlsx(output_list, "tables/SourceData_EDFigure3.xlsx", overwrite = TRUE, keepNA = TRUE, na.string = "NA")
