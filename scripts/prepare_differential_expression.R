# Prepare differential-expression and MAPD results

# 0. Setup ----
source("scripts/_project_setup.R")
use_packages(c("data.table", "tidyverse", "readxl", "limma"))
source("utils/imputation.R")
source("utils/differential_analysis.R")
stable_seed <- function(..., base_seed = 999L) {
    parts <- vapply(list(...), as.character, character(1))
    code <- utf8ToInt(enc2utf8(paste(parts, collapse = "|")))
    as.integer((base_seed + sum(as.double(code) * seq_along(code))) %% 2147483646 + 1)
}

median_finite <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) median(x) else NA_real_
}

# 1. Data loading and preprocessing ----
message(">>> Loading data...")
meta_batch <- read_xlsx("data/study_metadata.xlsx", sheet = "batch")
meta_sample <- read_xlsx("data/study_metadata.xlsx", sheet = "sample")

analysis_samples <- c("M", "Y", "P", "X", "F", "N")
meta_sample_analysis <- meta_sample %>% filter(Sample %in% analysis_samples)

feat_meta <- fread("data/feature_metadata.tsv.gz") %>% filter(!Is_Protein_Group, !Is_Unknown)
valid_features <- unique(feat_meta$UniqueID)

lod_status <- fread("results/detection_status.tsv.gz")
missing_detection_columns <- setdiff(analysis_samples, colnames(lod_status))
if (length(missing_detection_columns)) stop("Detection table is missing columns: ", paste(missing_detection_columns, collapse = ", "))

lod_join <- lod_status %>%
    select(Batch, UniqueID, all_of(analysis_samples)) %>%
    group_by(Batch, UniqueID) %>%
    summarize(across(all_of(analysis_samples), ~ any(.x %in% TRUE)), .groups = "drop")

long_df <- fread("data/protein_profiles_long.tsv.gz")
long_df_filter <- long_df %>%
    filter(UniqueID %in% valid_features, DataTier %in% c("Baseline", "Calibrated", "Reshaped"))

# 2. Stable technical-replicate selection ----
# Each Batch × Sample group receives an independent seed, so other sample groups cannot alter its selection.
replicate_selection <- meta_sample_analysis %>%
    arrange(Batch, Sample, ColName) %>%
    group_by(Batch, Sample) %>%
    filter(n() >= 3) %>%
    group_modify(~{
        set.seed(stable_seed(.y$Batch[[1]], .y$Sample[[1]], "DEA-replicates"))
        slice_sample(.x, n = 3, replace = FALSE)
    }) %>%
    ungroup()

selection_check <- replicate_selection %>% count(Batch, Sample)
if (any(selection_check$n != 3)) stop("Technical-replicate selection did not return exactly three replicates per eligible group.")

fwrite(replicate_selection %>% select(Batch, Sample, ColName),
       "results/dea_replicate_selection.tsv.gz", sep = "\t", na = "NA")

# 3. Multi-dimensional differential expression analysis ----
desired_contrasts <- c("M/F", "M/X", "M/P", "M/Y", "Y/F", "Y/X", "Y/P", "P/F", "P/X", "X/F",
                       "N/M", "N/Y", "N/P", "N/X", "N/F")
fc_th <- 1.2
p_th <- 0.05
limma_list <- list()

for (bat in sort(unique(replicate_selection$Batch))) {
    metaMat <- replicate_selection %>% filter(Batch == bat)
    levels_in_batch <- long_df_filter %>%
        filter(Batch == bat) %>%
        distinct(DataTier, ProcessLevel) %>%
        arrange(DataTier, ProcessLevel)

    for (level_index in seq_len(nrow(levels_in_batch))) {
        data_tier <- levels_in_batch$DataTier[level_index]
        proc_lvl <- levels_in_batch$ProcessLevel[level_index]
        cat(sprintf(">>> Processing: [%s] - Tier: [%s] - Level: [%s]\n", bat, data_tier, proc_lvl))

        expr_df <- long_df_filter %>%
            filter(Batch == bat, DataTier == data_tier, ProcessLevel == proc_lvl, ColName %in% metaMat$ColName) %>%
            select(UniqueID, ColName, Value) %>%
            pivot_wider(names_from = ColName, values_from = Value) %>%
            column_to_rownames("UniqueID")
        if (!nrow(expr_df)) next

        expr_bat_raw <- as.matrix(expr_df)
        valid_cols <- metaMat$ColName[metaMat$ColName %in% colnames(expr_bat_raw)]
        if (!length(valid_cols)) next

        expr_bat <- expr_bat_raw[, valid_cols, drop = FALSE]
        expr_bat <- expr_bat[rowSums(!is.na(expr_bat)) > 0, , drop = FALSE]
        if (!nrow(expr_bat)) next

        metaMat_available <- metaMat %>% filter(ColName %in% colnames(expr_bat))
        available_samples <- unique(metaMat_available$Sample)

        set.seed(stable_seed(bat, data_tier, proc_lvl, "LOD-imputation"))
        expr_bat_imputed <- impute_lod_noise(expr_bat, metaMat_available)

        for (contrast in desired_contrasts) {
            groups <- strsplit(contrast, "/", fixed = TRUE)[[1]]
            g1 <- groups[1]
            g2 <- groups[2]
            if (!all(groups %in% available_samples)) next

            cols_g1 <- metaMat_available %>% filter(Sample == g1) %>% pull(ColName)
            cols_g2 <- metaMat_available %>% filter(Sample == g2) %>% pull(ColName)
            if (length(cols_g1) != 3 || length(cols_g2) != 3) next

            valid_prots_local <- lod_join %>%
                filter(Batch == bat, .data[[g1]] %in% TRUE | .data[[g2]] %in% TRUE) %>%
                pull(UniqueID) %>%
                unique()

            cols_selected <- c(cols_g1, cols_g2)
            exprMat_contrast <- expr_bat_imputed[
                rownames(expr_bat_imputed) %in% valid_prots_local, cols_selected, drop = FALSE
            ]
            if (!nrow(exprMat_contrast)) next

            metaMat_sub <- metaMat_available %>% filter(ColName %in% cols_selected)
            dea_res <- dea_limma_flexible(exprMat_contrast, metaMat_sub, contrast_pair = contrast)

            if (!is.null(dea_res) && nrow(dea_res)) {
                dea_res <- dea_res %>% mutate(Batch = bat, DataTier = data_tier, ProcessLevel = proc_lvl)
                limma_list[[length(limma_list) + 1]] <- dea_res
            }
        }
    }
}

df_fc_pvalue <- bind_rows(limma_list) %>%
    select(-any_of(c("Platform", "Platform.x", "Platform.y"))) %>%
    left_join(meta_batch %>% select(Batch, Platform) %>% distinct(), by = "Batch")

duplicate_dea_keys <- df_fc_pvalue %>%
    count(UniqueID, Pair, Batch, DataTier, ProcessLevel) %>%
    filter(n > 1)
if (nrow(duplicate_dea_keys)) stop("Duplicated DEA keys were detected.")

fwrite(df_fc_pvalue, "results/dea_limma_results.tsv.gz", sep = "\t", na = "NA")
message(">>> Finished: results/dea_limma_results.tsv.gz")

# 4. Feature-level MAPD analysis ----
# MAPD uses only technical replicates of M, Y, P, X, F and N. BLK, CAL and QC samples are excluded.
task_combs <- long_df_filter %>%
    distinct(Batch, DataTier, ProcessLevel) %>%
    arrange(Batch, DataTier, ProcessLevel)

df_mad_protein <- lapply(seq_len(nrow(task_combs)), function(i) {
    bat <- task_combs$Batch[i]
    data_tier <- task_combs$DataTier[i]
    proc_lvl <- task_combs$ProcessLevel[i]

    df_wide <- long_df_filter %>%
        filter(Batch == bat, DataTier == data_tier, ProcessLevel == proc_lvl) %>%
        select(UniqueID, ColName, Value) %>%
        pivot_wider(names_from = ColName, values_from = Value)
    if (!nrow(df_wide)) return(NULL)

    expr_mat <- as.matrix(df_wide %>% select(-UniqueID))
    rownames(expr_mat) <- df_wide$UniqueID

    bm <- meta_sample_analysis %>%
        filter(Batch == bat, ColName %in% colnames(expr_mat))

    reps_groups <- bm %>%
        group_by(Sample) %>%
        summarize(cols = list(ColName), Replicates = n(), .groups = "drop") %>%
        filter(Replicates >= 2)
    if (!nrow(reps_groups)) return(NULL)

    all_diffs <- lapply(reps_groups$cols, function(col_names) {
        valid_cols <- intersect(col_names, colnames(expr_mat))
        if (length(valid_cols) < 2) return(NULL)
        combos <- combn(valid_cols, 2)
        apply(combos, 2, function(pair) abs(expr_mat[, pair[1]] - expr_mat[, pair[2]]))
    })

    all_diffs <- Filter(Negate(is.null), all_diffs)
    if (!length(all_diffs)) return(NULL)

    combined_diffs <- do.call(cbind, all_diffs)
    data.frame(
        UniqueID = rownames(expr_mat), Batch = bat, DataTier = data_tier, ProcessLevel = proc_lvl,
        MAPD = apply(combined_diffs, 1, median_finite), MAPD_SampleGroups = nrow(reps_groups),
        MAPD_PairwiseComparisons = ncol(combined_diffs), stringsAsFactors = FALSE
    )
}) %>%
    bind_rows()

duplicate_mapd_keys <- df_mad_protein %>%
    count(UniqueID, Batch, DataTier, ProcessLevel) %>%
    filter(n > 1)
if (nrow(duplicate_mapd_keys)) stop("Duplicated MAPD keys were detected.")

fwrite(df_mad_protein, "results/mapd_feature_level.tsv.gz", sep = "\t", na = "NA")
message(">>> Finished: results/mapd_feature_level.tsv.gz")

# 5. Merge and classify differential-expression results ----
dea_mapd_check <- df_fc_pvalue %>%
    left_join(df_mad_protein, by = c("UniqueID", "Batch", "DataTier", "ProcessLevel"))

dea_mapd_check %>%
    summarize(
        Total = n(),
        MAPD_not_evaluable = sum(!is.finite(MAPD)),
        Significant_MAPD_not_evaluable = sum(adj.P.Val < p_th & !is.finite(MAPD), na.rm = TRUE),
        Potential_false_verified = sum(adj.P.Val < p_th & abs(logFC) > log2(fc_th) & !is.finite(MAPD), na.rm = TRUE)
    ) %>%
    print()

dea_df_multi <- df_fc_pvalue %>%
    left_join(df_mad_protein, by = c("UniqueID", "Batch", "DataTier", "ProcessLevel")) %>%
    left_join(lod_join, by = c("UniqueID", "Batch")) %>%
    mutate(
        Classification = case_when(
            !is.finite(logFC) | !is.finite(adj.P.Val) | adj.P.Val >= p_th ~ "Non-significant",
            !is.finite(MAPD) ~ "MAPD-not-evaluable",
            abs(logFC) <= 3 * MAPD ~ "Precision-rejected",
            abs(logFC) <= log2(fc_th) ~ "Small-magnitude",
            TRUE ~ "Verified-DEP"
        )
    ) %>%
    relocate(Batch, Platform, DataTier, ProcessLevel, .after = Pair)

required_output_columns <- c("UniqueID", "logFC", "P.Value", "adj.P.Val", "Pair", "Batch", "Platform", "DataTier", "ProcessLevel", "MAPD", "Classification")
missing_output_columns <- setdiff(required_output_columns, colnames(dea_df_multi))
if (length(missing_output_columns)) stop("Final DEA output is missing columns: ", paste(missing_output_columns, collapse = ", "))

fwrite(dea_df_multi, "results/dea_df_multi.tsv.gz", sep = "\t", na = "NA")
message(">>> Finished: results/dea_df_multi.tsv.gz")
