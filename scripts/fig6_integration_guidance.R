# Figure 6 | Integration strategies and protein-specific analytical boundaries

# 0. Setup ----
source("scripts/_project_setup.R")
use_packages(c("data.table", "tidyverse", "readxl", "pbapply", "patchwork", "openxlsx", "showtext", "effsize"))
source("utils/figure_style.R")
source("utils/batch_correction.R")
source("utils/benchmark_metrics.R")
plasmix_theme <- if (is.function(theme_plasmix)) theme_plasmix() else theme_plasmix
showtext_auto(); showtext_opts(dpi = 600)
set.seed(2026)

# 1. Inputs ----
paths <- c(
    metadata = "data/study_metadata.xlsx",
    profiles = "data/protein_profiles_long.tsv.gz",
    feature_metadata = "data/feature_metadata.tsv.gz",
    mapd = "results/mapd_feature_level.tsv.gz",
    physchem = "data/physchem_matrix.tsv.gz",
    physchem_dictionary = "data/physchem_dictionary.tsv"
)

missing_inputs <- paths[!file.exists(paths)]
if (length(missing_inputs)) stop("Missing input files:\n", paste(missing_inputs, collapse = "\n"))

selected_batches <- c(
    "DIA_P1_B1", "DIA_P2_B1", "DIA_P3_B1", "DIA_P4_B1", "DIA_P5_B1", "DIA_P5_B2",
    "OLK_P2_B1", "OLK_P2_B2",
    "SOM_P1_B1", "SOM_P1_B2", "SOM_P2_B1", "SOM_P2_B2"
)

meta_batch <- read_xlsx(paths["metadata"], sheet = "batch") %>% as_tibble()
meta_sample <- read_xlsx(paths["metadata"], sheet = "sample") %>% as_tibble()
feature_meta <- fread(paths["feature_metadata"]) %>% as_tibble()
long_df <- fread(paths["profiles"]) %>% as_tibble()
mapd_df <- fread(paths["mapd"]) %>% as_tibble()
physchem_matrix <- fread(paths["physchem"]) %>% as_tibble()
physchem_dict <- fread(paths["physchem_dictionary"]) %>% as_tibble()
cat_colors <- c("Structure" = "#E64B35", "Surface" = "#4DBBD5", "Charge" = "#00A087", "Disorder" = "#F39B7F", "Secretory" = "#8491B4", "Abundance" = "#91D1C2")

valid_features <- feature_meta %>%
    filter(!Is_Protein_Group, !Is_Unknown) %>%
    pull(UniqueID) %>%
    unique()

meta_batch_ht <- meta_batch %>%
    filter(Batch %in% selected_batches, Platform %in% c("DIA", "OLK", "SOM")) %>%
    select(Batch, Platform, Protocol)

meta_sample_ht <- meta_sample %>%
    filter(Batch %in% selected_batches)

long_df <- long_df %>%
    filter(Batch %in% selected_batches, UniqueID %in% valid_features)

if (!"Platform" %in% colnames(long_df)) {
    long_df <- long_df %>% left_join(meta_batch_ht %>% select(Batch, Platform), by = "Batch")
}

mapd_df <- mapd_df %>%
    left_join(meta_batch_ht %>% select(Batch, Platform), by = "Batch")

# 2. Design definitions ----
TRC_DEVIATION_CUTOFF <- 0.25
SMAPE_CUTOFF <- 30

designs <- list(
    Balanced = list(
        b1_smp = c("M", "Y", "P", "X", "F", "N"),
        b2_smp = c("M", "Y", "P", "X", "F", "N"),
        eval_samples = c("M", "Y", "X", "F")
    ),
    Partial = list(
        b1_smp = c("M", "Y", "P", "X", "N"),
        b2_smp = c("Y", "P", "X", "F", "N"),
        eval_samples = c("Y", "X")
    ),
    Confounded = list(
        b1_smp = c("M", "Y", "P", "N"),
        b2_smp = c("X", "F", "P", "N"),
        eval_samples = character()
    )
)

df_bat_pairs <- as.data.frame(combn(meta_batch_ht$Batch, 2), stringsAsFactors = FALSE) %>%
    t() %>% as.data.frame() %>%
    setNames(c("Batch1", "Batch2")) %>%
    left_join(meta_batch_ht %>% select(Batch1 = Batch, Plat1 = Platform), by = "Batch1") %>%
    left_join(meta_batch_ht %>% select(Batch2 = Batch, Plat2 = Platform), by = "Batch2") %>%
    mutate(
        Detailed_Type = if_else(Plat1 == Plat2, paste0("Intra-", Plat1), paste0(pmin(Plat1, Plat2), "-", pmax(Plat1, Plat2))),
        Macro_Type = if_else(Plat1 == Plat2, "Intra-platform", "Cross-platform")
    )

# Strategy, comparison and boundary palettes are shared across panels b–f.
design_levels <- c("Balanced", "Partial", "Confounded")
detailed_type_levels <- c("Intra-DIA", "Intra-OLK", "Intra-SOM", "DIA-OLK", "DIA-SOM", "OLK-SOM")
strategy_levels <- c("Native", "RF-BECA", "RI-BECA", "SRR", "SRR+RF-BECA")
strategy_colors <- c("Native" = "#D4C8B8", "RF-BECA" = "#7E6148", "RI-BECA" = "#D9A05B",
                     "SRR" = "#CF4E9C", "SRR+RF-BECA" = "#6D2F7F")
comparison_colors <- c("RF-BECA - Native" = strategy_colors[["RF-BECA"]], "RI-BECA - Native" = strategy_colors[["RI-BECA"]],
                       "SRR - Native" = strategy_colors[["SRR"]], "SRR+RF-BECA - SRR" = strategy_colors[["SRR+RF-BECA"]])
boundary_colors <- c("Conventional success" = "#7E6148", "Reference-enabled rescue" = "#CF4E9C", "Unresolved" = "#D9D9D9")
category_order <- c("Structure", "Surface", "Charge", "Disorder", "Secretory", "Abundance")
category_colors <- c("Structure" = "#E64B35", "Surface" = "#4DBBD5", "Charge" = "#00A087",
                     "Disorder" = "#F39B7F", "Secretory" = "#8491B4", "Abundance" = "#91D1C2")

# 3. Helpers ----
select_branch_rows <- function(df, branch = c("baseline", "final")) {
    branch <- match.arg(branch)
    if (branch == "baseline") {
        df %>%
            filter(
                (Platform == "DIA" & DataTier == "Baseline" & ProcessLevel == "Intensity") |
                    (Platform %in% c("OLK", "SOM") & DataTier == "Baseline" & ProcessLevel == "HybNorm")
            )
    } else {
        df %>%
            filter(
                (Platform == "DIA" & DataTier == "Baseline" & ProcessLevel == "Intensity") |
                    (Platform == "OLK" & DataTier == "Calibrated" & str_detect(ProcessLevel, regex("NPX", ignore_case = TRUE))) |
                    (Platform == "SOM" & DataTier == "Calibrated" & ProcessLevel == "Calibrate")
            )
    }
}

safe_mean <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) mean(x) else NA_real_
}

calculate_smape_module <- function(df_long, batch1_name, batch2_name, eval_samples) {
    keys <- df_long %>% distinct(UniqueID)
    if (!length(eval_samples)) return(keys %>% mutate(sMAPE = NA_real_))

    df_long %>%
        filter(Sample %in% eval_samples) %>%
        group_by(UniqueID, Sample) %>%
        summarize(
            mean_b1 = safe_mean(Value[Batch == batch1_name]),
            mean_b2 = safe_mean(Value[Batch == batch2_name]),
            .groups = "drop"
        ) %>%
        mutate(
            denom = (abs(mean_b1) + abs(mean_b2)) / 2,
            sample_sMAPE = if_else(is.finite(denom) & denom > 0, abs(mean_b1 - mean_b2) / denom, NA_real_)
        ) %>%
        group_by(UniqueID) %>%
        summarize(
            sMAPE = if (sum(is.finite(sample_sMAPE)) == length(eval_samples)) mean(sample_sMAPE, na.rm = TRUE) * 100 else NA_real_,
            .groups = "drop"
        ) %>%
        right_join(keys, by = "UniqueID")
}

calculate_expected_response_module <- function(df_long, is_log2 = FALSE) {
    keys <- df_long %>% distinct(UniqueID)

    res <- calc_feature_titration_metrics(
        df_long %>% filter(Sample %in% c("M", "Y", "X", "F")),
        method = "Mean",
        is_log2 = is_log2,
        trc_deviation_cutoff = TRC_DEVIATION_CUTOFF,
        min_valid_replicates = 2
    )

    res %>%
        transmute(UniqueID, TitrationMonoValid, TRC_N_finite, MeanTRCDev, TRCDevValid, ExpectedResponseValid) %>%
        right_join(keys, by = "UniqueID")
}

prepare_branch_matrix <- function(batch1, batch2, samples1, samples2, branch = c("baseline", "final"), features_keep = NULL) {
    branch <- match.arg(branch)

    meta_1 <- meta_sample_ht %>%
        filter(Batch == batch1, Sample %in% samples1) %>%
        group_by(Sample) %>%
        slice_head(n = 3) %>%
        ungroup()

    meta_2 <- meta_sample_ht %>%
        filter(Batch == batch2, Sample %in% samples2) %>%
        group_by(Sample) %>%
        slice_head(n = 3) %>%
        ungroup()

    metadata <- bind_rows(meta_1, meta_2) %>%
        select(ColName, Sample, Batch)

    df_sub <- long_df %>%
        filter(Batch %in% c(batch1, batch2), ColName %in% metadata$ColName)

    df_sub <- select_branch_rows(df_sub, branch)

    if (!is.null(features_keep)) {
        df_sub <- df_sub %>% filter(UniqueID %in% features_keep)
    }

    if (!nrow(df_sub)) return(list(expr = NULL, metadata = metadata))

    df_wide <- df_sub %>%
        group_by(UniqueID, ColName) %>%
        summarize(Value = safe_mean(Value), .groups = "drop") %>%
        pivot_wider(names_from = ColName, values_from = Value)

    if (!nrow(df_wide)) return(list(expr = NULL, metadata = metadata))

    expr <- df_wide %>%
        column_to_rownames("UniqueID") %>%
        as.matrix()

    common_cols <- intersect(metadata$ColName, colnames(expr))
    metadata <- metadata %>% filter(ColName %in% common_cols)
    expr <- expr[, metadata$ColName, drop = FALSE]

    list(expr = expr, metadata = metadata)
}

map_method_to_strategy <- function(method) {
    case_when(
        method == "Native" ~ "Native",
        method %in% c("MAD", "Quantile", "ComBat", "RUV-III-C", "RUVg") ~ "RF-BECA",
        str_detect(method, "^RUVs \\(") ~ "RI-BECA",
        method == "SRR (P)" ~ "SRR(P)",
        method == "SRR (N)" ~ "SRR(N)",
        str_detect(method, "^SRR\\+") & str_detect(method, "\\(P\\)$") ~ "SRR(P)+RF-BECA",
        str_detect(method, "^SRR\\+") & str_detect(method, "\\(N\\)$") ~ "SRR(N)+RF-BECA",
        TRUE ~ "Other"
    )
}

rf_methods <- c("ComBat", "MAD", "Quantile", "RUV-III-C", "RUVg")

infer_method_reference <- function(method) {
    case_when(
        str_detect(method, "\\(P\\)$") ~ "P",
        str_detect(method, "\\(N\\)$") ~ "N",
        TRUE ~ NA_character_
    )
}

get_validation_plan <- function(method_name, design_name, study_samples) {
    if (design_name != "Confounded") {
        return(list(list(ValidationTrack = "Study", ValidationSamples = study_samples,
                              ValidationSample = paste(study_samples, collapse = "/"))))
    }

    if (method_name == "Native" || method_name %in% rf_methods) {
        return(list(
            list(ValidationTrack = "P-anchor", ValidationSamples = "N", ValidationSample = "N"),
            list(ValidationTrack = "N-anchor", ValidationSamples = "P", ValidationSample = "P")
        ))
    }

    method_reference <- infer_method_reference(method_name)
    if (identical(method_reference, "P")) {
        return(list(list(ValidationTrack = "P-anchor", ValidationSamples = "N", ValidationSample = "N")))
    }
    if (identical(method_reference, "N")) {
        return(list(list(ValidationTrack = "N-anchor", ValidationSamples = "P", ValidationSample = "P")))
    }

    stop("Could not assign a validation track to method: ", method_name)
}

collapse_paneld_strategy <- function(strategy) {
    case_when(
        strategy %in% c("SRR(P)", "SRR(N)") ~ "SRR",
        strategy %in% c("SRR(P)+RF-BECA", "SRR(N)+RF-BECA") ~ "SRR+RF-BECA",
        TRUE ~ strategy
    )
}

# 4. Native M/F signal table for MAPD-based eligibility ----
native_default_long <- select_branch_rows(long_df, "final")

mf_summary <- native_default_long %>%
    filter(Sample %in% c("M", "F")) %>%
    group_by(Batch, UniqueID, Sample) %>%
    summarize(
        NRep = sum(is.finite(Value)),
        MeanValue = safe_mean(Value),
        .groups = "drop"
    ) %>%
    pivot_wider(
        names_from = Sample,
        values_from = c(NRep, MeanValue),
        names_sep = "_"
    )

mapd_default <- select_branch_rows(mapd_df, "final") %>%
    select(Batch, UniqueID, MAPD)

native_signal_tbl <- mf_summary %>%
    left_join(mapd_default, by = c("Batch", "UniqueID")) %>%
    mutate(
        MF_Evaluable = coalesce(NRep_M, 0) >= 2 & coalesce(NRep_F, 0) >= 2 &
            is.finite(MeanValue_M) & is.finite(MeanValue_F) & is.finite(MAPD),
        MF_Signal = MF_Evaluable & abs(MeanValue_M - MeanValue_F) > 3 * MAPD
    ) %>%
    select(Batch, UniqueID, MAPD, MF_Evaluable, MF_Signal)

get_pair_eligibility <- function(batch1, batch2) {
    sig1 <- native_signal_tbl %>%
        filter(Batch == batch1) %>%
        select(UniqueID, MAPD_B1 = MAPD, MF_Evaluable_B1 = MF_Evaluable, MF_Signal_B1 = MF_Signal)

    sig2 <- native_signal_tbl %>%
        filter(Batch == batch2) %>%
        select(UniqueID, MAPD_B2 = MAPD, MF_Evaluable_B2 = MF_Evaluable, MF_Signal_B2 = MF_Signal)

    full_join(sig1, sig2, by = "UniqueID") %>%
        mutate(
            MF_Evaluable_B1 = replace_na(MF_Evaluable_B1, FALSE),
            MF_Evaluable_B2 = replace_na(MF_Evaluable_B2, FALSE),
            MF_Signal_B1 = replace_na(MF_Signal_B1, FALSE),
            MF_Signal_B2 = replace_na(MF_Signal_B2, FALSE),
            Eligible_Broad = MF_Evaluable_B1 & MF_Evaluable_B2,
            Eligible_Primary = Eligible_Broad & (MF_Signal_B1 | MF_Signal_B2)
        )
}

# 5. Method evaluation ----
evaluate_method_matrix <- function(expr_df, metadata, method_name, batch1, batch2, eval_samples,
                                   design_name, validation_track, validation_sample) {
    df_linear <- expr_df %>%
        rownames_to_column("UniqueID") %>%
        pivot_longer(-UniqueID, names_to = "ColName", values_to = "Corrected_log2") %>%
        left_join(metadata %>% select(ColName, Batch, Sample), by = "ColName") %>%
        mutate(Value = 2^Corrected_log2) %>%
        select(UniqueID, ColName, Batch, Sample, Value)

    smape_df <- calculate_smape_module(
        df_linear,
        batch1_name = batch1,
        batch2_name = batch2,
        eval_samples = eval_samples
    )

    response_df <- calculate_expected_response_module(df_linear, is_log2 = FALSE)

    full_join(smape_df, response_df, by = "UniqueID") %>%
        mutate(
            Method = method_name,
            Strategy = map_method_to_strategy(method_name),
            MethodReference = infer_method_reference(method_name),
            ValidationTrack = validation_track,
            ValidationSample = validation_sample,
            QuantAgreementValid = is.finite(sMAPE),
            OverallSuccess = coalesce(sMAPE <= SMAPE_CUTOFF, FALSE) & coalesce(ExpectedResponseValid, FALSE)
        )
}

process_one_task <- function(task_index) {
    i_pair <- task_grid$PairIndex[task_index]
    design_name <- task_grid$Design[task_index]
    design_info <- designs[[design_name]]

    pair_info <- df_bat_pairs[i_pair, ]
    batch1 <- pair_info$Batch1
    batch2 <- pair_info$Batch2

    pair_eligibility <- get_pair_eligibility(batch1, batch2)
    eligible_primary <- pair_eligibility %>% filter(Eligible_Primary) %>% pull(UniqueID)

    baseline_raw <- prepare_branch_matrix(batch1, batch2, design_info$b1_smp, design_info$b2_smp, branch = "baseline", features_keep = eligible_primary)
    final_raw <- prepare_branch_matrix(batch1, batch2, design_info$b1_smp, design_info$b2_smp, branch = "final", features_keep = eligible_primary)

    n_primary <- length(eligible_primary)
    n_complete_baseline <- if (is.null(baseline_raw$expr)) 0 else sum(complete.cases(baseline_raw$expr))
    n_complete_final <- if (is.null(final_raw$expr)) 0 else sum(complete.cases(final_raw$expr))

    features_analysis <- eligible_primary
    if (!is.null(baseline_raw$expr)) {
        features_analysis <- intersect(features_analysis, rownames(baseline_raw$expr)[complete.cases(baseline_raw$expr)])
    } else {
        features_analysis <- character()
    }
    if (!is.null(final_raw$expr)) {
        features_analysis <- intersect(features_analysis, rownames(final_raw$expr)[complete.cases(final_raw$expr)])
    } else {
        features_analysis <- character()
    }

    denominator_row <- tibble(
        Batch1 = batch1,
        Batch2 = batch2,
        Detailed_Type = pair_info$Detailed_Type,
        Macro_Type = pair_info$Macro_Type,
        Design = design_name,
        N_Eligible_Primary = n_primary,
        N_Complete_Baseline = n_complete_baseline,
        N_Complete_Final = n_complete_final,
        N_Analysis = length(features_analysis)
    )

    if (!length(features_analysis)) {
        return(list(
            denominator = denominator_row,
            feature_results = tibble(),
            status = tibble()
        ))
    }

    baseline_data <- list(
        expr = baseline_raw$expr[features_analysis, , drop = FALSE],
        metadata = baseline_raw$metadata
    )
    final_data <- list(
        expr = final_raw$expr[features_analysis, , drop = FALSE],
        metadata = final_raw$metadata
    )

    baseline_run <- apply_beca(baseline_data$expr, baseline_data$metadata, srr_ref_types = c("P", "N"), run_mode = "baseline")
    final_run <- apply_beca(final_data$expr, final_data$metadata, srr_ref_types = c("P", "N"), run_mode = "final")

    expr_list <- c(
        list(Native = final_run$expr$Native),
        final_run$expr[names(final_run$expr) %in% c("MAD", "Quantile", "ComBat", "RUV-III-C", "RUVg", "RUVs (P)", "RUVs (N)")],
        baseline_run$expr[names(baseline_run$expr) %in% c(
            "SRR (P)", "SRR (N)",
            "SRR+MAD (P)", "SRR+Quantile (P)", "SRR+ComBat (P)", "SRR+RUV-III-C (P)", "SRR+RUVg (P)",
            "SRR+MAD (N)", "SRR+Quantile (N)", "SRR+ComBat (N)", "SRR+RUV-III-C (N)", "SRR+RUVg (N)"
        )]
    )

    method_res <- bind_rows(lapply(names(expr_list), function(method_name) {
        expr_df <- expr_list[[method_name]]
        if (is.null(expr_df) || !nrow(expr_df)) return(NULL)

        method_metadata <- if (method_name == "Native" || method_name %in% c(rf_methods, "RUVs (P)", "RUVs (N)")) {
            final_data$metadata
        } else {
            baseline_data$metadata
        }

        validation_plan <- get_validation_plan(method_name, design_name, design_info$eval_samples)
        bind_rows(lapply(validation_plan, function(plan) {
            evaluate_method_matrix(
                expr_df = expr_df,
                metadata = method_metadata,
                method_name = method_name,
                batch1 = batch1,
                batch2 = batch2,
                eval_samples = plan$ValidationSamples,
                design_name = design_name,
                validation_track = plan$ValidationTrack,
                validation_sample = plan$ValidationSample
            )
        }))
    })) %>%
        mutate(
            Batch1 = batch1,
            Batch2 = batch2,
            Detailed_Type = pair_info$Detailed_Type,
            Macro_Type = pair_info$Macro_Type,
            Design = design_name
        )

    status_res <- bind_rows(
        final_run$status %>% mutate(Branch = "final"),
        baseline_run$status %>% mutate(Branch = "baseline")
    ) %>%
        mutate(
            Batch1 = batch1,
            Batch2 = batch2,
            Detailed_Type = pair_info$Detailed_Type,
            Macro_Type = pair_info$Macro_Type,
            Design = design_name
        )

    list(denominator = denominator_row, feature_results = method_res, status = status_res)
}

# 6. Run the complete analysis or load the frozen analytical cache ----
ANALYSIS_CACHE <- "cache/fig6_analysis_cache_v1_bridge_validation.rds"
REBUILD_ANALYSIS <- FALSE

if (!REBUILD_ANALYSIS && file.exists(ANALYSIS_CACHE)) {
    list2env(readRDS(ANALYSIS_CACHE), envir = .GlobalEnv)
    message("Loaded Figure 6 analytical cache: ", ANALYSIS_CACHE)
} else {
    # 6.1 Run and derive the reusable analytical mother tables ----

    task_grid <- expand.grid(
        PairIndex = seq_len(nrow(df_bat_pairs)),
        Design = names(designs),
        stringsAsFactors = FALSE
    ) %>% as_tibble()

    message("Running the complete Figure 6 integration benchmark ...")
    task_results <- pblapply(seq_len(nrow(task_grid)), process_one_task)

    denominator_audit <- bind_rows(lapply(task_results, `[[`, "denominator"))
    method_feature_results <- bind_rows(lapply(task_results, `[[`, "feature_results"))
    algorithm_status <- bind_rows(lapply(task_results, `[[`, "status"))

    # 6.2 Pair-level summaries and integrity checks ----
    method_track_pair_summary <- method_feature_results %>%
        group_by(Detailed_Type, Macro_Type, Batch1, Batch2, Design, ValidationTrack, ValidationSample,
                 MethodReference, Strategy, Method) %>%
        summarize(
            N_Features = n(),
            Rate_sMAPE = mean(coalesce(sMAPE <= SMAPE_CUTOFF, FALSE)),
            Rate_Response = mean(coalesce(ExpectedResponseValid, FALSE)),
            Rate_Harmonized = mean(coalesce(OverallSuccess, FALSE)),
            Median_sMAPE = median(sMAPE, na.rm = TRUE),
            Median_TRCDev = median(MeanTRCDev, na.rm = TRUE),
            .groups = "drop"
        )

    method_pair_summary <- method_track_pair_summary %>%
        group_by(Detailed_Type, Macro_Type, Batch1, Batch2, Design, MethodReference, Strategy, Method) %>%
        summarize(
            N_Features = max(N_Features),
            N_ValidationTracks = n(),
            Rate_sMAPE = mean(Rate_sMAPE),
            Rate_Response = mean(Rate_Response),
            Rate_Harmonized = mean(Rate_Harmonized),
            Median_sMAPE = mean(Median_sMAPE, na.rm = TRUE),
            Median_TRCDev = mean(Median_TRCDev, na.rm = TRUE),
            .groups = "drop"
        )

    validation_track_audit <- method_feature_results %>%
        distinct(Batch1, Batch2, Detailed_Type, Macro_Type, Design, Method, Strategy, MethodReference,
                 ValidationTrack, ValidationSample) %>%
        arrange(Design, Batch1, Batch2, Method, ValidationTrack)

    expected_track_n <- validation_track_audit %>%
        mutate(ExpectedTracks = if_else(Design == "Confounded" & (Method == "Native" | Method %in% rf_methods), 2L, 1L)) %>%
        count(Batch1, Batch2, Design, Method, ExpectedTracks, name = "ObservedTracks")

    if (nrow(df_bat_pairs) != 66L || nrow(denominator_audit) != 198L) stop("Incomplete batch-pair or design enumeration.")
    if (any(expected_track_n$ObservedTracks != expected_track_n$ExpectedTracks)) stop("Unexpected validation-track count.")
    if (any(method_pair_summary %>% count(Batch1, Batch2, Design) %>% pull(n) != 20L)) stop("A task does not contain all 20 methods.")

    # 6.3 Matched algorithm-level incremental effects ----
    result_keys <- c("UniqueID", "Batch1", "Batch2", "Detailed_Type", "Macro_Type", "Design")
    track_keys <- c(result_keys, "ValidationTrack", "ValidationSample")

    comparison_map <- bind_rows(
        tibble(TargetMethod = rf_methods, ComparatorMethod = "Native", Comparison = "RF-BECA - Native",
               Algorithm = rf_methods, Reference = NA_character_),
        tibble(TargetMethod = c("RUVs (P)", "RUVs (N)"), ComparatorMethod = "Native", Comparison = "RI-BECA - Native",
               Algorithm = "RUVs", Reference = c("P", "N")),
        tibble(TargetMethod = c("SRR (P)", "SRR (N)"), ComparatorMethod = "Native", Comparison = "SRR - Native",
               Algorithm = "SRR", Reference = c("P", "N")),
        expand_grid(Reference = c("P", "N"), Algorithm = rf_methods) %>%
            mutate(TargetMethod = paste0("SRR+", Algorithm, " (", Reference, ")"),
                   ComparatorMethod = paste0("SRR (", Reference, ")"), Comparison = "SRR+RF-BECA - SRR") %>%
            select(TargetMethod, ComparatorMethod, Comparison, Algorithm, Reference)
    )

    target_results <- method_feature_results %>%
        inner_join(comparison_map, by = c("Method" = "TargetMethod")) %>%
        transmute(across(all_of(track_keys)), Comparison, Algorithm, Reference, TargetMethod = Method, ComparatorMethod,
                  TargetSuccess = coalesce(OverallSuccess, FALSE), Target_sMAPE = sMAPE,
                  TargetExpectedResponse = coalesce(ExpectedResponseValid, FALSE))

    comparator_results <- method_feature_results %>%
        transmute(across(all_of(track_keys)), ComparatorMethod = Method,
                  ComparatorSuccess = coalesce(OverallSuccess, FALSE), Comparator_sMAPE = sMAPE,
                  ComparatorExpectedResponse = coalesce(ExpectedResponseValid, FALSE))

    matched_incremental_feature_results <- target_results %>%
        inner_join(comparator_results, by = c(track_keys, "ComparatorMethod")) %>%
        mutate(
            DeltaSuccess = as.integer(TargetSuccess) - as.integer(ComparatorSuccess),
            Transition = case_when(
                !ComparatorSuccess & TargetSuccess ~ "Gain",
                ComparatorSuccess & !TargetSuccess ~ "Loss",
                ComparatorSuccess & TargetSuccess ~ "Shared success",
                TRUE ~ "Shared failure"
            )
        )

    matched_incremental_track_pair_summary <- matched_incremental_feature_results %>%
        group_by(Detailed_Type, Macro_Type, Batch1, Batch2, Design, ValidationTrack, ValidationSample,
                 Comparison, Algorithm, Reference, TargetMethod, ComparatorMethod) %>%
        summarize(
            N_Features = n(),
            Comparator_Rate = mean(ComparatorSuccess),
            Target_Rate = mean(TargetSuccess),
            Delta_Rate = Target_Rate - Comparator_Rate,
            Gain_Rate = mean(Transition == "Gain"),
            Loss_Rate = mean(Transition == "Loss"),
            Shared_Success_Rate = mean(Transition == "Shared success"),
            Shared_Failure_Rate = mean(Transition == "Shared failure"),
            .groups = "drop"
        )

    matched_incremental_pair_summary <- matched_incremental_track_pair_summary %>%
        group_by(Detailed_Type, Macro_Type, Batch1, Batch2, Design, Comparison, Algorithm, Reference,
                 TargetMethod, ComparatorMethod) %>%
        summarize(
            N_Features = max(N_Features),
            N_ValidationTracks = n(),
            Comparator_Rate = mean(Comparator_Rate),
            Target_Rate = mean(Target_Rate),
            Delta_Rate = mean(Delta_Rate),
            Gain_Rate = mean(Gain_Rate),
            Loss_Rate = mean(Loss_Rate),
            Shared_Success_Rate = mean(Shared_Success_Rate),
            Shared_Failure_Rate = mean(Shared_Failure_Rate),
            .groups = "drop"
        )

    matched_incremental_overall_summary <- matched_incremental_pair_summary %>%
        group_by(Detailed_Type, Macro_Type, Design, Comparison, Algorithm, Reference, TargetMethod, ComparatorMethod) %>%
        summarize(
            N_Pairs = n(),
            Median_Comparator_Rate = median(Comparator_Rate),
            Median_Target_Rate = median(Target_Rate),
            Median_Delta_Rate = median(Delta_Rate),
            Q1_Delta_Rate = quantile(Delta_Rate, 0.25),
            Q3_Delta_Rate = quantile(Delta_Rate, 0.75),
            Median_Gain_Rate = median(Gain_Rate),
            Median_Loss_Rate = median(Loss_Rate),
            .groups = "drop"
        )

    # 6.4 Strategy-level coverage and Panel d composition ----
    strategy_levels <- c("Native", "RF-BECA", "RI-BECA", "SRR", "SRR+RF-BECA")

    strategy_success_feature <- method_feature_results %>%
        filter(Design == "Balanced") %>%
        mutate(
            StrategyCollapsed = case_when(
                Strategy == "Native" ~ "Native",
                Strategy == "RF-BECA" ~ "RF-BECA",
                Strategy == "RI-BECA" ~ "RI-BECA",
                Strategy %in% c("SRR(P)", "SRR(N)") ~ "SRR",
                Strategy %in% c("SRR(P)+RF-BECA", "SRR(N)+RF-BECA") ~ "SRR+RF-BECA",
                TRUE ~ NA_character_
            )
        ) %>%
        filter(!is.na(StrategyCollapsed)) %>%
        group_by(across(all_of(result_keys)), StrategyCollapsed) %>%
        summarize(StrategySuccess = any(coalesce(OverallSuccess, FALSE)), .groups = "drop") %>%
        mutate(StrategyCollapsed = factor(StrategyCollapsed, levels = strategy_levels)) %>%
        pivot_wider(names_from = StrategyCollapsed, values_from = StrategySuccess, values_fill = FALSE, names_expand = TRUE)

    panel_d_feature_classes <- strategy_success_feature %>%
        mutate(
            N_Success = as.integer(Native) + as.integer(`RF-BECA`) + as.integer(`RI-BECA`) +
                as.integer(SRR) + as.integer(`SRR+RF-BECA`),
            IntegrationOutcome = case_when(
                N_Success == 0 ~ "Unresolved",
                N_Success > 1 ~ "Shared success",
                Native ~ "Native",
                `RF-BECA` ~ "RF-BECA",
                `RI-BECA` ~ "RI-BECA",
                SRR ~ "SRR",
                `SRR+RF-BECA` ~ "SRR+RF-BECA"
            )
        )

    integration_outcome_levels <- c("Native", "RF-BECA", "RI-BECA", "SRR", "SRR+RF-BECA", "Shared success", "Unresolved")

    panel_d_pair_composition <- panel_d_feature_classes %>%
        count(Detailed_Type, Macro_Type, Batch1, Batch2, Design, IntegrationOutcome, name = "N_Features") %>%
        complete(nesting(Detailed_Type, Macro_Type, Batch1, Batch2, Design),
                 IntegrationOutcome = integration_outcome_levels, fill = list(N_Features = 0L)) %>%
        group_by(Detailed_Type, Macro_Type, Batch1, Batch2, Design) %>%
        mutate(Total_Features = sum(N_Features), Proportion = N_Features / Total_Features) %>%
        ungroup() %>%
        mutate(IntegrationOutcome = factor(IntegrationOutcome, levels = integration_outcome_levels))

    panel_d_summary <- panel_d_pair_composition %>%
        group_by(Detailed_Type, Macro_Type, Design, IntegrationOutcome) %>%
        summarize(N_Pairs = n(), Mean_Proportion = mean(Proportion), Median_Proportion = median(Proportion),
                  Q1_Proportion = quantile(Proportion, 0.25), Q3_Proportion = quantile(Proportion, 0.75), .groups = "drop")

    if (any(abs(panel_d_summary %>% group_by(Detailed_Type, Design) %>%
                summarize(Total = sum(Mean_Proportion), .groups = "drop") %>% pull(Total) - 1) > 1e-10)) {
        stop("Panel d mean proportions do not sum to 1.")
    }

    # 6.5 Direct integration-boundary classification ----
    panel_d_boundary_feature_classes <- strategy_success_feature %>%
        mutate(
            ConventionalSuccess = Native | `RF-BECA`,
            ReferenceEnabledSuccess = `RI-BECA` | SRR | `SRR+RF-BECA`,
            BoundaryOutcome = case_when(
                ConventionalSuccess ~ "Conventional success",
                !ConventionalSuccess & ReferenceEnabledSuccess ~ "Reference-enabled rescue",
                TRUE ~ "Unresolved"
            )
        )

    boundary_outcome_levels <- c("Conventional success", "Reference-enabled rescue", "Unresolved")

    panel_d_boundary_pair_composition <- panel_d_boundary_feature_classes %>%
        count(Detailed_Type, Macro_Type, Batch1, Batch2, Design, BoundaryOutcome, name = "N_Features") %>%
        complete(nesting(Detailed_Type, Macro_Type, Batch1, Batch2, Design),
                 BoundaryOutcome = boundary_outcome_levels, fill = list(N_Features = 0L)) %>%
        group_by(Detailed_Type, Macro_Type, Batch1, Batch2, Design) %>%
        mutate(Total_Features = sum(N_Features), Proportion = N_Features / Total_Features) %>%
        ungroup() %>%
        mutate(BoundaryOutcome = factor(BoundaryOutcome, levels = boundary_outcome_levels))

    panel_d_boundary_summary <- panel_d_boundary_pair_composition %>%
        group_by(Detailed_Type, Macro_Type, Design, BoundaryOutcome) %>%
        summarize(N_Pairs = n(), Mean_Proportion = mean(Proportion), Median_Proportion = median(Proportion),
                  Q1_Proportion = quantile(Proportion, 0.25), Q3_Proportion = quantile(Proportion, 0.75), .groups = "drop")

    if (any(abs(panel_d_boundary_summary %>% group_by(Detailed_Type, Design) %>%
                summarize(Total = sum(Mean_Proportion), .groups = "drop") %>% pull(Total) - 1) > 1e-10)) {
        stop("Boundary mean proportions do not sum to 1.")
    }

    cache_objects <- list(
        task_grid = task_grid, denominator_audit = denominator_audit, algorithm_status = algorithm_status,
        validation_track_audit = validation_track_audit, method_feature_results = method_feature_results,
        method_track_pair_summary = method_track_pair_summary, method_pair_summary = method_pair_summary,
        matched_incremental_feature_results = matched_incremental_feature_results,
        matched_incremental_track_pair_summary = matched_incremental_track_pair_summary,
        matched_incremental_pair_summary = matched_incremental_pair_summary,
        matched_incremental_overall_summary = matched_incremental_overall_summary,
        panel_d_feature_classes = panel_d_feature_classes, panel_d_pair_composition = panel_d_pair_composition,
        panel_d_summary = panel_d_summary, panel_d_boundary_feature_classes = panel_d_boundary_feature_classes,
        panel_d_boundary_pair_composition = panel_d_boundary_pair_composition,
        panel_d_boundary_summary = panel_d_boundary_summary, boundary_outcome_levels = boundary_outcome_levels
    )
    saveRDS(cache_objects, ANALYSIS_CACHE)
    message("Saved Figure 6 analytical cache: ", ANALYSIS_CACHE)
}

# Figure 6 panel b: harmonization performance across overlap designs ----
# Run after the analytical section has created method_pair_summary.
unified_res_annotated <- method_pair_summary %>%
    mutate(
        Design = factor(Design, levels = c("Balanced", "Partial", "Confounded")),
        Strategy = case_when(
            Method == "Native" ~ "Native",
            Method %in% c("ComBat", "MAD", "Quantile", "RUV-III-C", "RUVg") ~ "RF-BECA",
            Method %in% c("RUVs (P)", "RUVs (N)") ~ "RI-BECA",
            Method == "SRR (P)" ~ "SRR(P)",
            Method == "SRR (N)" ~ "SRR(N)",
            str_detect(Method, "SRR\\+.*\\(P\\)") ~ "SRR(P)+RF-BECA",
            str_detect(Method, "SRR\\+.*\\(N\\)") ~ "SRR(N)+RF-BECA"
        ),
        Strategy = factor(Strategy, levels = c("Native", "RF-BECA", "RI-BECA", "SRR(P)", "SRR(P)+RF-BECA", "SRR(N)", "SRR(N)+RF-BECA")),
        Algorithm = case_when(
            str_detect(Method, "ComBat") ~ "ComBat",
            str_detect(Method, "RUV-III-C") ~ "RUV-III-C",
            str_detect(Method, "RUVg") ~ "RUVg",
            str_detect(Method, "RUVs") ~ "RUVs",
            str_detect(Method, "MAD") ~ "MAD",
            str_detect(Method, "Quantile") ~ "Quantile",
            TRUE ~ "None"
        ),
        Algorithm = factor(Algorithm, levels = c("None", "ComBat", "RUV-III-C", "RUVg", "RUVs", "MAD", "Quantile"))
    )

plot_scatter_data <- unified_res_annotated %>%
    group_by(Detailed_Type, Design, Strategy, Algorithm, Method) %>%
    summarize(N_Pairs = n(), Rate_Harmonized = median(Rate_Harmonized, na.rm = TRUE), .groups = "drop") %>%
    mutate(Detailed_Type = factor(Detailed_Type, levels = c("Intra-DIA", "Intra-OLK", "Intra-SOM", "DIA-OLK", "DIA-SOM", "OLK-SOM")))

plot_median_data <- plot_scatter_data %>%
    group_by(Detailed_Type, Design, Strategy) %>%
    summarize(Median_Rate_Harmonized = median(Rate_Harmonized, na.rm = TRUE), .groups = "drop")

custom_strategy_colors <- c("Native" = "#D4C8B8FF", "RF-BECA" = "#7E6148FF", "RI-BECA" = "#D9A05B",
    "SRR(P)" = "#CF4E9CFF", "SRR(P)+RF-BECA" = "#6D2F7F", "SRR(N)" = "#00A087FF", "SRR(N)+RF-BECA" = "#226b6b")
custom_base_shapes <- c("None" = 16, "ComBat" = 15, "RUV-III-C" = 18, "RUVg" = 17, "RUVs" = 8, "MAD" = 3, "Quantile" = 4)
pd <- position_dodge(width = 0.5)
p_b <- ggplot() +
    geom_point(data = plot_scatter_data, aes(x = Design, y = Rate_Harmonized, color = Strategy, shape = Algorithm), position = pd, alpha = 0.3, size = 0.8) +
    geom_line(data = plot_median_data, aes(x = Design, y = Median_Rate_Harmonized, color = Strategy, group = Strategy), position = pd, linewidth = 0.5) +
    geom_point(data = plot_median_data, aes(x = Design, y = Median_Rate_Harmonized, color = Strategy), position = pd, shape = 21, fill = "white", size = 1.25, stroke = 1.25) +
    facet_wrap(~Detailed_Type, nrow = 1) +
    scale_x_discrete(labels = c("Balanced" = "Bal.", "Partial" = "Part.", "Confounded" = "Conf.")) +
    scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2), labels = c("0", "20", "40", "60", "80"), expand = c(0, 0)) +
    scale_color_manual(values = custom_strategy_colors) +
    scale_shape_manual(values = custom_base_shapes) +
    labs(x = NULL, y = "Harmonized-feature rate (%)") + plasmix_theme +
    theme(
        panel.grid.major.x = element_blank(),
        legend.position = "right", legend.box = "vertical", legend.box.just = "left",
        legend.key.size = unit(0.6, "lines"),
        legend.title = element_text(size = 7.5, face = "bold", vjust = 0.5, margin = margin(0, 0, 2, 0)),
        legend.text = element_text(size = 7.5, vjust = 0.5, margin = margin(0, 0, 0, 2)),
        legend.margin = margin(t = 5, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 0, b = 0, r = -6, l = 0)
    ) +
    guides(
        color = guide_legend(ncol = 1, order = 1, override.aes = list(size = 2, alpha = 1, shape = 21, fill = "white")),
        shape = guide_legend(ncol = 1, order = 2, byrow = TRUE, override.aes = list(size = 1.5, alpha = 1, color = "black"))
    )

# Panel c: exclusive integration outcomes across six settings and three designs ----
# Shared success should be at the bottom, Unresolved at the top, and labels must match stack order.
strategy_success_feature_all <- method_feature_results %>%
    mutate(
        StrategyCollapsed = case_when(
            Strategy == "Native" ~ "Native",
            Strategy == "RF-BECA" ~ "RF-BECA",
            Strategy == "RI-BECA" ~ "RI-BECA",
            Strategy %in% c("SRR(P)", "SRR(N)") ~ "SRR",
            Strategy %in% c("SRR(P)+RF-BECA", "SRR(N)+RF-BECA") ~ "SRR+RF-BECA",
            TRUE ~ NA_character_
        )
    ) %>%
    filter(!is.na(StrategyCollapsed)) %>%
    group_by(Detailed_Type, Batch1, Batch2, Design, UniqueID, StrategyCollapsed) %>%
    summarize(StrategySuccess = any(coalesce(OverallSuccess, FALSE)), .groups = "drop") %>%
    mutate(StrategyCollapsed = factor(StrategyCollapsed, levels = c("Native", "RF-BECA", "RI-BECA", "SRR", "SRR+RF-BECA"))) %>%
    pivot_wider(names_from = StrategyCollapsed, values_from = StrategySuccess, values_fill = FALSE, names_expand = TRUE)

stack_levels <- c("Shared success", "Native", "SRR", "SRR+RF-BECA", "RF-BECA", "RI-BECA", "Unresolved")

exclusive_stats_raw <- strategy_success_feature_all %>%
    mutate(
        Category = case_when(
            (SRR | `SRR+RF-BECA`) & (Native | `RF-BECA` | `RI-BECA`) ~ "Shared success",
            SRR ~ "SRR",
            `SRR+RF-BECA` ~ "SRR+RF-BECA",
            `RF-BECA` ~ "RF-BECA",
            `RI-BECA` ~ "RI-BECA",
            Native ~ "Native",
            TRUE ~ "Unresolved"
        ),
        Category = factor(Category, levels = stack_levels),
        Detailed_Type = factor(Detailed_Type, levels = detailed_type_levels),
        Design = factor(Design, levels = design_levels)
    )

plot_data_final <- exclusive_stats_raw %>%
    count(Detailed_Type, Design, Batch1, Batch2, Category, name = "Count") %>%
    complete(nesting(Detailed_Type, Design, Batch1, Batch2), Category, fill = list(Count = 0L)) %>%
    group_by(Detailed_Type, Design, Batch1, Batch2) %>%
    mutate(Prop = Count / sum(Count)) %>%
    ungroup() %>%
    group_by(Detailed_Type, Design, Category) %>%
    summarize(Mean_Prop = mean(Prop, na.rm = TRUE), .groups = "drop") %>%
    mutate(Category = factor(Category, levels = stack_levels))

label_data <- plot_data_final %>%
    group_by(Detailed_Type, Design) %>%
    arrange(Category, .by_group = TRUE) %>%
    mutate(ymin = lag(cumsum(Mean_Prop), default = 0), ymax = cumsum(Mean_Prop), ymid = (ymin + ymax) / 2) %>%
    ungroup() %>%
    filter(Category %in% c("Shared success", "SRR+RF-BECA", "Unresolved")) %>%
    mutate(
        Label = if_else(Mean_Prop >= 0.04, scales::percent(Mean_Prop, accuracy = 0.1), ""),
        TextColor = case_when(
            Category %in% c("Shared success", "SRR+RF-BECA") ~ "white",
            Category == "Unresolved" ~ "black",
            TRUE ~ "black"
        )
    )

category_colors <- c("Shared success" = "#2F5597", "Native" = "#D4C8B8FF", "SRR" = "#d6b8e6", "SRR+RF-BECA" = "#6D2F7F", "RF-BECA" = "#7E6148FF", "RI-BECA" = "#D9A05B", "Unresolved" = "#d4d4d4")
p_c <- ggplot(plot_data_final, aes(x = Design, y = Mean_Prop, fill = Category)) +
    geom_col(width = 0.75, color = "white", linewidth = 0.25, position = position_stack(reverse = TRUE)) +
    geom_text(data = label_data, aes(x = Design, y = ymid, label = Label, color = TextColor), inherit.aes = FALSE, size = 2.5, fontface = "bold", show.legend = FALSE) +
    facet_wrap(~Detailed_Type, nrow = 1) +
    scale_x_discrete(labels = c("Balanced" = "Bal.", "Partial" = "Part.", "Confounded" = "Conf.")) +
    scale_fill_manual(values = category_colors, breaks = c("Shared success", "Unresolved", "SRR", "SRR+RF-BECA", "RF-BECA", "RI-BECA", "Native")) +
    scale_color_manual(values = c("white" = "white", "black" = "black"), guide = "none") +
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = c("0", "25", "50", "75", "100"), expand = c(0, 0)) +
    labs(x = NULL, y = "Proportion of proteins (%)", fill = "Integration\noutcome") +
    plasmix_theme +
    theme(
        panel.grid.major.x = element_blank(),
        legend.position = "right", legend.box = "vertical", legend.box.just = "left",
        legend.key.size = unit(0.6, "lines"),
        legend.title = element_text(size = 7.5, face = "bold", vjust = 0.5, margin = margin(0, 0, 5, 0)),
        legend.text = element_text(size = 7.5, vjust = 0.5, margin = margin(0, 0, 0, 2)),
        legend.margin = margin(t = 5, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 0, b = 0, r = -5, l = 0),
        plot.margin = margin(t = 5, b = 5, r = 5, l = 5)
    ) +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE))

# 10. Panel d: consensus integration-success curve ----
hpa_sel <- physchem_matrix %>%
    transmute(
        UniProtID = Entry,
        BloodConc_log10_pgml = if ("BloodConc_log10_pgml" %in% colnames(physchem_matrix)) as.numeric(BloodConc_log10_pgml) else NA_real_,
        Abundance_Source = if ("Abundance_Source" %in% colnames(physchem_matrix)) as.character(Abundance_Source) else NA_character_
    ) %>%
    filter(!is.na(UniProtID)) %>%
    distinct(UniProtID, .keep_all = TRUE)

long_df_rank <- select_branch_rows(long_df, "final")
dict_rank_batch <- long_df_rank %>%
    filter(is.finite(Value)) %>%
    group_by(Batch, UniqueID) %>%
    summarize(Median_Val = median(Value, na.rm = TRUE), .groups = "drop") %>%
    group_by(Batch) %>%
    mutate(Batch_Rank_pct = percent_rank(Median_Val)) %>%
    ungroup() %>%
    select(Batch, UniqueID, Batch_Rank_pct)

entry_candidates <- c("UniProtID", "UniProt_ID", "UniProt", "Uniprot", "Entry", "ProteinID")
entry_column <- intersect(entry_candidates, colnames(feature_meta))[1]
feature_entry_map <- if (length(entry_column) && !is.na(entry_column)) feature_meta %>% transmute(UniqueID, UniProtID = as.character(.data[[entry_column]])) else feature_meta %>% transmute(UniqueID, UniProtID = as.character(UniqueID))
feature_entry_map <- feature_entry_map %>%
    mutate(UniProtID = str_trim(UniProtID), UniProtID = if_else(str_detect(UniProtID, "[;|, ]"), NA_character_, UniProtID)) %>%
    distinct(UniqueID, .keep_all = TRUE)

consensus_pair <- method_feature_results %>%
    filter(Design == "Balanced") %>%
    group_by(Detailed_Type, Batch1, Batch2, UniqueID, Method) %>%
    summarize(MethodSuccess = mean(coalesce(OverallSuccess, FALSE)), .groups = "drop") %>%
    group_by(Detailed_Type, Batch1, Batch2, UniqueID) %>%
    summarize(PairSuccessRate = mean(MethodSuccess), N_Methods = n(), .groups = "drop")
if (any(consensus_pair$N_Methods != 20)) stop("Panel d consensus does not contain exactly 20 methods for every feature–pair.")

consensus_voting <- consensus_pair %>%
    group_by(Detailed_Type, UniqueID) %>%
    summarize(Success_Rate = median(PairSuccessRate), .groups = "drop") %>%
    left_join(feature_entry_map, by = "UniqueID") %>%
    mutate(UniProtID = if_else(is.na(UniProtID) & as.character(UniqueID) %in% physchem_matrix$Entry, as.character(UniqueID), UniProtID)) %>%
    left_join(
        dict_rank_batch %>%
            group_by(UniqueID) %>% summarize(Native_Rank_pct = median(Batch_Rank_pct, na.rm = TRUE), .groups = "drop"),
        by = "UniqueID"
    ) %>%
    left_join(hpa_sel, by = "UniProtID") %>%
    mutate(
        Detailed_Type = factor(Detailed_Type, levels = detailed_type_levels),
        Rank_Source = case_when(
            is.finite(BloodConc_log10_pgml) ~ "HPA",
            is.finite(Native_Rank_pct) ~ "Native",
            TRUE ~ "Missing"
        ),
        Rank_Value = case_when(
            Rank_Source == "HPA" ~ percent_rank(BloodConc_log10_pgml),
            Rank_Source == "Native" ~ Native_Rank_pct,
            TRUE ~ NA_real_
        )
    ) %>%
    filter(!is.na(Rank_Value))

active_dt <- levels(droplevels(consensus_voting$Detailed_Type))
shade_and_text_data <- data.frame(Detailed_Type = factor(active_dt, levels = active_dt), Shade_Xmin = 0, Shade_Xmax = 0.25)

srr_status <- method_feature_results %>%
    filter(Design == "Balanced", Method == "SRR (P)") %>%
    group_by(Detailed_Type, UniqueID) %>%
    summarize(SRR_Rate = mean(coalesce(OverallSuccess, FALSE)), .groups = "drop") %>%
    mutate(SRR_Call = if_else(SRR_Rate > 0, "Pass", "Fail"), Detailed_Type = factor(Detailed_Type, levels = active_dt))

plot_combined_data <- consensus_voting %>%
    inner_join(srr_status, by = c("Detailed_Type", "UniqueID")) %>%
    group_by(Detailed_Type) %>%
    arrange(desc(Success_Rate), desc(SRR_Call), .by_group = TRUE) %>%
    mutate(Rank_Percentile = row_number() / n()) %>%
    ungroup() %>%
    mutate(Detailed_Type = factor(Detailed_Type, levels = active_dt))

cutoff_x <- 0.25
intersection_data <- plot_combined_data %>%
    group_by(Detailed_Type) %>%
    slice(which.min(abs(Rank_Percentile - cutoff_x))) %>%
    ungroup() %>%
    select(Detailed_Type, Success_Rate) %>%
    mutate(X_val = cutoff_x, Label_Text = scales::percent(Success_Rate, accuracy = 0.1))

p_cdf <- ggplot() +
    geom_rect(data = shade_and_text_data, aes(xmin = Shade_Xmin, xmax = Shade_Xmax, ymin = -Inf, ymax = Inf), fill = "#c0c0c0", alpha = 0.2, inherit.aes = FALSE) +
    geom_segment(data = intersection_data, aes(x = X_val, xend = X_val, y = 0, yend = Success_Rate), linetype = "dotted", color = "grey50", linewidth = 0.5) +
    geom_segment(data = intersection_data, aes(x = 0, xend = X_val, y = Success_Rate, yend = Success_Rate), linetype = "dotted", color = "grey50", linewidth = 0.5) +
    geom_rug(data = plot_combined_data, aes(x = pmin(pmax(Rank_Percentile, 0.015), 0.985), color = SRR_Call), sides = "b", length = unit(0.12, "npc"), alpha = 0.7, linewidth = 1) +
    geom_line(data = plot_combined_data, aes(x = Rank_Percentile, y = Success_Rate), linewidth = 0.5, color = "grey25") +
    geom_point(data = intersection_data, aes(x = X_val, y = Success_Rate), size = 1.5, color = "black") +
    geom_text(data = intersection_data, aes(x = X_val, y = Success_Rate, label = Label_Text), hjust = -0.15, vjust = -0.5, size = 2.5, fontface = "bold", color = "black") +
    facet_wrap(Detailed_Type ~ ., nrow = 2) +
    scale_x_continuous(breaks = seq(0, 1, 0.25), labels = c("0", "25", "50", "75", "100"), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.20), labels = c("0", "20", "40", "60", "80", "100"), expand = c(0, 0)) +
    scale_color_manual(values = c("Pass" = "#3171b8", "Fail" = "#DF6B6A")) +
    labs(x = "Ranked proteins (%)", y = "Consensus success rate (%)", color = "SRR diagnosis") + plasmix_theme +
    theme(
        panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = c(0.9, 0.5, 0.5, 0.5, 0.1)),
        legend.title = element_text(size = 7.5, face = "bold", vjust = 0.5), legend.text = element_text(size = 7.5, vjust = 0.5, margin = margin(0, 0, 0, 3)),
        legend.position = "top", legend.box.margin = margin(b = -10, t = -7)
    )
print(p_cdf)

# 11. Panel e: physicochemical differences between integration-success and failure groups ----
# 11.1 Prepare physicochemical comparison data ----
final_physchem <- physchem_dict %>% filter(Retained == "Yes") %>% pull(Feature)
consensus_voting_class <- consensus_voting %>%
    group_by(Detailed_Type) %>%
    mutate(
        Threshold_Q75 = quantile(Success_Rate, 0.75, na.rm = TRUE),
        Integration_Status = case_when(
            Success_Rate == 0 ~ "Failed",
            Success_Rate >= Threshold_Q75 ~ "Success",
            TRUE ~ "Intermediate"
        )
    ) %>%
    ungroup() %>%
    filter(Integration_Status %in% c("Success", "Failed"))

results_physchem <- map_dfr(unique(as.character(consensus_voting_class$Detailed_Type)), function(dt) {
    df_sub <- consensus_voting_class %>%
        filter(as.character(Detailed_Type) == dt) %>%
        inner_join(physchem_matrix, by = c("UniProtID" = "Entry"))
    map_dfr(final_physchem, function(feat) {
        v_success <- suppressWarnings(as.numeric(df_sub[[feat]][df_sub$Integration_Status == "Success"]))
        v_failed <- suppressWarnings(as.numeric(df_sub[[feat]][df_sub$Integration_Status == "Failed"]))
        v_success <- v_success[is.finite(v_success)]
        v_failed <- v_failed[is.finite(v_failed)]
        if (length(v_success) < 5 || length(v_failed) < 5) return(NULL)
        tibble(
            Detailed_Type = dt,
            Feature = feat,
            Cliff_Delta = as.numeric(effsize::cliff.delta(v_success, v_failed)$estimate),
            P_Value = suppressWarnings(wilcox.test(v_success, v_failed, exact = FALSE)$p.value)
        )
    })
})

df_physchem_summary <- results_physchem %>%
    left_join(
        physchem_dict %>%
            transmute(
                Feature,
                Property,
                Category = as.character(Category)
            ),
        by = "Feature"
    ) %>%
    mutate(
        Fill_Status = case_when(
            P_Value < 0.05 & Cliff_Delta > 0 ~ "Success-enriched",
            P_Value < 0.05 & Cliff_Delta < 0 ~ "Failure-enriched",
            TRUE ~ "Non-significant"
        ),
        Stars = case_when(
            P_Value < 0.001 ~ "***",
            P_Value < 0.01 ~ "**",
            P_Value < 0.05 ~ "*",
            TRUE ~ ""
        )
    )

sig_features_to_display <- df_physchem_summary %>%
    group_by(Feature) %>%
    summarize(
        Any_Significant = any(P_Value < 0.05, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    filter(Any_Significant) %>%
    pull(Feature)

df_physchem_display <- df_physchem_summary %>% filter(Feature %in% sig_features_to_display)

legend_category_order <- names(cat_colors)

# Preserve the original property order in physchem_dict
property_order_top_to_bottom <- physchem_dict %>%
    mutate(Dictionary_Order = row_number()) %>%
    transmute(Feature, Property, Category = as.character(Category), Dictionary_Order) %>%
    filter(Feature %in% unique(df_physchem_display$Feature),
           Property %in% unique(as.character(df_physchem_display$Property)),
           Category %in% legend_category_order) %>%
    distinct(Property, Category, .keep_all = TRUE) %>%
    mutate(Category = factor(Category, levels = legend_category_order)) %>%
    arrange(Category, Dictionary_Order) %>%
    pull(Property)

# The first factor level appears at the bottom of a ggplot discrete y-axis,
# so reverse the property levels to obtain the requested top-to-bottom order
df_physchem_display <- df_physchem_display %>%
    mutate(
        Category = factor(as.character(Category), levels = legend_category_order),
        Property = factor(as.character(Property), levels = rev(property_order_top_to_bottom)),
        Detailed_Type = factor(Detailed_Type, levels = active_dt)
    )

# Left physicochemical-category strip
df_cat <- df_physchem_display %>%
    distinct(Property, Category) %>%
    mutate(Ghost_Panel = " ")

p_cat <- ggplot(df_cat, aes(x = "1", y = Property, fill = Category)) +
    geom_tile(color = "white", linewidth = 1, width = 0.6) +
    facet_wrap(~Ghost_Panel) +
    scale_fill_manual(
        values = cat_colors,
        breaks = legend_category_order,
        limits = legend_category_order,
        name = "Category"
    ) +
    plasmix_theme +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color = "transparent", size = 8.5),
        plot.margin = margin(5, 0, 5, 5)
    ) +
    guides(
        fill = guide_legend(
            order = 2, ncol = 1,
            keywidth = unit(0.6, "lines"),
            keyheight = unit(0.6, "lines"),
            override.aes = list(color = "white", linewidth = 0.25)
        )
    )

# Main Cliff's-delta bar plot
cliff_colors <- c("Success-enriched" = "#3171b8", "Failure-enriched" = "#EECEB7", "Non-significant" = "#E0E0E0")

df_physchem_display$Fill_Status <- factor(df_physchem_display$Fill_Status, levels = names(cliff_colors))

star_data <- df_physchem_display %>%
    filter(P_Value < 0.05) %>%
    mutate(
        Required_Width = case_when(
            Stars == "*" ~ 0.055,
            Stars == "**" ~ 0.080,
            Stars == "***" ~ 0.105,
            TRUE ~ 0
        ),
        Star_Inside = Cliff_Delta > 0 & Cliff_Delta >= Required_Width,
        Star_x = case_when(
            Cliff_Delta > 0 & Star_Inside ~ 0.015,
            Cliff_Delta > 0 & !Star_Inside ~ Cliff_Delta + 0.015,
            TRUE ~ -0.015
        ),
        Star_hjust = case_when(
            Cliff_Delta > 0 ~ 0,
            TRUE ~ 1
        ),
        Star_color = case_when(
            Cliff_Delta > 0 & Star_Inside ~ "white",
            Cliff_Delta > 0 & !Star_Inside ~ "#3171b8",
            TRUE ~ "black"
        )
    )

axis_guard <- df_physchem_display %>%
    group_by(Detailed_Type) %>%
    summarize(Property = first(Property), .groups = "drop") %>%
    crossing(Cliff_Delta = c(-0.15, 0.15))

p_main <- ggplot(df_physchem_display, aes(x = Cliff_Delta, y = Property, fill = Fill_Status)) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.35) +
    geom_blank(data = axis_guard, aes(x = Cliff_Delta, y = Property), inherit.aes = FALSE) +
    geom_bar(stat = "identity", color = "black", width = 0.75, alpha = 0.9, linewidth = 0.1) +
    geom_text(
        data = star_data,
        aes(x = Star_x, y = Property, label = Stars, hjust = Star_hjust, color = Star_color),
        inherit.aes = FALSE,
        size = 3.5,
        vjust = 0.75
    ) +
    scale_color_identity() +
    facet_wrap(~Detailed_Type, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = cliff_colors, name = "Direction") +
    scale_x_continuous(breaks = c(-0.47, -0.33, -0.15, 0, 0.15, 0.33, 0.47), expand = c(0, 0)) +
    labs(x = expression("Effect size (Cliff's " * delta * ")"), y = NULL) +
    plasmix_theme +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_text(margin = margin(0, 0, 0, 0)),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(0, 5, 5, 0)
    ) +
    guides(
        fill = guide_legend(
            order = 1, ncol = 1,
            keywidth = unit(0.6, "lines"),
            keyheight = unit(0.6, "lines"),
            override.aes = list(color = "white", linewidth = 0.25)
        )
    )

# Join the left strip and the main bar plot; the two parts share exactly the same Property factor
sp_cliff <- p_main + p_cat +
    plot_layout(design = "BA", widths = c(0.03, 1), guides = "collect") &
    theme(
        legend.position = "right",
        legend.box = "vertical",
        legend.justification = "center",
        legend.margin = margin(l = -8, r = 0),
        plot.margin = margin(2, -2, 1.5, 2)
    )

print(sp_cliff)

# 12. Main Figure 6 assembly and export ----
label_style <- list(size = 12, face = "bold")
row1 <- ggarrange(p_b, p_c, nrow = 2, heights = c(1, 1), font.label = label_style, labels = c("b", "c"), label.x = 0, label.y = 1, hjust = -0.2, vjust = 1.2)
row2 <- ggarrange(p_cdf, sp_cliff, nrow = 1, widths = c(1, 2.05), font.label = label_style, labels = c("d", "e"), label.x = 0, label.y = 1, hjust = -0.2, vjust = 1.2)
merge_fig6 <- ggarrange(row1, row2, ncol = 1, heights = c(1.4, 1))
ggsave("figures/fig6_integration_guidance.pdf", merge_fig6, width = 10, height = 8)
ggsave("figures/fig6_integration_guidance.png", merge_fig6, width = 10, height = 8, dpi = 600, bg = "white")

# Save reusable Figure 6 results for Extended Data Figures 8–9 ----
fig6_ed_inputs <- list(
    analysis_version = "fig6_final_20260724",
    thresholds = list(sMAPE = SMAPE_CUTOFF, TRC_deviation = TRC_DEVIATION_CUTOFF),
    denominator_audit = denominator_audit,
    method_feature_results = method_feature_results,
    method_pair_summary = method_pair_summary,
    srr_reference_feature_results = method_feature_results %>% filter(Method %in% c("SRR (P)", "SRR (N)")),
    srr_reference_pair_summary = method_pair_summary %>% filter(Method %in% c("SRR (P)", "SRR (N)")),
    consensus_voting = consensus_voting,
    consensus_voting_class = consensus_voting_class,
    physicochemical_results = df_physchem_summary
)

saveRDS(fig6_ed_inputs, "results/fig6_extended_data_inputs.rds")
message("Reusable Figure 6 results were saved to results/fig6_extended_data_inputs.rds")

# 13. Source data ----
# Run after method_pair_summary, exclusive_stats_raw, consensus_voting_class and df_physchem_summary have been created.

if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Please install required package: openxlsx")
required_objects <- c("method_pair_summary", "exclusive_stats_raw", "consensus_voting_class", "df_physchem_summary", "physchem_matrix", "physchem_dict", "final_physchem", "cat_colors")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1), inherits = TRUE)]
if (length(missing_objects)) stop("Missing objects: ", paste(missing_objects, collapse = ", "))

fmt_int <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_range <- function(median_n, min_n, max_n) paste0(fmt_int(median_n), " [", fmt_int(min_n), "–", fmt_int(max_n), "]")
scenario_levels <- c("Intra-DIA", "Intra-OLK", "Intra-SOM", "DIA-OLK", "DIA-SOM", "OLK-SOM")
design_levels_export <- c("Balanced", "Partial", "Confounded")

# 13.1 Integration performance ----
integration_performance <- method_pair_summary %>%
    mutate(Scenario = factor(Detailed_Type, levels = scenario_levels), Design = factor(Design, levels = design_levels_export),
           Reference = coalesce(as.character(MethodReference), "None"),
           `Validation scheme` = case_when(
               Design != "Confounded" ~ "Overlapping study samples",
               Reference == "P" ~ "P anchor; N validation",
               Reference == "N" ~ "N anchor; P validation",
               TRUE ~ "Dual-track P-anchor/N-validation and N-anchor/P-validation; equally weighted"
           )) %>%
    group_by(Scenario, Design, Strategy, Method, Reference, `Validation scheme`) %>%
    summarize(Pairs = n(), Protein_Median = median(N_Features), Protein_Min = min(N_Features), Protein_Max = max(N_Features),
              `Quantitative agreement (%)` = median(Rate_sMAPE, na.rm = TRUE), `Expected response (%)` = median(Rate_Response, na.rm = TRUE),
              `Harmonized (%)` = median(Rate_Harmonized, na.rm = TRUE), .groups = "drop") %>%
    mutate(`Proteins per pair` = fmt_range(Protein_Median, Protein_Min, Protein_Max)) %>%
    select(Scenario, Design, Strategy, Method, Reference, `Validation scheme`, Pairs, `Proteins per pair`,
           `Quantitative agreement (%)`, `Expected response (%)`, `Harmonized (%)`) %>%
    arrange(Scenario, Design, Strategy, Method)

# 13.2 Exclusive harmonization outcomes ----
exclusive_outcome_levels <- c("Shared success", "Native", "SRR", "SRR+RF-BECA", "RF-BECA", "RI-BECA", "Unresolved")
exclusive_outcome_pair <- exclusive_stats_raw %>%
    mutate(Scenario = factor(Detailed_Type, levels = scenario_levels), Design = factor(Design, levels = design_levels_export),
           Category = factor(Category, levels = exclusive_outcome_levels)) %>%
    count(Scenario, Design, Batch1, Batch2, Category, name = "N_Features") %>%
    complete(nesting(Scenario, Design, Batch1, Batch2), Category = factor(exclusive_outcome_levels, levels = exclusive_outcome_levels), fill = list(N_Features = 0L)) %>%
    group_by(Scenario, Design, Batch1, Batch2) %>% mutate(Proteins = sum(N_Features), Proportion = N_Features / Proteins) %>% ungroup()

exclusive_outcome_sizes <- exclusive_outcome_pair %>% distinct(Scenario, Design, Batch1, Batch2, Proteins) %>%
    group_by(Scenario, Design) %>%
    summarize(Pairs = n(), Protein_Median = median(Proteins), Protein_Min = min(Proteins), Protein_Max = max(Proteins), .groups = "drop") %>%
    mutate(`Proteins per pair` = fmt_range(Protein_Median, Protein_Min, Protein_Max)) %>% select(-Protein_Median, -Protein_Min, -Protein_Max)

exclusive_outcome_proportions <- exclusive_outcome_pair %>% group_by(Scenario, Design, Category) %>% summarize(Value = mean(Proportion), .groups = "drop") %>%
    pivot_wider(names_from = Category, values_from = Value, names_glue = "{Category} (%)")
exclusive_outcomes <- exclusive_outcome_sizes %>% left_join(exclusive_outcome_proportions, by = c("Scenario", "Design")) %>%
    select(Scenario, Design, Pairs, `Proteins per pair`, all_of(paste0(exclusive_outcome_levels, " (%)"))) %>% arrange(Scenario, Design)
exclusive_outcome_total <- rowSums(exclusive_outcomes %>% select(all_of(paste0(exclusive_outcome_levels, " (%)"))), na.rm = TRUE)
if (any(abs(exclusive_outcome_total - 1) > 1e-10)) stop("Exclusive-outcome proportions do not sum to 100%.")

# 13.3 Physicochemical associations ----
physchem_association_counts <- map_dfr(unique(as.character(consensus_voting_class$Detailed_Type)), function(dt) {
    df_sub <- consensus_voting_class %>% filter(as.character(Detailed_Type) == dt) %>% inner_join(physchem_matrix, by = c("UniProtID" = "Entry"))
    map_dfr(final_physchem, function(feat) tibble(
        Detailed_Type = dt, Feature = feat,
        `Success n` = sum(is.finite(suppressWarnings(as.numeric(df_sub[[feat]][df_sub$Integration_Status == "Success"])))),
        `Failure n` = sum(is.finite(suppressWarnings(as.numeric(df_sub[[feat]][df_sub$Integration_Status == "Failed"]))))
    ))
})

physchem_association_thresholds <- consensus_voting_class %>% group_by(Detailed_Type) %>%
    summarize(`Success threshold (%)` = first(Threshold_Q75), .groups = "drop") %>% mutate(Detailed_Type = as.character(Detailed_Type))
dictionary_order <- physchem_dict %>% mutate(Dictionary_Order = row_number()) %>% select(Feature, Dictionary_Order)

physicochemical_associations <- df_physchem_summary %>%
    left_join(physchem_association_counts, by = c("Detailed_Type", "Feature")) %>% left_join(physchem_association_thresholds, by = "Detailed_Type") %>%
    left_join(dictionary_order, by = "Feature") %>% group_by(Detailed_Type) %>% mutate(FDR = p.adjust(P_Value, method = "BH")) %>% ungroup() %>%
    transmute(Scenario = factor(Detailed_Type, levels = scenario_levels), Category = factor(Category, levels = names(cat_colors)), Property,
              `Success threshold (%)`, `Success n`, `Failure n`, `Cliff's delta (δ)` = Cliff_Delta, `P-value` = P_Value, FDR,
              `Enrichment status` = Fill_Status, Significance = Stars, Dictionary_Order) %>%
    arrange(Scenario, Category, Dictionary_Order, Property) %>% select(-Dictionary_Order)

# Titles, descriptions and column definitions ----
table_metadata <- tribble(
    ~Table, ~Title, ~Description,
    "Integration performance", "Integration performance.", "Pair-level quantitative agreement, expected-response retention and overall harmonization rates across methods, scenarios and study designs.",
    "Exclusive outcomes", "Exclusive harmonization outcomes.", "Pair-weighted proportions of proteins assigned to mutually exclusive harmonization outcomes under the overall success criterion.",
    "Physicochemical associations", "Physicochemical associations.", "Cliff’s δ effect sizes and two-sided Wilcoxon rank-sum P values comparing physicochemical properties between high-consensus and zero-success proteins."
)

definitions <- bind_rows(
    tribble(
        ~Table, ~Column, ~Definition,
        "Integration performance", "Scenario", "Within- or cross-platform harmonization scenario.",
        "Integration performance", "Design", "Study-sample overlap design: Balanced, Partial or Confounded.",
        "Integration performance", "Strategy", "Methodological strategy class, retaining P- and N-anchored variants where applicable.",
        "Integration performance", "Method", "Specific harmonization algorithm.",
        "Integration performance", "Reference", "Reference sample used by the method (P or N); None denotes methods without an explicit reference anchor.",
        "Integration performance", "Validation scheme", "Samples used to assess quantitative agreement. In Confounded designs, P-anchored methods are validated on N, N-anchored methods on P, and Native/RF-BECA results are equally averaged across both tracks.",
        "Integration performance", "Pairs", "Number of evaluable batch pairs.",
        "Integration performance", "Proteins per pair", "Median number of evaluated proteins across batch pairs, followed by the minimum–maximum range in brackets.",
        "Integration performance", "Quantitative agreement (%)", "Median pair-level proportion of proteins with sMAPE ≤ 30%.",
        "Integration performance", "Expected response (%)", "Median pair-level proportion satisfying monotonic titration response, at least two finite TRC estimates and mean TRC deviation < 0.25.",
        "Integration performance", "Harmonized (%)", "Median pair-level proportion satisfying both quantitative agreement and expected-response criteria."
    ),
    tribble(
        ~Table, ~Column, ~Definition,
        "Exclusive outcomes", "Scenario", "Within- or cross-platform harmonization scenario.",
        "Exclusive outcomes", "Design", "Study-sample overlap design: Balanced, Partial or Confounded.",
        "Exclusive outcomes", "Pairs", "Number of evaluable batch pairs.",
        "Exclusive outcomes", "Proteins per pair", "Median number of evaluated proteins across batch pairs, followed by the minimum–maximum range in brackets.",
        "Exclusive outcomes", "Shared success (%)", "Mean pair-level proportion successful with at least one SRR-enabled strategy and at least one non-SRR strategy.",
        "Exclusive outcomes", "Native (%)", "Mean pair-level proportion satisfying the overall success criterion only without additional harmonization.",
        "Exclusive outcomes", "SRR (%)", "Mean pair-level proportion successful with at least one base SRR method, with no non-SRR strategy succeeding.",
        "Exclusive outcomes", "SRR+RF-BECA (%)", "Mean pair-level proportion unsuccessful with base SRR and non-SRR strategies but successful with at least one SRR+RF-BECA combination.",
        "Exclusive outcomes", "RF-BECA (%)", "Mean pair-level proportion successful with at least one reference-free batch-effect correction method after preceding categories are excluded.",
        "Exclusive outcomes", "RI-BECA (%)", "Mean pair-level proportion successful with at least one reference-informed RUVs method after preceding categories are excluded.",
        "Exclusive outcomes", "Unresolved (%)", "Mean pair-level proportion not satisfying the overall success criterion under any evaluated strategy."
    ),
    tribble(
        ~Table, ~Column, ~Definition,
        "Physicochemical associations", "Scenario", "Within- or cross-platform harmonization scenario.",
        "Physicochemical associations", "Category", "Physicochemical property class from the curated annotation dictionary.",
        "Physicochemical associations", "Property", "Display name of the tested physicochemical property.",
        "Physicochemical associations", "Success threshold (%)", "Scenario-specific 75th percentile of consensus success rate used to define the high-consensus group.",
        "Physicochemical associations", "Success n", "Number of high-consensus proteins with a finite value for the tested property.",
        "Physicochemical associations", "Failure n", "Number of zero-success proteins with a finite value for the tested property.",
        "Physicochemical associations", "Cliff's delta (δ)", "Non-parametric effect size; positive values indicate higher values in high-consensus proteins and negative values indicate higher values in zero-success proteins.",
        "Physicochemical associations", "P-value", "Nominal two-sided Wilcoxon rank-sum P value.",
        "Physicochemical associations", "FDR", "Benjamini–Hochberg adjusted P value within each scenario.",
        "Physicochemical associations", "Enrichment status", "Direction based on the sign of Cliff’s δ for nominally significant associations; otherwise Non-significant.",
        "Physicochemical associations", "Significance", "Nominal significance annotation: *P < 0.05, **P < 0.01 and ***P < 0.001."
    )
) %>% left_join(table_metadata, by = "Table") %>% select(Table, Title, Description, Column, Definition)

# Excel workbook ----
wb <- createWorkbook()
title_style <- createStyle(fontSize = 12, textDecoration = "bold", fontColour = "#FFFFFF", fgFill = "#1F4E78", halign = "left", valign = "center")
desc_style <- createStyle(fontSize = 10, fontColour = "#404040", fgFill = "#D9EAF7", wrapText = TRUE, valign = "center")
wrap_style <- createStyle(wrapText = TRUE, valign = "top")
percent_style <- createStyle(numFmt = "0.0%")
integer_style <- createStyle(numFmt = "#,##0")
delta_style <- createStyle(numFmt = "0.000")
pvalue_style <- createStyle(numFmt = "0.00E+00")

write_source_sheet <- function(sheet, title, description, data, widths, percent_cols = character(), integer_cols = character(), delta_cols = character(), pvalue_cols = character()) {
    addWorksheet(wb, sheet, gridLines = FALSE); n_col <- ncol(data); n_row <- nrow(data)
    mergeCells(wb, sheet, cols = 1:n_col, rows = 1); writeData(wb, sheet, title, startRow = 1, startCol = 1)
    addStyle(wb, sheet, title_style, rows = 1, cols = 1:n_col, gridExpand = TRUE)
    mergeCells(wb, sheet, cols = 1:n_col, rows = 2); writeData(wb, sheet, description, startRow = 2, startCol = 1)
    addStyle(wb, sheet, desc_style, rows = 2, cols = 1:n_col, gridExpand = TRUE)
    writeDataTable(wb, sheet, data, startRow = 4, tableStyle = "TableStyleMedium2", withFilter = TRUE)
    addStyle(wb, sheet, wrap_style, rows = 5:(n_row + 4), cols = 1:n_col, gridExpand = TRUE, stack = TRUE)
    if (length(percent_cols)) addStyle(wb, sheet, percent_style, rows = 5:(n_row + 4), cols = match(percent_cols, names(data)), gridExpand = TRUE, stack = TRUE)
    if (length(integer_cols)) addStyle(wb, sheet, integer_style, rows = 5:(n_row + 4), cols = match(integer_cols, names(data)), gridExpand = TRUE, stack = TRUE)
    if (length(delta_cols)) addStyle(wb, sheet, delta_style, rows = 5:(n_row + 4), cols = match(delta_cols, names(data)), gridExpand = TRUE, stack = TRUE)
    if (length(pvalue_cols)) addStyle(wb, sheet, pvalue_style, rows = 5:(n_row + 4), cols = match(pvalue_cols, names(data)), gridExpand = TRUE, stack = TRUE)
    setColWidths(wb, sheet, cols = 1:n_col, widths = widths); setRowHeights(wb, sheet, rows = 1, heights = 24)
    setRowHeights(wb, sheet, rows = 2, heights = 42); freezePane(wb, sheet, firstActiveRow = 5)
}

write_source_sheet("Integration_performance", table_metadata$Title[1], table_metadata$Description[1], integration_performance,
               c(14, 12, 22, 24, 10, 42, 9, 22, 18, 18, 16),
               percent_cols = c("Quantitative agreement (%)", "Expected response (%)", "Harmonized (%)"), integer_cols = "Pairs")
write_source_sheet("Exclusive_outcomes", table_metadata$Title[2], table_metadata$Description[2], exclusive_outcomes,
               c(14, 12, 9, 22, rep(18, 7)), percent_cols = paste0(exclusive_outcome_levels, " (%)"), integer_cols = "Pairs")
write_source_sheet("Physicochemical_associations", table_metadata$Title[3], table_metadata$Description[3], physicochemical_associations,
               c(14, 14, 34, 18, 11, 11, 17, 14, 14, 20, 12), percent_cols = "Success threshold (%)",
               integer_cols = c("Success n", "Failure n"), delta_cols = "Cliff's delta (δ)", pvalue_cols = c("P-value", "FDR"))
write_source_sheet("Definitions", "Source-data definitions", "Titles, descriptions and column definitions for the Figure 6 source-data sheets.", definitions,
               c(10, 30, 58, 28, 86))

output_file <- "tables/SourceData_Figure6.xlsx"
saveWorkbook(wb, output_file, overwrite = TRUE)
message("Figure 6 source data were exported to ", output_file)
