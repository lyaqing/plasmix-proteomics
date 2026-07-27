# Build equal-depth replicate plans while preserving acquisition-round structure.
# Structurally unavailable wells may be replaced by the nearest unused valid replicate;
# feature-level missing values are never used to alter the replicate plan.
make_replicate_plan <- function(sample_metadata, available_columns, n_replicates = 3,
                                required_samples = character()) {
    required <- c("Sample", "Replicate", "ColName")
    missing_columns <- setdiff(required, colnames(sample_metadata))
    if (length(missing_columns)) {
        columns <- paste(missing_columns, collapse = ", ")
        stop(sprintf("sample_metadata is missing required columns: %s.", columns), call. = FALSE)
    }
    if (length(n_replicates) != 1 || !is.finite(n_replicates) || n_replicates < 2) {
        stop("n_replicates must be one integer of at least 2.", call. = FALSE)
    }

    metadata <- sample_metadata |>
        dplyr::filter(!is.na(.data$Sample), !is.na(.data$Replicate), !is.na(.data$ColName)) |>
        dplyr::mutate(Replicate = suppressWarnings(as.numeric(as.character(.data$Replicate)))) |>
        dplyr::filter(is.finite(.data$Replicate)) |>
        dplyr::distinct(.data$Sample, .data$Replicate, .data$ColName, .keep_all = TRUE)
    if (nrow(metadata) == 0) return(list())

    included <- rep(TRUE, nrow(metadata))
    if ("Include" %in% colnames(metadata)) {
        include_value <- as.character(metadata$Include)
        included <- is.na(metadata$Include) | include_value %in% c("TRUE", "True", "true", "1", "YES", "Yes", "yes")
    }
    metadata$StructurallyAvailable <- included & metadata$ColName %in% available_columns

    eligible_samples <- metadata |>
        dplyr::group_by(.data$Sample) |>
        dplyr::summarize(
            Available_replicates = dplyr::n_distinct(.data$Replicate[.data$StructurallyAvailable]),
            .groups = "drop"
        ) |>
        dplyr::filter(.data$Available_replicates >= n_replicates) |>
        dplyr::pull(.data$Sample)
    if (!length(eligible_samples) || !all(required_samples %in% eligible_samples)) return(list())

    metadata <- metadata |>
        dplyr::filter(.data$Sample %in% eligible_samples)

    designed_sets <- metadata |>
        dplyr::group_by(.data$Sample) |>
        dplyr::summarize(Replicates = list(sort(unique(.data$Replicate))), .groups = "drop") |>
        dplyr::pull(.data$Replicates)
    if (!length(designed_sets)) return(list())

    design_replicates <- Reduce(intersect, designed_sets)
    if (length(design_replicates) < n_replicates) return(list())
    target_combinations <- utils::combn(sort(design_replicates), n_replicates, simplify = FALSE)

    resolve_sample <- function(sample_name, target_replicates) {
        sample_rows <- metadata |>
            dplyr::filter(.data$Sample == sample_name, .data$StructurallyAvailable) |>
            dplyr::arrange(.data$Replicate, .data$ColName) |>
            dplyr::distinct(.data$Replicate, .keep_all = TRUE)
        if (nrow(sample_rows) < n_replicates) return(NULL)

        direct_replicates <- intersect(target_replicates, sample_rows$Replicate)
        used_replicates <- direct_replicates
        chosen_rows <- vector("list", length(target_replicates))

        for (position in seq_along(target_replicates)) {
            target_replicate <- target_replicates[position]
            direct_row <- sample_rows |>
                dplyr::filter(.data$Replicate == target_replicate)

            if (nrow(direct_row) > 0) {
                chosen <- direct_row[1, , drop = FALSE]
            } else {
                candidates <- sample_rows |>
                    dplyr::filter(!.data$Replicate %in% used_replicates) |>
                    dplyr::mutate(
                        Distance = abs(.data$Replicate - target_replicate),
                        Forward = .data$Replicate >= target_replicate
                    ) |>
                    dplyr::arrange(.data$Distance, dplyr::desc(.data$Forward), .data$Replicate)
                if (nrow(candidates) == 0) return(NULL)
                chosen <- candidates[1, setdiff(colnames(candidates), c("Distance", "Forward")), drop = FALSE]
                used_replicates <- c(used_replicates, chosen$Replicate)
            }

            chosen_rows[[position]] <- chosen |>
                dplyr::mutate(
                    TargetReplicate = target_replicate,
                    SelectedReplicate = .data$Replicate,
                    WasSubstituted = .data$Replicate != target_replicate
                )
        }
        dplyr::bind_rows(chosen_rows)
    }

    plans <- lapply(target_combinations, function(target_replicates) {
        resolved <- lapply(eligible_samples, resolve_sample, target_replicates = target_replicates)
        if (any(vapply(resolved, is.null, logical(1)))) return(NULL)
        resolved <- dplyr::bind_rows(resolved)
        if (nrow(resolved) != length(eligible_samples) * n_replicates) return(NULL)

        substitutions <- resolved |>
            dplyr::filter(.data$WasSubstituted) |>
            dplyr::transmute(
                Label = paste0(.data$Sample, ":", .data$TargetReplicate, "->", .data$SelectedReplicate)
            ) |>
            dplyr::pull(.data$Label)

        list(
            TargetReplicates = target_replicates,
            Metadata = resolved,
            SubstitutionMap = if (length(substitutions)) paste(substitutions, collapse = ";") else "None"
        )
    })
    Filter(Negate(is.null), plans)
}


# Calculate the original PCA-based SNR used in the Plasmix Figure 2 workflow.
calc_pca_snr <- function(expr_mat, group, center = TRUE, scale. = FALSE) {
    expr_mat <- as.matrix(expr_mat)
    if (!is.numeric(expr_mat)) stop("expr_mat must contain only numeric values.", call. = FALSE)
    if (is.null(colnames(expr_mat))) stop("expr_mat must have sample names in its column names.", call. = FALSE)
    if (length(group) != ncol(expr_mat)) stop("group must contain one value per sample column.", call. = FALSE)
    if (anyNA(group)) stop("group contains missing values.", call. = FALSE)
    if (anyNA(expr_mat) || any(!is.finite(expr_mat))) stop("expr_mat must be complete and finite.", call. = FALSE)

    group <- as.character(group)
    if (length(unique(group)) < 2) stop("At least two sample groups are required.", call. = FALSE)
    if (!any(table(group) >= 2)) stop("At least one sample group must contain replicates.", call. = FALSE)

    feature_var <- apply(expr_mat, 1, stats::var)
    expr_mat <- expr_mat[is.finite(feature_var) & feature_var > 0, , drop = FALSE]
    if (nrow(expr_mat) < 3) stop("Fewer than three variable features remain for PCA.", call. = FALSE)

    pca <- stats::prcomp(t(expr_mat), retx = TRUE, center = center, scale. = scale.)
    if (ncol(pca$x) < 2) stop("PCA returned fewer than two components.", call. = FALSE)

    weights <- summary(pca)$importance[2, 1:2]
    pairs <- utils::combn(seq_len(nrow(pca$x)), 2)
    delta <- pca$x[pairs[1, ], 1:2, drop = FALSE] - pca$x[pairs[2, ], 1:2, drop = FALSE]
    distance <- rowSums(sweep(delta^2, 2, weights, `*`))
    same_group <- group[pairs[1, ]] == group[pairs[2, ]]

    intra <- distance[same_group]
    inter <- distance[!same_group]
    if (!length(intra) || !length(inter)) {
        stop("Both within-group and between-group pairs are required.", call. = FALSE)
    }

    snr <- mean(inter) / mean(intra)
    if (!is.finite(snr) || snr <= 0) return(NA_real_)
    10 * log10(snr)
}


# Calculate feature-level replicate CV in linear space and summarize across available groups.
calc_feature_cv <- function(expr_mat, group, log2_input = TRUE, min_replicates = 2) {
    expr_mat <- as.matrix(expr_mat)
    if (!is.numeric(expr_mat)) stop("expr_mat must contain only numeric values.", call. = FALSE)
    if (length(group) != ncol(expr_mat)) stop("group must contain one value per sample column.", call. = FALSE)
    if (anyNA(group)) stop("group contains missing values.", call. = FALSE)

    linear_mat <- if (log2_input) 2^expr_mat else expr_mat
    group <- as.character(group)
    eligible_groups <- names(which(table(group) >= min_replicates))
    if (!length(eligible_groups)) {
        result <- rep(NA_real_, nrow(linear_mat))
        names(result) <- rownames(linear_mat)
        return(result)
    }

    cv_by_group <- vapply(eligible_groups, function(current_group) {
        values <- linear_mat[, group == current_group, drop = FALSE]
        valid_n <- rowSums(is.finite(values))
        cv <- apply(values, 1, stats::sd, na.rm = TRUE) / rowMeans(values, na.rm = TRUE)
        cv[valid_n < min_replicates | !is.finite(cv)] <- NA_real_
        cv
    }, numeric(nrow(linear_mat)))

    result <- if (is.null(dim(cv_by_group))) cv_by_group else rowMeans(cv_by_group, na.rm = TRUE)
    result[!is.finite(result)] <- NA_real_
    names(result) <- rownames(linear_mat)
    result
}


# Calculate feature-level titration centers and TRCs within one replicate subsample.
calc_feature_titration_metrics <- function(df_long, method = c("Mean", "Median"), is_log2 = TRUE,
                                           trc_deviation_cutoff = 0.25, min_valid_replicates = 2,
                                           trc_error_cutoff = NULL) {
    required <- c("UniqueID", "Sample", "Value")
    missing_columns <- setdiff(required, colnames(df_long))
    if (length(missing_columns)) {
        columns <- paste(missing_columns, collapse = ", ")
        stop(sprintf("df_long is missing required columns: %s.", columns), call. = FALSE)
    }
    if (!is.null(trc_error_cutoff)) trc_deviation_cutoff <- trc_error_cutoff
    if (length(trc_deviation_cutoff) != 1 || !is.finite(trc_deviation_cutoff) || trc_deviation_cutoff <= 0) {
        stop("trc_deviation_cutoff must be one positive finite value.", call. = FALSE)
    }
    if (length(min_valid_replicates) != 1 || !is.finite(min_valid_replicates) || min_valid_replicates < 1) {
        stop("min_valid_replicates must be one positive integer.", call. = FALSE)
    }

    method <- match.arg(method)
    target_order <- c("M", "Y", "P", "X", "F")
    df_agg <- df_long |>
        dplyr::filter(.data$Sample %in% target_order) |>
        dplyr::group_by(.data$UniqueID, .data$Sample) |>
        dplyr::summarize(
            val = {
                values <- .data$Value[is.finite(.data$Value)]
                if (length(values) < min_valid_replicates) {
                    NA_real_
                } else if (method == "Mean") {
                    mean(values)
                } else {
                    stats::median(values)
                }
            },
            .groups = "drop"
        ) |>
        tidyr::pivot_wider(names_from = "Sample", values_from = "val")

    for (sample_name in target_order) {
        if (!sample_name %in% colnames(df_agg)) df_agg[[sample_name]] <- NA_real_
    }

    if (is_log2) {
        m_value <- 2^df_agg$M
        y_value <- 2^df_agg$Y
        p_value <- 2^df_agg$P
        x_value <- 2^df_agg$X
        f_value <- 2^df_agg$F
    } else {
        m_value <- df_agg$M
        y_value <- df_agg$Y
        p_value <- df_agg$P
        x_value <- df_agg$X
        f_value <- df_agg$F
    }

    denominator <- m_value - f_value
    denominator[is.finite(denominator) & denominator == 0] <- NA_real_
    df_agg$TRCAnchorValid <- is.finite(denominator)
    df_agg$TRC_Y <- (m_value - y_value) / denominator
    df_agg$TRC_P <- (m_value - p_value) / denominator
    df_agg$TRC_X <- (m_value - x_value) / denominator
    df_agg$TRC_err_Y <- df_agg$TRC_Y - 0.25
    df_agg$TRC_err_P <- df_agg$TRC_P - 0.50
    df_agg$TRC_err_X <- df_agg$TRC_X - 0.75
    df_agg$TRC_abs_err_Y <- abs(df_agg$TRC_err_Y)
    df_agg$TRC_abs_err_P <- abs(df_agg$TRC_err_P)
    df_agg$TRC_abs_err_X <- abs(df_agg$TRC_err_X)

    trc_coordinates <- cbind(
        M = ifelse(df_agg$TRCAnchorValid, 0, NA_real_),
        Y = df_agg$TRC_Y,
        P = df_agg$TRC_P,
        X = df_agg$TRC_X,
        F = ifelse(df_agg$TRCAnchorValid, 1, NA_real_)
    )
    df_agg$TitrationMonoValid <- apply(trc_coordinates, 1, function(values) {
        values <- values[is.finite(values)]
        length(values) >= 3 && all(diff(values) > 0)
    })

    deviation_mat <- as.matrix(df_agg[, c("TRC_abs_err_Y", "TRC_abs_err_P", "TRC_abs_err_X")])
    df_agg$TRC_N_finite <- rowSums(is.finite(deviation_mat))
    df_agg$MeanTRCDev <- apply(deviation_mat, 1, function(values) {
        if (any(is.finite(values))) mean(values[is.finite(values)]) else NA_real_
    })
    df_agg$TRCDevValid <- df_agg$TRC_N_finite >= 2 & is.finite(df_agg$MeanTRCDev) &
        df_agg$MeanTRCDev < trc_deviation_cutoff
    df_agg$ExpectedResponseValid <- df_agg$TitrationMonoValid & df_agg$TRCDevValid

    # Deprecated compatibility aliases for scripts that have not yet been updated.
    df_agg$MCValid <- df_agg$TitrationMonoValid
    df_agg$TRC_AE_mean <- df_agg$MeanTRCDev
    df_agg$TRCValid <- df_agg$TRCDevValid
    df_agg$TitrationValid <- df_agg$ExpectedResponseValid
    df_agg
}


# Summarize equal-depth replicate subsamples using majority support for every adjacent TRC relation.
summarize_titration_subsamples <- function(df_iterations,
                                           group_cols = c("Batch", "Platform", "ProcessLevel", "DataTier"),
                                           majority_cutoff = 0.5, trc_deviation_cutoff = 0.25) {
    required <- c("UniqueID", "Iteration", "M", "Y", "P", "X", "F",
                  "TRC_Y", "TRC_P", "TRC_X",
                  "TRC_abs_err_Y", "TRC_abs_err_P", "TRC_abs_err_X")
    missing_columns <- setdiff(c(required, group_cols), colnames(df_iterations))
    if (length(missing_columns)) {
        columns <- paste(missing_columns, collapse = ", ")
        stop(sprintf("df_iterations is missing required columns: %s.", columns), call. = FALSE)
    }

    safe_mean <- function(values) {
        values <- values[is.finite(values)]
        if (!length(values)) NA_real_ else mean(values)
    }
    nominal_positions <- c(M = 0, Y = 0.25, P = 0.50, X = 0.75, F = 1)
    feature_keys <- c(group_cols, "UniqueID")

    feature_summary <- df_iterations |>
        dplyr::group_by(dplyr::across(dplyr::all_of(feature_keys))) |>
        dplyr::summarize(
            M = safe_mean(.data$M), Y = safe_mean(.data$Y), P = safe_mean(.data$P),
            X = safe_mean(.data$X), F = safe_mean(.data$F),
            TRC_Y = safe_mean(.data$TRC_Y), TRC_P = safe_mean(.data$TRC_P), TRC_X = safe_mean(.data$TRC_X),
            TRC_abs_err_Y = safe_mean(.data$TRC_abs_err_Y),
            TRC_abs_err_P = safe_mean(.data$TRC_abs_err_P),
            TRC_abs_err_X = safe_mean(.data$TRC_abs_err_X),
            CompleteChainSupport = safe_mean(as.numeric(.data$TitrationMonoValid)),
            Iterations = dplyr::n_distinct(.data$Iteration),
            .groups = "drop"
        ) |>
        dplyr::mutate(
            TRC_err_Y = .data$TRC_Y - 0.25,
            TRC_err_P = .data$TRC_P - 0.50,
            TRC_err_X = .data$TRC_X - 0.75,
            TRC_N_finite = rowSums(cbind(
                is.finite(.data$TRC_abs_err_Y),
                is.finite(.data$TRC_abs_err_P),
                is.finite(.data$TRC_abs_err_X)
            )),
            MeanTRCDev = rowMeans(cbind(
                .data$TRC_abs_err_Y,
                .data$TRC_abs_err_P,
                .data$TRC_abs_err_X
            ), na.rm = TRUE),
            MeanTRCDev = dplyr::if_else(.data$TRC_N_finite > 0, .data$MeanTRCDev, NA_real_),
            TRCDevValid = .data$TRC_N_finite >= 2 & is.finite(.data$MeanTRCDev) &
                .data$MeanTRCDev < trc_deviation_cutoff,
            CompleteChainMajorityValid = is.finite(.data$CompleteChainSupport) &
                .data$CompleteChainSupport > majority_cutoff
        )

    anchor_valid <- if ("TRCAnchorValid" %in% colnames(df_iterations)) {
        df_iterations$TRCAnchorValid %in% TRUE
    } else {
        is.finite(df_iterations$M) & is.finite(df_iterations$F) & df_iterations$M != df_iterations$F
    }

    trc_long <- df_iterations |>
        dplyr::mutate(
            .TRC_M = dplyr::if_else(anchor_valid, 0, NA_real_),
            .TRC_F = dplyr::if_else(anchor_valid, 1, NA_real_)
        ) |>
        dplyr::transmute(
            dplyr::across(dplyr::all_of(feature_keys)),
            Iteration = .data$Iteration,
            M = .data$.TRC_M, Y = .data$TRC_Y, P = .data$TRC_P,
            X = .data$TRC_X, F = .data$.TRC_F
        ) |>
        tidyr::pivot_longer(c("M", "Y", "P", "X", "F"), names_to = "Sample", values_to = "TRC")

    available_samples <- trc_long |>
        dplyr::group_by(dplyr::across(dplyr::all_of(c(feature_keys, "Sample")))) |>
        dplyr::summarize(Available = any(is.finite(.data$TRC)), .groups = "drop") |>
        dplyr::filter(.data$Available) |>
        dplyr::mutate(NominalTRC = unname(nominal_positions[as.character(.data$Sample)]))

    feature_relations <- available_samples |>
        dplyr::group_by(dplyr::across(dplyr::all_of(feature_keys))) |>
        dplyr::group_modify(function(data, key) {
            data <- data |>
                dplyr::filter(is.finite(.data$NominalTRC)) |>
                dplyr::arrange(.data$NominalTRC)
            samples <- data$Sample
            if (length(samples) < 3 || !all(c("M", "F") %in% samples)) {
                return(tibble::tibble(LeftSample = character(), RightSample = character(),
                                      Relation = character(), LeftNominal = numeric(), RightNominal = numeric()))
            }
            tibble::tibble(
                LeftSample = head(samples, -1),
                RightSample = tail(samples, -1),
                Relation = paste0(head(samples, -1), "-", tail(samples, -1)),
                LeftNominal = head(data$NominalTRC, -1),
                RightNominal = tail(data$NominalTRC, -1)
            )
        }) |>
        dplyr::ungroup()

    if (nrow(feature_relations) > 0) {
        left_values <- trc_long |>
            dplyr::rename(LeftSample = .data$Sample, LeftTRC = .data$TRC)
        right_values <- trc_long |>
            dplyr::rename(RightSample = .data$Sample, RightTRC = .data$TRC)

        relation_iteration <- feature_relations |>
            dplyr::inner_join(left_values, by = c(feature_keys, "LeftSample")) |>
            dplyr::left_join(right_values,
                             by = c(feature_keys, "Iteration", "RightSample")) |>
            dplyr::mutate(
                Evaluable = is.finite(.data$LeftTRC) & is.finite(.data$RightTRC),
                Supported = .data$Evaluable & .data$LeftTRC < .data$RightTRC
            )

        relation_summary <- relation_iteration |>
            dplyr::group_by(dplyr::across(dplyr::all_of(c(
                feature_keys, "Relation", "LeftSample", "RightSample", "LeftNominal", "RightNominal"
            )))) |>
            dplyr::summarize(
                EvaluableIterations = sum(.data$Evaluable),
                SupportedIterations = sum(.data$Supported),
                .groups = "drop"
            ) |>
            dplyr::mutate(
                SupportRate = dplyr::if_else(
                    .data$EvaluableIterations > 0,
                    .data$SupportedIterations / .data$EvaluableIterations,
                    NA_real_
                )
            )

        monotonicity_summary <- relation_summary |>
            dplyr::group_by(dplyr::across(dplyr::all_of(feature_keys))) |>
            dplyr::summarize(
                MonoRelationN = dplyr::n(),
                MonoRelationEvaluableN = sum(.data$EvaluableIterations > 0),
                MeanMonoSupport = if (all(.data$EvaluableIterations > 0)) mean(.data$SupportRate) else NA_real_,
                MinMonoSupport = if (all(.data$EvaluableIterations > 0)) min(.data$SupportRate) else NA_real_,
                .groups = "drop"
            ) |>
            dplyr::mutate(
                TitrationMonoValid = .data$MonoRelationN >= 2 &
                    .data$MonoRelationEvaluableN == .data$MonoRelationN &
                    is.finite(.data$MinMonoSupport) & .data$MinMonoSupport > majority_cutoff
            )
    } else {
        relation_summary <- feature_summary |> dplyr::select(dplyr::all_of(feature_keys)) |> dplyr::slice(0) |>
            dplyr::mutate(
                Relation = character(), LeftSample = character(), RightSample = character(),
                LeftNominal = numeric(), RightNominal = numeric(),
                EvaluableIterations = integer(), SupportedIterations = integer(), SupportRate = numeric()
            )
        monotonicity_summary <- feature_summary |> dplyr::select(dplyr::all_of(feature_keys)) |> dplyr::slice(0) |>
            dplyr::mutate(
                MonoRelationN = integer(), MonoRelationEvaluableN = integer(),
                MeanMonoSupport = numeric(), MinMonoSupport = numeric(),
                TitrationMonoValid = logical()
            )
    }

    feature_summary <- feature_summary |>
        dplyr::left_join(monotonicity_summary, by = feature_keys) |>
        dplyr::mutate(
            MonoRelationN = tidyr::replace_na(.data$MonoRelationN, 0L),
            MonoRelationEvaluableN = tidyr::replace_na(.data$MonoRelationEvaluableN, 0L),
            TitrationMonoValid = tidyr::replace_na(.data$TitrationMonoValid, FALSE),
            ExpectedResponseValid = .data$TitrationMonoValid & .data$TRCDevValid,
            MC_frequency = .data$CompleteChainSupport,
            MCValid = .data$TitrationMonoValid,
            TRC_AE_mean = .data$MeanTRCDev,
            TRCValid = .data$TRCDevValid,
            TitrationValid = .data$ExpectedResponseValid
        )

    list(feature = feature_summary, relation = relation_summary)
}


# Calculate batch-level gradient-fit R² for each Y, P and X mixture independently.
calc_gradient_fit <- function(df_features, group_cols = c("Batch", "Platform", "ProcessLevel", "DataTier")) {
    required <- c("UniqueID", "Is_Detected", "M", "Y", "P", "X", "F")
    missing_columns <- setdiff(required, colnames(df_features))
    if (length(missing_columns)) {
        columns <- paste(missing_columns, collapse = ", ")
        stop(sprintf("df_features is missing required columns: %s.", columns), call. = FALSE)
    }

    missing_group_cols <- setdiff(group_cols, colnames(df_features))
    if (length(missing_group_cols)) {
        columns <- paste(missing_group_cols, collapse = ", ")
        stop(sprintf("df_features is missing grouping columns: %s.", columns), call. = FALSE)
    }

    calc_fixed_prediction_r2 <- function(observed, expected) {
        valid <- is.finite(observed) & is.finite(expected)
        if (sum(valid) < 3) return(NA_real_)

        observed <- observed[valid]
        expected <- expected[valid]
        total_ss <- sum((observed - mean(observed))^2)
        if (!is.finite(total_ss) || total_ss == 0) return(NA_real_)

        1 - sum((observed - expected)^2) / total_ss
    }

    gradients <- c("Y", "P", "X")
    m_fraction <- c(Y = 0.75, P = 0.50, X = 0.25)
    f_fraction <- 1 - m_fraction
    nominal_trc <- c(Y = 0.25, P = 0.50, X = 0.75)

    df_features |>
        dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
        dplyr::group_modify(function(data, key) {
            m_minus_f <- data$M - data$F

            purrr::map_dfr(gradients, function(gradient) {
                observed_log2_ratio <- data[[gradient]] - data$F
                expected_log2_ratio <- log2(
                    m_fraction[[gradient]] * 2^m_minus_f +
                        f_fraction[[gradient]]
                )

                detected <- data$Is_Detected %in% TRUE
                all_valid <- is.finite(observed_log2_ratio) & is.finite(expected_log2_ratio)
                detected_valid <- all_valid & detected

                tibble::tibble(
                    Gradient = gradient,
                    Nominal_TRC = nominal_trc[[gradient]],
                    All_R2 = calc_fixed_prediction_r2(observed_log2_ratio, expected_log2_ratio),
                    Detected_R2 = calc_fixed_prediction_r2(
                        observed_log2_ratio[detected],
                        expected_log2_ratio[detected]
                    ),
                    All_N = sum(all_valid),
                    Detected_N = sum(detected_valid)
                )
            })
        }) |>
        dplyr::ungroup()
}
