# Batch-correction and sample-to-reference workflows used in Figure 6.
#
# M, Y, X and F are study samples; P and N are co-measured bridge references. Every method returns
# all design-available samples so non-anchor bridge agreement can be evaluated. Reference-free methods
# process P and N as ordinary samples and do not link their identities across batches. RUV-III-C uses
# technical replicate groups defined within each batch only.
#
# expr: numeric matrix or data frame with proteins in rows and samples in columns.
# metadata: sample-level data frame containing ColName, Sample and Batch.
# The main function returns corrected matrices in $expr and a transparent run log in $status.

apply_beca <- function(expr, metadata, srr_ref_types = c("P", "N"), run_mode = c("all", "final", "baseline")) {
    run_mode <- match.arg(run_mode)
    expr <- as.matrix(expr)
    required_metadata <- c("ColName", "Sample", "Batch")
    missing_metadata <- setdiff(required_metadata, names(metadata))

    if (length(missing_metadata)) {
        stop(sprintf("metadata is missing: %s.", paste(missing_metadata, collapse = ", ")), call. = FALSE)
    }
    if (!is.numeric(expr) || is.null(rownames(expr)) || is.null(colnames(expr))) {
        stop("expr must be numeric and have protein and sample names.", call. = FALSE)
    }
    if (anyDuplicated(rownames(expr)) || anyDuplicated(colnames(expr))) {
        stop("expr contains duplicated protein or sample names.", call. = FALSE)
    }
    if (anyDuplicated(metadata$ColName)) stop("metadata$ColName contains duplicated sample names.", call. = FALSE)
    if (anyNA(metadata[, required_metadata])) stop("ColName, Sample and Batch cannot contain missing values.", call. = FALSE)

    missing_samples <- setdiff(metadata$ColName, colnames(expr))
    if (length(missing_samples)) {
        stop(sprintf("metadata samples are absent from expr: %s.", paste(missing_samples, collapse = ", ")), call. = FALSE)
    }

    analysis_cols <- metadata$ColName
    analysis_expr <- expr[, analysis_cols, drop = FALSE]
    study_metadata <- metadata[metadata$Sample %in% c("M", "Y", "X", "F"), , drop = FALSE]
    if (!nrow(study_metadata)) stop("No M, Y, X or F study samples were found.", call. = FALSE)

    corrected_expr <- list(Native = as.data.frame(analysis_expr, check.names = FALSE))
    status_rows <- list(data.frame(
        Method = "Native", Reference = NA_character_, Status = "success", Error = NA_character_,
        FeaturesInput = nrow(analysis_expr), FeaturesOutput = nrow(analysis_expr), stringsAsFactors = FALSE
    ))

    run_method <- function(method, reference, input_expr, input_metadata, method_function) {
        tryCatch({
            value <- as.data.frame(method_function(input_expr, input_metadata), check.names = FALSE)
            missing_output_cols <- setdiff(input_metadata$ColName, colnames(value))
            if (length(missing_output_cols)) {
                stop(sprintf("Corrected output is missing samples: %s.", paste(missing_output_cols, collapse = ", ")), call. = FALSE)
            }
            value <- value[, input_metadata$ColName, drop = FALSE]
            status <- data.frame(
                Method = method, Reference = reference, Status = "success", Error = NA_character_,
                FeaturesInput = nrow(input_expr), FeaturesOutput = nrow(value), stringsAsFactors = FALSE
            )
            list(value = value, status = status)
        }, error = function(error) {
            status <- data.frame(
                Method = method, Reference = reference, Status = "failed", Error = conditionMessage(error),
                FeaturesInput = nrow(input_expr), FeaturesOutput = NA_integer_, stringsAsFactors = FALSE
            )
            list(value = NULL, status = status)
        })
    }

    # P and N are processed as ordinary samples; no cross-batch reference identity is supplied here.
    standard_methods <- list(
        MAD = function(expr, metadata) mad_scaling(expr, metadata, log2_input = TRUE),
        Quantile = function(expr, metadata) quantile_correction(expr, metadata),
        ComBat = function(expr, metadata) combat_correction(expr, metadata, use_model_matrix = FALSE),
        `RUV-III-C` = function(expr, metadata) ruv_correction(expr, metadata, k = 2, n_controls = 50),
        RUVg = function(expr, metadata) ruvg_correction(expr, metadata, k = 2)
    )

    if (run_mode %in% c("all", "final")) {
        for (method_name in names(standard_methods)) {
            run <- run_method(method_name, NA_character_, analysis_expr, metadata, standard_methods[[method_name]])
            status_rows[[length(status_rows) + 1]] <- run$status
            if (!is.null(run$value)) corrected_expr[[method_name]] <- run$value
        }

        # Only the selected anchor is linked across batches; the non-anchor reference remains an ordinary sample.
        for (reference in unique(srr_ref_types)) {
            reference_metadata <- metadata[metadata$Sample == reference, , drop = FALSE]
            if (!nrow(reference_metadata)) {
                status_rows[[length(status_rows) + 1]] <- data.frame(
                    Method = "RUVs", Reference = reference, Status = "skipped",
                    Error = "No reference samples were available.", FeaturesInput = nrow(analysis_expr),
                    FeaturesOutput = NA_integer_, stringsAsFactors = FALSE
                )
                next
            }

            run <- run_method("RUVs", reference, analysis_expr, metadata, function(expr, metadata) {
                ruvs_correction(expr, metadata, k = 2, ref_type = reference)
            })
            status_rows[[length(status_rows) + 1]] <- run$status
            if (!is.null(run$value)) corrected_expr[[sprintf("RUVs (%s)", reference)]] <- run$value
        }
    }

    if (run_mode %in% c("all", "baseline")) {
        # SRR uses one anchor but retains the non-anchor bridge sample for independent validation.
        for (reference in unique(srr_ref_types)) {
            reference_metadata <- metadata[metadata$Sample == reference, , drop = FALSE]
            if (!nrow(reference_metadata)) {
                status_rows[[length(status_rows) + 1]] <- data.frame(
                    Method = "SRR", Reference = reference, Status = "skipped",
                    Error = "No reference samples were available.", FeaturesInput = nrow(analysis_expr),
                    FeaturesOutput = NA_integer_, stringsAsFactors = FALSE
                )
                next
            }

            run <- run_method("SRR", reference, analysis_expr, metadata, function(expr, metadata) {
                srr_correction(expr, metadata, ref_type = reference)
            })
            status_rows[[length(status_rows) + 1]] <- run$status
            if (is.null(run$value)) next

            srr_all <- run$value
            corrected_expr[[sprintf("SRR (%s)", reference)]] <- srr_all
            for (method_name in names(standard_methods)) {
                combined_name <- sprintf("SRR+%s", method_name)
                run <- run_method(combined_name, reference, srr_all, metadata, standard_methods[[method_name]])
                status_rows[[length(status_rows) + 1]] <- run$status
                if (!is.null(run$value)) corrected_expr[[sprintf("%s (%s)", combined_name, reference)]] <- run$value
            }
        }
    }

    method_status <- do.call(rbind, status_rows)
    rownames(method_status) <- NULL
    if (any(method_status$Status == "failed")) {
        warning("Some batch-correction methods failed; inspect result$status.", call. = FALSE)
    }
    list(expr = corrected_expr, status = method_status)
}

# Subtract the arithmetic mean of log2 reference intensities within each batch from all samples in that batch.
srr_correction <- function(expr, metadata, ref_type = "P") {
    expr <- as.matrix(expr)
    cols <- metadata$ColName
    if (!all(cols %in% colnames(expr))) stop("Some metadata samples are absent from expr.", call. = FALSE)
    metadata <- metadata[match(cols, metadata$ColName), , drop = FALSE]
    if (anyNA(metadata$Batch)) stop("Batch labels cannot be missing.", call. = FALSE)

    corrected <- matrix(NA_real_, nrow(expr), length(cols), dimnames = list(rownames(expr), cols))
    batch_labels <- as.character(metadata$Batch)
    for (batch in unique(batch_labels)) {
        batch_metadata <- metadata[batch_labels == batch, , drop = FALSE]
        batch_cols <- batch_metadata$ColName
        reference_cols <- batch_metadata$ColName[batch_metadata$Sample == ref_type]
        if (!length(reference_cols)) {
            stop(sprintf("Batch %s has no %s reference samples.", batch, ref_type), call. = FALSE)
        }

        reference_mean <- rowMeans(expr[, reference_cols, drop = FALSE], na.rm = TRUE)
        reference_mean[!is.finite(reference_mean)] <- NA_real_
        corrected[, batch_cols] <- sweep(expr[, batch_cols, drop = FALSE], 1, reference_mean, `-`)
    }
    as.data.frame(corrected, check.names = FALSE)
}

# Estimate and apply a common batch MAD on the linear scale, then return to the original scale.
mad_scaling <- function(expr, metadata, log2_input = TRUE) {
    expr <- as.matrix(expr)
    cols <- metadata$ColName
    if (!all(cols %in% colnames(expr))) stop("Some metadata samples are absent from expr.", call. = FALSE)
    metadata <- metadata[match(cols, metadata$ColName), , drop = FALSE]
    if (anyNA(metadata$Batch)) stop("Batch labels cannot be missing.", call. = FALSE)

    # This is a generic batch-level comparator, not a reproduction of SomaScan normalization.
    linear_expr <- if (log2_input) 2^expr[, cols, drop = FALSE] else expr[, cols, drop = FALSE]
    batch_labels <- as.character(metadata$Batch)
    batches <- unique(batch_labels)
    batch_mads <- vapply(batches, function(batch) {
        batch_cols <- metadata$ColName[batch_labels == batch]
        stats::mad(as.vector(linear_expr[, batch_cols, drop = FALSE]), na.rm = TRUE)
    }, numeric(1))
    if (any(!is.finite(batch_mads) | batch_mads <= 0)) {
        stop("Every batch must have a positive finite MAD.", call. = FALSE)
    }

    target_mad <- stats::median(batch_mads)
    for (batch in batches) {
        batch_cols <- metadata$ColName[batch_labels == batch]
        linear_expr[, batch_cols] <- linear_expr[, batch_cols, drop = FALSE] * (target_mad / batch_mads[batch])
    }
    corrected <- if (log2_input) log2(linear_expr) else linear_expr
    as.data.frame(corrected, check.names = FALSE)
}

# Force all sample distributions to share the same empirical quantiles.
quantile_correction <- function(expr, metadata) {
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("Package 'limma' is required for quantile correction.", call. = FALSE)
    }
    cols <- metadata$ColName
    if (!all(cols %in% colnames(expr))) stop("Some metadata samples are absent from expr.", call. = FALSE)

    corrected <- limma::normalizeQuantiles(as.matrix(expr[, cols, drop = FALSE]))
    rownames(corrected) <- rownames(expr)
    colnames(corrected) <- cols
    as.data.frame(corrected, check.names = FALSE)
}

# Estimate batch shifts with empirical Bayes; constant proteins are removed because ComBat cannot fit them.
combat_correction <- function(expr, metadata, use_model_matrix = FALSE) {
    if (!requireNamespace("sva", quietly = TRUE)) stop("Package 'sva' is required for ComBat.", call. = FALSE)
    if (!requireNamespace("matrixStats", quietly = TRUE)) {
        stop("Package 'matrixStats' is required for ComBat.", call. = FALSE)
    }

    cols <- metadata$ColName
    if (!all(cols %in% colnames(expr))) stop("Some metadata samples are absent from expr.", call. = FALSE)
    metadata <- metadata[match(cols, metadata$ColName), , drop = FALSE]
    expr <- as.matrix(expr[, cols, drop = FALSE])
    if (anyNA(metadata$Batch)) stop("Batch labels cannot be missing.", call. = FALSE)
    if (anyNA(expr) || any(!is.finite(expr))) stop("ComBat requires a complete finite matrix.", call. = FALSE)
    if (length(unique(metadata$Batch)) < 2) stop("ComBat requires at least two batches.", call. = FALSE)

    row_variance <- matrixStats::rowVars(expr)
    expr <- expr[is.finite(row_variance) & row_variance > 0, , drop = FALSE]
    if (!nrow(expr)) stop("No variable proteins remain for ComBat.", call. = FALSE)
    model <- if (use_model_matrix) stats::model.matrix(~factor(Sample), data = metadata) else NULL

    corrected <- sva::ComBat(
        dat = expr, batch = factor(metadata$Batch), mod = model,
        par.prior = TRUE, prior.plots = FALSE
    )
    as.data.frame(corrected, check.names = FALSE)
}

# Select stable proteins on the linear scale, then fit RUV-III-C to the log2 expression matrix.
# Study-sample identities are not linked across batches; only within-batch technical replicates are supplied.
ruv_correction <- function(expr, metadata, k = 2, n_controls = 50) {
    if (!requireNamespace("RUVIIIC", quietly = TRUE)) {
        stop("Package 'RUVIIIC' is required for RUV-III-C.", call. = FALSE)
    }
    if (!requireNamespace("matrixStats", quietly = TRUE)) {
        stop("Package 'matrixStats' is required for RUV-III-C.", call. = FALSE)
    }

    cols <- metadata$ColName
    if (!all(cols %in% colnames(expr))) stop("Some metadata samples are absent from expr.", call. = FALSE)
    metadata <- metadata[match(cols, metadata$ColName), , drop = FALSE]
    expr <- as.matrix(expr[, cols, drop = FALSE])
    if (anyNA(metadata$Sample)) stop("Sample labels cannot be missing.", call. = FALSE)
    if (anyNA(metadata$Batch)) stop("Batch labels cannot be missing.", call. = FALSE)

    linear_expr <- 2^expr
    feature_cvs <- matrixStats::rowSds(linear_expr, na.rm = TRUE) / rowMeans(linear_expr, na.rm = TRUE)
    complete_controls <- rowSums(is.finite(linear_expr)) == ncol(linear_expr)
    candidates <- which(is.finite(feature_cvs) & complete_controls)
    if (!length(candidates)) stop("No complete finite control proteins are available for RUV-III-C.", call. = FALSE)
    control_indices <- candidates[order(feature_cvs[candidates])][seq_len(min(n_controls, length(candidates)))]
    control_proteins <- rownames(expr)[control_indices]

    # Technical replicates are linked only within the same batch; cross-batch study-sample identity is withheld.
    replicate_factor <- interaction(metadata$Batch, metadata$Sample, drop = TRUE, lex.order = TRUE)
    mapping_matrix <- stats::model.matrix(~replicate_factor - 1)
    colnames(mapping_matrix) <- levels(replicate_factor)
    replication_rank <- nrow(mapping_matrix) - ncol(mapping_matrix)
    actual_k <- min(k, length(control_proteins), replication_rank)
    if (actual_k < 1) stop("RUV-III-C requires at least one replicated sample group.", call. = FALSE)

    corrected <- RUVIIIC::RUVIII_C(
        k = actual_k, Y = t(expr), M = mapping_matrix, toCorrect = rownames(expr),
        controls = control_proteins, progress = FALSE
    )
    corrected <- t(corrected)
    rownames(corrected) <- rownames(expr)
    colnames(corrected) <- cols
    as.data.frame(corrected, check.names = FALSE)
}

# Use the most stable 25% of proteins as empirical negative controls for RUVg.
ruvg_correction <- function(expr, metadata, k = 2) {
    if (!requireNamespace("RUVSeq", quietly = TRUE)) stop("Package 'RUVSeq' is required for RUVg.", call. = FALSE)
    if (!requireNamespace("matrixStats", quietly = TRUE)) {
        stop("Package 'matrixStats' is required for RUVg.", call. = FALSE)
    }

    cols <- metadata$ColName
    if (!all(cols %in% colnames(expr))) stop("Some metadata samples are absent from expr.", call. = FALSE)
    expr <- as.matrix(expr[, cols, drop = FALSE])
    if (anyNA(expr) || any(!is.finite(expr))) stop("RUVg requires a complete finite matrix.", call. = FALSE)

    linear_expr <- 2^expr
    feature_cvs <- matrixStats::rowSds(linear_expr) / rowMeans(linear_expr)
    candidates <- which(is.finite(feature_cvs))
    if (!length(candidates)) stop("No finite control proteins are available for RUVg.", call. = FALSE)
    n_controls <- max(1, round(length(candidates) * 0.25))
    control_indices <- candidates[order(feature_cvs[candidates])][seq_len(n_controls)]
    actual_k <- min(k, length(control_indices) - 1, ncol(expr) - 1)
    if (actual_k < 1) stop("Insufficient dimensions for RUVg.", call. = FALSE)

    result <- RUVSeq::RUVg(x = expr, cIdx = control_indices, k = actual_k, isLog = TRUE)
    corrected <- result$normalizedCounts
    rownames(corrected) <- rownames(expr)
    colnames(corrected) <- cols
    as.data.frame(corrected, check.names = FALSE)
}

# Use repeated P or N measurements as the replicate group for RUVs; study samples remain independent.
ruvs_correction <- function(expr, metadata, k = 2, ref_type = c("P", "N")) {
    if (!requireNamespace("RUVSeq", quietly = TRUE)) stop("Package 'RUVSeq' is required for RUVs.", call. = FALSE)

    cols <- metadata$ColName
    if (!all(cols %in% colnames(expr))) stop("Some metadata samples are absent from expr.", call. = FALSE)
    metadata <- metadata[match(cols, metadata$ColName), , drop = FALSE]
    expr <- as.matrix(expr[, cols, drop = FALSE])
    if (anyNA(metadata$Sample)) stop("Sample labels cannot be missing.", call. = FALSE)
    if (anyNA(expr) || any(!is.finite(expr))) stop("RUVs requires a complete finite matrix.", call. = FALSE)

    group_labels <- as.character(metadata$Sample)
    is_reference <- group_labels %in% ref_type
    if (sum(is_reference) < 2) stop("RUVs requires at least two reference samples.", call. = FALSE)
    group_labels[!is_reference] <- paste0("Study_Sample_", seq_len(sum(!is_reference)))
    replicate_groups <- RUVSeq::makeGroups(group_labels)
    control_indices <- seq_len(nrow(expr))
    actual_k <- min(k, length(control_indices) - 1, sum(is_reference) - 1)
    if (actual_k < 1) stop("Insufficient dimensions for RUVs.", call. = FALSE)

    result <- RUVSeq::RUVs(
        x = expr, cIdx = control_indices, k = actual_k,
        scIdx = replicate_groups, isLog = TRUE
    )
    corrected <- result$normalizedCounts
    rownames(corrected) <- rownames(expr)
    colnames(corrected) <- cols
    as.data.frame(corrected, check.names = FALSE)
}
