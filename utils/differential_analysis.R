# Differential protein analysis using limma.

dea_limma_flexible <- function(expr_mat, meta_mat, contrast_pair,
                               logfc_th = log2(1.2), fdr_th = 0.05,
                               min_samples_per_group = 3) {
    if (!requireNamespace("limma", quietly = TRUE)) stop("Package 'limma' is required.", call. = FALSE)

    required_meta <- c("ColName", "Sample")
    missing_meta <- setdiff(required_meta, colnames(meta_mat))
    if (length(missing_meta)) stop(sprintf("meta_mat is missing: %s.", paste(missing_meta, collapse = ", ")), call. = FALSE)

    groups <- strsplit(contrast_pair, "/", fixed = TRUE)[[1]]
    if (length(groups) != 2 || any(groups == "") || groups[[1]] == groups[[2]]) {
        stop("contrast_pair must contain two different groups in the format 'Group1/Group2'.", call. = FALSE)
    }
    group1 <- groups[[1]]
    group2 <- groups[[2]]
    if (!all(groups %in% meta_mat$Sample)) stop("One or both contrast groups are absent from meta_mat.", call. = FALSE)

    # Match samples strictly so that metadata errors cannot silently remove columns.
    expr_mat <- as.matrix(expr_mat)
    if (!is.numeric(expr_mat)) stop("expr_mat must contain only numeric values.", call. = FALSE)
    if (is.null(colnames(expr_mat)) || is.null(rownames(expr_mat))) {
        stop("expr_mat must have sample column names and feature row names.", call. = FALSE)
    }
    if (anyDuplicated(colnames(expr_mat)) || anyDuplicated(rownames(expr_mat))) {
        stop("expr_mat contains duplicated sample or feature names.", call. = FALSE)
    }
    if (anyDuplicated(meta_mat$ColName)) stop("meta_mat$ColName contains duplicated sample names.", call. = FALSE)

    selected_meta <- meta_mat[meta_mat$Sample %in% groups, , drop = FALSE]
    missing_expr <- setdiff(selected_meta$ColName, colnames(expr_mat))
    if (length(missing_expr)) {
        stop(sprintf("Samples in meta_mat are absent from expr_mat: %s.", paste(missing_expr, collapse = ", ")), call. = FALSE)
    }
    selected_expr <- expr_mat[, selected_meta$ColName, drop = FALSE]

    group1_cols <- selected_meta$ColName[selected_meta$Sample == group1]
    group2_cols <- selected_meta$ColName[selected_meta$Sample == group2]
    if (length(group1_cols) < min_samples_per_group || length(group2_cols) < min_samples_per_group) {
        stop("Both groups must contain at least min_samples_per_group samples.", call. = FALSE)
    }

    # Retain features with enough finite measurements in both comparison groups.
    keep <- apply(selected_expr, 1, function(values) {
        sum(is.finite(values[group1_cols])) >= min_samples_per_group &&
            sum(is.finite(values[group2_cols])) >= min_samples_per_group
    })
    selected_expr <- selected_expr[keep, , drop = FALSE]
    if (!nrow(selected_expr)) stop("No features satisfy the finite-value requirement.", call. = FALSE)
    selected_expr[!is.finite(selected_expr)] <- NA_real_

    # Fit Group1 - Group2 on the log2 measurement scale.
    group_factor <- factor(selected_meta$Sample, levels = groups)
    design <- stats::model.matrix(~0 + group_factor)
    safe_names <- make.names(groups, unique = TRUE)
    colnames(design) <- safe_names
    rownames(design) <- selected_meta$ColName
    contrast_matrix <- limma::makeContrasts(contrasts = sprintf("%s-%s", safe_names[[1]], safe_names[[2]]), levels = design)

    fit <- limma::lmFit(selected_expr, design)
    fit <- limma::contrasts.fit(fit, contrast_matrix)
    fit <- limma::eBayes(fit)
    result <- limma::topTable(fit, coef = 1, adjust.method = "BH", n = Inf, p.value = 1)
    result <- data.frame(UniqueID = rownames(result), result, row.names = NULL, check.names = FALSE)

    result$Pair <- contrast_pair
    result$Label <- "NotSig"
    result$Label[which(result$logFC >= logfc_th & result$adj.P.Val < fdr_th)] <- "Up"
    result$Label[which(result$logFC <= -logfc_th & result$adj.P.Val < fdr_th)] <- "Down"
    result$Significance <- "ns"
    result$Significance[which(result$adj.P.Val < 0.05)] <- "*"
    result$Significance[which(result$adj.P.Val < 0.01)] <- "**"
    result$Significance[which(result$adj.P.Val < 0.001)] <- "***"
    result
}
