check_imputation_input <- function(expr, meta, require_batch = TRUE) {
    if (is.null(colnames(expr)) || anyNA(colnames(expr)) || any(!nzchar(colnames(expr)))) {
        stop("expr must have non-missing sample names in colnames(expr).", call. = FALSE)
    }

    required <- c("ColName", if (require_batch) "Batch")
    missing_fields <- setdiff(required, colnames(meta))
    if (length(missing_fields)) {
        fields <- paste(missing_fields, collapse = ", ")
        stop(sprintf("meta is missing required columns: %s.", fields), call. = FALSE)
    }

    expr_names <- as.character(colnames(expr))
    meta_names <- as.character(meta$ColName)
    if (anyDuplicated(expr_names)) stop("expr contains duplicated sample names.", call. = FALSE)
    if (anyDuplicated(meta_names)) stop("meta$ColName contains duplicated sample names.", call. = FALSE)

    unmatched <- setdiff(expr_names, meta_names)
    if (length(unmatched)) {
        shown <- paste(utils::head(unmatched, 8), collapse = ", ")
        suffix <- if (length(unmatched) > 8) ", ..." else ""
        stop(sprintf("Expression samples are missing from meta$ColName: %s%s.", shown, suffix), call. = FALSE)
    }

    dat <- as.matrix(expr[, expr_names, drop = FALSE])
    if (!is.numeric(dat)) stop("expr must contain only numeric values.", call. = FALSE)
    storage.mode(dat) <- "double"
    if (any(is.infinite(dat))) stop("expr contains Inf or -Inf values.", call. = FALSE)

    aligned_meta <- meta[match(expr_names, meta_names), , drop = FALSE]
    if (require_batch && anyNA(aligned_meta$Batch)) {
        stop("meta$Batch contains missing values for expression samples.", call. = FALSE)
    }
    list(data = dat, meta = aligned_meta)
}

with_local_seed <- function(seed, code) {
    if (is.null(seed)) return(force(code))

    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
        if (had_seed) {
            assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
        }
    }, add = TRUE)

    set.seed(seed)
    force(code)
}

calc_horizon_parameters <- function(values, quantile_prob = 0.001, noise_width = 0.05) {
    values <- values[is.finite(values)]
    if (length(values) < 6) return(NULL)

    horizon <- as.numeric(stats::quantile(values, quantile_prob, names = FALSE, type = 7))
    noise_sd <- abs(horizon) * noise_width

    # Retain the original half-minimum bound for non-negative intensities.
    lower_bound <- if (min(values) >= 0) min(values) * 0.5 else horizon - 3 * max(noise_sd, .Machine$double.eps)
    list(center = horizon, sd = noise_sd, lower = lower_bound)
}

impute_knn_core <- function(dat, k = 10, rowmax = 0.9, seed = 362436069) {
    if (!requireNamespace("impute", quietly = TRUE)) {
        stop("Package 'impute' is required but not installed.", call. = FALSE)
    }
    if (!is.numeric(k) || length(k) != 1 || is.na(k) || k < 1) stop("k must be positive.", call. = FALSE)
    if (!is.numeric(rowmax) || length(rowmax) != 1 || is.na(rowmax) || rowmax <= 0 || rowmax > 1) {
        stop("rowmax must be in the interval (0, 1].", call. = FALSE)
    }

    missing_before <- is.na(dat)
    if (!any(missing_before)) return(list(data = dat, imputed = missing_before))
    if (any(colSums(!missing_before) == 0)) {
        stop("KNN cannot impute a sample column with no observed values.", call. = FALSE)
    }

    # impute.knn searches for neighbours among feature rows, not sample columns.
    k_value <- min(as.integer(k), nrow(dat) - 1L)
    if (k_value < 1) {
        warning("Too few feature rows for KNN imputation; residual NAs were retained.", call. = FALSE)
        empty_mask <- matrix(FALSE, nrow(dat), ncol(dat), dimnames = dimnames(dat))
        return(list(data = dat, imputed = empty_mask))
    }

    result <- impute::impute.knn(dat, k = k_value, rowmax = rowmax, colmax = 1,
                                 rng.seed = as.integer(seed))$data
    list(data = result, imputed = missing_before & is.finite(result))
}

# Impute all missing values at the batch-specific 0.1% intensity horizon.
impute_lod_noise <- function(expr, meta, noise_width = 0.05, quantile_prob = 0.001,
                             seed = 42, return_details = FALSE) {
    if (noise_width < 0) stop("noise_width must be non-negative.", call. = FALSE)
    if (quantile_prob < 0 || quantile_prob > 1) stop("quantile_prob must be between 0 and 1.", call. = FALSE)

    checked <- check_imputation_input(expr, meta)
    imputed <- checked$data
    batch_info <- as.character(checked$meta$Batch)
    mask <- matrix(FALSE, nrow(imputed), ncol(imputed), dimnames = dimnames(imputed))
    skipped <- character()

    with_local_seed(seed, {
        for (batch in unique(batch_info)) {
            cols <- which(batch_info == batch)
            batch_data <- imputed[, cols, drop = FALSE]
            params <- calc_horizon_parameters(as.vector(batch_data), quantile_prob, noise_width)
            if (is.null(params)) {
                skipped <- c(skipped, batch)
                next
            }

            missing <- is.na(batch_data)
            if (!any(missing)) next
            values <- stats::rnorm(sum(missing), mean = params$center, sd = params$sd)
            batch_data[missing] <- pmax(values, params$lower)
            imputed[, cols] <- batch_data

            batch_mask <- mask[, cols, drop = FALSE]
            batch_mask[missing] <- TRUE
            mask[, cols] <- batch_mask
        }
    })

    if (length(skipped)) {
        batches <- paste(skipped, collapse = ", ")
        warning(sprintf("Insufficient observed values in batch(es): %s.", batches), call. = FALSE)
    }
    result <- as.data.frame(imputed, check.names = FALSE)
    if (return_details) return(list(data = result, imputed = mask))
    result
}

# Impute sporadic missing values with global K-nearest neighbours.
impute_knn <- function(expr, meta, k = 10, rowmax = 0.9, seed = 362436069,
                       return_details = FALSE) {
    checked <- check_imputation_input(expr, meta, require_batch = FALSE)
    knn <- with_local_seed(seed, impute_knn_core(checked$data, k, rowmax, seed))
    result <- as.data.frame(knn$data, check.names = FALSE)
    if (return_details) return(list(data = result, imputed = knn$imputed))
    result
}

# Fill batch-wide missing rows at the intensity horizon, then apply global KNN.
impute_hybrid_global <- function(expr, meta, k = 10, noise_width = 0.05, quantile_prob = 0.001,
                                 rowmax = 0.9, seed = 42, knn_seed = 362436069,
                                 return_details = FALSE) {
    if (noise_width < 0) stop("noise_width must be non-negative.", call. = FALSE)
    if (quantile_prob < 0 || quantile_prob > 1) stop("quantile_prob must be between 0 and 1.", call. = FALSE)

    checked <- check_imputation_input(expr, meta)
    imputed <- checked$data
    batch_info <- as.character(checked$meta$Batch)
    horizon_mask <- matrix(FALSE, nrow(imputed), ncol(imputed), dimnames = dimnames(imputed))
    skipped <- character()
    knn <- NULL

    with_local_seed(seed, {
        for (batch in unique(batch_info)) {
            cols <- which(batch_info == batch)
            batch_data <- imputed[, cols, drop = FALSE]
            params <- calc_horizon_parameters(as.vector(batch_data), quantile_prob, noise_width)
            if (is.null(params)) {
                skipped <- c(skipped, batch)
                next
            }

            full_missing <- rowSums(is.na(batch_data)) == length(cols)
            if (!any(full_missing)) next
            values <- stats::rnorm(sum(full_missing) * length(cols), mean = params$center, sd = params$sd)
            values <- pmax(values, params$lower)
            imputed[full_missing, cols] <- matrix(values, nrow = sum(full_missing), ncol = length(cols))
            horizon_mask[full_missing, cols] <- TRUE
        }
        knn <- impute_knn_core(imputed, k, rowmax, knn_seed)
    })

    if (length(skipped)) {
        batches <- paste(skipped, collapse = ", ")
        warning(sprintf("Insufficient observed values in batch(es): %s.", batches), call. = FALSE)
    }
    result <- as.data.frame(knn$data, check.names = FALSE)
    if (return_details) {
        return(list(data = result, imputed = horizon_mask | knn$imputed,
                    horizon_imputed = horizon_mask, knn_imputed = knn$imputed))
    }
    result
}
