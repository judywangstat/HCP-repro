#' Compare density-score and residual-score HCP regions
#'
#' @description
#' Constructs two HCP prediction regions under the same data splits,
#' calibration subsamples, propensity weights, and candidate response grid.
#' The first region uses the density score
#' \eqn{R(x,y) = -\hat f(y \mid x)}, while the second uses the residual score
#' \eqn{R(x,y) = |y - \hat\mu(x)|}.
#'
#' This function is intended for reproducing Table S.5 and Table S.6, where the
#' goal is to compare the proposed density-based nonconformity score with a
#' residual-based score within the same HCP framework.
#'
#' @param dat A data frame containing clustered observations.
#' @param id_col Name of the subject or cluster identifier column.
#' @param y_col Name of the outcome column.
#' @param delta_col Name of the missingness indicator column, where 1 indicates
#'   an observed outcome and 0 indicates a missing outcome.
#' @param x_cols Character vector giving the covariate column names.
#' @param x_test New covariate value(s), either a numeric vector, matrix, or data frame.
#' @param y_grid Numeric vector of candidate response values.
#' @param alpha Miscoverage level.
#' @param train_frac Fraction of subjects assigned to training in each split.
#' @param S Number of repeated subject-level splits.
#' @param B Number of calibration subsamples per split.
#' @param combine_B Method for combining p-values over calibration subsamples.
#' @param combine_S Method for combining p-values over subject-level splits.
#' @param seed Optional random seed.
#' @param return_details Logical; if TRUE, return split-level p-values.
#' @param dens_method Method used by fit_cond_density_qp().
#' @param dens_taus Quantile grid used by fit_cond_density_qp().
#' @param dens_h Optional bandwidth passed to fit_cond_density_qp().
#' @param dens_h_scenario Optional scenario name passed to quantile_levels()
#'   when dens_h is NULL.
#' @param enforce_monotone Logical argument passed to fit_cond_density_qp().
#' @param tail_decay Logical argument passed to fit_cond_density_qp().
#' @param reg_method Regression engine for the residual score.
#' @param prop_method Method used by fit_propensity_model().
#' @param prop_eps Clipping constant for propensity scores.
#' @param weight_cap Upper bound for inverse-propensity weights.
#' @param ... Additional arguments passed to fit_propensity_model().
#'
#' @return A list with two main components:
#' \describe{
#'   \item{density}{HCP result using the density score.}
#'   \item{residual}{HCP result using the residual score.}
#' }
#'
#' Each component contains region, lo_hi, p_final, and y_grid.
#'
#' @export
hcp_score_compare <- function(
    dat,
    id_col,
    y_col = "Y",
    delta_col = "delta",
    x_cols,
    x_test,
    y_grid,
    alpha = 0.1,
    train_frac = 0.5,
    S = 5,
    B = 5,
    combine_B = c("cct", "mean"),
    combine_S = c("cct", "mean"),
    seed = NULL,
    return_details = FALSE,
    dens_method = c("rq", "qrf"),
    dens_taus = seq(0.05, 0.95, by = 0.02),
    dens_h = NULL,
    dens_h_scenario = NULL,
    enforce_monotone = FALSE,
    tail_decay = TRUE,
    reg_method = c("linear", "grf"),
    prop_method = c("logistic", "grf", "boosting"),
    prop_eps = 1e-6,
    weight_cap = Inf,
    ...
) {
  
  # ---------------------------------------------------------------------------
  # Step 0: Match arguments and perform basic input checks
  # ---------------------------------------------------------------------------
  combine_B <- match.arg(combine_B)
  combine_S <- match.arg(combine_S)
  dens_method <- match.arg(dens_method)
  reg_method <- match.arg(reg_method)
  prop_method <- match.arg(prop_method)
  
  if (!is.data.frame(dat)) {
    stop("dat must be a data.frame.")
  }
  
  if (!all(c(id_col, y_col, delta_col) %in% names(dat))) {
    stop("id_col, y_col, and delta_col must all be columns in dat.")
  }
  
  if (!is.character(x_cols) || length(x_cols) < 1L) {
    stop("x_cols must be a non-empty character vector.")
  }
  
  if (anyDuplicated(x_cols)) {
    stop("x_cols contains duplicated column names.")
  }
  
  if (!all(x_cols %in% names(dat))) {
    stop("Some columns specified in x_cols are missing from dat.")
  }
  
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number in (0,1).")
  }
  
  if (!is.numeric(train_frac) || length(train_frac) != 1L ||
      !is.finite(train_frac) || train_frac <= 0 || train_frac >= 1) {
    stop("train_frac must be a single number in (0,1).")
  }
  
  if (!is.numeric(S) || length(S) != 1L || !is.finite(S) ||
      S < 1 || S != floor(S)) {
    stop("S must be a positive integer.")
  }
  
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) ||
      B < 1 || B != floor(B)) {
    stop("B must be a positive integer.")
  }
  
  S <- as.integer(S)
  B <- as.integer(B)
  
  if (!is.logical(return_details) || length(return_details) != 1L) {
    stop("return_details must be TRUE or FALSE.")
  }
  
  if (!is.numeric(prop_eps) || length(prop_eps) != 1L || !is.finite(prop_eps) ||
      prop_eps <= 0 || prop_eps >= 0.5) {
    stop("prop_eps must be a single finite number in (0,0.5).")
  }
  
  if (!is.numeric(weight_cap) || length(weight_cap) != 1L ||
      is.na(weight_cap) || weight_cap < 1) {
    stop("weight_cap must be a single numeric value >= 1, or Inf.")
  }
  
  dens_taus <- sort(unique(as.numeric(dens_taus)))
  if (anyNA(dens_taus) || any(dens_taus <= 0 | dens_taus >= 1)) {
    stop("dens_taus must be a numeric vector with values strictly inside (0,1).")
  }
  
  if (!is.null(dens_h)) {
    dens_h <- as.numeric(dens_h)
    if (anyNA(dens_h) || any(!is.finite(dens_h)) || any(dens_h <= 0)) {
      stop("If provided, dens_h must contain only positive finite numeric values.")
    }
    if (!(length(dens_h) %in% c(1L, length(dens_taus)))) {
      stop("If provided, dens_h must have length 1 or length(dens_taus).")
    }
  }
  
  y_grid <- as.numeric(y_grid)
  if (anyNA(y_grid) || any(!is.finite(y_grid))) {
    stop("y_grid must contain only finite numeric values.")
  }
  y_grid <- sort(unique(y_grid))
  if (length(y_grid) < 2L) {
    stop("y_grid must contain at least two distinct values.")
  }
  Ny <- length(y_grid)
  
  # ---------------------------------------------------------------------------
  # Step 1: Standardize x_test input
  # ---------------------------------------------------------------------------
  if (is.vector(x_test) && !is.list(x_test)) {
    x_test <- matrix(as.numeric(x_test), nrow = 1L)
  } else if (is.data.frame(x_test)) {
    x_test <- as.matrix(x_test)
  } else {
    x_test <- as.matrix(x_test)
  }
  
  if (ncol(x_test) != length(x_cols)) {
    stop("x_test must have length(x_cols) columns.")
  }
  
  storage.mode(x_test) <- "numeric"
  colnames(x_test) <- x_cols
  
  if (anyNA(x_test) || any(!is.finite(x_test))) {
    stop("x_test contains NA or non-finite values after coercion.")
  }
  
  K <- nrow(x_test)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # ---------------------------------------------------------------------------
  # Helper: Robust binary checks for missingness indicators
  # ---------------------------------------------------------------------------
  is_delta1 <- function(d_raw) {
    d_chr <- toupper(trimws(as.character(d_raw)))
    d_num <- suppressWarnings(as.numeric(d_chr))
    (!is.na(d_num) & d_num == 1) | d_chr %in% c("1", "TRUE", "T")
  }
  
  is_delta0 <- function(d_raw) {
    d_chr <- toupper(trimws(as.character(d_raw)))
    d_num <- suppressWarnings(as.numeric(d_chr))
    (!is.na(d_num) & d_num == 0) | d_chr %in% c("0", "FALSE", "F")
  }
  
  # ---------------------------------------------------------------------------
  # Helper: Combine p-values across the first dimension of an array
  # ---------------------------------------------------------------------------
  combine_array_firstdim <- function(p_arr, method) {
    if (method == "mean") {
      return(apply(p_arr, c(2, 3), function(v) {
        m <- mean(v, na.rm = TRUE)
        if (!is.finite(m)) return(1)
        pmin(pmax(m, 0), 1)
      }))
    }
    
    p_arr <- pmin(pmax(p_arr, 1e-15), 1 - 1e-15)
    t_arr <- tan((0.5 - p_arr) * pi)
    tbar <- apply(t_arr, c(2, 3), mean, na.rm = TRUE)
    
    if (any(!is.finite(tbar))) {
      tbar[is.nan(tbar)] <- 0
      idx_inf <- is.infinite(tbar)
      tbar[idx_inf] <- sign(tbar[idx_inf]) * .Machine$double.xmax
    }
    
    p <- 0.5 - atan(tbar) / pi
    pmin(pmax(p, 0), 1)
  }
  
  # ---------------------------------------------------------------------------
  # Helper: Split data by subject into training and calibration sets
  # ---------------------------------------------------------------------------
  split_by_subject <- function(dat0, seed0 = NULL) {
    if (!is.null(seed0)) {
      set.seed(seed0)
    }
    
    ids <- unique(dat0[[id_col]])
    n_id <- length(ids)
    
    if (n_id < 2L) {
      stop("Need at least 2 unique subjects for splitting.")
    }
    
    n_tr <- floor(train_frac * n_id)
    n_tr <- max(1L, min(n_tr, n_id - 1L))
    
    tr_ids <- sample(ids, size = n_tr, replace = FALSE)
    idx_tr <- dat0[[id_col]] %in% tr_ids
    
    list(
      train = dat0[idx_tr, , drop = FALSE],
      calib = dat0[!idx_tr, , drop = FALSE]
    )
  }
  
  # ---------------------------------------------------------------------------
  # Helper: Draw one observation per calibration subject
  # ---------------------------------------------------------------------------
  subsample_one_per_subject <- function(dat_cal, seed0 = NULL) {
    if (!is.null(seed0)) {
      set.seed(seed0)
    }
    
    ids <- unique(dat_cal[[id_col]])
    pick <- lapply(ids, function(id) {
      rows <- which(dat_cal[[id_col]] == id)
      sample(rows, size = 1L)
    })
    
    dat_cal[unlist(pick), , drop = FALSE]
  }
  
  # ---------------------------------------------------------------------------
  # Helper: Determine split-specific bandwidth for density estimation
  # ---------------------------------------------------------------------------
  resolve_dens_h_split <- function(dat_tr_split) {
    if (!is.null(dens_h)) {
      return(dens_h)
    }
    
    obs_idx <- is_delta1(dat_tr_split[[delta_col]]) &
      !is.na(dat_tr_split[[y_col]])
    n_obs_tr <- sum(obs_idx)
    
    if (n_obs_tr < 2L) {
      stop("Too few observed training outcomes in this split to recompute dens_h.")
    }
    
    quantile_levels(
      taus = dens_taus,
      n = n_obs_tr,
      scenario = dens_h_scenario
    )$h
  }
  
  # ---------------------------------------------------------------------------
  # Helper: Fit regression model for residual score
  # ---------------------------------------------------------------------------
  fit_regression_model <- function(dat_tr_obs) {
    X_obs <- as.matrix(dat_tr_obs[, x_cols, drop = FALSE])
    Y_obs <- as.numeric(dat_tr_obs[[y_col]])
    
    if (reg_method == "linear") {
      lm_dat <- data.frame(Y_obs = Y_obs, as.data.frame(X_obs))
      fit <- stats::lm(Y_obs ~ ., data = lm_dat)
      
      predict_fun <- function(x_new) {
        x_new <- as.data.frame(as.matrix(x_new))
        colnames(x_new) <- x_cols
        as.numeric(stats::predict(fit, newdata = x_new))
      }
      
      return(list(fit_object = fit, predict = predict_fun))
    }
    
    if (reg_method == "grf") {
      if (!requireNamespace("grf", quietly = TRUE)) {
        stop("Package 'grf' is required for reg_method = 'grf'.")
      }
      
      fit <- grf::regression_forest(X = X_obs, Y = Y_obs)
      
      predict_fun <- function(x_new) {
        x_new <- as.matrix(x_new)
        storage.mode(x_new) <- "numeric"
        as.numeric(stats::predict(fit, newdata = x_new)$predictions)
      }
      
      return(list(fit_object = fit, predict = predict_fun))
    }
    
    stop("Unsupported reg_method.")
  }
  
  # ---------------------------------------------------------------------------
  # Step 2: Initialize storage for split-level p-values
  # ---------------------------------------------------------------------------
  p_density_split <- array(NA_real_, dim = c(S, K, Ny))
  p_residual_split <- array(NA_real_, dim = c(S, K, Ny))
  
  # ---------------------------------------------------------------------------
  # Step 3: Repeat over subject-level splits
  # ---------------------------------------------------------------------------
  for (s in seq_len(S)) {
    
    seed_s <- if (!is.null(seed)) seed + 1000L * s else NULL
    sp <- split_by_subject(dat, seed0 = seed_s)
    
    dat_tr <- sp$train
    dat_ca <- sp$calib
    
    obs_tr <- is_delta1(dat_tr[[delta_col]]) & !is.na(dat_tr[[y_col]])
    dat_tr_obs <- dat_tr[obs_tr, , drop = FALSE]
    
    if (nrow(dat_tr_obs) < 10L) {
      stop("Too few observed training outcomes in this split.")
    }
    
    # -------------------------------------------------------------------------
    # Step 3A: Fit conditional density model
    # -------------------------------------------------------------------------
    dens_h_use <- resolve_dens_h_split(dat_tr)
    
    dens_fit <- fit_cond_density_qp(
      dat = dat_tr,
      y_col = y_col,
      delta_col = delta_col,
      x_cols = x_cols,
      taus = dens_taus,
      h = dens_h_use,
      method = dens_method,
      enforce_monotone = enforce_monotone,
      tail_decay = tail_decay
    )
    
    # -------------------------------------------------------------------------
    # Step 3B: Fit regression model for residual score
    # -------------------------------------------------------------------------
    reg_fit <- fit_regression_model(dat_tr_obs)
    
    # -------------------------------------------------------------------------
    # Step 3C: Fit propensity model unless all training outcomes are observed
    # -------------------------------------------------------------------------
    dtr <- dat_tr[[delta_col]]
    fully_observed_split <- all(is_delta1(dtr)) && !any(is_delta0(dtr))
    
    prop_fit <- NULL
    if (!fully_observed_split) {
      prop_fit <- fit_propensity_model(
        dat = dat_tr,
        delta_col = delta_col,
        x_cols = x_cols,
        method = prop_method,
        eps = prop_eps,
        ...
      )
    }
    
    # -------------------------------------------------------------------------
    # Step 3D: Compute test-side weights
    # -------------------------------------------------------------------------
    if (fully_observed_split) {
      w_test <- rep(1, K)
    } else {
      x_test_df <- as.data.frame(x_test)
      colnames(x_test_df) <- x_cols
      
      p_test <- as.numeric(prop_fit$predict(x_test_df))
      if (any(!is.finite(p_test))) {
        stop("predict returned non-finite propensity values for x_test.")
      }
      
      p_test <- pmax(pmin(p_test, 1 - prop_eps), prop_eps)
      w_test <- pmin(1 / p_test, weight_cap)
    }
    
    # -------------------------------------------------------------------------
    # Step 3E: Compute test-side scores over the full y-grid
    # -------------------------------------------------------------------------
    dens_te_mat <- dens_fit$predict_density_grid(x_test, y_grid)
    if (any(!is.finite(dens_te_mat))) {
      stop("Density prediction returned non-finite values for x_test.")
    }
    R_density_test <- -dens_te_mat
    
    mu_test <- reg_fit$predict(x_test)
    R_residual_test <- matrix(NA_real_, nrow = K, ncol = Ny)
    for (k in seq_len(K)) {
      R_residual_test[k, ] <- abs(y_grid - mu_test[k])
    }
    
    # -------------------------------------------------------------------------
    # Step 3F: Repeated calibration subsampling
    # -------------------------------------------------------------------------
    p_density_bky <- array(NA_real_, dim = c(B, K, Ny))
    p_residual_bky <- array(NA_real_, dim = c(B, K, Ny))
    
    for (b in seq_len(B)) {
      seed_b <- if (!is.null(seed_s)) seed_s + 10000L + b else NULL
      dat_sub <- subsample_one_per_subject(dat_ca, seed0 = seed_b)
      
      obs_sub <- is_delta1(dat_sub[[delta_col]]) & !is.na(dat_sub[[y_col]])
      dat_obs <- dat_sub[obs_sub, , drop = FALSE]
      
      if (nrow(dat_obs) == 0L) {
        next
      }
      
      Xc_df <- as.data.frame(as.matrix(dat_obs[, x_cols, drop = FALSE]))
      colnames(Xc_df) <- x_cols
      Yc <- as.numeric(dat_obs[[y_col]])
      
      # Calibration weights
      if (fully_observed_split) {
        w_c <- rep(1, nrow(Xc_df))
      } else {
        p_c <- as.numeric(prop_fit$predict(Xc_df))
        if (any(!is.finite(p_c))) {
          stop("predict returned non-finite propensity values in calibration.")
        }
        
        p_c <- pmax(pmin(p_c, 1 - prop_eps), prop_eps)
        w_c <- pmin(1 / p_c, weight_cap)
      }
      
      # Density score on calibration observations
      dens_ca_mat <- dens_fit$predict_density_grid(Xc_df, Yc)
      R_density_c <- -diag(dens_ca_mat)
      
      if (any(!is.finite(R_density_c))) {
        stop("Density prediction returned non-finite calibration scores.")
      }
      
      # Residual score on calibration observations
      mu_c <- reg_fit$predict(Xc_df)
      R_residual_c <- abs(Yc - mu_c)
      
      if (any(!is.finite(R_residual_c))) {
        stop("Regression prediction returned non-finite calibration scores.")
      }
      
      # Sort scores once for fast weighted p-value calculation
      ord_d <- order(R_density_c)
      R_density_sorted <- R_density_c[ord_d]
      w_density_sorted <- w_c[ord_d]
      W_density_ge <- rev(cumsum(rev(w_density_sorted)))
      sum_w_density <- sum(w_c)
      
      ord_r <- order(R_residual_c)
      R_residual_sorted <- R_residual_c[ord_r]
      w_residual_sorted <- w_c[ord_r]
      W_residual_ge <- rev(cumsum(rev(w_residual_sorted)))
      sum_w_residual <- sum(w_c)
      
      for (k in seq_len(K)) {
        # Density-score p-values
        denom_d <- sum_w_density + w_test[k]
        rt_d <- R_density_test[k, ]
        j_d <- findInterval(rt_d, R_density_sorted, left.open = TRUE) + 1L
        w_ge_d <- ifelse(j_d <= length(W_density_ge), W_density_ge[j_d], 0)
        p_density_bky[b, k, ] <- (w_ge_d + w_test[k]) / denom_d
        
        # Residual-score p-values
        denom_r <- sum_w_residual + w_test[k]
        rt_r <- R_residual_test[k, ]
        j_r <- findInterval(rt_r, R_residual_sorted, left.open = TRUE) + 1L
        w_ge_r <- ifelse(j_r <= length(W_residual_ge), W_residual_ge[j_r], 0)
        p_residual_bky[b, k, ] <- (w_ge_r + w_test[k]) / denom_r
      }
    }
    
    # -------------------------------------------------------------------------
    # Step 3G: Combine p-values across B calibration subsamples
    # -------------------------------------------------------------------------
    p_density_split[s, , ] <- combine_array_firstdim(p_density_bky, combine_B)
    p_residual_split[s, , ] <- combine_array_firstdim(p_residual_bky, combine_B)
  }
  
  # ---------------------------------------------------------------------------
  # Step 4: Combine p-values across S subject-level splits
  # ---------------------------------------------------------------------------
  p_density_final <- combine_array_firstdim(p_density_split, combine_S)
  p_residual_final <- combine_array_firstdim(p_residual_split, combine_S)
  
  # ---------------------------------------------------------------------------
  # Step 5: Convert final p-values into grid-based prediction regions
  # ---------------------------------------------------------------------------
  build_region_output <- function(p_final) {
    regions <- vector("list", K)
    lo_hi <- matrix(NA_real_, nrow = K, ncol = 2)
    colnames(lo_hi) <- c("lo", "hi")
    
    for (k in seq_len(K)) {
      reg <- y_grid[p_final[k, ] >= alpha]
      regions[[k]] <- reg
      
      if (length(reg) == 0L) {
        lo_hi[k, ] <- c(NA_real_, NA_real_)
      } else {
        lo_hi[k, ] <- c(min(reg), max(reg))
      }
    }
    
    list(
      region = regions,
      lo_hi = lo_hi,
      p_final = p_final,
      y_grid = y_grid
    )
  }
  
  out_density <- build_region_output(p_density_final)
  out_residual <- build_region_output(p_residual_final)
  
  out <- list(
    density = out_density,
    residual = out_residual
  )
  
  if (isTRUE(return_details)) {
    out$density$p_split <- p_density_split
    out$residual$p_split <- p_residual_split
  }
  
  out
}