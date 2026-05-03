#' Construct a cross-sectional LC prediction region with per-subject subsampling
#'
#' @description
#' This function implements an LC-style conformal prediction region for clustered
#' data with missing outcomes by first drawing one observation per subject to form
#' a cross-sectional sample, and then applying weighted conformal prediction on
#' the resulting subsampled data.
#'
#' To adapt the LC method to clustered data, the function first performs a single
#' subject-level subsampling step, producing one observation per subject. The
#' resulting cross-sectional data are then split into training and calibration
#' subsets. A conditional quantile model is fitted on the observed training data,
#' and a missingness propensity model is fitted on the full training data. Weighted
#' conformal \eqn{p}-values are computed on the observed calibration subset.
#'
#' @param dat A data frame containing clustered observations.
#' @param id_col Name of the subject or cluster identifier column.
#' @param y_col Name of the outcome column.
#' @param delta_col Name of the missingness indicator column, where 1 indicates
#'   an observed outcome and 0 indicates a missing outcome.
#' @param x_cols Character vector giving the names of the covariate columns.
#' @param x_test New covariate value(s). This can be a numeric vector
#'   (treated as one test point), or a matrix/data.frame with one row per test point
#'   and \code{length(x_cols)} columns.
#' @param y_grid Numeric vector of candidate response values at which conformal
#'   \eqn{p}-values are evaluated.
#' @param alpha Miscoverage level in \eqn{(0,1)}.
#' @param train_frac Fraction of subsampled subjects assigned to training.
#' @param seed Optional random seed.
#' @param quant_method Quantile model used to estimate lower and upper conditional
#'   quantiles. Must be \code{"linear"} or \code{"grf"}.
#' @param prop_method Method used in \code{\link{fit_propensity_model}}.
#' @param prop_eps Symmetric clipping constant used in
#'   \code{\link{fit_propensity_model}}.
#' @param lower_tau Lower quantile level.
#' @param upper_tau Upper quantile level.
#' @param ... Additional arguments passed to \code{\link{fit_propensity_model}}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{region}}{A list of prediction regions, one for each row of
#'   \code{x_test}.}
#'   \item{\code{lo_hi}}{A matrix with columns \code{"lo"} and \code{"hi"}
#'   giving the outer interval envelope on \code{y_grid}.}
#'   \item{\code{p_final}}{A matrix of conformal \eqn{p}-values on \code{y_grid}.}
#'   \item{\code{y_grid}}{The candidate grid used in the conformal search.}
#' }
#'
#' @export
lc_region <- function(
    dat,
    id_col,
    y_col = "Y",
    delta_col = "delta",
    x_cols,
    x_test,
    y_grid,
    alpha = 0.1,
    train_frac = 0.5,
    seed = NULL,
    quant_method = c("linear", "grf"),
    prop_method = c("logistic", "grf", "boosting"),
    prop_eps = 1e-6,
    lower_tau = 0.05,
    upper_tau = 0.95,
    ...
) {
  
  # ---------------------------------------------------------------------------
  # Step 0: Match arguments and perform basic input checks
  # ---------------------------------------------------------------------------
  quant_method <- match.arg(quant_method)
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
  
  if (!all(x_cols %in% names(dat))) {
    stop("Some columns specified in x_cols are missing from dat.")
  }
  
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number in (0,1).")
  }
  
  if (!is.numeric(train_frac) || length(train_frac) != 1L || !is.finite(train_frac) ||
      train_frac <= 0 || train_frac >= 1) {
    stop("train_frac must be a single number in (0,1).")
  }
  
  if (!is.numeric(prop_eps) || length(prop_eps) != 1L || !is.finite(prop_eps) ||
      prop_eps <= 0 || prop_eps >= 0.5) {
    stop("prop_eps must be a single finite number in (0, 0.5).")
  }
  
  if (!is.numeric(lower_tau) || !is.numeric(upper_tau) ||
      length(lower_tau) != 1L || length(upper_tau) != 1L ||
      lower_tau <= 0 || upper_tau >= 1 || lower_tau >= upper_tau) {
    stop("lower_tau and upper_tau must satisfy 0 < lower_tau < upper_tau < 1.")
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
    x_test <- matrix(as.numeric(x_test), nrow = 1)
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
  # Helper: Robust binary checks for the missingness indicator
  # ---------------------------------------------------------------------------
  is_delta1 <- function(d_raw) {
    d_chr <- toupper(trimws(as.character(d_raw)))
    d_num <- suppressWarnings(as.numeric(d_chr))
    (!is.na(d_num) & d_num == 1) | (d_chr %in% c("1", "TRUE", "T"))
  }
  
  # ---------------------------------------------------------------------------
  # Helper: Draw one observation per subject
  # ---------------------------------------------------------------------------
  subsample_one_per_subject <- function(dat0, seed0 = NULL) {
    if (!is.null(seed0)) {
      set.seed(seed0)
    }
    
    ids <- unique(dat0[[id_col]])
    pick <- lapply(ids, function(id) {
      rows <- which(dat0[[id_col]] == id)
      sample(rows, size = 1L)
    })
    
    dat0[unlist(pick), , drop = FALSE]
  }
  
  # ---------------------------------------------------------------------------
  # Helper: Split cross-sectional data into training and calibration subsets
  # ---------------------------------------------------------------------------
  split_rows <- function(n_total, train_frac, seed0 = NULL) {
    if (!is.null(seed0)) {
      set.seed(seed0)
    }
    
    n_tr <- floor(train_frac * n_total)
    n_tr <- max(1L, min(n_tr, n_total - 1L))
    
    idx_tr <- sample(seq_len(n_total), size = n_tr, replace = FALSE)
    idx_ca <- setdiff(seq_len(n_total), idx_tr)
    
    list(train = idx_tr, calib = idx_ca)
  }
  
  # ---------------------------------------------------------------------------
  # Step 2-8: Test-point-specific subsampling, fitting, and weighted p-values
  # ---------------------------------------------------------------------------
  p_final <- matrix(NA_real_, nrow = K, ncol = Ny)
  
  for (k in seq_len(K)) {
    
    seed_k_sub <- if (!is.null(seed)) seed + 1000L * k else NULL
    seed_k_sp  <- if (!is.null(seed)) seed + 1000L * k + 1L else NULL
    
    # -------------------------------------------------------------------------
    # Step 2: Subsample one observation per subject for this test point
    # -------------------------------------------------------------------------
    dat_sub <- subsample_one_per_subject(dat, seed0 = seed_k_sub)
    
    n_sub <- nrow(dat_sub)
    if (n_sub < 2L) {
      stop("Need at least two subsampled observations.")
    }
    
    # -------------------------------------------------------------------------
    # Step 3: Split subsampled data into training and calibration subsets
    # -------------------------------------------------------------------------
    sp <- split_rows(n_sub, train_frac = train_frac, seed0 = seed_k_sp)
    
    dat_tr <- dat_sub[sp$train, , drop = FALSE]
    dat_ca <- dat_sub[sp$calib, , drop = FALSE]
    
    # -------------------------------------------------------------------------
    # Step 4: Fit the missingness propensity model on the training subset
    # -------------------------------------------------------------------------
    prop_fit <- fit_propensity_model(
      dat = dat_tr,
      delta_col = delta_col,
      x_cols = x_cols,
      method = prop_method,
      eps = prop_eps,
      ...
    )
    
    # -------------------------------------------------------------------------
    # Step 5: Fit the conditional quantile model on observed training data
    # -------------------------------------------------------------------------
    dtr_obs <- is_delta1(dat_tr[[delta_col]]) & !is.na(dat_tr[[y_col]])
    dat_tr_obs <- dat_tr[dtr_obs, , drop = FALSE]
    
    if (nrow(dat_tr_obs) < 10L) {
      stop("Too few observed training outcomes for quantile fitting.")
    }
    
    X_tr_obs <- dat_tr_obs[, x_cols, drop = FALSE]
    Y_tr_obs <- as.numeric(dat_tr_obs[[y_col]])
    
    if (quant_method == "linear") {
      if (!requireNamespace("quantreg", quietly = TRUE)) {
        stop("Package 'quantreg' is required for quant_method = 'linear'.")
      }
      
      rq_fit <- quantreg::rq(
        formula = stats::as.formula("Y_tr_obs ~ ."),
        data = data.frame(Y_tr_obs = Y_tr_obs, X_tr_obs),
        tau = c(lower_tau, upper_tau)
      )
      
      predict_quantiles <- function(x_new) {
        qhat <- stats::predict(rq_fit, newdata = as.data.frame(x_new))
        qhat <- as.matrix(qhat)
        colnames(qhat) <- c("q_lo", "q_hi")
        qhat
      }
      
      quant_fit_object <- rq_fit
      
    } else {
      if (!requireNamespace("grf", quietly = TRUE)) {
        stop("Package 'grf' is required for quant_method = 'grf'.")
      }
      
      qrf_fit <- grf::quantile_forest(
        X = as.matrix(X_tr_obs),
        Y = Y_tr_obs,
        quantiles = c(lower_tau, upper_tau)
      )
      
      predict_quantiles <- function(x_new) {
        qhat <- predict(
          qrf_fit,
          newdata = as.matrix(x_new),
          quantiles = c(lower_tau, upper_tau)
        )$predictions
        qhat <- as.matrix(qhat)
        colnames(qhat) <- c("q_lo", "q_hi")
        qhat
      }
      
      quant_fit_object <- qrf_fit
    }
    
    # -------------------------------------------------------------------------
    # Step 6: Compute test-side weight and conformity scores for this test point
    # -------------------------------------------------------------------------
    x_test_k <- as.data.frame(matrix(x_test[k, ], nrow = 1))
    colnames(x_test_k) <- x_cols
    
    p_test <- as.numeric(prop_fit$predict(x_test_k))
    if (any(!is.finite(p_test))) {
      stop("predict returned non-finite propensity values for x_test.")
    }
    p_test <- pmax(pmin(p_test, 1 - prop_eps), prop_eps)
    w_test <- 1 / p_test
    
    q_test <- predict_quantiles(x_test_k)
    
    q_lo <- q_test[1, 1]
    q_hi <- q_test[1, 2]
    score_test <- pmax(q_lo - y_grid, y_grid - q_hi)
    
    # -------------------------------------------------------------------------
    # Step 7: Compute calibration-side weights and conformity scores
    # -------------------------------------------------------------------------
    dca_obs <- is_delta1(dat_ca[[delta_col]]) & !is.na(dat_ca[[y_col]])
    dat_ca_obs <- dat_ca[dca_obs, , drop = FALSE]
    
    if (nrow(dat_ca_obs) == 0L) {
      stop("No observed calibration outcomes are available.")
    }
    
    X_ca_obs <- dat_ca_obs[, x_cols, drop = FALSE]
    Y_ca_obs <- as.numeric(dat_ca_obs[[y_col]])
    
    p_cali <- as.numeric(prop_fit$predict(X_ca_obs))
    if (any(!is.finite(p_cali))) {
      stop("predict returned non-finite propensity values in calibration.")
    }
    p_cali <- pmax(pmin(p_cali, 1 - prop_eps), prop_eps)
    w_cali <- 1 / p_cali
    
    q_cali <- predict_quantiles(X_ca_obs)
    score_cali <- pmax(q_cali[, 1] - Y_ca_obs, Y_ca_obs - q_cali[, 2])
    
    # -------------------------------------------------------------------------
    # Step 8: Compute weighted conformal p-values for this test point
    # -------------------------------------------------------------------------
    ord <- order(score_cali)
    score_sorted <- score_cali[ord]
    w_sorted <- w_cali[ord]
    W_ge <- rev(cumsum(rev(w_sorted)))
    sum_w <- sum(w_cali)
    
    denom <- sum_w + w_test
    st <- score_test
    
    j_ge <- findInterval(st, score_sorted, left.open = TRUE) + 1L
    w_ge <- ifelse(j_ge <= length(W_ge), W_ge[j_ge], 0)
    
    p_final[k, ] <- (w_ge + w_test) / denom
  }
  
  # ---------------------------------------------------------------------------
  # Step 9: Construct prediction regions on the y-grid
  # ---------------------------------------------------------------------------
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
  
  # ---------------------------------------------------------------------------
  # Return final output
  # ---------------------------------------------------------------------------
  list(
    region = regions,
    lo_hi = lo_hi,
    p_final = p_final,
    y_grid = y_grid
  )
}