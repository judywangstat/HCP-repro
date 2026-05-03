# ============================================================
# dwr_region_legacy.R
# Legacy DWR-style prediction region
# ============================================================

#' Construct a legacy DWR-style prediction region
#'
#' @description
#' Legacy-compatible DWR implementation for gallstones reproduction. For each
#' test point, this function independently draws one observation per subject,
#' splits the resulting cross-sectional sample into training and calibration
#' sets, fits a score model on observed training outcomes, and computes
#' unweighted conformal p-values.
#'
#' @param score_type Nonconformity score type. Use \code{"residual"} for the
#'   original residual-score DWR, or \code{"quantile"} for a quantile-score
#'   DWR under the same legacy per-test-point resampling framework.
#'
#' @export
dwr_region_legacy <- function(
    dat,
    id_col,
    y_col = "Y",
    delta_col = "delta",
    x_cols,
    x_test,
    y_grid,
    alpha = 0.1,
    train_frac = 0.5,
    B = 1L,
    seed = NULL,
    score_type = c("residual", "quantile"),
    reg_method = c("linear", "grf"),
    quant_method = c("linear", "grf"),
    lower_tau = 0.05,
    upper_tau = 0.95
) {
  score_type <- match.arg(score_type)
  reg_method <- match.arg(reg_method)
  quant_method <- match.arg(quant_method)
  
  if (!is.data.frame(dat)) stop("dat must be a data.frame.")
  if (!all(c(id_col, y_col, delta_col, x_cols) %in% names(dat))) {
    stop("dat is missing required columns.")
  }
  
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number in (0,1).")
  }
  
  if (!is.numeric(train_frac) || length(train_frac) != 1L ||
      !is.finite(train_frac) || train_frac <= 0 || train_frac >= 1) {
    stop("train_frac must be a single number in (0,1).")
  }
  
  B <- as.integer(B)
  if (!is.finite(B) || B < 1L) {
    stop("B must be a positive integer.")
  }
  
  if (!is.numeric(lower_tau) || !is.numeric(upper_tau) ||
      length(lower_tau) != 1L || length(upper_tau) != 1L ||
      lower_tau <= 0 || upper_tau >= 1 || lower_tau >= upper_tau) {
    stop("lower_tau and upper_tau must satisfy 0 < lower_tau < upper_tau < 1.")
  }
  
  y_grid <- sort(unique(as.numeric(y_grid)))
  if (anyNA(y_grid) || any(!is.finite(y_grid)) || length(y_grid) < 2L) {
    stop("y_grid must contain at least two finite numeric values.")
  }
  Ny <- length(y_grid)
  
  if (is.vector(x_test) && !is.list(x_test)) {
    x_test <- matrix(as.numeric(x_test), nrow = 1L)
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
  
  is_delta1 <- function(d_raw) {
    d_chr <- toupper(trimws(as.character(d_raw)))
    d_num <- suppressWarnings(as.numeric(d_chr))
    (!is.na(d_num) & d_num == 1) | (d_chr %in% c("1", "TRUE", "T"))
  }
  
  subsample_one_per_subject <- function(dat0) {
    ids <- unique(dat0[[id_col]])
    pick <- lapply(ids, function(id) {
      rows <- which(dat0[[id_col]] == id)
      sample(rows, size = 1L)
    })
    dat0[unlist(pick), , drop = FALSE]
  }
  
  split_rows <- function(n_total) {
    n_tr <- floor(train_frac * n_total)
    n_tr <- max(1L, min(n_tr, n_total - 1L))
    
    idx_tr <- sample(seq_len(n_total), size = n_tr, replace = FALSE)
    idx_ca <- setdiff(seq_len(n_total), idx_tr)
    
    list(train = idx_tr, calib = idx_ca)
  }
  
  fit_regression_model <- function(dat_tr_obs) {
    X_tr_obs <- dat_tr_obs[, x_cols, drop = FALSE]
    Y_tr_obs <- as.numeric(dat_tr_obs[[y_col]])
    
    if (reg_method == "linear") {
      reg_dat <- data.frame(Y_tr_obs = Y_tr_obs, X_tr_obs)
      
      lm_fit <- stats::lm(
        formula = stats::as.formula("Y_tr_obs ~ ."),
        data = reg_dat
      )
      
      return(function(x_new) {
        as.numeric(stats::predict(lm_fit, newdata = as.data.frame(x_new)))
      })
    }
    
    if (!requireNamespace("grf", quietly = TRUE)) {
      stop("Package 'grf' is required for reg_method = 'grf'.")
    }
    
    rf_fit <- grf::regression_forest(
      X = as.matrix(X_tr_obs),
      Y = Y_tr_obs
    )
    
    function(x_new) {
      as.numeric(stats::predict(
        rf_fit,
        newdata = as.matrix(x_new)
      )$predictions)
    }
  }
  
  fit_quantile_model <- function(dat_tr_obs) {
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
      
      return(function(x_new) {
        qhat <- stats::predict(rq_fit, newdata = as.data.frame(x_new))
        qhat <- as.matrix(qhat)
        if (ncol(qhat) != 2L) {
          qhat <- matrix(qhat, ncol = 2L)
        }
        colnames(qhat) <- c("q_lo", "q_hi")
        qhat
      })
    }
    
    if (!requireNamespace("grf", quietly = TRUE)) {
      stop("Package 'grf' is required for quant_method = 'grf'.")
    }
    
    qrf_fit <- grf::quantile_forest(
      X = as.matrix(X_tr_obs),
      Y = Y_tr_obs,
      quantiles = c(lower_tau, upper_tau)
    )
    
    function(x_new) {
      qhat <- stats::predict(
        qrf_fit,
        newdata = as.matrix(x_new),
        quantiles = c(lower_tau, upper_tau)
      )$predictions
      qhat <- as.matrix(qhat)
      colnames(qhat) <- c("q_lo", "q_hi")
      qhat
    }
  }
  
  p_final <- matrix(NA_real_, nrow = K, ncol = Ny)
  
  for (k in seq_len(K)) {
    p_by <- matrix(NA_real_, nrow = B, ncol = Ny)
    x_test_one <- x_test[k, , drop = FALSE]
    
    for (b in seq_len(B)) {
      dat_sub <- subsample_one_per_subject(dat)
      
      if (nrow(dat_sub) < 2L) {
        stop("Need at least two subsampled observations.")
      }
      
      sp <- split_rows(nrow(dat_sub))
      dat_tr <- dat_sub[sp$train, , drop = FALSE]
      dat_ca <- dat_sub[sp$calib, , drop = FALSE]
      
      dtr_obs <- is_delta1(dat_tr[[delta_col]]) & !is.na(dat_tr[[y_col]])
      dat_tr_obs <- dat_tr[dtr_obs, , drop = FALSE]
      
      if (nrow(dat_tr_obs) < 10L) {
        next
      }
      
      dca_obs <- is_delta1(dat_ca[[delta_col]]) & !is.na(dat_ca[[y_col]])
      dat_ca_obs <- dat_ca[dca_obs, , drop = FALSE]
      
      if (nrow(dat_ca_obs) == 0L) {
        next
      }
      
      X_ca_obs <- dat_ca_obs[, x_cols, drop = FALSE]
      Y_ca_obs <- as.numeric(dat_ca_obs[[y_col]])
      
      if (score_type == "residual") {
        predict_regression <- fit_regression_model(dat_tr_obs)
        
        reg_test <- predict_regression(x_test_one)
        score_test <- abs(reg_test - y_grid)
        
        reg_cali <- predict_regression(X_ca_obs)
        score_cali <- abs(Y_ca_obs - reg_cali)
        
      } else {
        predict_quantiles <- fit_quantile_model(dat_tr_obs)
        
        q_test <- predict_quantiles(x_test_one)
        score_test <- pmax(q_test[, 1] - y_grid, y_grid - q_test[, 2])
        
        q_cali <- predict_quantiles(X_ca_obs)
        score_cali <- pmax(q_cali[, 1] - Y_ca_obs, Y_ca_obs - q_cali[, 2])
      }
      
      ord <- order(score_cali)
      score_sorted <- score_cali[ord]
      n_cali <- length(score_sorted)
      
      j_ge <- findInterval(score_test, score_sorted, left.open = TRUE) + 1L
      n_ge <- ifelse(j_ge <= n_cali, n_cali - j_ge + 1L, 0L)
      
      p_by[b, ] <- (n_ge + 1) / (n_cali + 1)
    }
    
    p_final[k, ] <- colMeans(p_by, na.rm = TRUE)
  }
  
  regions <- vector("list", K)
  lo_hi <- matrix(NA_real_, nrow = K, ncol = 2)
  colnames(lo_hi) <- c("lo", "hi")
  
  for (k in seq_len(K)) {
    reg <- y_grid[p_final[k, ] >= alpha]
    regions[[k]] <- reg
    
    if (length(reg) == 0L || all(is.na(reg))) {
      lo_hi[k, ] <- c(min(y_grid), max(y_grid))
      regions[[k]] <- y_grid
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