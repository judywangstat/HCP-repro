#' Construct a legacy LC prediction region for real-data reproduction
#'
#' @description
#' This function implements the legacy LC procedure used for reproducing the
#' CD4 and gallstones real-data analyses. It adapts an LC-style conformal
#' prediction method to clustered data with missing outcomes by drawing one
#' observation per subject, splitting the resulting cross-sectional sample into
#' training and calibration subsets, and computing weighted conformal
#' \eqn{p}-values on a candidate response grid.
#'
#' Compared with \code{\link{lc_region}}, this legacy version keeps the
#' real-data reproduction workflow used in the original analysis. In particular,
#' it supports optional time-matched subsampling for longitudinal real-data
#' examples and uses an internal legacy implementation of the quantile and
#' propensity fitting steps.
#'
#' The procedure consists of the following main steps:
#' \enumerate{
#'   \item Optionally restrict the data to observations matching the test time.
#'   \item Draw one observation per subject to form a cross-sectional sample.
#'   \item Split the subsampled data into training and calibration subsets.
#'   \item Fit a missingness propensity model on the training subset.
#'   \item Fit lower and upper conditional quantile models on observed training
#'   outcomes.
#'   \item Compute weighted conformal \eqn{p}-values using observed calibration
#'   outcomes.
#'   \item Return the prediction region on \code{y_grid}.
#' }
#'
#' @param dat A data frame containing clustered or longitudinal observations.
#' @param id_col Name of the subject or cluster identifier column.
#' @param y_col Name of the outcome column.
#' @param delta_col Name of the missingness indicator column, where 1 indicates
#'   an observed outcome and 0 indicates a missing outcome.
#' @param x_cols Character vector giving the names of the covariate columns.
#' @param x_test New covariate value(s). This can be a numeric vector treated as
#'   one test point, or a matrix/data.frame with one row per test point and
#'   \code{length(x_cols)} columns.
#' @param y_grid Numeric vector of candidate response values at which conformal
#'   \eqn{p}-values are evaluated.
#' @param alpha Miscoverage level in \eqn{(0,1)}.
#' @param train_frac Fraction of subsampled observations assigned to training.
#' @param seed Optional random seed.
#' @param quant_method Quantile model used to estimate lower and upper
#'   conditional quantiles. Must be \code{"linear"} or \code{"grf"}.
#' @param prop_method Propensity model used to estimate the probability of an
#'   observed response. Must be \code{"logistic"}, \code{"grf"}, or
#'   \code{"boosting"}.
#' @param lower_tau Lower quantile level.
#' @param upper_tau Upper quantile level.
#' @param weight_cap Upper bound used to truncate inverse-propensity weights
#'   for numerical stability.
#' @param time_matched Logical; if \code{TRUE}, the calibration pool is
#'   restricted to observations with the same time value as the corresponding
#'   test point.
#' @param time_col Name of the time variable used when
#'   \code{time_matched = TRUE}.
#' @param test_time Test-time value(s) used when \code{time_matched = TRUE}.
#'   Its length must match the number of rows in \code{x_test}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{region}}{A list of prediction regions, one for each row of
#'   \code{x_test}.}
#'   \item{\code{lo_hi}}{A matrix with columns \code{"lo"} and \code{"hi"}
#'   giving the outer interval envelope on \code{y_grid}.}
#'   \item{\code{p_final}}{A matrix of conformal \eqn{p}-values on
#'   \code{y_grid}.}
#'   \item{\code{y_grid}}{The candidate grid used in the conformal search.}
#' }
lc_region_legacy <- function(
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
    lower_tau = 0.05,
    upper_tau = 0.95,
    weight_cap = 50,
    time_matched = FALSE,
    time_col = "time",
    test_time = NULL
) {
  quant_method <- match.arg(quant_method)
  prop_method <- match.arg(prop_method)
  
  if (!is.data.frame(dat)) stop("dat must be a data.frame.")
  if (!all(c(id_col, y_col, delta_col, x_cols) %in% names(dat))) {
    stop("dat is missing required columns.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0,1).")
  }
  if (!is.numeric(train_frac) || length(train_frac) != 1L ||
      train_frac <= 0 || train_frac >= 1) {
    stop("train_frac must be in (0,1).")
  }
  if (!is.numeric(weight_cap) || length(weight_cap) != 1L ||
      weight_cap < 1) {
    stop("weight_cap must be >= 1.")
  }
  if (!is.logical(time_matched) || length(time_matched) != 1L) {
    stop("time_matched must be TRUE or FALSE.")
  }
  
  y_grid <- sort(unique(as.numeric(y_grid)))
  if (length(y_grid) < 2L || anyNA(y_grid) || any(!is.finite(y_grid))) {
    stop("y_grid must contain at least two finite numeric values.")
  }
  Ny <- length(y_grid)
  
  if (is.vector(x_test) && !is.list(x_test)) {
    x_test <- matrix(as.numeric(x_test), nrow = 1)
  } else {
    x_test <- as.matrix(x_test)
  }
  storage.mode(x_test) <- "numeric"
  if (ncol(x_test) != length(x_cols)) {
    stop("x_test must have length(x_cols) columns.")
  }
  colnames(x_test) <- x_cols
  K <- nrow(x_test)
  
  if (isTRUE(time_matched)) {
    if (!time_col %in% names(dat)) {
      stop("time_col is not found in dat.")
    }
    if (is.null(test_time)) {
      stop("test_time must be provided when time_matched = TRUE.")
    }
    if (length(test_time) != K) {
      stop("test_time must have length equal to nrow(x_test).")
    }
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  is_delta1 <- function(d_raw) {
    d_chr <- toupper(trimws(as.character(d_raw)))
    d_num <- suppressWarnings(as.numeric(d_chr))
    (!is.na(d_num) & d_num == 1) | (d_chr %in% c("1", "TRUE", "T"))
  }
  
  drop_constant_x <- function(dat0, x_cols0) {
    keep <- vapply(x_cols0, function(nm) {
      x <- dat0[[nm]]
      length(unique(x[!is.na(x)])) > 1L
    }, logical(1))
    
    x_use <- x_cols0[keep]
    if (length(x_use) == 0L) return(character(0))
    x_use
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
  
  fit_propensity_legacy <- function(dat_tr) {
    delta_num <- as.integer(is_delta1(dat_tr[[delta_col]]))
    
    if (length(unique(delta_num)) == 1L) {
      p_const <- ifelse(unique(delta_num) == 1L, 1, 1e-6)
      
      return(function(x_new) {
        n_new <- nrow(as.data.frame(x_new))
        rep(p_const, n_new)
      })
    }
    
    x_use <- drop_constant_x(dat_tr, x_cols)
    
    if (length(x_use) == 0L) {
      p_const <- mean(delta_num)
      p_const <- pmax(pmin(p_const, 1 - 1e-6), 1e-6)
      
      return(function(x_new) {
        n_new <- nrow(as.data.frame(x_new))
        rep(p_const, n_new)
      })
    }
    
    X_tr <- dat_tr[, x_use, drop = FALSE]
    Delta_tr <- as.factor(delta_num)
    
    if (prop_method == "logistic") {
      fit <- stats::glm(
        Delta_tr ~ .,
        data = as.data.frame(X_tr),
        family = stats::binomial()
      )
      
      return(function(x_new) {
        x_new <- as.data.frame(x_new)[, x_use, drop = FALSE]
        as.numeric(stats::predict(
          fit,
          newdata = x_new,
          type = "response"
        ))
      })
    }
    
    if (prop_method == "grf") {
      if (!requireNamespace("grf", quietly = TRUE)) {
        stop("Package 'grf' is required for prop_method = 'grf'.")
      }
      
      fit <- grf::probability_forest(
        X = as.matrix(X_tr),
        Y = Delta_tr
      )
      
      return(function(x_new) {
        x_new <- as.data.frame(x_new)[, x_use, drop = FALSE]
        pred <- stats::predict(fit, newdata = as.matrix(x_new))$predictions
        
        if (is.matrix(pred)) {
          if ("1" %in% colnames(pred)) {
            as.numeric(pred[, "1"])
          } else {
            as.numeric(pred[, ncol(pred)])
          }
        } else {
          as.numeric(pred)
        }
      })
    }
    
    if (prop_method == "boosting") {
      if (!requireNamespace("gbm", quietly = TRUE)) {
        stop("Package 'gbm' is required for prop_method = 'boosting'.")
      }
      
      boost_dat <- data.frame(delta_num = delta_num, X_tr)
      
      fit <- gbm::gbm(
        formula = delta_num ~ .,
        data = boost_dat,
        distribution = "bernoulli",
        n.trees = 2000,
        interaction.depth = 2,
        shrinkage = 0.01,
        bag.fraction = 0.8,
        train.fraction = 1.0,
        verbose = FALSE
      )
      
      return(function(x_new) {
        x_new <- as.data.frame(x_new)[, x_use, drop = FALSE]
        as.numeric(stats::predict(
          fit,
          newdata = x_new,
          n.trees = 2000,
          type = "response"
        ))
      })
    }
  }
  
  fit_quantile_legacy <- function(dat_tr_obs) {
    Y_obs <- as.numeric(dat_tr_obs[[y_col]])
    x_use <- drop_constant_x(dat_tr_obs, x_cols)
    
    if (quant_method == "linear") {
      if (!requireNamespace("quantreg", quietly = TRUE)) {
        stop("Package 'quantreg' is required for quant_method = 'linear'.")
      }
      
      if (length(x_use) == 0L) {
        return(function(x_new) {
          n_new <- nrow(as.data.frame(x_new))
          qhat_one <- stats::quantile(
            Y_obs,
            probs = c(lower_tau, upper_tau),
            names = FALSE,
            na.rm = TRUE
          )
          matrix(rep(qhat_one, times = n_new), nrow = n_new, byrow = TRUE)
        })
      }
      
      X_obs <- dat_tr_obs[, x_use, drop = FALSE]
      
      rq_fit <- quantreg::rq(
        Y_obs ~ .,
        data = as.data.frame(X_obs),
        tau = c(lower_tau, upper_tau)
      )
      
      return(function(x_new) {
        x_new <- as.data.frame(x_new)[, x_use, drop = FALSE]
        qhat <- stats::predict(rq_fit, newdata = x_new)
        qhat <- as.matrix(qhat)
        
        if (ncol(qhat) == 1L) {
          qhat <- matrix(qhat, ncol = 2L)
        }
        
        qhat
      })
    }
    
    if (quant_method == "grf") {
      if (!requireNamespace("grf", quietly = TRUE)) {
        stop("Package 'grf' is required for quant_method = 'grf'.")
      }
      
      if (length(x_use) == 0L) {
        qhat_one <- stats::quantile(
          Y_obs,
          probs = c(lower_tau, upper_tau),
          names = FALSE,
          na.rm = TRUE
        )
        
        return(function(x_new) {
          n_new <- nrow(as.data.frame(x_new))
          matrix(rep(qhat_one, times = n_new), nrow = n_new, byrow = TRUE)
        })
      }
      
      X_obs <- dat_tr_obs[, x_use, drop = FALSE]
      
      qrf_fit <- grf::quantile_forest(
        X = as.matrix(X_obs),
        Y = Y_obs,
        quantiles = c(lower_tau, upper_tau)
      )
      
      return(function(x_new) {
        x_new <- as.data.frame(x_new)[, x_use, drop = FALSE]
        as.matrix(stats::predict(
          qrf_fit,
          newdata = as.matrix(x_new),
          quantiles = c(lower_tau, upper_tau)
        )$predictions)
      })
    }
  }
  
  p_final <- matrix(NA_real_, nrow = K, ncol = Ny)
  
  for (k in seq_len(K)) {
    x_test_one <- x_test[k, , drop = FALSE]
    
    dat_pool <- dat
    if (isTRUE(time_matched)) {
      dat_pool <- dat_pool[dat_pool[[time_col]] == test_time[k], , drop = FALSE]
      
      if (length(unique(dat_pool[[id_col]])) < 2L) {
        next
      }
    }
    
    dat_sub <- subsample_one_per_subject(dat_pool)
    sp <- split_rows(nrow(dat_sub))
    
    dat_tr <- dat_sub[sp$train, , drop = FALSE]
    dat_ca <- dat_sub[sp$calib, , drop = FALSE]
    
    prop_predict <- fit_propensity_legacy(dat_tr)
    
    dtr_obs <- is_delta1(dat_tr[[delta_col]]) & !is.na(dat_tr[[y_col]])
    dat_tr_obs <- dat_tr[dtr_obs, , drop = FALSE]
    
    if (nrow(dat_tr_obs) < 5L) {
      next
    }
    
    quant_predict <- fit_quantile_legacy(dat_tr_obs)
    
    p_test <- prop_predict(x_test_one)
    p_test <- pmax(pmin(p_test, 1 - 1e-6), 1e-6)
    w_test <- pmin(1 / p_test, weight_cap)
    
    q_test <- quant_predict(x_test_one)
    score_test <- pmax(q_test[, 1] - y_grid, y_grid - q_test[, 2])
    
    dca_obs <- is_delta1(dat_ca[[delta_col]]) & !is.na(dat_ca[[y_col]])
    dat_ca_obs <- dat_ca[dca_obs, , drop = FALSE]
    
    if (nrow(dat_ca_obs) == 0L) {
      next
    }
    
    X_ca_obs <- dat_ca_obs[, x_cols, drop = FALSE]
    Y_ca_obs <- as.numeric(dat_ca_obs[[y_col]])
    
    p_cali <- prop_predict(X_ca_obs)
    p_cali <- pmax(pmin(p_cali, 1 - 1e-6), 1e-6)
    w_cali <- pmin(1 / p_cali, weight_cap)
    
    q_cali <- quant_predict(X_ca_obs)
    score_cali <- pmax(q_cali[, 1] - Y_ca_obs, Y_ca_obs - q_cali[, 2])
    
    ord <- order(score_cali)
    score_sorted <- score_cali[ord]
    w_sorted <- w_cali[ord]
    
    W_ge <- rev(cumsum(rev(w_sorted)))
    sum_w <- sum(w_cali)
    
    j_ge <- findInterval(score_test, score_sorted, left.open = TRUE) + 1L
    w_ge <- ifelse(j_ge <= length(W_ge), W_ge[j_ge], 0)
    
    p_final[k, ] <- (w_ge + w_test) / (sum_w + w_test)
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