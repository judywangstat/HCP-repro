#' Construct a localized HCP conformal prediction region
#'
#' @description
#' This function constructs a localized conformal prediction region under the
#' hierarchical conformal prediction (HCP) framework for clustered data with
#' missing outcomes.
#'
#' Relative to \code{\link{hcp_region}}, this function modifies the calibration
#' step by using local calibration groups. For each subject-level split, local
#' groups are learned from the training split using
#' \code{\link{build_local_groups}}. For each test point, calibration
#' observations are restricted to the local group assigned to that test point.
#'
#' The procedure consists of five main steps:
#' \enumerate{
#'   \item Fit a pooled conditional density model and, if needed, a pooled
#'   missingness propensity model on a subject-level training split.
#'   \item Learn local groups from the training split and assign both calibration
#'   observations and test points to these groups.
#'   \item Repeatedly construct calibration samples by drawing one observation
#'   per subject from the calibration split.
#'   \item Compute localized weighted conformal \eqn{p}-values on a candidate
#'   response grid using only calibration observations in the same local group
#'   as the test point.
#'   \item Aggregate dependent \eqn{p}-values across repeated subsamples and
#'   repeated data splits.
#' }
#'
#' The final localized prediction region is returned as the subset of
#' \code{y_grid} satisfying \eqn{p_{\mathrm{final}}(y) > \alpha}.
#'
#' @param dat A data frame containing clustered observations.
#' @param id_col Name of the subject or cluster identifier column.
#' @param y_col Name of the outcome column.
#' @param delta_col Name of the missingness indicator column, where 1 indicates
#'   an observed outcome and 0 indicates a missing outcome.
#' @param x_cols Character vector giving the names of the covariate columns.
#' @param x_test New covariate value(s). This can be a numeric vector
#'   treated as one test point, or a matrix/data.frame with one row per test
#'   point and \code{length(x_cols)} columns.
#' @param y_grid Numeric vector of candidate response values at which conformal
#'   \eqn{p}-values are evaluated.
#' @param alpha Miscoverage level in \eqn{(0,1)}.
#' @param train_frac Fraction of subjects assigned to training in each split.
#' @param S Number of repeated subject-level data splits.
#' @param B Number of repeated calibration subsamples per split.
#' @param combine_B Method used to combine \eqn{p}-values across the \code{B}
#'   subsamples. Must be \code{"cct"} or \code{"mean"}.
#' @param combine_S Method used to combine \eqn{p}-values across the \code{S}
#'   splits. Must be \code{"cct"} or \code{"mean"}.
#' @param seed Optional random seed.
#' @param return_details Logical; if \code{TRUE}, also return split-level
#'   \eqn{p}-values and split-level metadata.
#'
#' @param dens_method Method used in \code{\link{fit_cond_density_qp}}.
#'   Must be one of \code{"rq"} or \code{"qrf"}.
#' @param dens_taus Tau grid passed to \code{\link{fit_cond_density_qp}}.
#' @param dens_h Bandwidth argument passed to \code{\link{fit_cond_density_qp}}.
#'   If \code{NULL}, then in each split the bandwidth is recomputed using the
#'   observed outcomes in that split's training data via \code{quantile_levels()}.
#'   If supplied, the same user-provided bandwidth is used in every split.
#' @param dens_h_scenario Optional scenario name passed to
#'   \code{quantile_levels()} when \code{dens_h = NULL}. This is mainly useful
#'   in simulation settings such as \code{"Bimo"} and \code{"Bimo-fix"}.
#'   For real data, this can usually be left as \code{NULL}.
#' @param enforce_monotone Logical argument passed to
#'   \code{\link{fit_cond_density_qp}}.
#' @param tail_decay Logical argument passed to
#'   \code{\link{fit_cond_density_qp}}.
#'
#' @param prop_method Method used in \code{\link{fit_propensity_model}}.
#' @param prop_eps Symmetric clipping constant used in
#'   \code{\link{fit_propensity_model}}.
#' @param weight_cap Optional upper bound used to truncate inverse-propensity
#'   weights for numerical stability. Use \code{Inf} for no truncation.
#'
#' @param n_local_groups Number of local groups used in the localized
#'   calibration step.
#' @param t_num Number of candidate threshold values used when constructing
#'   local groups in \code{\link{build_local_groups}}.
#' @param min_local_size Minimum number of local calibration observations
#'   required for a test point. If the local group is smaller than this value,
#'   the behavior is controlled by \code{fallback_to_global}.
#' @param fallback_to_global Logical; if \code{TRUE}, use the full calibration
#'   sample when the assigned local group has fewer than
#'   \code{min_local_size} observations. If \code{FALSE}, the corresponding
#'   localized \eqn{p}-values are set to missing for that subsample.
#' @param ... Additional arguments passed to \code{\link{fit_propensity_model}}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{region}}{A list of localized prediction regions, one for each
#'   row of \code{x_test}.}
#'   \item{\code{lo_hi}}{A matrix with columns \code{"lo"} and \code{"hi"}
#'   giving the outer interval envelope on \code{y_grid}.}
#'   \item{\code{p_final}}{A matrix of final localized conformal \eqn{p}-values
#'   on \code{y_grid}.}
#'   \item{\code{y_grid}}{The candidate grid used in the conformal search.}
#' }
#'
#' If \code{return_details = TRUE}, the output also includes:
#' \describe{
#'   \item{\code{p_split}}{An array of split-level \eqn{p}-values.}
#'   \item{\code{split_meta}}{A list containing split-level metadata, including
#'   training subject IDs, whether the split is fully observed, and the learned
#'   local grouping information.}
#' }
#'
#' @export
hcp_localized_region <- function(
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
    prop_method = c("logistic", "grf", "boosting"),
    prop_eps = 1e-6,
    weight_cap = Inf,
    n_local_groups = 5L,
    t_num = 200L,
    min_local_size = 5L,
    fallback_to_global = FALSE,
    ...
) {
  # ---------------------------------------------------------------------------
  # Step 0: Match arguments and perform basic input checks
  # ---------------------------------------------------------------------------
  combine_B <- match.arg(combine_B)
  combine_S <- match.arg(combine_S)
  dens_method <- match.arg(dens_method)
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
  
  if (!is.numeric(prop_eps) || length(prop_eps) != 1L ||
      !is.finite(prop_eps) || prop_eps <= 0 || prop_eps >= 0.5) {
    stop("prop_eps must be a single finite number in (0, 0.5).")
  }
  
  if (!is.numeric(weight_cap) || length(weight_cap) != 1L ||
      is.na(weight_cap) || weight_cap < 1) {
    stop("weight_cap must be a single numeric value >= 1, or Inf.")
  }
  
  n_local_groups <- as.integer(n_local_groups)
  if (!is.finite(n_local_groups) || n_local_groups < 1L) {
    stop("n_local_groups must be a positive integer.")
  }
  
  t_num <- as.integer(t_num)
  if (!is.finite(t_num) || t_num < 10L) {
    stop("t_num must be an integer greater than or equal to 10.")
  }
  
  min_local_size <- as.integer(min_local_size)
  if (!is.finite(min_local_size) || min_local_size < 1L) {
    stop("min_local_size must be a positive integer.")
  }
  
  if (!is.logical(fallback_to_global) || length(fallback_to_global) != 1L) {
    stop("fallback_to_global must be TRUE or FALSE.")
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
  # Helper: Combine p-values across the first array dimension
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
      inf <- is.infinite(tbar)
      tbar[inf] <- sign(tbar[inf]) * .Machine$double.xmax
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
      calib = dat0[!idx_tr, , drop = FALSE],
      train_ids = tr_ids
    )
  }
  
  # ---------------------------------------------------------------------------
  # Helper: Draw one observation per subject from a calibration split
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
  # Helper: Robust binary checks for the missingness indicator
  # ---------------------------------------------------------------------------
  is_delta1 <- function(d_raw) {
    d_chr <- toupper(trimws(as.character(d_raw)))
    d_num <- suppressWarnings(as.numeric(d_chr))
    (!is.na(d_num) & d_num == 1) | (d_chr %in% c("1", "TRUE", "T"))
  }
  
  is_delta0 <- function(d_raw) {
    d_chr <- toupper(trimws(as.character(d_raw)))
    d_num <- suppressWarnings(as.numeric(d_chr))
    (!is.na(d_num) & d_num == 0) | (d_chr %in% c("0", "FALSE", "F"))
  }
  
  # ---------------------------------------------------------------------------
  # Helper: determine split-specific bandwidth for density estimation
  # ---------------------------------------------------------------------------
  resolve_dens_h_split <- function(dat_tr_split) {
    if (!is.null(dens_h)) {
      return(dens_h)
    }
    
    delta_tr <- dat_tr_split[[delta_col]]
    obs_idx <- is_delta1(delta_tr) & !is.na(dat_tr_split[[y_col]])
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
  # Step 2: Initialize storage for split-level p-values
  # ---------------------------------------------------------------------------
  p_split_arr <- array(NA_real_, dim = c(S, K, Ny))
  split_meta <- if (isTRUE(return_details)) vector("list", S) else NULL
  
  # ---------------------------------------------------------------------------
  # Step 3: Repeat over subject-level splits
  # ---------------------------------------------------------------------------
  for (s in seq_len(S)) {
    seed_s <- if (!is.null(seed)) seed + 1000L * s else NULL
    sp <- split_by_subject(dat, seed0 = seed_s)
    
    dat_tr <- sp$train
    dat_ca <- sp$calib
    
    # -------------------------------------------------------------------------
    # Step 3A: Resolve bandwidth h for this split
    # -------------------------------------------------------------------------
    dens_h_use <- resolve_dens_h_split(dat_tr)
    
    # -------------------------------------------------------------------------
    # Step 3B: Fit the pooled conditional density model on the training split
    # -------------------------------------------------------------------------
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
    # Step 3C: Fit the pooled missingness propensity model, if needed
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
      w_test_vec <- rep(1, K)
    } else {
      x_test_df <- as.data.frame(x_test)
      colnames(x_test_df) <- x_cols
      
      p_test <- as.numeric(prop_fit$predict(x_test_df))
      if (any(!is.finite(p_test))) {
        stop("predict returned non-finite propensity values for x_test.")
      }
      
      p_test <- pmax(pmin(p_test, 1 - prop_eps), prop_eps)
      w_test_vec <- 1 / p_test
      w_test_vec <- pmin(w_test_vec, weight_cap)
    }
    
    # -------------------------------------------------------------------------
    # Step 3E: Compute test-side nonconformity scores on the full y_grid
    # -------------------------------------------------------------------------
    dens_te_mat <- dens_fit$predict_density_grid(x_test, y_grid)
    
    if (any(!is.finite(dens_te_mat))) {
      stop("predict_density_grid returned non-finite values for the test grid.")
    }
    
    R_test_mat <- -dens_te_mat
    
    # -------------------------------------------------------------------------
    # Step 3F: Build local groups once for this split
    # -------------------------------------------------------------------------
    local_groups_split <- build_local_groups(
      X_train = dat_tr[, x_cols, drop = FALSE],
      X_cal = dat_ca[, x_cols, drop = FALSE],
      X_test = x_test,
      dens_fit = dens_fit,
      y_grid = y_grid,
      n_local_groups = n_local_groups,
      t_num = t_num,
      seed = seed_s
    )
    
    group_cal_full <- local_groups_split$group_cal
    group_test <- local_groups_split$group_test
    
    # -------------------------------------------------------------------------
    # Step 3G: Build repeated calibration subsamples
    # -------------------------------------------------------------------------
    sub_cache <- vector("list", B)
    
    for (b in seq_len(B)) {
      seed_b <- if (!is.null(seed_s)) seed_s + 10000L + b else NULL
      dat_sub <- subsample_one_per_subject(dat_ca, seed0 = seed_b)
      
      sampled_row_idx <- match(rownames(dat_sub), rownames(dat_ca))
      if (anyNA(sampled_row_idx)) {
        stop("Internal error: calibration subsample row mapping failed.")
      }
      
      d_raw <- dat_sub[[delta_col]]
      delta1 <- is_delta1(d_raw)
      
      d_obs <- dat_sub[delta1 & !is.na(dat_sub[[y_col]]), , drop = FALSE]
      obs_row_idx <- sampled_row_idx[delta1 & !is.na(dat_sub[[y_col]])]
      
      if (nrow(d_obs) == 0L) {
        sub_cache[[b]] <- list(empty = TRUE)
        next
      }
      
      Xc_df <- as.data.frame(as.matrix(d_obs[, x_cols, drop = FALSE]))
      colnames(Xc_df) <- x_cols
      Yc <- as.numeric(d_obs[[y_col]])
      
      dens_ca_mat <- dens_fit$predict_density_grid(Xc_df, Yc)
      R_c <- -diag(dens_ca_mat)
      
      if (any(!is.finite(R_c))) {
        stop("predict_density_grid returned non-finite values in calibration.")
      }
      
      if (fully_observed_split) {
        w_c <- rep(1, nrow(Xc_df))
      } else {
        p_c <- as.numeric(prop_fit$predict(Xc_df))
        if (any(!is.finite(p_c))) {
          stop("predict returned non-finite propensity values in calibration.")
        }
        
        p_c <- pmax(pmin(p_c, 1 - prop_eps), prop_eps)
        w_c <- 1 / p_c
        w_c <- pmin(w_c, weight_cap)
      }
      
      group_cal_obs <- group_cal_full[obs_row_idx]
      
      sub_cache[[b]] <- list(
        empty = FALSE,
        R_c = R_c,
        w_c = w_c,
        group_cal_obs = group_cal_obs
      )
    }
    
    # -------------------------------------------------------------------------
    # Step 3H: Compute localized conformal p-values for each subsample
    # -------------------------------------------------------------------------
    p_bky <- array(NA_real_, dim = c(B, K, Ny))
    
    for (b in seq_len(B)) {
      cb <- sub_cache[[b]]
      
      if (isTRUE(cb$empty)) {
        p_bky[b, , ] <- NA_real_
        next
      }
      
      for (k in seq_len(K)) {
        idx_local <- which(cb$group_cal_obs == group_test[k])
        
        if (length(idx_local) < min_local_size) {
          if (isTRUE(fallback_to_global)) {
            idx_local <- seq_along(cb$R_c)
          } else {
            p_bky[b, k, ] <- NA_real_
            next
          }
        }
        
        R_local <- cb$R_c[idx_local]
        w_local <- cb$w_c[idx_local]
        
        ord <- order(R_local)
        R_sorted <- R_local[ord]
        w_sorted <- w_local[ord]
        W_ge <- rev(cumsum(rev(w_sorted)))
        sum_w <- sum(w_local)
        
        denom <- sum_w + w_test_vec[k]
        rt <- R_test_mat[k, ]
        
        j_ge <- findInterval(rt, R_sorted, left.open = TRUE) + 1L
        w_ge <- ifelse(j_ge <= length(W_ge), W_ge[j_ge], 0)
        
        p_bky[b, k, ] <- (w_ge + w_test_vec[k]) / denom
      }
    }
    
    # -------------------------------------------------------------------------
    # Step 3I: Aggregate over the B subsamples within this split
    # -------------------------------------------------------------------------
    p_split_arr[s, , ] <- combine_array_firstdim(p_bky, method = combine_B)
    
    if (isTRUE(return_details)) {
      split_meta[[s]] <- list(
        train_ids = sp$train_ids,
        fully_observed = fully_observed_split,
        local_groups = local_groups_split
      )
    }
  }
  
  # ---------------------------------------------------------------------------
  # Step 4: Aggregate over the S repeated splits
  # ---------------------------------------------------------------------------
  p_final <- combine_array_firstdim(p_split_arr, method = combine_S)
  
  # ---------------------------------------------------------------------------
  # Step 5: Construct localized prediction regions on the y_grid
  # ---------------------------------------------------------------------------
  regions <- vector("list", K)
  lo_hi <- matrix(NA_real_, nrow = K, ncol = 2)
  colnames(lo_hi) <- c("lo", "hi")
  
  for (k in seq_len(K)) {
    reg <- y_grid[p_final[k, ] > alpha]
    regions[[k]] <- reg
    
    if (length(reg) == 0L) {
      lo_hi[k, ] <- c(NA_real_, NA_real_)
    } else {
      lo_hi[k, ] <- c(min(reg), max(reg))
    }
  }
  
  out <- list(
    region = regions,
    lo_hi = lo_hi,
    p_final = p_final,
    y_grid = y_grid
  )
  
  if (isTRUE(return_details)) {
    out$p_split <- p_split_arr
    out$split_meta <- split_meta
  }
  
  out
}