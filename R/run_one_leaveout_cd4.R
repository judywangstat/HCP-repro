# ============================================================
# run_one_leaveout_cd4.R
# Leave-one-subject-out wrapper for the CD4 real-data analysis
# ============================================================

#' Run one leave-one-subject-out CD4 analysis
#'
#' @description
#' Holds out one CD4 subject, fits one selected method on the remaining subjects,
#' and summarizes coverage and prediction-region length for the held-out subject.
#'
#' For pointwise prediction, the method uses alpha directly. For simultaneous
#' prediction, the method uses alpha divided by the number of observations for
#' the held-out subject.
#'
#' @param dat CD4 data frame in long format.
#' @param test_id Subject ID to hold out.
#' @param method Method name: "HCP", "DWR", "LC", or "LMEM".
#' @param prediction_type "pointwise" or "simultaneous".
#' @param missing_rate Missingness setting: "20" or "50".
#' @param alpha Miscoverage level.
#' @param n_grid Number of grid points.
#' @param seed Optional random seed.
#'
#' @return A named numeric vector with entries coverage and length.
#'
#' @export
run_one_leaveout_cd4 <- function(
    dat,
    test_id,
    method = c("HCP", "DWR", "LC", "LMEM"),
    prediction_type = c("pointwise", "simultaneous"),
    missing_rate = c("20", "50"),
    alpha = 0.1,
    n_grid = 200,
    seed = NULL
) {
  method <- match.arg(method)
  prediction_type <- match.arg(prediction_type)
  missing_rate <- match.arg(as.character(missing_rate), choices = c("20", "50"))
  
  if (!is.data.frame(dat)) {
    stop("dat must be a data.frame.")
  }
  
  x_cols <- get_cd4_x_cols()
  required_cols <- c("id", "Y", "delta", x_cols)
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0L) {
    stop("dat is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (!test_id %in% unique(dat$id)) {
    stop("test_id is not found in dat$id.")
  }
  
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number in (0,1).")
  }
  
  n_grid <- as.integer(n_grid)
  if (!is.finite(n_grid) || n_grid < 2L) {
    stop("n_grid must be an integer >= 2.")
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  dat_test <- dat[dat$id == test_id, , drop = FALSE]
  dat_sample <- dat[dat$id != test_id, , drop = FALSE]
  
  y_true <- as.numeric(dat_test$Y)
  x_test <- as.matrix(dat_test[, x_cols, drop = FALSE])
  storage.mode(x_test) <- "numeric"
  colnames(x_test) <- x_cols
  
  y_grid <- make_cd4_y_grid(
    dat_sample = dat_sample,
    y_col = "Y",
    n_grid = n_grid
  )
  
  alpha_use <- if (prediction_type == "simultaneous") {
    alpha / nrow(dat_test)
  } else {
    alpha
  }
  
  weight_cap_use <- get_cd4_weight_cap(
    method = method,
    missing_rate = missing_rate
  )
  
  if (prediction_type == "simultaneous" && method %in% c("HCP", "LC")) {
    weight_cap_use <- 3
  }
  
  res <- switch(
    method,
    
    HCP = hcp_region(
      dat = dat_sample,
      id_col = "id",
      y_col = "Y",
      delta_col = "delta",
      x_cols = x_cols,
      x_test = x_test,
      y_grid = y_grid,
      alpha = alpha_use,
      train_frac = 0.3,
      S = 5,
      B = 5,
      combine_B = "cct",
      combine_S = "cct",
      dens_method = "grf",
      dens_taus = (1:(2^4 - 1)) / (2^4),
      dens_h = NULL,
      dens_h_scenario = NULL,
      enforce_monotone = FALSE,
      tail_decay = TRUE,
      prop_method = "logistic",
      weight_cap = weight_cap_use,
      seed = seed
    ),
    
    DWR = dwr_region_realdata(
      dat = dat_sample,
      id_col = "id",
      y_col = "Y",
      delta_col = "delta",
      x_cols = x_cols,
      x_test = x_test,
      y_grid = y_grid,
      alpha = alpha_use,
      train_frac = 0.5,
      B = 1,
      seed = seed,
      score_type = "residual",
      reg_method = "grf",
      quant_method = "grf",
      lower_tau = 0.05,
      upper_tau = 0.95
    ),
    
    LC = lc_region_realdata(
      dat = dat_sample,
      id_col = "id",
      y_col = "Y",
      delta_col = "delta",
      x_cols = x_cols,
      x_test = x_test,
      y_grid = y_grid,
      alpha = alpha_use,
      train_frac = 0.5,
      seed = seed,
      quant_method = "grf",
      prop_method = "logistic",
      lower_tau = 0.05,
      upper_tau = 0.95,
      weight_cap = weight_cap_use,
      time_matched = FALSE
    ),
    
    LMEM = lmem_region(
      dat = dat_sample, x_test = dat_test, setting = "cd4",
      alpha = alpha_use, n_sims = 10000L, seed = seed
    )
  )
  
  eval_mat <- evaluate_interval_region(res, y_true)
  
  summarize_cd4_leaveout(
    eval_mat = eval_mat,
    prediction_type = prediction_type
  )
}