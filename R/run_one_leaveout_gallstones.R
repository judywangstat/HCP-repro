# ============================================================
# run_one_leaveout_gallstones.R
# Leave-one-subject-out evaluation for the gallstones data
# ============================================================

#' Run one leave-one-subject-out evaluation for the gallstones data
#'
#' @description
#' Holds out one subject from the gallstones data, constructs prediction regions
#' using HCP, DWR, LC, and LMEM on the remaining subjects, and evaluates either
#' pointwise or simultaneous coverage on the held-out subject.
#'
#' The outcome column \code{Y} is used for model fitting with its observed
#' missingness pattern. The column \code{Y_eval} is used only as the evaluation
#' target for the held-out subject.
#'
#' @param test_id Subject identifier to hold out.
#' @param dat Gallstones data frame returned by \code{impute_gallstones_outcomes()}.
#' @param prediction_type Either \code{"pointwise"} or \code{"simultaneous"}.
#' @param alpha Miscoverage level. For simultaneous prediction, this should
#'   usually be \code{0.1 / 4}.
#' @param n_grid Number of response-grid points.
#' @param x_cols Covariate column names.
#' @param seed Optional random seed.
#'
#' @return A named numeric vector containing coverage and average length for
#'   each method.
#'
#' @export
run_one_leaveout_gallstones <- function(
    test_id,
    dat,
    prediction_type = c("pointwise", "simultaneous"),
    alpha = 0.1,
    n_grid = 200,
    x_cols = c("T1", "T2", "T3", "Treat"),
    seed = NULL
) {
  # ---------------------------------------------------------------------------
  # Step 0: Basic input checks
  # ---------------------------------------------------------------------------
  prediction_type <- match.arg(prediction_type)
  
  if (!is.data.frame(dat)) {
    stop("dat must be a data.frame.")
  }
  
  required_cols <- c("id", "Y", "Y_eval", "delta", x_cols)
  if (!all(required_cols %in% names(dat))) {
    stop(
      "dat must contain columns: ",
      paste(required_cols, collapse = ", ")
    )
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
  
  # ---------------------------------------------------------------------------
  # Step 1: Split held-out subject and training/calibration subjects
  # ---------------------------------------------------------------------------
  dat_test <- dat[dat$id == test_id, , drop = FALSE]
  dat_sample <- dat[dat$id != test_id, , drop = FALSE]
  
  y_true <- as.numeric(dat_test$Y_eval)
  x_test <- as.matrix(dat_test[, x_cols, drop = FALSE])
  storage.mode(x_test) <- "numeric"
  colnames(x_test) <- x_cols
  
  y_grid <- make_gallstones_y_grid(
    dat_sample = dat_sample,
    y_col = "Y",
    n_grid = n_grid
  )
  
  # ---------------------------------------------------------------------------
  # Helper: evaluate pointwise or simultaneous coverage
  # ---------------------------------------------------------------------------
  summarize_eval <- function(eval_mat) {
    covered <- as.numeric(eval_mat[, "covered"])
    lengths <- as.numeric(eval_mat[, "length"])
    
    if (prediction_type == "simultaneous") {
      c(
        cov = as.numeric(all(covered == 1)),
        len = mean(lengths, na.rm = TRUE)
      )
    } else {
      c(
        cov = mean(covered, na.rm = TRUE),
        len = mean(lengths, na.rm = TRUE)
      )
    }
  }
  
  # ---------------------------------------------------------------------------
  # Step 2: HCP
  # ---------------------------------------------------------------------------
  res_hcp <- hcp_region(
    dat = dat_sample,
    id_col = "id",
    y_col = "Y",
    delta_col = "delta",
    x_cols = x_cols,
    x_test = x_test,
    y_grid = y_grid,
    alpha = alpha,
    train_frac = if (prediction_type == "pointwise") 0.3 else 0.4,
    S = 12,
    B = 5,
    dens_method = "rq",
    dens_taus = (1:(2^4 - 1)) / (2^4),
    dens_h = NULL,
    dens_h_scenario = NULL,
    prop_method = "logistic",
    weight_cap = 50,
    seed = seed
  )
  
  eval_hcp <- evaluate_interval_region(res_hcp, y_true)
  out_hcp <- summarize_eval(eval_hcp)
  
  # ---------------------------------------------------------------------------
  # Step 3: DWR
  # ---------------------------------------------------------------------------
  res_dwr <- dwr_region_realdata(
    dat = dat_sample,
    id_col = "id",
    y_col = "Y",
    delta_col = "delta",
    x_cols = x_cols,
    x_test = x_test,
    y_grid = y_grid,
    alpha = alpha,
    train_frac = 0.5,
    B = 1,
    score_type = "residual",
    reg_method = "linear",
    seed = seed
  )
  
  eval_dwr <- evaluate_interval_region(res_dwr, y_true)
  out_dwr <- summarize_eval(eval_dwr)
  
  # ---------------------------------------------------------------------------
  # Step 4: LC
  # ---------------------------------------------------------------------------
  res_lc <- lc_region_realdata(
    dat = dat_sample,
    id_col = "id",
    y_col = "Y",
    delta_col = "delta",
    x_cols = x_cols,
    x_test = x_test,
    y_grid = y_grid,
    alpha = alpha,
    train_frac = 0.5,
    quant_method = "linear",
    prop_method = "logistic",
    weight_cap = 50,
    time_matched = FALSE,
#    time_col = "time",
#    test_time = dat_test$time,
    seed = seed
  )
  
  eval_lc <- evaluate_interval_region(res_lc, y_true)
  out_lc <- summarize_eval(eval_lc)
  
  # ---------------------------------------------------------------------------
  # Step 5: LMEM
  # ---------------------------------------------------------------------------
  res_lmem <- lmem_region(
    dat = dat_sample, x_test = dat_test, setting = "gallstones",
    alpha = alpha, n_sims = 10000L, seed = seed
  )
  
  eval_lmem <- evaluate_interval_region(res_lmem, y_true)
  out_lmem <- summarize_eval(eval_lmem)
  
  # ---------------------------------------------------------------------------
  # Step 6: Return one-row summary
  # ---------------------------------------------------------------------------
  c(
    HCP_cov = out_hcp["cov"],
    HCP_len = out_hcp["len"],
    DWR_cov = out_dwr["cov"],
    DWR_len = out_dwr["len"],
    LC_cov = out_lc["cov"],
    LC_len = out_lc["len"],
    LMEM_cov = out_lmem["cov"],
    LMEM_len = out_lmem["len"]
  )
}