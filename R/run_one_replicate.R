#' Run one simulation replicate for the benchmark comparison
#'
#' @description
#' This function runs one full simulation replicate under a specified sample size
#' and scenario. It generates one clustered dataset, holds out one new subject
#' as the test subject, applies all benchmark methods to the remaining sample,
#' evaluates coverage and prediction-region length, and returns one row of
#' replicate-level summary results.
#'
#' The function is intended for simulation studies such as the reproduction of
#' Table 1. It assumes that the required method functions and helper functions
#' have already been sourced into the current R session.
#'
#' @param n Number of subjects in the generated dataset.
#' @param scenario Scenario name. Must be one of \code{"Homo"},
#'   \code{"Heter"}, \code{"Asym"}, \code{"Bimo"}, or \code{"Bimo-fix"}.
#' @param alpha Miscoverage level in \eqn{(0,1)}.
#' @param n_test Number of held-out test subjects. The current implementation
#'   assumes \code{n_test = 1}.
#' @param n_grid Number of grid points used to construct the candidate response grid.
#' @param theta_mean Mean vector of the random-effects distribution.
#' @param theta_sigma Covariance matrix of the random-effects distribution.
#' @param x_mean Mean vector of the covariate distribution.
#' @param x_sigma Covariance matrix of the covariate distribution.
#' @param beta Logistic regression coefficients for the missingness mechanism.
#' @param weight_cap Optional truncation level for inverse-propensity weights.
#'   If \code{NULL}, the default rule \code{get_max_wt(n)} is used.
#' @param seed Optional random seed for reproducibility.
#'
#' @return A named numeric vector containing:
#' \describe{
#'   \item{\code{HCP_cov}}{Coverage of HCP in this replicate.}
#'   \item{\code{HCP_len}}{Prediction-region length of HCP in this replicate.}
#'   \item{\code{DWR_cov}}{Coverage of DWR in this replicate.}
#'   \item{\code{DWR_len}}{Prediction-region length of DWR in this replicate.}
#'   \item{\code{LC_cov}}{Coverage of LC in this replicate.}
#'   \item{\code{LC_len}}{Prediction-region length of LC in this replicate.}
#'   \item{\code{LMEM_cov}}{Coverage of LMEM in this replicate.}
#'   \item{\code{LMEM_len}}{Prediction-region length of LMEM in this replicate.}
#'   \item{\code{Oracle_cov}}{Coverage of Oracle in this replicate.}
#'   \item{\code{Oracle_len}}{Prediction-region length of Oracle in this replicate.}
#' }
#'
#' @export
run_one_replicate <- function(
    n,
    scenario,
    alpha = 0.1,
    n_test = 1,
    n_grid = 200,
    theta_mean = rep(2, 4),
    theta_sigma = diag(1, 4),
    x_mean = rep(0, 4),
    x_sigma = diag(1, 4),
    beta = c(3, 0, 2, 2, 2),
    weight_cap = NULL,
    seed = NULL
) {
  # ---------------------------------------------------------------------------
  # Step 0: Basic input checks
  # ---------------------------------------------------------------------------
  valid_scenarios <- c("Homo", "Heter", "Asym", "Bimo", "Bimo-fix")
  if (!scenario %in% valid_scenarios) {
    stop("scenario must be one of: ", paste(valid_scenarios, collapse = ", "))
  }
  
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 2 || n != floor(n)) {
    stop("n must be an integer >= 2.")
  }
  
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number in (0,1).")
  }
  
  if (!is.numeric(n_test) || length(n_test) != 1L || !is.finite(n_test) ||
      n_test < 1 || n_test != floor(n_test)) {
    stop("n_test must be a positive integer.")
  }
  
  if (n_test != 1L) {
    stop("The current implementation assumes n_test = 1.")
  }
  
  if (!is.numeric(n_grid) || length(n_grid) != 1L || !is.finite(n_grid) ||
      n_grid < 2 || n_grid != floor(n_grid)) {
    stop("n_grid must be an integer >= 2.")
  }
  
  if (!is.null(weight_cap)) {
    if (!is.numeric(weight_cap) || length(weight_cap) != 1L ||
        !is.finite(weight_cap) || weight_cap < 1) {
      stop("weight_cap must be NULL or a single numeric value >= 1.")
    }
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # ---------------------------------------------------------------------------
  # Step 1: Generate one clustered dataset under the specified DGP
  # ---------------------------------------------------------------------------
  dat <- generate_simulation_data(
    n = n,
    m = 5,
    scenario = scenario,
    theta_mean = theta_mean,
    theta_sigma = theta_sigma,
    x_mean = x_mean,
    x_sigma = x_sigma,
    beta = beta
  )
  
  # ---------------------------------------------------------------------------
  # Step 2: Hold out one new subject as the test subject
  # ---------------------------------------------------------------------------
  test_id <- sample(unique(dat$id), size = n_test, replace = FALSE)
  
  dat_test <- dat[dat$id %in% test_id, , drop = FALSE]
  dat_sample <- dat[!dat$id %in% test_id, , drop = FALSE]
  
  y_true <- dat_test$Y
  x_test <- as.matrix(dat_test[, c("X1", "X2", "X3", "X4"), drop = FALSE])
  
  # ---------------------------------------------------------------------------
  # Step 3: Construct the candidate response grid from the remaining sample
  # ---------------------------------------------------------------------------
  y_grid <- seq(
    from = min(dat_sample$Y, na.rm = TRUE),
    to   = max(dat_sample$Y, na.rm = TRUE),
    length.out = n_grid
  )
  
  # ---------------------------------------------------------------------------
  # Step 4: Simulation-specific tuning choices for HCP
  # ---------------------------------------------------------------------------
  # These tuning rules are used in the simulation study to preserve the intended
  # density structure, especially in the bimodal settings.
  S_use <- get_S(n, scenario)
  b_use <- get_b(n, scenario)
  dens_taus_use <- (1:(2^b_use - 1)) / (2^b_use)
  
  # Weight truncation used for numerical stability in inverse-propensity weighting.
  # If weight_cap is supplied by the caller, use it; otherwise use the default rule.
  weight_cap_use <- if (is.null(weight_cap)) get_max_wt(n) else weight_cap
  
  # ---------------------------------------------------------------------------
  # Step 5: Apply HCP
  # ---------------------------------------------------------------------------
  # Important:
  # dens_h = NULL means hcp_region() will recompute the bandwidth inside each
  # subject-level training split using the observed outcomes in that split.
  # dens_h_scenario = scenario preserves the scenario-specific bandwidth rule
  # used in the simulation study, especially for Bimo / Bimo-fix.
  res_hcp <- hcp_region(
    dat = dat_sample,
    id_col = "id",
    y_col = "Y",
    delta_col = "delta",
    x_cols = c("X1", "X2", "X3", "X4"),
    x_test = x_test,
    y_grid = y_grid,
    alpha = alpha,
    train_frac = 0.3,
    S = S_use,
    B = 5,
    dens_method = "rq",
    dens_taus = dens_taus_use,
    dens_h = NULL,
    dens_h_scenario = scenario,
    prop_method = "logistic",
    weight_cap = weight_cap_use,
    seed = seed
  )
  
  if (scenario %in% c("Bimo", "Bimo-fix")) {
    eval_hcp <- evaluate_bimodal_region(res_hcp, y_true)
  } else {
    eval_hcp <- evaluate_interval_region(res_hcp, y_true)
  }
  
  # ---------------------------------------------------------------------------
  # Step 6: Apply DWR
  # ---------------------------------------------------------------------------
  res_dwr <- dwr_region(
    dat = dat_sample,
    id_col = "id",
    y_col = "Y",
    delta_col = "delta",
    x_cols = c("X1", "X2", "X3", "X4"),
    x_test = x_test,
    y_grid = y_grid,
    alpha = alpha,
    train_frac = 0.5,
    quant_method = "linear",
    seed = seed
  )
  eval_dwr <- evaluate_interval_region(res_dwr, y_true)
  
  # ---------------------------------------------------------------------------
  # Step 7: Apply LC
  # ---------------------------------------------------------------------------
  res_lc <- lc_region(
    dat = dat_sample,
    id_col = "id",
    y_col = "Y",
    delta_col = "delta",
    x_cols = c("X1", "X2", "X3", "X4"),
    x_test = x_test,
    y_grid = y_grid,
    alpha = alpha,
    train_frac = 0.5,
    quant_method = "linear",
    prop_method = "logistic",
    seed = seed
  )
  eval_lc <- evaluate_interval_region(res_lc, y_true)
  
  # ---------------------------------------------------------------------------
  # Step 8: Apply LMEM
  # ---------------------------------------------------------------------------
  res_lmem <- lmem_region(
    dat = dat_sample,
    id_col = "id",
    y_col = "Y",
    delta_col = "delta",
    x_cols = c("X1", "X2", "X3", "X4"),
    x_test = x_test,
    y_grid = y_grid,
    alpha = alpha,
    fixed_formula = "X1 + X2 + X3 + X4",
    random_formula = "(1 | id)",
    n_sims = 500,
    seed = seed
  )
  eval_lmem <- evaluate_interval_region(res_lmem, y_true)
  
  # ---------------------------------------------------------------------------
  # Step 9: Apply Oracle
  # ---------------------------------------------------------------------------
  res_oracle <- oracle_region(
    x_test = x_test,
    y_grid = y_grid,
    scenario = scenario,
    alpha = alpha,
    B_mc = 5000,
    theta_mean = 2,
    theta_sd = 1,
    seed = seed
  )
  eval_oracle <- evaluate_interval_region(res_oracle, y_true)
  
  # ---------------------------------------------------------------------------
  # Step 10: Return replicate-level summary results
  # ---------------------------------------------------------------------------
  c(
    HCP_cov    = mean(eval_hcp[, "covered"], na.rm = TRUE),
    HCP_len    = mean(eval_hcp[, "length"],  na.rm = TRUE),
    DWR_cov    = mean(eval_dwr[, "covered"], na.rm = TRUE),
    DWR_len    = mean(eval_dwr[, "length"],  na.rm = TRUE),
    LC_cov     = mean(eval_lc[, "covered"],  na.rm = TRUE),
    LC_len     = mean(eval_lc[, "length"],   na.rm = TRUE),
    LMEM_cov   = mean(eval_lmem[, "covered"], na.rm = TRUE),
    LMEM_len   = mean(eval_lmem[, "length"],  na.rm = TRUE),
    Oracle_cov = mean(eval_oracle[, "covered"], na.rm = TRUE),
    Oracle_len = mean(eval_oracle[, "length"],  na.rm = TRUE)
  )
}