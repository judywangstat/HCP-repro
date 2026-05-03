# ============================================================
# run_one_replicate_sensitivity.R
# HCP sensitivity simulation helper
# ============================================================

#' Run one HCP sensitivity-analysis replicate
#'
#' @description
#' This function runs one simulation replicate for the sensitivity analyses
#' in Figures S.8 and S.9. It generates one clustered dataset, holds out one
#' subject as the test subject, applies HCP with user-specified values of
#' \code{B} and \code{S}, and returns replicate-level marginal coverage and
#' prediction-region length.
#'
#' This helper is used to study the sensitivity of HCP to the number of
#' calibration subsamples \code{B} and the number of repeated data splits
#' \code{S}. Other tuning choices are kept consistent with the main simulation
#' study.
#'
#' @param n Number of subjects in the generated dataset.
#' @param scenario Simulation scenario. One of \code{"Homo"}, \code{"Heter"},
#'   \code{"Asym"}, \code{"Bimo"}, or \code{"Bimo-fix"}.
#' @param B Number of repeated calibration subsamples within each split.
#' @param S Number of repeated subject-level training/calibration splits.
#' @param alpha Miscoverage level.
#' @param n_test Number of held-out test subjects. The current implementation
#'   assumes \code{n_test = 1}.
#' @param n_grid Number of grid points used to construct the candidate response
#'   grid.
#' @param theta_mean Mean vector of the random-effects distribution.
#' @param theta_sigma Covariance matrix of the random-effects distribution.
#' @param x_mean Mean vector of the covariate distribution.
#' @param x_sigma Covariance matrix of the covariate distribution.
#' @param beta Logistic regression coefficients for the missingness mechanism.
#' @param weight_cap Upper truncation level for inverse-propensity weights.
#' @param seed Optional random seed.
#'
#' @return A named numeric vector with entries:
#' \describe{
#'   \item{\code{coverage}}{Marginal coverage in this replicate.}
#'   \item{\code{length}}{Average prediction-region length in this replicate.}
#' }
#'
#' @export
run_one_replicate_sensitivity <- function(
    n,
    scenario = "Heter",
    B = 5,
    S = 1,
    alpha = 0.1,
    n_test = 1,
    n_grid = 200,
    theta_mean = rep(2, 4),
    theta_sigma = diag(1, 4),
    x_mean = rep(0, 4),
    x_sigma = diag(1, 4),
    beta = c(3, 0, 2, 2, 2),
    weight_cap = 30,
    seed = NULL
) {
  # ----------------------------------------------------------
  # Step 0: Basic input checks
  # ----------------------------------------------------------
  valid_scenarios <- c("Homo", "Heter", "Asym", "Bimo", "Bimo-fix")
  if (!scenario %in% valid_scenarios) {
    stop("scenario must be one of: ", paste(valid_scenarios, collapse = ", "))
  }
  
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 2 || n != floor(n)) {
    stop("n must be an integer >= 2.")
  }
  
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) ||
      B < 1 || B != floor(B)) {
    stop("B must be a positive integer.")
  }
  
  if (!is.numeric(S) || length(S) != 1L || !is.finite(S) ||
      S < 1 || S != floor(S)) {
    stop("S must be a positive integer.")
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
  
  if (!is.numeric(weight_cap) || length(weight_cap) != 1L ||
      !is.finite(weight_cap) || weight_cap < 1) {
    stop("weight_cap must be a single finite numeric value >= 1.")
  }
  
  B <- as.integer(B)
  S <- as.integer(S)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # ----------------------------------------------------------
  # Step 1: Generate clustered simulation data
  # ----------------------------------------------------------
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
  
  # ----------------------------------------------------------
  # Step 2: Hold out one subject as the test subject
  # ----------------------------------------------------------
  test_id <- sample(unique(dat$id), size = n_test, replace = FALSE)
  
  dat_test <- dat[dat$id %in% test_id, , drop = FALSE]
  dat_sample <- dat[!dat$id %in% test_id, , drop = FALSE]
  
  y_true <- dat_test$Y
  x_test <- as.matrix(dat_test[, c("X1", "X2", "X3", "X4"), drop = FALSE])
  
  # ----------------------------------------------------------
  # Step 3: Construct candidate response grid
  # ----------------------------------------------------------
  y_grid <- seq(
    from = min(dat_sample$Y, na.rm = TRUE),
    to   = max(dat_sample$Y, na.rm = TRUE),
    length.out = n_grid
  )
  
  # ----------------------------------------------------------
  # Step 4: Set density tau grid
  # ----------------------------------------------------------
  b_use <- get_b(n, scenario)
  dens_taus_use <- (1:(2^b_use - 1)) / (2^b_use)
  
  # ----------------------------------------------------------
  # Step 5: Apply HCP with user-specified B and S
  # ----------------------------------------------------------
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
    S = S,
    B = B,
    dens_method = "rq",
    dens_taus = dens_taus_use,
    dens_h = NULL,
    dens_h_scenario = scenario,
    prop_method = "logistic",
    weight_cap = weight_cap,
    seed = seed
  )
  
  # ----------------------------------------------------------
  # Step 6: Evaluate coverage and prediction-region length
  # ----------------------------------------------------------
  if (scenario %in% c("Bimo", "Bimo-fix")) {
    eval_hcp <- evaluate_bimodal_region(res_hcp, y_true)
  } else {
    eval_hcp <- evaluate_interval_region(res_hcp, y_true)
  }
  
  # ----------------------------------------------------------
  # Step 7: Return replicate-level summary
  # ----------------------------------------------------------
  c(
    coverage = mean(eval_hcp[, "covered"], na.rm = TRUE),
    length   = mean(eval_hcp[, "length"], na.rm = TRUE)
  )
}