#' Run one simulation replicate for conditional coverage evaluation
#'
#' @description
#' This function runs one simulation replicate for the conditional-coverage
#' experiments used in Figures S.1--S.3. It generates one clustered dataset,
#' constructs a fixed grid of test covariates with varying X2, applies localized
#' HCP, and evaluates coverage and prediction-region length at each test point.
#'
#' The test covariates are fixed as
#' \deqn{X = (0.6, x_2, 0, 0),}
#' where \code{x_2} ranges over \code{x2_grid}. The true test outcomes are
#' sampled from the data-generating model conditional on these covariates.
#'
#' @param n Number of subjects in the generated dataset.
#' @param scenario Scenario name. Must be one of \code{"Homo"},
#'   \code{"Heter"}, \code{"Asym"}, \code{"Bimo"}, or \code{"Bimo-fix"}.
#' @param n_local_groups Number of localized groups. Use \code{1} for the
#'   non-localized version.
#' @param alpha Miscoverage level.
#' @param n_grid Number of grid points used for candidate response values.
#' @param x2_grid Numeric grid of X2 values used for conditional evaluation.
#' @param theta_mean Mean vector of the random-effects distribution.
#' @param theta_sigma Covariance matrix of the random-effects distribution.
#' @param x_mean Mean vector of the covariate distribution.
#' @param x_sigma Covariance matrix of the covariate distribution.
#' @param beta Logistic regression coefficients for the missingness mechanism.
#' @param train_frac Fraction of subjects assigned to training in each split.
#' @param S Number of repeated subject-level splits.
#' @param B Number of repeated calibration subsamples per split.
#' @param n_test Dummy argument kept for consistency with other replicate runners.
#' @param weight_cap Optional upper bound for inverse-propensity weights.
#' @param seed Optional random seed.
#'
#' @return A numeric matrix with one row per value in \code{x2_grid} and two
#' columns:
#' \describe{
#'   \item{\code{covered}}{Coverage indicator at the corresponding test point.}
#'   \item{\code{length}}{Prediction-region length at the corresponding test point.}
#' }
#'
#' @export
run_one_replicate_conditional <- function(
    n,
    scenario,
    n_local_groups,
    alpha = 0.1,
    n_grid = 200,
    x2_grid = seq(-3, 3, length.out = 30),
    theta_mean = rep(2, 4),
    theta_sigma = diag(1, 4),
    x_mean = rep(0, 4),
    x_sigma = diag(1, 4),
    beta = c(3, 0, 2, 2, 2),
    train_frac = 0.3,
    S = 2,
    B = 2,
    n_test = 0,
    weight_cap = 30,
    seed = NULL
) {
  # ---------------------------------------------------------------------------
  # Step 0: Basic input checks
  # ---------------------------------------------------------------------------
  valid_scenarios <- c("Homo", "Heter", "Asym", "Bimo", "Bimo-fix")
  if (!scenario %in% valid_scenarios) {
    stop("scenario must be one of: ", paste(valid_scenarios, collapse = ", "))
  }
  
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 2 || n != floor(n)) {
    stop("n must be an integer >= 2.")
  }
  
  if (!is.numeric(n_local_groups) || length(n_local_groups) != 1L ||
      !is.finite(n_local_groups) || n_local_groups < 1 ||
      n_local_groups != floor(n_local_groups)) {
    stop("n_local_groups must be a positive integer.")
  }
  
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number in (0,1).")
  }
  
  if (!is.numeric(n_grid) || length(n_grid) != 1L || !is.finite(n_grid) ||
      n_grid < 2 || n_grid != floor(n_grid)) {
    stop("n_grid must be an integer >= 2.")
  }
  
  if (!is.numeric(x2_grid) || length(x2_grid) < 2L ||
      anyNA(x2_grid) || any(!is.finite(x2_grid))) {
    stop("x2_grid must be a numeric vector with at least two finite values.")
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
  
  if (!is.numeric(weight_cap) || length(weight_cap) != 1L ||
      is.na(weight_cap) || weight_cap < 1) {
    stop("weight_cap must be a single numeric value >= 1, or Inf.")
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # ---------------------------------------------------------------------------
  # Step 1: Generate one clustered training dataset
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
  
  x_cols <- c("X1", "X2", "X3", "X4")
  
  # ---------------------------------------------------------------------------
  # Step 2: Construct fixed test covariates for conditional evaluation
  # ---------------------------------------------------------------------------
  x_test <- as.matrix(data.frame(
    X1 = rep(0.6, length(x2_grid)),
    X2 = x2_grid,
    X3 = rep(0, length(x2_grid)),
    X4 = rep(0, length(x2_grid))
  ))
  
  colnames(x_test) <- x_cols
  
  # ---------------------------------------------------------------------------
  # Step 3: Generate true test outcomes conditional on fixed X_test
  # ---------------------------------------------------------------------------
  y_true <- generate_Y_given_X(
    X = x_test,
    scenario = scenario,
    theta_mean = theta_mean,
    theta_sigma = theta_sigma
  )
  
  y_true <- as.numeric(y_true)
  
  # ---------------------------------------------------------------------------
  # Step 4: Construct candidate response grid from the training sample
  # ---------------------------------------------------------------------------
  y_grid <- seq(
    from = min(dat$Y, na.rm = TRUE),
    to   = max(dat$Y, na.rm = TRUE),
    length.out = n_grid
  )
  
  # ---------------------------------------------------------------------------
  # Step 5: Simulation-specific tuning choices
  # ---------------------------------------------------------------------------
  b_use <- get_b(n, scenario)
  dens_taus_use <- (1:(2^b_use - 1)) / (2^b_use)
  
  # ---------------------------------------------------------------------------
  # Step 6: Apply localized HCP
  # ---------------------------------------------------------------------------
  res_hcp <- hcp_localized_region(
    dat = dat,
    id_col = "id",
    y_col = "Y",
    delta_col = "delta",
    x_cols = x_cols,
    x_test = x_test,
    y_grid = y_grid,
    alpha = alpha,
    train_frac = train_frac,
    S = S,
    B = B,
    dens_method = "rq",
    dens_taus = dens_taus_use,
    dens_h = NULL,
    dens_h_scenario = scenario,
    enforce_monotone = FALSE,
    tail_decay = TRUE,
    prop_method = "logistic",
    weight_cap = weight_cap,
    n_local_groups = n_local_groups,
    t_num = 500L,
    min_local_size = 2L,
    fallback_to_global = FALSE,
    seed = seed
  )
  
  # ---------------------------------------------------------------------------
  # Step 7: Evaluate pointwise conditional coverage and length
  # ---------------------------------------------------------------------------
  if (scenario %in% c("Bimo", "Bimo-fix")) {
    eval_hcp <- evaluate_bimodal_region(res_hcp, y_true)
  } else {
    eval_hcp <- evaluate_interval_region(res_hcp, y_true)
  }
  
  colnames(eval_hcp) <- c("covered", "length")
  eval_hcp
}