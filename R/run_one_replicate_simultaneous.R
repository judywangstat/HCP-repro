#' Run one simulation replicate for simultaneous prediction
#'
#' @description
#' This function runs one full simulation replicate under a specified sample size,
#' scenario, and missingness level. It generates one clustered dataset, holds out
#' one new subject as the test subject, applies the HCP method to the remaining
#' sample, evaluates simultaneous coverage over all observations from the held-out
#' subject, and returns replicate-level simultaneous coverage and average
#' prediction-region length.
#'
#' The function is intended for simulation studies such as the reproduction of
#' Table S.1 and Figure S.7. In contrast to \code{run_one_replicate()}, which
#' reports average pointwise coverage over the held-out subject, this function
#' reports whether all test observations from the held-out subject are covered
#' simultaneously.
#'
#' It assumes that the required method functions and helper functions have
#' already been sourced into the current R session.
#'
#' @param n Number of subjects in the generated dataset.
#' @param scenario Scenario name. Must be one of \code{"Homo"},
#'   \code{"Heter"}, \code{"Asym"}, \code{"Bimo"}, or \code{"Bimo-fix"}.
#' @param missing Missingness level. Must be \code{"20"} or \code{"50"}.
#' @param alpha Target subject-level miscoverage level. The internal pointwise
#'   level is set to \code{alpha / m}, where \code{m = 5} in the simulation.
#' @param n_test Number of held-out test subjects. The current implementation
#'   assumes \code{n_test = 1}.
#' @param n_grid Number of grid points used to construct the candidate response
#'   grid.
#' @param theta_mean Mean vector of the random-effects distribution.
#' @param theta_sigma Covariance matrix of the random-effects distribution.
#' @param x_mean Mean vector of the covariate distribution.
#' @param x_sigma Covariance matrix of the covariate distribution.
#' @param beta Optional logistic regression coefficients for the missingness
#'   mechanism. If \code{NULL}, the default value is chosen according to
#'   \code{missing}.
#' @param weight_cap Optional truncation level for inverse-propensity weights.
#'   If \code{NULL}, the default value is chosen according to \code{missing}.
#' @param train_frac Optional fraction of subjects assigned to training in each
#'   split. If \code{NULL}, the default value is \code{0.3}.
#' @param seed Optional random seed for reproducibility.
#'
#' @return A named numeric vector containing:
#' \describe{
#'   \item{\code{HCP_cov}}{Simultaneous coverage indicator for the held-out
#'   subject. This equals 1 if all test observations are covered and 0 otherwise.}
#'   \item{\code{HCP_len}}{Average prediction-region length across the held-out
#'   subject's test observations.}
#' }
#'
#' @export
run_one_replicate_simultaneous <- function(
    n,
    scenario,
    missing = c("20", "50"),
    alpha = 0.1,
    n_test = 1,
    n_grid = 200,
    theta_mean = rep(2, 4),
    theta_sigma = diag(1, 4),
    x_mean = rep(0, 4),
    x_sigma = diag(1, 4),
    beta = NULL,
    weight_cap = NULL,
    train_frac = NULL,
    seed = NULL
) {
  # ---------------------------------------------------------------------------
  # Step 0: Basic input checks
  # ---------------------------------------------------------------------------
  valid_scenarios <- c("Homo", "Heter", "Asym", "Bimo", "Bimo-fix")
  if (!scenario %in% valid_scenarios) {
    stop("scenario must be one of: ", paste(valid_scenarios, collapse = ", "))
  }
  
  missing <- match.arg(missing)
  
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
  
  if (!is.null(train_frac)) {
    if (!is.numeric(train_frac) || length(train_frac) != 1L ||
        !is.finite(train_frac) || train_frac <= 0 || train_frac >= 1) {
      stop("train_frac must be NULL or a single number in (0,1).")
    }
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # ---------------------------------------------------------------------------
  # Step 1: Missingness-specific simulation settings
  # ---------------------------------------------------------------------------
  if (is.null(beta)) {
    beta <- switch(
      missing,
      "20" = c(3, 0, 2, 2, 2),
      "50" = c(0, 0, 2, 2, 2)
    )
  }
  
  if (is.null(train_frac)) {
    train_frac <- 0.3
  }
  
  if (is.null(weight_cap)) {
    weight_cap <- 30
  }
  
  # ---------------------------------------------------------------------------
  # Step 2: Generate one clustered dataset under the specified DGP
  # ---------------------------------------------------------------------------
  m <- 5L
  
  dat <- generate_simulation_data(
    n = n,
    m = m,
    scenario = scenario,
    theta_mean = theta_mean,
    theta_sigma = theta_sigma,
    x_mean = x_mean,
    x_sigma = x_sigma,
    beta = beta
  )
  
  # ---------------------------------------------------------------------------
  # Step 3: Hold out one new subject as the test subject
  # ---------------------------------------------------------------------------
  test_id <- sample(unique(dat$id), size = n_test, replace = FALSE)
  
  dat_test <- dat[dat$id %in% test_id, , drop = FALSE]
  dat_sample <- dat[!dat$id %in% test_id, , drop = FALSE]
  
  y_true <- dat_test$Y
  x_cols <- c("X1", "X2", "X3", "X4")
  x_test <- as.matrix(dat_test[, x_cols, drop = FALSE])
  
  # ---------------------------------------------------------------------------
  # Step 4: Construct the candidate response grid from the remaining sample
  # ---------------------------------------------------------------------------
  y_grid <- seq(
    from = min(dat_sample$Y, na.rm = TRUE),
    to   = max(dat_sample$Y, na.rm = TRUE),
    length.out = n_grid
  )
  
  # ---------------------------------------------------------------------------
  # Step 5: Simulation-specific tuning choices for HCP
  # ---------------------------------------------------------------------------
  alpha_use <- alpha / m
  
  # Use the shared simulation tuning rules for the split number and tau-grid size.
  S_use <- get_S(n, scenario)
  b_use <- get_b(n, scenario)
  dens_taus_use <- (1:(2^b_use - 1)) / (2^b_use)
  
  # ---------------------------------------------------------------------------
  # Step 6: Apply HCP
  # ---------------------------------------------------------------------------
  res_hcp <- hcp_region(
    dat = dat_sample,
    id_col = "id",
    y_col = "Y",
    delta_col = "delta",
    x_cols = x_cols,
    x_test = x_test,
    y_grid = y_grid,
    alpha = alpha_use,
    train_frac = train_frac,
    S = S_use,
    B = 5,
    dens_method = "rq",
    dens_taus = dens_taus_use,
    dens_h = NULL,
    dens_h_scenario = scenario,
    prop_method = "logistic",
    weight_cap = weight_cap,
    seed = seed
  )
  
  # ---------------------------------------------------------------------------
  # Step 7: Evaluate simultaneous coverage and average length
  # ---------------------------------------------------------------------------
  if (scenario %in% c("Bimo", "Bimo-fix")) {
    eval_hcp <- evaluate_bimodal_region(res_hcp, y_true)
  } else {
    eval_hcp <- evaluate_interval_region(res_hcp, y_true)
  }
  
  c(
    HCP_cov = as.numeric(all(eval_hcp[, "covered"] == 1)),
    HCP_len = mean(eval_hcp[, "length"], na.rm = TRUE)
  )
}