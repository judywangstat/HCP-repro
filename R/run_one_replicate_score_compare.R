#' Run one score-comparison simulation replicate
#'
#' @description
#' Runs one simulation replicate for Table S.5 or Table S.6. The function
#' generates one clustered dataset, holds out one subject as the test subject,
#' applies HCP with both density and residual nonconformity scores, and returns
#' replicate-level coverage and region length.
#'
#' @param n Number of subjects in the generated dataset.
#' @param scenario Simulation scenario. One of "Homo", "Heter", "Asym",
#'   "Bimo", or "Bimo-fix".
#' @param alpha Miscoverage level.
#' @param n_test Number of held-out test subjects. Currently must be 1.
#' @param n_grid Number of grid points for the candidate response grid.
#' @param train_frac Fraction of subjects assigned to training in each split.
#' @param theta_mean Mean vector of the random-effects distribution.
#' @param theta_sigma Covariance matrix of the random-effects distribution.
#' @param x_mean Mean vector of the covariate distribution.
#' @param x_sigma Covariance matrix of the covariate distribution.
#' @param beta Logistic regression coefficients for the missingness mechanism.
#' @param weight_cap Optional inverse-propensity weight cap. If NULL, use
#'   get_max_wt(n).
#' @param seed Optional random seed.
#'
#' @return A named numeric vector with density and residual coverage/length.
#'
#' @export
run_one_replicate_score_compare <- function(
    n,
    scenario,
    alpha = 0.1,
    n_test = 1,
    n_grid = 200,
    train_frac = 0.5,
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
  
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 2 || n != floor(n)) {
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
  
  if (!is.numeric(train_frac) || length(train_frac) != 1L ||
      !is.finite(train_frac) || train_frac <= 0 || train_frac >= 1) {
    stop("train_frac must be a single number in (0,1).")
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
  # Step 2: Hold out one subject as the test subject
  # ---------------------------------------------------------------------------
  test_id <- sample(unique(dat$id), size = n_test, replace = FALSE)
  
  dat_test <- dat[dat$id %in% test_id, , drop = FALSE]
  dat_sample <- dat[!dat$id %in% test_id, , drop = FALSE]
  
  y_true <- dat_test$Y
  x_cols <- c("X1", "X2", "X3", "X4")
  x_test <- as.matrix(dat_test[, x_cols, drop = FALSE])
  
  # ---------------------------------------------------------------------------
  # Step 3: Construct the candidate response grid
  # ---------------------------------------------------------------------------
  y_grid <- seq(
    from = min(dat_sample$Y, na.rm = TRUE),
    to = max(dat_sample$Y, na.rm = TRUE),
    length.out = n_grid
  )
  
  # ---------------------------------------------------------------------------
  # Step 4: Simulation-specific HCP tuning
  # ---------------------------------------------------------------------------
  S_use <- get_S(n, scenario)
  b_use <- get_b(n, scenario)
  dens_taus_use <- (1:(2^b_use - 1)) / (2^b_use)
  
  weight_cap_use <- if (is.null(weight_cap)) get_max_wt(n) else weight_cap
  
  # ---------------------------------------------------------------------------
  # Step 5: Run HCP score comparison
  # ---------------------------------------------------------------------------
  res <- hcp_score_compare(
    dat = dat_sample,
    id_col = "id",
    y_col = "Y",
    delta_col = "delta",
    x_cols = x_cols,
    x_test = x_test,
    y_grid = y_grid,
    alpha = alpha,
    train_frac = train_frac,
    S = S_use,
    B = 5,
    combine_B = "cct",
    combine_S = "cct",
    dens_method = "rq",
    dens_taus = dens_taus_use,
    dens_h = NULL,
    dens_h_scenario = scenario,
    enforce_monotone = FALSE,
    tail_decay = TRUE,
    reg_method = "linear",
    prop_method = "logistic",
    weight_cap = weight_cap_use,
    seed = seed
  )
  
  # ---------------------------------------------------------------------------
  # Step 6: Evaluate coverage and region length
  # ---------------------------------------------------------------------------
  if (scenario %in% c("Bimo", "Bimo-fix")) {
    eval_density <- evaluate_bimodal_region(res$density, y_true)
  } else {
    eval_density <- evaluate_interval_region(res$density, y_true)
  }
  
  eval_residual <- evaluate_interval_region(res$residual, y_true)
  
  # ---------------------------------------------------------------------------
  # Step 7: Return replicate-level summary
  # ---------------------------------------------------------------------------
  c(
    density_cov = mean(eval_density[, "covered"], na.rm = TRUE),
    density_len = mean(eval_density[, "length"], na.rm = TRUE),
    residual_cov = mean(eval_residual[, "covered"], na.rm = TRUE),
    residual_len = mean(eval_residual[, "length"], na.rm = TRUE)
  )
}