# ============================================================
# run_one_replicate_local.R
# One-replicate runner for local coverage heatmap experiments
# ============================================================

#' Run one simulation replicate for local coverage evaluation
#'
#' @description
#' This function runs one simulation replicate for the local-coverage heatmap
#' experiments used in Figures S.4--S.6. It generates one clustered dataset,
#' constructs fixed test covariates on a two-dimensional grid
#' \deqn{X = (0.6, x_2, x_3, 0),}
#' applies localized HCP, evaluates pointwise coverage and prediction-region
#' length, and finally averages these quantities within rectangular cells in the
#' \eqn{(x_2, x_3)} plane.
#'
#' The function is intended for simulation studies. It assumes that the required
#' method functions and helper functions have already been sourced.
#'
#' @param n Number of subjects in the generated clustered dataset.
#' @param scenario Simulation scenario. Must be one of \code{"Homo"},
#'   \code{"Heter"}, \code{"Asym"}, \code{"Bimo"}, or \code{"Bimo-fix"}.
#' @param n_local_groups Number of localized groups. Use \code{1} for the
#'   non-localized version.
#' @param alpha Miscoverage level.
#' @param n_grid Number of grid points used for the candidate response grid.
#' @param x2_grid Numeric grid of \eqn{x_2} values for test covariates.
#' @param x3_grid Numeric grid of \eqn{x_3} values for test covariates.
#' @param x2_breaks Breakpoints used to define local cells along \eqn{x_2}.
#' @param x3_breaks Breakpoints used to define local cells along \eqn{x_3}.
#' @param theta_mean Mean vector of the random-effects distribution.
#' @param theta_sigma Covariance matrix of the random-effects distribution.
#' @param x_mean Mean vector of the covariate distribution.
#' @param x_sigma Covariance matrix of the covariate distribution.
#' @param beta Logistic regression coefficients for the missingness mechanism.
#' @param train_frac Fraction of subjects assigned to training.
#' @param S Number of repeated subject-level splits.
#' @param B Number of repeated calibration subsamples per split.
#' @param weight_cap Upper truncation level for inverse-propensity weights.
#' @param seed Optional random seed.
#'
#' @return A data frame with one row per local cell and columns
#' \code{cell_id}, \code{x2_mid}, \code{x3_mid}, \code{covered}, and
#' \code{length}.
#'
#' @export
run_one_replicate_local <- function(
    n,
    scenario = "Homo",
    n_local_groups,
    alpha = 0.1,
    n_grid = 200,
    x2_grid = seq(-3, 3, length.out = 42),
    x3_grid = seq(-3, 3, length.out = 42),
    x2_breaks = seq(-3, 3, length.out = 7),
    x3_breaks = seq(-3, 3, length.out = 7),
    theta_mean = rep(2, 4),
    theta_sigma = diag(1, 4),
    x_mean = rep(0, 4),
    x_sigma = diag(1, 4),
    beta = c(3, 0, 2, 2, 2),
    train_frac = 0.3,
    S = 1,
    B = 1,
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
  
  if (!is.numeric(x3_grid) || length(x3_grid) < 2L ||
      anyNA(x3_grid) || any(!is.finite(x3_grid))) {
    stop("x3_grid must be a numeric vector with at least two finite values.")
  }
  
  if (!is.numeric(x2_breaks) || length(x2_breaks) < 2L ||
      anyNA(x2_breaks) || any(!is.finite(x2_breaks)) ||
      is.unsorted(x2_breaks, strictly = TRUE)) {
    stop("x2_breaks must be a strictly increasing numeric vector.")
  }
  
  if (!is.numeric(x3_breaks) || length(x3_breaks) < 2L ||
      anyNA(x3_breaks) || any(!is.finite(x3_breaks)) ||
      is.unsorted(x3_breaks, strictly = TRUE)) {
    stop("x3_breaks must be a strictly increasing numeric vector.")
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
  
  x_cols <- c("X1", "X2", "X3", "X4")
  
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
  
  # ---------------------------------------------------------------------------
  # Step 2: Construct fixed two-dimensional test covariates
  # ---------------------------------------------------------------------------
  x_test_df <- expand.grid(
    X1 = 0.6,
    X2 = x2_grid,
    X3 = x3_grid,
    X4 = 0
  )
  
  x_test <- as.matrix(x_test_df[, x_cols, drop = FALSE])
  colnames(x_test) <- x_cols
  
  # ---------------------------------------------------------------------------
  # Step 3: Generate true outcomes conditional on fixed X_test
  # ---------------------------------------------------------------------------
  y_true <- generate_Y_given_X(
    X = x_test,
    scenario = scenario,
    theta_mean = theta_mean,
    theta_sigma = theta_sigma
  )
  y_true <- as.numeric(y_true)
  
  # ---------------------------------------------------------------------------
  # Step 4: Construct candidate response grid from the generated sample
  # ---------------------------------------------------------------------------
  y_grid <- seq(
    from = min(dat$Y, na.rm = TRUE),
    to = max(dat$Y, na.rm = TRUE),
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
  # Step 7: Evaluate pointwise coverage and prediction-region length
  # ---------------------------------------------------------------------------
  if (scenario %in% c("Bimo", "Bimo-fix")) {
    eval_hcp <- evaluate_bimodal_region(res_hcp, y_true)
  } else {
    eval_hcp <- evaluate_interval_region(res_hcp, y_true)
  }
  
  eval_df <- data.frame(
    X2 = x_test_df$X2,
    X3 = x_test_df$X3,
    covered = as.numeric(eval_hcp[, "covered"]),
    length = as.numeric(eval_hcp[, "length"])
  )
  
  # ---------------------------------------------------------------------------
  # Step 8: Average coverage and length within rectangular local cells
  # ---------------------------------------------------------------------------
  out_list <- vector(
    "list",
    (length(x2_breaks) - 1L) * (length(x3_breaks) - 1L)
  )
  
  cell_id <- 1L
  
  for (j in seq_len(length(x3_breaks) - 1L)) {
    for (i in seq_len(length(x2_breaks) - 1L)) {
      
      # Match the old implementation: both endpoints are included.
      idx <- which(
        eval_df$X2 >= x2_breaks[i] &
          eval_df$X2 <= x2_breaks[i + 1L] &
          eval_df$X3 >= x3_breaks[j] &
          eval_df$X3 <= x3_breaks[j + 1L]
      )
      
      out_list[[cell_id]] <- data.frame(
        cell_id = cell_id,
        x2_mid = (x2_breaks[i] + x2_breaks[i + 1L]) / 2,
        x3_mid = (x3_breaks[j] + x3_breaks[j + 1L]) / 2,
        covered = mean(eval_df$covered[idx], na.rm = TRUE),
        length = mean(eval_df$length[idx], na.rm = TRUE)
      )
      
      cell_id <- cell_id + 1L
    }
  }
  
  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}