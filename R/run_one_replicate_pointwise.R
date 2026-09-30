#' Run one simulation replicate for pointwise benchmark results
#'
#' @description
#' This function runs one full simulation replicate under a specified sample size
#' and scenario. It generates one clustered dataset, holds out one new subject
#' as the test subject, applies all benchmark methods to the remaining sample,
#' evaluates pointwise coverage and prediction-region length, and returns one row
#' per test observation.
#'
#' The function is intended for simulation studies such as the reproduction of
#' Figure S.10. In contrast to \code{run_one_replicate()}, which returns
#' replicate-level averages, this function returns pointwise results. When
#' \code{n_test = 1} and each test subject has 5 observations, each replicate
#' returns a \code{5 x 10} matrix.
#'
#' It assumes that the required method functions and helper functions have
#' already been sourced into the current R session.
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
#' @param lmem_n_sims Number of new-subject predictive draws for LMEM.
#' @param cluster_size_design Balanced size 5 or independent sizes 5/50.
#' @param method_names Methods to compute; skipped methods consume no RNG.
#' @param seed Optional random seed for reproducibility.
#'
#' @return A numeric matrix with one row per test observation and columns:
#' \describe{
#'   \item{\code{HCP_cov}}{Pointwise coverage indicator of HCP.}
#'   \item{\code{HCP_len}}{Pointwise prediction-region length of HCP.}
#'   \item{\code{DWR_cov}}{Pointwise coverage indicator of DWR.}
#'   \item{\code{DWR_len}}{Pointwise prediction-region length of DWR.}
#'   \item{\code{LC_cov}}{Pointwise coverage indicator of LC.}
#'   \item{\code{LC_len}}{Pointwise prediction-region length of LC.}
#'   \item{\code{LMEM_cov}}{Pointwise coverage indicator of LMEM.}
#'   \item{\code{LMEM_len}}{Pointwise prediction-region length of LMEM.}
#'   \item{\code{Oracle_cov}}{Pointwise coverage indicator of Oracle.}
#'   \item{\code{Oracle_len}}{Pointwise prediction-region length of Oracle.}
#' }
#'
#' @export
run_one_replicate_pointwise <- function(
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
    lmem_n_sims = 500L,
    seed = NULL,
    return_details = FALSE,
    cluster_size_design = c("balanced", "imbalanced"),
    method_names = simulation_methods()
) {
  # ---------------------------------------------------------------------------
  # Step 0: Basic input checks
  # ---------------------------------------------------------------------------
  cluster_size_design <- match.arg(cluster_size_design)
  if (!length(method_names) || anyDuplicated(method_names) || !all(method_names %in% simulation_methods())) stop("Invalid method_names.")
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
  
  if (!isTRUE(all.equal(as.numeric(theta_mean), rep(2,4))) ||
      !isTRUE(all.equal(as.matrix(theta_sigma), diag(1,4))))
    stop("Oracle requires theta_mean=rep(2,4), theta_sigma=I4 DGP.")

  lmem_n_sims <- as.integer(lmem_n_sims)
  if (!is.finite(lmem_n_sims) || lmem_n_sims < 1L) {
    stop("lmem_n_sims must be a positive integer.")
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # ---------------------------------------------------------------------------
  # Step 1: Generate one clustered dataset under the specified DGP
  # ---------------------------------------------------------------------------
  cluster_sizes <- if (cluster_size_design == "imbalanced") sample(c(5L, 50L), n, replace = TRUE) else rep(5L, n)
  dat <- generate_simulation_data(
    n = n,
    m = cluster_sizes,
    scenario = scenario,
    theta_mean = theta_mean,
    theta_sigma = theta_sigma,
    x_mean = x_mean,
    x_sigma = x_sigma,
    beta = beta,
    covariate_generator = if (cluster_size_design == "imbalanced") generate_imbalanced_covariates else generate_subject_covariates
  )
  
  # ---------------------------------------------------------------------------
  # Step 2: Hold out one new subject as the test subject
  # ---------------------------------------------------------------------------
  test_id <- sample(unique(dat$id), size = n_test, replace = FALSE)
  
  dat_test <- dat[dat$id %in% test_id, , drop = FALSE]
  dat_sample <- dat[!dat$id %in% test_id, , drop = FALSE]
  
  y_true <- dat_test$Y
  x_cols <- c("X1", "X2", "X3", "X4")
  x_test <- as.matrix(dat_test[, x_cols, drop = FALSE])
  
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
  S_use <- get_S(n, scenario)
  b_use <- get_b(n, scenario)
  dens_taus_use <- (1:(2^b_use - 1)) / (2^b_use)
  
  weight_cap_use <- if (is.null(weight_cap)) get_max_wt(n) else weight_cap
  
  # ---------------------------------------------------------------------------
  # Evaluate each method separately; successful peers survive any method failure.
  methods <- list()
  if ("HCP" %in% method_names) methods$HCP <- run_simulation_method("HCP", function(state) hcp_finite_region(
    dat = dat_sample,
    id_col = "id",
    y_col = "Y",
    delta_col = "delta",
    x_cols = x_cols,
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
  ),
    evaluate_hcp_finite, y_true, keep_prediction=return_details)

  if ("DWR" %in% method_names) methods$DWR <- run_simulation_method("DWR", function(state) dwr_region(
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
    seed = seed
  ),
    evaluate_interval_region, y_true, keep_prediction=return_details)

  if ("LC" %in% method_names) methods$LC <- run_simulation_method("LC", function(state) lc_region(
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
    seed = seed
  ),
    evaluate_interval_region, y_true, keep_prediction=return_details)

  if ("LMEM" %in% method_names) methods$LMEM <- run_simulation_method("LMEM", function(state)
    lmem_region(dat_sample, x_test, setting = "simulation", alpha = alpha,
      n_sims = lmem_n_sims, seed = seed, state = state),
    evaluate_interval_region, y_true, keep_prediction = return_details)
  if ("Oracle" %in% method_names) methods$Oracle <- run_simulation_method("Oracle", function(state) {
    state$stage <- "prediction"
    oracle_region(x_test=x_test, scenario=scenario, alpha=alpha)
  }, evaluate_oracle_union, y_true, keep_prediction=return_details)

  out <- do.call(cbind, lapply(methods, `[[`, "evaluation"))
  colnames(out) <- as.vector(rbind(paste0(names(methods), "_cov"), paste0(names(methods), "_len"))); rownames(out) <- NULL
  attr(out, "method_diagnostics") <- lapply(methods, `[[`, "diagnostics")
  attr(out, "oracle_diagnostics") <- methods$Oracle$oracle_diagnostics
  attr(out, "input_metadata") <- list(test_id=test_id, test_rows=nrow(dat_test),
    seed=seed, n=n, scenario=scenario, rng=RNGkind(), cluster_size_design=cluster_size_design, cluster_sizes=cluster_sizes)
  if (return_details) attr(out, "details") <- list(train=dat_sample, test=dat_test,
    y_grid=y_grid, predictions=lapply(methods, `[[`, "prediction"))
  out
}
