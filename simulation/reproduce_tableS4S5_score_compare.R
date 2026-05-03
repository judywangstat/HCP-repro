# ============================================================
# Reproduce Table S.4 and Table S.5
# Score-comparison simulation study
# HCPclust-repro
#
# This script reproduces the score-comparison tables in the
# supplementary materials.
#
# Table S.4: 20% missing responses.
# Table S.5: 50% missing responses.
#
# The comparison is between:
#   1. HCP with density score
#   2. HCP with residual score
#
# Output:
#   results/tableS4_score_compare_missing20.csv
#   results/tableS5_score_compare_missing50.csv
#
# Note:
#   Update repo_dir below to the local path of this repository
#   before running the script.
# ============================================================

suppressPackageStartupMessages({
  library(doParallel)
  library(doRNG)
  library(foreach)
})

# ------------------------------------------------------------
# User-adjustable paths
# ------------------------------------------------------------
repo_dir <- normalizePath("~/Desktop/HCPclust-repro", mustWork = TRUE)
r_dir <- file.path(repo_dir, "R")
results_dir <- file.path(repo_dir, "results")

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# ------------------------------------------------------------
# Load required function files
# ------------------------------------------------------------
source(file.path(r_dir, "simulation_dgp.R"), chdir = TRUE)
source(file.path(r_dir, "tuning_helpers.R"), chdir = TRUE)
source(file.path(r_dir, "evaluation_helpers.R"), chdir = TRUE)
source(file.path(r_dir, "fit_propensity_model.R"), chdir = TRUE)
source(file.path(r_dir, "fit_cond_density_qp.R"), chdir = TRUE)
source(file.path(r_dir, "hcp_score_compare.R"), chdir = TRUE)
source(file.path(r_dir, "run_one_replicate_score_compare.R"), chdir = TRUE)

#' Run the score-comparison simulation study in parallel
#'
#' @description
#' This helper runs the HCP score-comparison simulation study for one missingness
#' setting. For each sample size and simulation scenario, it compares HCP using
#' the density score with HCP using the residual score, then aggregates
#' replicate-level coverage and prediction-region length.
#'
#' @param table_name Character label used in progress messages.
#' @param num_simulations_run Number of simulation replicates per setting.
#' @param num_cores Number of parallel workers.
#' @param sample_sizes Sample sizes to evaluate.
#' @param scenarios Simulation scenarios to evaluate.
#' @param alpha Miscoverage level.
#' @param n_test Number of held-out test subjects per replicate.
#' @param n_grid Number of response-grid points.
#' @param theta_mean Mean vector of subject-specific random effects.
#' @param theta_sigma Covariance matrix of subject-specific random effects.
#' @param x_mean Mean vector for covariate generation.
#' @param x_sigma Covariance matrix for covariate generation.
#' @param beta Coefficient vector for the logistic missingness model.
#' @param train_frac_rule Function returning the training fraction for a given
#'   sample size and scenario.
#' @param weight_cap_rule Function returning the weight truncation level for a
#'   given sample size and scenario.
#' @param output_file Optional path to save the final CSV summary.
#' @param run_sanity_check Logical; if `TRUE`, run one test replicate before
#'   the full simulation.
#'
#' @return Invisibly returns a data frame containing coverage and
#' prediction-region length summaries for the two HCP score choices.

run_score_compare_parallel <- function(
    table_name,
    num_simulations_run = 1000,
    num_cores = 2L,
    sample_sizes = c(100, 300, 500),
    scenarios = c("Homo", "Heter", "Asym", "Bimo", "Bimo-fix"),
    alpha = 0.1,
    n_test = 1,
    n_grid = 200,
    theta_mean = rep(2, 4),
    theta_sigma = diag(1, 4),
    x_mean = rep(0, 4),
    x_sigma = diag(1, 4),
    beta,
    train_frac_rule,
    weight_cap_rule,
    output_file = NULL,
    run_sanity_check = FALSE
) {
  
  # ----------------------------------------------------------
  # Input checks
  # ----------------------------------------------------------
  if (!is.character(table_name) || length(table_name) != 1L) {
    stop("table_name must be a single character string.")
  }
  
  if (!is.numeric(num_simulations_run) ||
      length(num_simulations_run) != 1L ||
      !is.finite(num_simulations_run) ||
      num_simulations_run < 1 ||
      num_simulations_run != floor(num_simulations_run)) {
    stop("num_simulations_run must be a positive integer.")
  }
  
  if (!is.numeric(num_cores) ||
      length(num_cores) != 1L ||
      !is.finite(num_cores)) {
    stop("num_cores must be a single finite number.")
  }
  
  if (!is.function(train_frac_rule)) {
    stop("train_frac_rule must be a function.")
  }
  
  if (!is.function(weight_cap_rule)) {
    stop("weight_cap_rule must be a function.")
  }
  
  num_cores <- as.integer(max(1L, num_cores))
  
  # ----------------------------------------------------------
  # Optional sanity check before starting the full simulation
  # ----------------------------------------------------------
  if (isTRUE(run_sanity_check)) {
    cat(sprintf("Running one sanity-check replicate for %s...\n", table_name))
    
    n0 <- sample_sizes[1]
    scen0 <- scenarios[1]
    
    test_run <- run_one_replicate_score_compare(
      n = n0,
      scenario = scen0,
      alpha = alpha,
      n_test = n_test,
      n_grid = n_grid,
      train_frac = train_frac_rule(n0, scen0),
      theta_mean = theta_mean,
      theta_sigma = theta_sigma,
      x_mean = x_mean,
      x_sigma = x_sigma,
      beta = beta,
      weight_cap = weight_cap_rule(n0, scen0),
      seed = 1
    )
    
    print(test_run)
  }
  
  # ----------------------------------------------------------
  # Build the sample-size/scenario grid
  # ----------------------------------------------------------
  param_grid <- expand.grid(
    sample_size = sample_sizes,
    scenario = scenarios,
    stringsAsFactors = FALSE
  )
  
  param_grid$sample_size <- factor(param_grid$sample_size, levels = sample_sizes)
  param_grid$scenario <- factor(param_grid$scenario, levels = scenarios)
  param_grid <- param_grid[order(param_grid$sample_size, param_grid$scenario), , drop = FALSE]
  rownames(param_grid) <- NULL
  
  results_all <- matrix(NA_real_, nrow = nrow(param_grid), ncol = 4)
  colnames(results_all) <- c(
    "density_cov",
    "density_len",
    "residual_cov",
    "residual_len"
  )
  
  # ----------------------------------------------------------
  # Export all user-defined functions needed by workers
  # ----------------------------------------------------------
  custom_exports <- c(
    "run_one_replicate_score_compare",
    "hcp_score_compare",
    "generate_theta",
    "generate_subject_covariates",
    "generate_errors",
    "generate_outcomes",
    "generate_missing_indicator",
    "generate_simulation_data",
    "evaluate_interval_region",
    "evaluate_bimodal_region",
    "fit_cond_density_qp",
    "fit_propensity_model",
    "adjust_bandwidth_bimodal",
    "quantile_levels",
    "get_S",
    "get_b",
    "get_max_wt",
    "train_frac_rule",
    "weight_cap_rule"
  )
  
  # ----------------------------------------------------------
  # Start parallel backend
  # ----------------------------------------------------------
  cat(sprintf("Running %s using %d cores\n", table_name, num_cores))
  foreach::registerDoSEQ()
  
  cl <- parallel::makeCluster(num_cores, outfile = "")
  doParallel::registerDoParallel(cl)
  parallel::clusterSetRNGStream(cl, iseed = 123)
  
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    try(foreach::registerDoSEQ(), silent = TRUE)
  }, add = TRUE)
  
  # ----------------------------------------------------------
  # Run simulations for each setting
  # ----------------------------------------------------------
  for (row_index in seq_len(nrow(param_grid))) {
    n_val <- as.numeric(as.character(param_grid$sample_size[row_index]))
    scen <- as.character(param_grid$scenario[row_index])
    
    train_frac_use <- train_frac_rule(n_val, scen)
    weight_cap_use <- weight_cap_rule(n_val, scen)
    
    cat(sprintf(
      "Running %s: n = %d, scenario = %s, train_frac = %.2f, weight_cap = %.4f\n",
      table_name, n_val, scen, train_frac_use, weight_cap_use
    ))
    flush.console()
    
    res_mat <- foreach(
      i = seq_len(num_simulations_run),
      .combine = rbind,
      .errorhandling = "pass",
      .export = custom_exports
    ) %dopar% {
      tryCatch(
        run_one_replicate_score_compare(
          n = n_val,
          scenario = scen,
          alpha = alpha,
          n_test = n_test,
          n_grid = n_grid,
          train_frac = train_frac_rule(n_val, scen),
          theta_mean = theta_mean,
          theta_sigma = theta_sigma,
          x_mean = x_mean,
          x_sigma = x_sigma,
          beta = beta,
          weight_cap = weight_cap_rule(n_val, scen),
          seed = i
        ),
        error = function(e) {
          message(sprintf(
            "Replicate %d failed for %s, n=%d, scenario=%s: %s",
            i, table_name, n_val, scen, conditionMessage(e)
          ))
          rep(NA_real_, 4)
        }
      )
    }
    
    num_failed <- sum(!stats::complete.cases(res_mat))
    cat(sprintf(
      "Completed %s: n = %d, scenario = %s with %d failed replicates\n",
      table_name, n_val, scen, num_failed
    ))
    
    results_all[row_index, ] <- colMeans(res_mat, na.rm = TRUE)
  }
  
  # ----------------------------------------------------------
  # Format and save final results
  # ----------------------------------------------------------
  final_results <- cbind(
    data.frame(
      sample_size = as.numeric(as.character(param_grid$sample_size)),
      scenario = as.character(param_grid$scenario),
      stringsAsFactors = FALSE
    ),
    results_all
  )
  
  final_results$scenario[final_results$scenario == "Bimo-fix"] <- "Bimo-Fix"
  
  print(final_results)
  
  if (!is.null(output_file)) {
    write.csv(final_results, file = output_file, row.names = FALSE)
  }
  
  invisible(final_results)
}

# ------------------------------------------------------------
# Table S.4 settings: 20% missing responses
# ------------------------------------------------------------
train_frac_rule_s4 <- function(n, scenario) {
  0.3
}

weight_cap_rule_s4 <- function(n, scenario) {
  30
}

# ------------------------------------------------------------
# Table S.5 settings: 50% missing responses
# ------------------------------------------------------------
train_frac_rule_s5 <- function(n, scenario) {
  0.3
}

weight_cap_rule_s5 <- function(n, scenario) {
  pmin(28, 2.2 * (n / 100)^1.72)
}

# ------------------------------------------------------------
# Run Table S.4
# ------------------------------------------------------------
results_s4 <- run_score_compare_parallel(
  table_name = "Table S.4",
  num_simulations_run = 1000,
  num_cores = 7,
  beta = c(3, 0, 2, 2, 2),
  train_frac_rule = train_frac_rule_s4,
  weight_cap_rule = weight_cap_rule_s4,
  output_file = file.path(results_dir, "tableS4_score_compare_missing20.csv"),
  run_sanity_check = FALSE
)

# ------------------------------------------------------------
# Run Table S.5
# ------------------------------------------------------------
results_s5 <- run_score_compare_parallel(
  table_name = "Table S.5",
  num_simulations_run = 1000,
  num_cores = 7,
  beta = c(0, 0, 2, 2, 2),
  train_frac_rule = train_frac_rule_s5,
  weight_cap_rule = weight_cap_rule_s5,
  output_file = file.path(results_dir, "tableS5_score_compare_missing50.csv"),
  run_sanity_check = FALSE
)