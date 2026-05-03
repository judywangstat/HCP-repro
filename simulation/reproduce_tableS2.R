# ============================================================
# Reproduce Table S.2
# Main benchmark simulation study with 50% missingness
# HCPclust-repro
#
# This script reproduces Table S.2 in the supplementary materials.
# The simulation setting is the same as Table 1, except that the
# missingness level is set to 50%.
#
# Output:
#   results/tableS2_final_results_missing50.csv
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

repo_dir <- normalizePath("~/Desktop/HCPclust-repro", mustWork = TRUE)
r_dir <- file.path(repo_dir, "R")
results_dir <- file.path(repo_dir, "results")

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

source(file.path(r_dir, "simulation_dgp.R"), chdir = TRUE)
source(file.path(r_dir, "tuning_helpers.R"), chdir = TRUE)
source(file.path(r_dir, "evaluation_helpers.R"), chdir = TRUE)
source(file.path(r_dir, "hcp_region.R"), chdir = TRUE)
source(file.path(r_dir, "dwr_region.R"), chdir = TRUE)
source(file.path(r_dir, "lc_region.R"), chdir = TRUE)
source(file.path(r_dir, "lmem_region.R"), chdir = TRUE)
source(file.path(r_dir, "oracle_region.R"), chdir = TRUE)
source(file.path(r_dir, "run_one_replicate.R"), chdir = TRUE)
source(file.path(r_dir, "fit_propensity_model.R"), chdir = TRUE)
source(file.path(r_dir, "fit_cond_density_qp.R"), chdir = TRUE)


#' Run the Table S.2 simulation study in parallel
#'
#' @description
#' This helper runs the benchmark simulation study under the 50% missingness
#' setting. For each sample size and simulation scenario, it runs multiple
#' simulation replicates, evaluates all benchmark methods, and aggregates
#' replicate-level coverage and prediction-region length.
#'
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
#' @param weight_cap_rule Function returning the weight truncation level as a
#'   function of sample size.
#' @param output_file Optional path to save the final CSV summary.
#' @param run_sanity_check Logical; if `TRUE`, run one test replicate before
#'   the full simulation.
#'
#' @return Invisibly returns a data frame containing average coverage and
#' prediction-region length for each method, sample size, and scenario.

run_tableS2_parallel <- function(
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
    beta = c(0, 0, 2, 2, 2),
    weight_cap_rule = function(n) {
      pmin(28, 2.2 * (n / 100)^1.72)
    },
    output_file = NULL,
    run_sanity_check = FALSE
) {
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
  
  num_cores <- as.integer(max(1L, num_cores))
  
  if (isTRUE(run_sanity_check)) {
    cat("Running one sanity-check replicate for Table S.2...\n")
    test_run <- run_one_replicate(
      n = sample_sizes[1],
      scenario = scenarios[1],
      alpha = alpha,
      n_test = n_test,
      n_grid = n_grid,
      theta_mean = theta_mean,
      theta_sigma = theta_sigma,
      x_mean = x_mean,
      x_sigma = x_sigma,
      beta = beta,
      weight_cap = weight_cap_rule(sample_sizes[1]),
      seed = 1
    )
    print(test_run)
  }
  
  param_grid <- expand.grid(
    sample_size = sample_sizes,
    scenario = scenarios,
    stringsAsFactors = FALSE
  )
  
  param_grid$sample_size <- factor(param_grid$sample_size, levels = sample_sizes)
  param_grid$scenario <- factor(param_grid$scenario, levels = scenarios)
  param_grid <- param_grid[order(param_grid$sample_size, param_grid$scenario), , drop = FALSE]
  rownames(param_grid) <- NULL
  
  results_all <- matrix(NA_real_, nrow = nrow(param_grid), ncol = 10)
  colnames(results_all) <- c(
    "HCP_cov", "HCP_len",
    "DWR_cov", "DWR_len",
    "LC_cov", "LC_len",
    "LMEM_cov", "LMEM_len",
    "Oracle_cov", "Oracle_len"
  )
  
  custom_exports <- c(
    "run_one_replicate",
    "generate_theta",
    "generate_subject_covariates",
    "generate_errors",
    "generate_outcomes",
    "generate_missing_indicator",
    "generate_simulation_data",
    "evaluate_interval_region",
    "evaluate_bimodal_region",
    "hcp_region",
    "dwr_region",
    "lc_region",
    "lmem_region",
    "oracle_interval_from_dgp",
    "oracle_region",
    "fit_cond_density_qp",
    "fit_propensity_model",
    "adjust_bandwidth_bimodal",
    "quantile_levels",
    "get_S",
    "get_b",
    "get_max_wt",
    "weight_cap_rule"
  )
  
  cat("Using", num_cores, "cores\n")
  foreach::registerDoSEQ()
  
  cl <- parallel::makeCluster(num_cores, outfile = "")
  doParallel::registerDoParallel(cl)
  parallel::clusterSetRNGStream(cl, iseed = 123)
  
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    try(foreach::registerDoSEQ(), silent = TRUE)
  }, add = TRUE)
  
  for (row_index in seq_len(nrow(param_grid))) {
    n_val <- as.numeric(as.character(param_grid$sample_size[row_index]))
    scen <- as.character(param_grid$scenario[row_index])
    
    cat(sprintf("Running n = %d, scenario = %s\n", n_val, scen))
    flush.console()
    
    res_mat <- foreach(
      i = seq_len(num_simulations_run),
      .combine = rbind,
      .errorhandling = "pass",
      .export = custom_exports
    ) %dopar% {
      tryCatch(
        run_one_replicate(
          n = n_val,
          scenario = scen,
          alpha = alpha,
          n_test = n_test,
          n_grid = n_grid,
          theta_mean = theta_mean,
          theta_sigma = theta_sigma,
          x_mean = x_mean,
          x_sigma = x_sigma,
          beta = beta,
          weight_cap = weight_cap_rule(n_val),
          seed = i
        ),
        error = function(e) {
          message(sprintf(
            "Replicate %d failed for n=%d, scenario=%s: %s",
            i, n_val, scen, conditionMessage(e)
          ))
          rep(NA_real_, 10)
        }
      )
    }
    
    num_failed <- sum(!stats::complete.cases(res_mat))
    cat(sprintf(
      "Completed n = %d, scenario = %s with %d failed replicates\n",
      n_val, scen, num_failed
    ))
    
    results_all[row_index, ] <- colMeans(res_mat, na.rm = TRUE)
  }
  
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

results <- run_tableS2_parallel(
  num_simulations_run = 1000,
  num_cores = 7,
  beta = c(0, 0, 2, 2, 2),
  output_file = file.path(results_dir, "tableS2_final_results_missing50.csv")
)