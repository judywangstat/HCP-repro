# ============================================================
# Reproduce Figure S.8 and Figure S.9 Data
# HCPclust-repro
#
# This script generates simulation data for Figure S.8 and
# Figure S.9. The two figures evaluate the sensitivity of HCP
# to the number of calibration subsamples B and the number of
# repeated data splits S.
#
# Outputs:
#   results/figureS8_sensitivity_B_data.csv
#   results/figureS9_sensitivity_S_data.csv
#
# Note:
#   Update repo_dir below to the local path of this repository
#   before running the script.
# ============================================================

suppressPackageStartupMessages({
  library(doParallel)
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
source(file.path(r_dir, "hcp_region.R"), chdir = TRUE)
source(file.path(r_dir, "run_one_replicate_sensitivity.R"), chdir = TRUE)

#' Run one sensitivity setting in parallel
#'
#' @description
#' Runs repeated HCP simulation replicates for a fixed sample size and fixed
#' values of B and S. The function returns the replicate mean and standard error
#' for marginal coverage and prediction-region length.
#'
#' @param n Sample size.
#' @param scenario Simulation scenario.
#' @param B Number of calibration subsamples.
#' @param S Number of repeated data splits.
#' @param num_simulations Number of Monte Carlo replicates.
#' @param num_cores Number of parallel workers.
#' @param alpha Miscoverage level.
#' @param n_grid Number of response-grid points.
#' @param weight_cap Upper truncation level for inverse-propensity weights.
#'
#' @return A one-row data frame with summary statistics.
run_one_sensitivity_setting <- function(
    n,
    scenario,
    B,
    S,
    num_simulations = 1000,
    num_cores = 2L,
    alpha = 0.1,
    n_grid = 200,
    weight_cap = 30
) {
  if (!is.numeric(num_simulations) ||
      length(num_simulations) != 1L ||
      !is.finite(num_simulations) ||
      num_simulations < 1 ||
      num_simulations != floor(num_simulations)) {
    stop("num_simulations must be a positive integer.")
  }
  
  num_cores <- as.integer(max(1L, num_cores))
  
  custom_exports <- c(
    "run_one_replicate_sensitivity",
    "generate_theta",
    "generate_subject_covariates",
    "generate_errors",
    "generate_outcomes",
    "generate_missing_indicator",
    "generate_simulation_data",
    "evaluate_interval_region",
    "evaluate_bimodal_region",
    "hcp_region",
    "fit_cond_density_qp",
    "fit_propensity_model",
    "adjust_bandwidth_bimodal",
    "quantile_levels",
    "get_S",
    "get_b",
    "get_max_wt"
  )
  
  cl <- parallel::makeCluster(num_cores, outfile = "")
  doParallel::registerDoParallel(cl)
  
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    try(foreach::registerDoSEQ(), silent = TRUE)
  }, add = TRUE)
  
  res_mat <- foreach(
    i = seq_len(num_simulations),
    .combine = rbind,
    .errorhandling = "pass",
    .export = custom_exports
  ) %dopar% {
    tryCatch(
      run_one_replicate_sensitivity(
        n = n,
        scenario = scenario,
        B = B,
        S = S,
        alpha = alpha,
        n_grid = n_grid,
        weight_cap = weight_cap,
        seed = i
      ),
      error = function(e) {
        message(sprintf(
          "Replicate %d failed for n=%d, scenario=%s, B=%d, S=%d: %s",
          i, n, scenario, B, S, conditionMessage(e)
        ))
        c(coverage = NA_real_, length = NA_real_)
      }
    )
  }
  
  num_failed <- sum(!stats::complete.cases(res_mat))
  
  data.frame(
    n = n,
    scenario = scenario,
    B = B,
    S = S,
    coverage = mean(res_mat[, "coverage"], na.rm = TRUE),
    length = mean(res_mat[, "length"], na.rm = TRUE),
    coverage_se = stats::sd(res_mat[, "coverage"], na.rm = TRUE) / sqrt(sum(!is.na(res_mat[, "coverage"]))),
    length_se = stats::sd(res_mat[, "length"], na.rm = TRUE) / sqrt(sum(!is.na(res_mat[, "length"]))),
    num_failed = num_failed,
    stringsAsFactors = FALSE
  )
}

#' Run the Figure S.8 sensitivity study
#'
#' @description
#' Figure S.8 varies B while fixing S = 1.
#'
#' @param num_simulations Number of simulation replicates per setting.
#' @param num_cores Number of parallel workers.
#' @param n_values Sample sizes to evaluate.
#' @param B_values Values of B to evaluate.
#' @param scenario Simulation scenario.
#' @param S_fixed Fixed value of S.
#' @param alpha Miscoverage level.
#' @param n_grid Number of response-grid points.
#' @param weight_cap Upper bound for inverse-propensity weights.
#' @param output_file Path to save the output CSV file.
#'
#' @return Invisibly returns a data frame containing Figure S.8 sensitivity
#' results.
run_figureS8_data <- function(
    num_simulations = 1000,
    num_cores = 7L,
    n_values = c(100, 300, 500),
    B_values = c(1, 3, 5, 9, 11, 13, 15, 17, 19, 21),
    scenario = "Heter",
    S_fixed = 1L,
    alpha = 0.1,
    n_grid = 200,
    weight_cap = 30,
    output_file = file.path(results_dir, "figureS8_sensitivity_B_data.csv")
) {
  out <- list()
  counter <- 1L
  
  for (n_val in n_values) {
    for (B_val in B_values) {
      cat(sprintf("Figure S.8: running n = %d, B = %d, S = %d\n",
                  n_val, B_val, S_fixed))
      flush.console()
      
      out[[counter]] <- run_one_sensitivity_setting(
        n = n_val,
        scenario = scenario,
        B = B_val,
        S = S_fixed,
        num_simulations = num_simulations,
        num_cores = num_cores,
        alpha = alpha,
        n_grid = n_grid,
        weight_cap = weight_cap
      )
      
      counter <- counter + 1L
    }
  }
  
  dat <- do.call(rbind, out)
  
  if (!is.null(output_file)) {
    write.csv(dat, file = output_file, row.names = FALSE)
  }
  
  invisible(dat)
}

#' Run the Figure S.9 sensitivity study
#'
#' @description
#' Figure S.9 varies S while fixing B = 2.
#'
#' @param num_simulations Number of simulation replicates per setting.
#' @param num_cores Number of parallel workers.
#' @param n_values Sample sizes to evaluate.
#' @param S_values Values of S to evaluate.
#' @param scenario Simulation scenario.
#' @param B_fixed Fixed value of B.
#' @param alpha Miscoverage level.
#' @param n_grid Number of response-grid points.
#' @param weight_cap Upper bound for inverse-propensity weights.
#' @param output_file Path to save the output CSV file.
#'
#' @return Invisibly returns a data frame containing Figure S.9 sensitivity
#' results.
run_figureS9_data <- function(
    num_simulations = 1000,
    num_cores = 7L,
    n_values = c(100, 300, 500),
    S_values = c(1, 3, 5, 9, 11, 13, 15, 17, 19, 21),
    scenario = "Heter",
    B_fixed = 2L,
    alpha = 0.1,
    n_grid = 200,
    weight_cap = 30,
    output_file = file.path(results_dir, "figureS9_sensitivity_S_data.csv")
) {
  out <- list()
  counter <- 1L
  
  for (n_val in n_values) {
    for (S_val in S_values) {
      cat(sprintf("Figure S.9: running n = %d, B = %d, S = %d\n",
                  n_val, B_fixed, S_val))
      flush.console()
      
      out[[counter]] <- run_one_sensitivity_setting(
        n = n_val,
        scenario = scenario,
        B = B_fixed,
        S = S_val,
        num_simulations = num_simulations,
        num_cores = num_cores,
        alpha = alpha,
        n_grid = n_grid,
        weight_cap = weight_cap
      )
      
      counter <- counter + 1L
    }
  }
  
  dat <- do.call(rbind, out)
  
  if (!is.null(output_file)) {
    write.csv(dat, file = output_file, row.names = FALSE)
  }
  
  invisible(dat)
}

# ------------------------------------------------------------
# Run experiments
# ------------------------------------------------------------
figureS8_data <- run_figureS8_data(
  num_simulations = 1000,
  num_cores = 7L,
  weight_cap = 30
)

figureS9_data <- run_figureS9_data(
  num_simulations = 1000,
  num_cores = 7L,
  weight_cap = 30
)