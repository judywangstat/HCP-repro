# ============================================================
# Reproduce Figure S.7 Data
# HCP-repro
#
# This script generates simulation data for Figure S.7.
# It evaluates simultaneous prediction performance of HCP as
# the sample size n varies.
#
# For each scenario, it saves one CSV file to results/:
#   figureS7_Homo_results.csv
#   figureS7_Heter_results.csv
#   figureS7_Asym_results.csv
#   figureS7_Bimo_results.csv
#   figureS7_Bimo-Fix_results.csv
#
# Each file has columns:
#   cov_n100, len_n100, cov_n200, len_n200, ...
# where each row corresponds to one simulation replicate.
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
repo_dir <- normalizePath(Sys.getenv("HCP_REPO_DIR", "."), mustWork = TRUE)
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

source(file.path(r_dir, "run_one_replicate_simultaneous.R"), chdir = TRUE)



#' Generate sample-size sensitivity data for Figure S.7
#'
#' @description
#' This helper generates simulation data for Figure S.7 by evaluating
#' simultaneous prediction performance across sample sizes and simulation
#' scenarios. For each scenario, it saves one CSV file containing replicate-level
#' coverage and prediction-region length for all sample sizes.
#'
#' @param num_simulations_run Number of simulation replicates per setting.
#' @param num_cores Number of parallel workers.
#' @param n_values Sample sizes to evaluate.
#' @param scenarios Simulation scenarios to evaluate.
#' @param missing Missingness level. Supported values are `"20"` and `"50"`.
#' @param alpha Miscoverage level.
#' @param n_test Number of held-out test subjects per replicate.
#' @param n_grid Number of response-grid points.
#' @param train_frac Fraction of subjects assigned to training.
#' @param theta_mean Mean vector of subject-specific random effects.
#' @param theta_sigma Covariance matrix of subject-specific random effects.
#' @param x_mean Mean vector for covariate generation.
#' @param x_sigma Covariance matrix for covariate generation.
#' @param weight_cap_rule Function returning the weight truncation level as a
#'   function of sample size.
#' @param output_prefix Prefix used in the output CSV file names.
#' @param run_sanity_check Logical; if `TRUE`, run one test replicate before
#'   the full simulation.
#'
#' @return Invisibly returns `TRUE` after saving all scenario-specific CSV files.
run_figureS7_data_parallel <- function(
    num_simulations_run = 1000,
    num_cores = 7L,
    n_values = seq(100, 2000, by = 100),
    scenarios = c("Homo", "Heter", "Asym", "Bimo", "Bimo-fix"),
    missing = "20",
    alpha = 0.1,
    n_test = 1,
    n_grid = 200,
    train_frac = 0.3,
    theta_mean = rep(2, 4),
    theta_sigma = diag(1, 4),
    x_mean = rep(0, 4),
    x_sigma = diag(1, 4),
    weight_cap_rule = function(n) {
      30
    },
    output_prefix = "figureS7",
    run_sanity_check = FALSE
) {
  if (!is.numeric(num_simulations_run) ||
      length(num_simulations_run) != 1L ||
      !is.finite(num_simulations_run) ||
      num_simulations_run < 1 ||
      num_simulations_run != floor(num_simulations_run)) {
    stop("num_simulations_run must be a positive integer.")
  }
  
  num_cores <- as.integer(max(1L, num_cores))
  
  if (!missing %in% c("20", "50")) {
    stop("missing must be either '20' or '50'.")
  }
  
  if (!is.numeric(train_frac) ||
      length(train_frac) != 1L ||
      !is.finite(train_frac) ||
      train_frac <= 0 ||
      train_frac >= 1) {
    stop("train_frac must be a single number in (0,1).")
  }
  
  custom_exports <- c(
    "run_one_replicate_simultaneous",
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
  
  beta_use <- switch(
    missing,
    "20" = c(3, 0, 2, 2, 2),
    "50" = c(0, 0, 2, 2, 2)
  )
  
  if (isTRUE(run_sanity_check)) {
    cat("Running one sanity-check replicate...\n")
    print(
      run_one_replicate_simultaneous(
        n = n_values[1],
        scenario = scenarios[1],
        missing = missing,
        alpha = alpha,
        n_test = n_test,
        n_grid = n_grid,
        theta_mean = theta_mean,
        theta_sigma = theta_sigma,
        x_mean = x_mean,
        x_sigma = x_sigma,
        beta = beta_use,
        weight_cap = weight_cap_rule(n_values[1]),
        train_frac = train_frac,
        seed = 1
      )
    )
  }
  
  cat("Using", num_cores, "cores\n")
  
  foreach::registerDoSEQ()
  cl <- parallel::makeCluster(num_cores, outfile = "")
  doParallel::registerDoParallel(cl)
  
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    try(foreach::registerDoSEQ(), silent = TRUE)
  }, add = TRUE)
  
  for (scen in scenarios) {
    cat(sprintf("\n===== Running scenario = %s =====\n", scen))
    flush.console()
    
    scenario_mat <- matrix(
      NA_real_,
      nrow = num_simulations_run,
      ncol = 2 * length(n_values)
    )
    
    colnames(scenario_mat) <- as.vector(rbind(
      paste0("cov_n", n_values),
      paste0("len_n", n_values)
    ))
    
    for (j in seq_along(n_values)) {
      n_val <- n_values[j]
      weight_cap_use <- weight_cap_rule(n_val)
      
      cat(sprintf(
        "Running scenario = %s, n = %d, missing = %s%%, weight_cap = %.4f\n",
        scen, n_val, missing, weight_cap_use
      ))
      flush.console()
      
      res_mat <- foreach(
        i = seq_len(num_simulations_run),
        .combine = rbind,
        .export = custom_exports,
        .errorhandling = "pass"
      ) %dopar% {
        tryCatch(
          run_one_replicate_simultaneous(
            n = n_val,
            scenario = scen,
            missing = missing,
            alpha = alpha,
            n_test = n_test,
            n_grid = n_grid,
            theta_mean = theta_mean,
            theta_sigma = theta_sigma,
            x_mean = x_mean,
            x_sigma = x_sigma,
            beta = beta_use,
            weight_cap = weight_cap_use,
            train_frac = train_frac,
            seed = i
          ),
          error = function(e) {
            message(sprintf(
              "Replicate %d failed for scenario=%s, n=%d, missing=%s: %s",
              i, scen, n_val, missing, conditionMessage(e)
            ))
            c(HCP_cov = NA_real_, HCP_len = NA_real_)
          }
        )
      }
      
      num_failed <- sum(!stats::complete.cases(res_mat))
      cat(sprintf(
        "Completed scenario = %s, n = %d with %d failed replicates\n",
        scen, n_val, num_failed
      ))
      
      scenario_mat[, paste0("cov_n", n_val)] <- res_mat[, "HCP_cov"]
      scenario_mat[, paste0("len_n", n_val)] <- res_mat[, "HCP_len"]
      
      gc()
    }
    
    scen_file <- ifelse(scen == "Bimo-fix", "Bimo-Fix", scen)
    output_file <- file.path(
      results_dir,
      paste0(output_prefix, "_", scen_file, "_results.csv")
    )
    
    write.csv(scenario_mat, file = output_file, row.names = FALSE)
    cat("Saved:", output_file, "\n")
  }
  
  invisible(TRUE)
}

# ------------------------------------------------------------
# Run Figure S.7 simulation
# ------------------------------------------------------------
figureS7_done <- run_figureS7_data_parallel(
  num_simulations_run = 1000,
  num_cores = 7,
  n_values = seq(100, 2000, by = 100),
  scenarios = c("Homo", "Heter", "Asym", "Bimo", "Bimo-fix"),
  missing = "20",
  train_frac = 0.3,
  weight_cap_rule = function(n) {
    30
  },
  output_prefix = "figureS7",
  run_sanity_check = FALSE
)