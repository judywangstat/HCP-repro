# ============================================================
# Reproduce Figure S.10 Data
# Pointwise benchmark comparison
# HCPclust-repro
#
# This script generates the pointwise simulation results used
# for Figure S.10 in the supplementary materials. It uses
# run_one_replicate_pointwise(), where each replicate returns one
# row per test observation.
#
# With n = 500, n_test = 1, m = 5, and 1000 replicates, each
# scenario produces 5000 pointwise results. The final outputs are
# method-specific CSV files saved in results/.
#
# Output:
#   results/figureS10_HCP_final_mat_n500.csv
#   results/figureS10_DWR_final_mat_n500.csv
#   results/figureS10_LC_final_mat_n500.csv
#   results/figureS10_LMEM_final_mat_n500.csv
#   results/figureS10_Oracle_final_mat_n500.csv
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
source(file.path(r_dir, "dwr_region.R"), chdir = TRUE)
source(file.path(r_dir, "lc_region.R"), chdir = TRUE)
source(file.path(r_dir, "lmem_region.R"), chdir = TRUE)
source(file.path(r_dir, "oracle_region.R"), chdir = TRUE)

source(file.path(r_dir, "run_one_replicate_pointwise.R"), chdir = TRUE)

#' Generate pointwise simulation data for Figure S.10
#'
#' @description
#' This helper runs pointwise benchmark simulations across multiple scenarios
#' and saves method-specific result matrices. Each replicate holds out one test
#' subject and returns pointwise coverage and prediction-region length for all
#' methods.
#'
#' @param num_simulations_run Number of simulation replicates.
#' @param num_cores Number of parallel workers.
#' @param n Number of training subjects.
#' @param scenarios Simulation scenarios to evaluate.
#' @param alpha Miscoverage level.
#' @param n_test Number of held-out test subjects per replicate.
#' @param n_grid Number of response-grid points.
#' @param theta_mean Mean vector of subject-specific random effects.
#' @param theta_sigma Covariance matrix of subject-specific random effects.
#' @param x_mean Mean vector for covariate generation.
#' @param x_sigma Covariance matrix for covariate generation.
#' @param beta Coefficient vector for the logistic missingness model.
#' @param weight_cap Upper bound for inverse-propensity weights.
#' @param oracle_B_mc Monte Carlo sample size used by the oracle method.
#' @param lmem_n_sims Number of simulations used by the LMEM benchmark.
#' @param output_dir Directory where CSV outputs are saved.
#'
#' @return Invisibly returns a list of scenario-specific pointwise result
#' matrices.
run_figureS10_data <- function(
    num_simulations_run = 1000,
    num_cores = 7L,
    n = 500,
    scenarios = c("Homo", "Heter", "Asym", "Bimo", "Bimo-fix"),
    alpha = 0.1,
    n_test = 1,
    n_grid = 200,
    theta_mean = rep(2, 4),
    theta_sigma = diag(1, 4),
    x_mean = rep(0, 4),
    x_sigma = diag(1, 4),
    beta = c(3, 0, 2, 2, 2),
    weight_cap = 30,
    oracle_B_mc = 5000L,
    lmem_n_sims = 500L,
    output_dir = results_dir
) {
  if (!is.numeric(num_simulations_run) ||
      length(num_simulations_run) != 1L ||
      !is.finite(num_simulations_run) ||
      num_simulations_run < 1 ||
      num_simulations_run != floor(num_simulations_run)) {
    stop("num_simulations_run must be a positive integer.")
  }
  
  num_cores <- as.integer(max(1L, num_cores))
  
  methods <- c("HCP", "DWR", "LC", "LMEM", "Oracle")
  
  custom_exports <- c(
    "run_one_replicate_pointwise",
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
    "get_max_wt"
  )
  
  cat("Using", num_cores, "cores\n")
  
  foreach::registerDoSEQ()
  cl <- parallel::makeCluster(num_cores, outfile = "")
  doParallel::registerDoParallel(cl)
  
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    try(foreach::registerDoSEQ(), silent = TRUE)
  }, add = TRUE)
  
  all_by_scenario <- list()
  
  for (scen in scenarios) {
    cat(sprintf("Running Figure S.10 data: n = %d, scenario = %s\n", n, scen))
    flush.console()
    
    res_scen <- foreach(
      i = seq_len(num_simulations_run),
      .combine = rbind,
      .export = custom_exports,
      .errorhandling = "pass"
    ) %dopar% {
      tryCatch(
        run_one_replicate_pointwise(
          n = n,
          scenario = scen,
          alpha = alpha,
          n_test = n_test,
          n_grid = n_grid,
          theta_mean = theta_mean,
          theta_sigma = theta_sigma,
          x_mean = x_mean,
          x_sigma = x_sigma,
          beta = beta,
          weight_cap = weight_cap,
          oracle_B_mc = oracle_B_mc,
          lmem_n_sims = lmem_n_sims,
          seed = i
        ),
        error = function(e) {
          message(sprintf(
            "Replicate %d failed for scenario=%s: %s",
            i, scen, conditionMessage(e)
          ))
          
          matrix(
            NA_real_,
            nrow = 5,
            ncol = 10,
            dimnames = list(
              NULL,
              c(
                "HCP_cov", "HCP_len",
                "DWR_cov", "DWR_len",
                "LC_cov", "LC_len",
                "LMEM_cov", "LMEM_len",
                "Oracle_cov", "Oracle_len"
              )
            )
          )
        }
      )
    }
    
    num_failed_rows <- sum(!stats::complete.cases(res_scen))
    cat(sprintf(
      "Completed scenario = %s with %d failed pointwise rows\n",
      scen, num_failed_rows
    ))
    
    all_by_scenario[[scen]] <- res_scen
  }
  
  # ----------------------------------------------------------
  # Save method-specific result matrices used by the plotting script
  # ----------------------------------------------------------
  for (mtd in methods) {
    mat_mtd <- do.call(cbind, lapply(names(all_by_scenario), function(scen) {
      tmp <- all_by_scenario[[scen]][
        ,
        c(paste0(mtd, "_cov"), paste0(mtd, "_len")),
        drop = FALSE
      ]
      
      colnames(tmp) <- c(
        paste0(scen, "_cov"),
        paste0(scen, "_len")
      )
      
      tmp
    }))
    
    # read.csv() converts Bimo-fix to Bimo.fix, so make.names keeps naming
    # consistent with the plotting script.
    colnames(mat_mtd) <- make.names(colnames(mat_mtd))
    
    output_file <- file.path(output_dir, paste0("figureS10_", mtd, "_final_mat_n500.csv"))
    write.csv(mat_mtd, file = output_file, row.names = FALSE)
    
    cat("Saved:", output_file, "\n")
  }
  
  invisible(all_by_scenario)
}

# ------------------------------------------------------------
# Run Figure S.10 data generation
# ------------------------------------------------------------
figureS10_results <- run_figureS10_data(
  num_simulations_run = 1000,
  num_cores = 7,
  n = 500,
  weight_cap = 30,
  oracle_B_mc = 5000L,
  lmem_n_sims = 500L,
  output_dir = results_dir
)