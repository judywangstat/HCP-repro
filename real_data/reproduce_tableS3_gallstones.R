# ============================================================
# Reproduce Table S.3
# Gallstones real-data analysis
# HCPclust-repro
#
# This script reproduces the gallstones real-data results using
# leave-one-subject-out evaluation. It reports pointwise and
# simultaneous coverage and prediction-region length for all methods.
#
# Input:
#   data/gallstones.txt
#
# Output:
#   results/tableS3_gallstones_results.csv
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
data_dir <- file.path(repo_dir, "data")
results_dir <- file.path(repo_dir, "results")

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# ------------------------------------------------------------
# Load required function files
# ------------------------------------------------------------
source(file.path(r_dir, "fit_cond_density_qp.R"), chdir = TRUE)
source(file.path(r_dir, "fit_propensity_model.R"), chdir = TRUE)
source(file.path(r_dir, "tuning_helpers.R"), chdir = TRUE)
source(file.path(r_dir, "evaluation_helpers.R"), chdir = TRUE)

source(file.path(r_dir, "hcp_region.R"), chdir = TRUE)
source(file.path(r_dir, "dwr_region_legacy.R"), chdir = TRUE)
source(file.path(r_dir, "lc_region_legacy.R"), chdir = TRUE)
source(file.path(r_dir, "lmem_region.R"), chdir = TRUE)

source(file.path(r_dir, "gallstones_data_helpers.R"), chdir = TRUE)
source(file.path(r_dir, "run_one_leaveout_gallstones.R"), chdir = TRUE)

# ------------------------------------------------------------
# Main runner
# ------------------------------------------------------------

#' Run gallstones leave-one-subject-out analysis for Table S.3
#'
#' @description
#' This helper runs the gallstones leave-one-subject-out analysis for pointwise
#' and simultaneous prediction settings. It evaluates HCP, DWR, LC, and LMEM,
#' then returns and saves the summary table.
#'
#' @param num_cores Number of parallel workers.
#' @param n_grid Number of response-grid points used to construct prediction
#'   regions.
#' @param seed Random seed.
#' @param output_file Path to the output CSV file.
#'
#' @return A matrix containing coverage and prediction-region length summaries
#'   for all methods.
run_tableS3_gallstones <- function(
    num_cores = 4L,
    n_grid = 200,
    seed = 123,
    output_file = file.path(results_dir, "tableS3_gallstones_results.csv")
) {
  data_file <- file.path(data_dir, "gallstones.txt")
  
  dat <- load_gallstones_data(data_file)
  dat <- impute_gallstones_outcomes(dat, seed = 12345)
  
  test_ids <- sort(unique(dat$id))
  
  run_one_type <- function(prediction_type) {
    alpha_use <- if (prediction_type == "simultaneous") 0.1 / 4 else 0.1
    
    cat("Running", prediction_type, "leave-one-subject-out evaluation\n")
    flush.console()
    
    cl <- parallel::makeCluster(num_cores, outfile = "")
    doParallel::registerDoParallel(cl)
    parallel::clusterSetRNGStream(cl, iseed = seed)
    
    on.exit({
      try(parallel::stopCluster(cl), silent = TRUE)
      try(foreach::registerDoSEQ(), silent = TRUE)
    }, add = TRUE)
    
    res_mat <- foreach(
      ii = seq_along(test_ids),
      .combine = rbind,
      .errorhandling = "pass",
      .export = c(
        "dat",
        "test_ids",
        "prediction_type",
        "alpha_use",
        "n_grid",
        "seed",
        "run_one_leaveout_gallstones",
        "make_gallstones_y_grid",
        "impute_gallstones_outcomes",
        "hcp_region",
        "dwr_region_legacy",
        "lc_region_legacy",
        "lmem_region",
        "fit_cond_density_qp",
        "fit_propensity_model",
        "quantile_levels",
        "evaluate_interval_region"
      ),
      .packages = c("quantreg", "lme4", "merTools")
    ) %dopar% {
      test_id <- test_ids[ii]
      
      tryCatch(
        run_one_leaveout_gallstones(
          test_id = test_id,
          dat = dat,
          prediction_type = prediction_type,
          alpha = alpha_use,
          n_grid = n_grid,
          seed = seed + ii
        ),
        error = function(e) {
          message(sprintf(
            "Held-out subject %s failed under %s: %s",
            as.character(test_id),
            prediction_type,
            conditionMessage(e)
          ))
          rep(NA_real_, 8)
        }
      )
    }
    
    colnames(res_mat) <- c(
      "HCP_cov", "HCP_len",
      "DWR_cov", "DWR_len",
      "LC_cov", "LC_len",
      "LMEM_cov", "LMEM_len"
    )
    
    num_failed <- sum(!stats::complete.cases(res_mat))
    cat(sprintf(
      "Completed %s evaluation with %d failed leave-one-out runs\n",
      prediction_type,
      num_failed
    ))
    
    colMeans(res_mat, na.rm = TRUE)
  }
  
  pointwise_res <- run_one_type("pointwise")
  simultaneous_res <- run_one_type("simultaneous")
  
  final_table <- rbind(
    pointwise_coverage = pointwise_res[c("HCP_cov", "DWR_cov", "LC_cov", "LMEM_cov")],
    pointwise_length = pointwise_res[c("HCP_len", "DWR_len", "LC_len", "LMEM_len")],
    simultaneous_coverage = simultaneous_res[c("HCP_cov", "DWR_cov", "LC_cov", "LMEM_cov")],
    simultaneous_length = simultaneous_res[c("HCP_len", "DWR_len", "LC_len", "LMEM_len")]
  )
  
  colnames(final_table) <- c("HCP", "DWR", "LC", "LMEM")
  
  print(final_table)
  write.csv(final_table, file = output_file, row.names = TRUE)
  
  invisible(final_table)
}

# ------------------------------------------------------------
# Run Table S.3 reproduction
# ------------------------------------------------------------
tableS3_results <- run_tableS3_gallstones(
  num_cores = 7,
  n_grid = 200,
  seed = 123,
  output_file = file.path(results_dir, "tableS3_gallstones_results.csv")
)