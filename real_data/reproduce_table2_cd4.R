# ============================================================
# Reproduce Table 2
# CD4 real-data analysis
# HCPclust-repro
#
# This script reproduces the CD4 real-data results under two
# artificial missingness settings: 20% and 50%. It runs
# leave-one-subject-out analyses for pointwise and simultaneous
# prediction settings.
#
# Input:
#   data/CD4_data.txt
#
# Output:
#   results/table2_cd4_results.csv
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
# Load required functions
# ------------------------------------------------------------
source(file.path(r_dir, "cd4_data_helpers.R"), chdir = TRUE)
source(file.path(r_dir, "fit_propensity_model.R"), chdir = TRUE)
source(file.path(r_dir, "fit_cond_density_qp.R"), chdir = TRUE)
source(file.path(r_dir, "tuning_helpers.R"), chdir = TRUE)
source(file.path(r_dir, "evaluation_helpers.R"), chdir = TRUE)

source(file.path(r_dir, "hcp_region.R"), chdir = TRUE)
source(file.path(r_dir, "dwr_region_legacy.R"), chdir = TRUE)
source(file.path(r_dir, "lc_region_legacy.R"), chdir = TRUE)
source(file.path(r_dir, "lmem_region.R"), chdir = TRUE)

source(file.path(r_dir, "run_one_leaveout_cd4.R"), chdir = TRUE)

#' Run CD4 leave-one-subject-out analysis
#'
#' @description
#' This helper runs the CD4 leave-one-subject-out analysis for one prediction
#' type and one missingness setting. It applies each method to every held-out
#' subject and returns method-level average coverage and prediction-region
#' length.
#'
#' @param prediction_type Type of prediction summary. Supported values are
#'   `"pointwise"` and `"simultaneous"`.
#' @param missing_rate Artificial missingness setting. Supported values are
#'   `"20"` and `"50"`.
#' @param methods Character vector of method names to evaluate.
#' @param alpha Miscoverage level.
#' @param n_grid Number of response-grid points used to construct prediction
#'   regions.
#' @param num_cores Number of parallel workers.
#' @param seed Random seed.
#'
#' @return A data frame with method-level coverage, prediction-region length,
#'   number of subjects, and number of failed leave-one-subject-out fits.
run_cd4_table <- function(
    prediction_type = c("pointwise", "simultaneous"),
    missing_rate = c("20", "50"),
    methods = c("HCP", "DWR", "LC", "LMEM"),
    alpha = 0.1,
    n_grid = 200,
    num_cores = 2L,
    seed = 123
) {
  prediction_type <- match.arg(prediction_type)
  missing_rate <- match.arg(as.character(missing_rate), choices = c("20", "50"))
  
  cd4_file <- file.path(data_dir, "CD4_data.txt")
  
  cd4_dat <- read_cd4_data(cd4_file)
  cd4_dat <- add_cd4_missingness(
    dat = cd4_dat,
    missing_rate = missing_rate,
    seed = seed
  )
  
  test_ids <- sort(unique(cd4_dat$id))
  
  custom_exports <- c(
    "run_one_leaveout_cd4",
    "read_cd4_data",
    "add_cd4_missingness",
    "get_cd4_x_cols",
    "get_cd4_missing_beta",
    "get_cd4_weight_cap",
    "generate_cd4_missing_indicator",
    "make_cd4_y_grid",
    "summarize_cd4_leaveout",
    "evaluate_interval_region",
    "hcp_region",
    "dwr_region_legacy",
    "lc_region_legacy",
    "lmem_region",
    "fit_cond_density_qp",
    "fit_propensity_model",
    "quantile_levels"
  )
  
  num_cores <- as.integer(max(1L, num_cores))
  
  cl <- parallel::makeCluster(num_cores, outfile = "")
  doParallel::registerDoParallel(cl)
  parallel::clusterSetRNGStream(cl, iseed = seed)
  
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    try(foreach::registerDoSEQ(), silent = TRUE)
  }, add = TRUE)
  
  out_list <- vector("list", length(methods))
  names(out_list) <- methods
  
  for (method in methods) {
    cat(sprintf(
      "Running CD4: prediction_type = %s, missing = %s, method = %s\n",
      prediction_type, missing_rate, method
    ))
    flush.console()
    
    res_mat <- foreach(
      ii = seq_along(test_ids),
      .combine = rbind,
      .errorhandling = "pass",
      .export = custom_exports
    ) %dorng% {
      test_id <- test_ids[ii]
      
      tryCatch(
        run_one_leaveout_cd4(
          dat = cd4_dat,
          test_id = test_id,
          method = method,
          prediction_type = prediction_type,
          missing_rate = missing_rate,
          alpha = alpha,
          n_grid = n_grid,
          seed = seed + ii
        ),
        error = function(e) {
          message(sprintf(
            "Failed: prediction_type=%s, missing=%s, method=%s, test_id=%s: %s",
            prediction_type, missing_rate, method, as.character(test_id),
            conditionMessage(e)
          ))
          c(coverage = NA_real_, length = NA_real_)
        }
      )
    }
    
    out_list[[method]] <- data.frame(
      prediction_type = prediction_type,
      missing_rate = missing_rate,
      method = method,
      coverage = mean(res_mat[, "coverage"], na.rm = TRUE),
      length = mean(res_mat[, "length"], na.rm = TRUE),
      n_subjects = length(test_ids),
      n_failed = sum(!stats::complete.cases(res_mat)),
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, out_list)
}

# ------------------------------------------------------------
# Run Table 2
# ------------------------------------------------------------
methods <- c("HCP", "DWR", "LC", "LMEM")
missing_settings <- c("20", "50")

table2_pointwise <- do.call(
  rbind,
  lapply(missing_settings, function(miss) {
    run_cd4_table(
      prediction_type = "pointwise",
      missing_rate = miss,
      methods = methods,
      alpha = 0.1,
      n_grid = 200,
      num_cores = 7,
      seed = 123
    )
  })
)

table2_simultaneous <- do.call(
  rbind,
  lapply(missing_settings, function(miss) {
    run_cd4_table(
      prediction_type = "simultaneous",
      missing_rate = miss,
      methods = methods,
      alpha = 0.1,
      n_grid = 200,
      num_cores = 7,
      seed = 123
    )
  })
)

table2_results <- rbind(table2_pointwise, table2_simultaneous)

# ------------------------------------------------------------
# Save output
# ------------------------------------------------------------
write.csv(
  table2_results,
  file = file.path(results_dir, "table2_cd4_results.csv"),
  row.names = FALSE
)

print(table2_results)