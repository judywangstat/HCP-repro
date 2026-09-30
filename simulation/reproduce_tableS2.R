# ============================================================
# Reproduce Table S.2
# Simultaneous prediction simulation study
# HCP-repro
#
# This script reproduces Table S.2 in the supplementary materials.
# It evaluates simultaneous prediction using run_one_replicate_simultaneous(),
# where each replicate returns one subject-level simultaneous coverage
# indicator and the average prediction-region length over the held-out subject.
#
# Output:
#   results/tableS2_final_results.csv
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


#' Run the Table S.2 simultaneous-prediction simulation study
#'
#' @description
#' This helper runs the simultaneous-prediction simulation study across
#' missingness levels, sample sizes, and simulation scenarios. For each setting,
#' it runs multiple simulation replicates, evaluates HCP, and aggregates
#' subject-level simultaneous coverage and average prediction-region length.
#'
#' @param num_simulations_run Number of simulation replicates per setting.
#' @param num_cores Number of parallel workers.
#' @param sample_sizes Sample sizes to evaluate.
#' @param scenarios Simulation scenarios to evaluate.
#' @param missing_levels Missingness levels to evaluate.
#' @param alpha Miscoverage level.
#' @param n_test Number of held-out test subjects per replicate.
#' @param n_grid Number of response-grid points.
#' @param train_frac Fraction of subjects assigned to training in each split.
#' @param theta_mean Mean vector of subject-specific random effects.
#' @param theta_sigma Covariance matrix of subject-specific random effects.
#' @param x_mean Mean vector for covariate generation.
#' @param x_sigma Covariance matrix for covariate generation.
#' @param weight_cap_rule Function returning the weight truncation level as a
#'   function of sample size for the 50% missingness setting.
#' @param output_file Path to save the final CSV summary.
#' @param run_sanity_check Logical; if `TRUE`, run one test replicate before
#'   the full simulation.
#'
#' @return Invisibly returns a data frame containing simultaneous coverage and
#' prediction-region length summaries for each setting.

run_tableS2_parallel <- function(
    num_simulations_run = 1000,
    num_cores = 7L,
    sample_sizes = c(100, 300, 500),
    scenarios = c("Homo", "Heter", "Asym", "Bimo", "Bimo-fix"),
    missing_levels = c("20", "50"),
    alpha = 0.1,
    n_test = 1,
    n_grid = 200,
    train_frac = 0.3,
    theta_mean = rep(2, 4),
    theta_sigma = diag(1, 4),
    x_mean = rep(0, 4),
    x_sigma = diag(1, 4),
    weight_cap_rule = function(n) {
      pmin(28, 2.2 * (n / 100)^1.72)
    },
    output_file = file.path(results_dir, "tableS2_final_results.csv"),
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

  if (isTRUE(run_sanity_check)) {
    cat("Running one sanity-check replicate...\n")
    print(
      run_one_replicate_simultaneous(
        n = sample_sizes[1],
        scenario = scenarios[1],
        missing = missing_levels[1],
        alpha = alpha,
        n_test = n_test,
        n_grid = n_grid,
        train_frac = train_frac,
        weight_cap = 30,
        seed = 1
      )
    )
  }

  results_info <- expand.grid(
    missing = missing_levels,
    sample_size = sample_sizes,
    stringsAsFactors = FALSE
  )

  results_info$missing <- factor(results_info$missing, levels = missing_levels)
  results_info$sample_size <- factor(results_info$sample_size, levels = sample_sizes)
  results_info <- results_info[order(results_info$missing, results_info$sample_size), , drop = FALSE]
  rownames(results_info) <- NULL

  results_all <- matrix(
    NA_real_,
    nrow = nrow(results_info),
    ncol = 2 * length(scenarios)
  )

  colnames(results_all) <- as.vector(rbind(
    paste0(scenarios, "_cov"),
    paste0(scenarios, "_len")
  ))

  cat("Using", num_cores, "cores\n")

  foreach::registerDoSEQ()
  cl <- parallel::makeCluster(num_cores, outfile = "")
  doParallel::registerDoParallel(cl)

  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    try(foreach::registerDoSEQ(), silent = TRUE)
  }, add = TRUE)

  for (row_index in seq_len(nrow(results_info))) {
    n_val <- as.numeric(as.character(results_info$sample_size[row_index]))
    miss <- as.character(results_info$missing[row_index])

    scenario_results <- numeric(2 * length(scenarios))
    col_index <- 1L

    for (scen in scenarios) {
      beta_use <- switch(
        miss,
        "20" = c(3, 0, 2, 2, 2),
        "50" = c(0, 0, 2, 2, 2)
      )

      weight_cap_use <- switch(
        miss,
        "20" = 30,
        "50" = weight_cap_rule(n_val)
      )

      cat(sprintf(
        "Running missing = %s%%, n = %d, scenario = %s, weight_cap = %.4f\n",
        miss, n_val, scen, weight_cap_use
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
            missing = miss,
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
              "Replicate %d failed for missing=%s, n=%d, scenario=%s: %s",
              i, miss, n_val, scen, conditionMessage(e)
            ))
            c(HCP_cov = NA_real_, HCP_len = NA_real_)
          }
        )
      }

      if (is.null(dim(res_mat))) {
        res_mat <- matrix(res_mat, nrow = 1L, dimnames = list(NULL, names(res_mat)))
      }
      num_failed <- sum(!stats::complete.cases(res_mat))
      cat(sprintf(
        "Completed missing = %s%%, n = %d, scenario = %s with %d failed replicates\n",
        miss, n_val, scen, num_failed
      ))

      vals <- colMeans(res_mat, na.rm = TRUE)
      scenario_results[col_index] <- vals["HCP_cov"]
      scenario_results[col_index + 1L] <- vals["HCP_len"]

      col_index <- col_index + 2L
    }

    results_all[row_index, ] <- scenario_results
  }

  final_results <- cbind(
    data.frame(
      missing = paste0(as.character(results_info$missing), "%"),
      sample_size = as.numeric(as.character(results_info$sample_size)),
      stringsAsFactors = FALSE
    ),
    results_all
  )

  colnames(final_results) <- gsub("Bimo-fix", "Bimo-Fix", colnames(final_results), fixed = TRUE)

  print(final_results)

  if (!is.null(output_file)) {
    write.csv(final_results, file = output_file, row.names = FALSE)
    cat("Saved:", output_file, "\n")
  }

  invisible(final_results)
}

# ------------------------------------------------------------
# Run Table S.2 simulation
# ------------------------------------------------------------
if (sys.nframe() == 0L && "--run" %in% commandArgs(TRUE)) {
tableS2_results <- run_tableS2_parallel(
  num_simulations_run = 1000,
  num_cores = 7,
  sample_sizes = c(100, 300, 500),
  scenarios = c("Homo", "Heter", "Asym", "Bimo", "Bimo-fix"),
  missing_levels = c("20", "50"),
  train_frac = 0.3,
  output_file = file.path(results_dir, "tableS2_final_results.csv")
)
}
