# ============================================================
# Reproduce Figures S.1--S.3 Data
# HCP-repro
#
# This script generates conditional-coverage simulation data
# for Figures S.1--S.3.
#
# Figure S.1: n = 300, n_local_groups = 1 vs 10
# Figure S.2: n = 100, n_local_groups = 1 vs 5
# Figure S.3: n = 500, n_local_groups = 1 vs 20
#
# Each replicate evaluates conditional coverage and length over
# fixed test covariates:
#   X = (0.6, x2, 0, 0), x2 in [-3, 3].
#
# Outputs are saved to results/:
#   figureS1_conditional_data.csv
#   figureS2_conditional_data.csv
#   figureS3_conditional_data.csv
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

source(file.path(r_dir, "build_local_groups.R"), chdir = TRUE)
source(file.path(r_dir, "hcp_localized_region.R"), chdir = TRUE)

source(file.path(r_dir, "run_one_replicate_conditional.R"), chdir = TRUE)



#' Generate conditional-coverage data for one supplementary figure
#'
#' @description
#' This helper generates conditional-coverage simulation data for one of
#' Figures S.1--S.3. It evaluates localized HCP at fixed test covariates over
#' a grid of x2 values and compares two choices of the number of local groups.
#'
#' @param figure_id Figure identifier used in the output file name.
#' @param n Number of training subjects.
#' @param local_group_values Two choices of the number of local groups.
#' @param num_simulations_run Number of simulation replicates.
#' @param num_cores Number of parallel workers.
#' @param scenario Simulation scenario.
#' @param alpha Miscoverage level.
#' @param n_grid Number of response-grid points.
#' @param x2_grid Grid of fixed x2 values used for conditional evaluation.
#' @param train_frac Fraction of subjects assigned to training.
#' @param S Number of repeated subject-level splits.
#' @param B Number of repeated calibration subsamples per split.
#' @param theta_mean Mean vector of subject-specific random effects.
#' @param theta_sigma Covariance matrix of subject-specific random effects.
#' @param x_mean Mean vector for covariate generation.
#' @param x_sigma Covariance matrix for covariate generation.
#' @param beta Coefficient vector for the logistic missingness model.
#' @param weight_cap Upper bound for inverse-propensity weights.
#' @param output_dir Directory where the CSV output is saved.
#' @param run_sanity_check Logical; if `TRUE`, run one test replicate before
#'   the full simulation.
#'
#' @return Invisibly returns a data frame containing conditional coverage and
#' length summaries over the fixed x2 grid.
run_one_conditional_figure <- function(
    figure_id,
    n,
    local_group_values,
    num_simulations_run = 1000,
    num_cores = 7L,
    scenario = "Homo",
    alpha = 0.1,
    n_grid = 200,
    x2_grid = seq(-3, 3, length.out = 30),
    train_frac = 0.3,
    S = 2,
    B = 2,
    theta_mean = rep(2, 4),
    theta_sigma = diag(1, 4),
    x_mean = rep(0, 4),
    x_sigma = diag(1, 4),
    beta = c(3, 0, 2, 2, 2),
    weight_cap = 30,
    output_dir = results_dir,
    run_sanity_check = FALSE
) {
  if (!is.character(figure_id) || length(figure_id) != 1L) {
    stop("figure_id must be a single character string.")
  }
  
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 2 || n != floor(n)) {
    stop("n must be an integer >= 2.")
  }
  
  if (!is.numeric(local_group_values) ||
      anyNA(local_group_values) ||
      any(!is.finite(local_group_values)) ||
      any(local_group_values < 1) ||
      any(local_group_values != floor(local_group_values))) {
    stop("local_group_values must contain positive integers.")
  }
  
  if (!is.numeric(num_simulations_run) ||
      length(num_simulations_run) != 1L ||
      !is.finite(num_simulations_run) ||
      num_simulations_run < 1 ||
      num_simulations_run != floor(num_simulations_run)) {
    stop("num_simulations_run must be a positive integer.")
  }
  
  num_cores <- as.integer(max(1L, num_cores))
  
  custom_exports <- c(
    "run_one_replicate_conditional",
    "generate_theta",
    "generate_subject_covariates",
    "generate_errors",
    "generate_outcomes",
    "generate_missing_indicator",
    "generate_simulation_data",
    "generate_Y_given_X",
    "evaluate_interval_region",
    "evaluate_bimodal_region",
    "hcp_localized_region",
    "build_local_groups",
    "fit_cond_density_qp",
    "fit_propensity_model",
    "adjust_bandwidth_bimodal",
    "quantile_levels",
    "get_S",
    "get_b",
    "get_max_wt"
  )
  
  if (isTRUE(run_sanity_check)) {
    cat("Running one sanity-check replicate for ", figure_id, "...\n", sep = "")
    print(
      run_one_replicate_conditional(
        n = n,
        scenario = scenario,
        n_local_groups = local_group_values[1],
        alpha = alpha,
        n_grid = n_grid,
        x2_grid = x2_grid,
        theta_mean = theta_mean,
        theta_sigma = theta_sigma,
        x_mean = x_mean,
        x_sigma = x_sigma,
        beta = beta,
        train_frac = train_frac,
        S = S,
        B = B,
        weight_cap = weight_cap,
        seed = 1
      )
    )
  }
  
  cat(sprintf(
    "\n===== Running %s: n = %d, local groups = %s =====\n",
    figure_id,
    n,
    paste(local_group_values, collapse = ", ")
  ))
  cat("Using", num_cores, "cores\n")
  flush.console()
  
  foreach::registerDoSEQ()
  cl <- parallel::makeCluster(num_cores, outfile = "")
  doParallel::registerDoParallel(cl)
  
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    try(foreach::registerDoSEQ(), silent = TRUE)
  }, add = TRUE)
  
  all_outputs <- list()
  
  for (K_local in local_group_values) {
    cat(sprintf(
      "Running %s, n = %d, n_local_groups = %d\n",
      figure_id, n, K_local
    ))
    flush.console()
    
    res_list <- foreach(
      i = seq_len(num_simulations_run),
      .combine = "c",
      .export = custom_exports,
      .errorhandling = "pass"
    ) %dopar% {
      tryCatch(
        list(
          run_one_replicate_conditional(
            n = n,
            scenario = scenario,
            n_local_groups = K_local,
            alpha = alpha,
            n_grid = n_grid,
            x2_grid = x2_grid,
            theta_mean = theta_mean,
            theta_sigma = theta_sigma,
            x_mean = x_mean,
            x_sigma = x_sigma,
            beta = beta,
            train_frac = train_frac,
            S = S,
            B = B,
            weight_cap = weight_cap,
            seed = i
          )
        ),
        error = function(e) {
          message(sprintf(
            "Replicate %d failed for %s, n=%d, n_local_groups=%d: %s",
            i, figure_id, n, K_local, conditionMessage(e)
          ))
          
          list(
            matrix(
              NA_real_,
              nrow = length(x2_grid),
              ncol = 2,
              dimnames = list(NULL, c("covered", "length"))
            )
          )
        }
      )
    }
    
    # --------------------------------------------------------
    # Aggregate over replicates
    # --------------------------------------------------------
    cover_mat <- do.call(cbind, lapply(res_list, function(z) z[, "covered"]))
    len_mat   <- do.call(cbind, lapply(res_list, function(z) z[, "length"]))
    
    mean_cover <- rowMeans(cover_mat, na.rm = TRUE)
    mean_len   <- rowMeans(len_mat, na.rm = TRUE)
    
    sd_cover <- sqrt(mean_cover * (1 - mean_cover) / num_simulations_run)
    sd_len   <- apply(len_mat, 1, stats::sd, na.rm = TRUE) /
      sqrt(rowSums(!is.na(len_mat)))
    
    out_K <- data.frame(
      figure = figure_id,
      n = n,
      scenario = scenario,
      n_local_groups = K_local,
      x2 = x2_grid,
      Mean_Cover = mean_cover,
      SD_Cover = sd_cover,
      Mean_Leng = mean_len,
      SD_Leng = sd_len,
      stringsAsFactors = FALSE
    )
    
    all_outputs[[paste0("K", K_local)]] <- out_K
    
    cat(sprintf(
      "Completed %s, n = %d, n_local_groups = %d\n",
      figure_id, n, K_local
    ))
    flush.console()
  }
  
  final_df <- do.call(rbind, all_outputs)
  rownames(final_df) <- NULL
  
  output_file <- file.path(
    output_dir,
    paste0(figure_id, "_conditional_data.csv")
  )
  
  write.csv(final_df, file = output_file, row.names = FALSE)
  cat("Saved:", output_file, "\n")
  
  invisible(final_df)
}

# ------------------------------------------------------------
# Run Figures S.1--S.3
# ------------------------------------------------------------
figureS1_data <- run_one_conditional_figure(
  figure_id = "figureS1",
  n = 300,
  local_group_values = c(1, 10),
  num_simulations_run = 1000,
  num_cores = 7,
  scenario = "Homo",
  train_frac = 0.3,
  S = 2,
  B = 2,
  weight_cap = 30,
  output_dir = results_dir
)

figureS2_data <- run_one_conditional_figure(
  figure_id = "figureS2",
  n = 100,
  local_group_values = c(1, 5),
  num_simulations_run = 1000,
  num_cores = 7,
  scenario = "Homo",
  train_frac = 0.3,
  S = 2,
  B = 2,
  weight_cap = 30,
  output_dir = results_dir
)

figureS3_data <- run_one_conditional_figure(
  figure_id = "figureS3",
  n = 500,
  local_group_values = c(1, 20),
  num_simulations_run = 1000,
  num_cores = 7,
  scenario = "Homo",
  train_frac = 0.3,
  S = 2,
  B = 2,
  weight_cap = 30,
  output_dir = results_dir
)