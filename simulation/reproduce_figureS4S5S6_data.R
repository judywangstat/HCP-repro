# ============================================================
# Reproduce Figures S.4--S.6 Data
# HCPclust-repro
#
# This script generates simulation data for Figures S.4--S.6.
# These figures evaluate local coverage and local length over a
# two-dimensional test grid:
#
#   X = (0.6, x2, x3, 0),  x2, x3 in [-3, 3].
#
# For each figure, the script compares the non-localized version
# with a localized version using profile-based groups.
#
# Figure S.4: n = 300, K = 1 vs K = 10
# Figure S.5: n = 100, K = 1 vs K = 5
# Figure S.6: n = 500, K = 1 vs K = 20
#
# Outputs are saved to results/:
#   figureS4_local_data.csv
#   figureS5_local_data.csv
#   figureS6_local_data.csv
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

source(file.path(r_dir, "build_local_groups.R"), chdir = TRUE)
source(file.path(r_dir, "hcp_localized_region.R"), chdir = TRUE)

source(file.path(r_dir, "run_one_replicate_local.R"), chdir = TRUE)



#' Generate local-grid simulation data for one supplementary figure
#'
#' @description
#' This helper generates local-grid simulation data for one of Figures S.4--S.6.
#' It evaluates local coverage and local prediction-region length over a
#' two-dimensional test grid and compares different choices of the number of
#' local groups.
#'
#' @param figure_id Figure identifier used in the output file name.
#' @param n Number of training subjects.
#' @param local_group_values Choices of the number of local groups.
#' @param num_simulations_run Number of simulation replicates.
#' @param num_cores Number of parallel workers.
#' @param scenario Simulation scenario.
#' @param alpha Miscoverage level.
#' @param n_grid Number of response-grid points.
#' @param x2_grid Grid of fixed x2 values.
#' @param x3_grid Grid of fixed x3 values.
#' @param x2_breaks Breakpoints used to define x2 local cells.
#' @param x3_breaks Breakpoints used to define x3 local cells.
#' @param theta_mean Mean vector of subject-specific random effects.
#' @param theta_sigma Covariance matrix of subject-specific random effects.
#' @param x_mean Mean vector for covariate generation.
#' @param x_sigma Covariance matrix for covariate generation.
#' @param beta Coefficient vector for the logistic missingness model.
#' @param train_frac Fraction of subjects assigned to training.
#' @param S Number of repeated subject-level splits.
#' @param B Number of repeated calibration subsamples per split.
#' @param weight_cap Upper bound for inverse-propensity weights.
#' @param output_dir Directory where the CSV output is saved.
#' @param run_sanity_check Logical; if `TRUE`, run one test replicate before
#'   the full simulation.
#'
#' @return Invisibly returns a data frame containing local coverage and
#' prediction-region length summaries over the two-dimensional test grid.
run_one_local_figure <- function(
    figure_id,
    n,
    local_group_values,
    num_simulations_run = 1000,
    num_cores = 7L,
    scenario = "Homo",
    alpha = 0.1,
    n_grid = 200,
    x2_grid = seq(-3, 3, length.out = 42),
    x3_grid = seq(-3, 3, length.out = 42),
    x2_breaks = seq(-3, 3, length.out = 7),
    x3_breaks = seq(-3, 3, length.out = 7),
    theta_mean = rep(2, 4),
    theta_sigma = diag(1, 4),
    x_mean = rep(0, 4),
    x_sigma = diag(1, 4),
    beta = c(3, 0, 2, 2, 2),
    train_frac = 0.3,
    S = 1,
    B = 1,
    weight_cap = 30,
    output_dir = results_dir,
    run_sanity_check = FALSE
) {
  # ----------------------------------------------------------
  # Step 0: Basic input checks
  # ----------------------------------------------------------
  if (!is.character(figure_id) || length(figure_id) != 1L) {
    stop("figure_id must be a single character string.")
  }
  
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 2 || n != floor(n)) {
    stop("n must be an integer >= 2.")
  }
  
  if (!is.numeric(local_group_values) ||
      length(local_group_values) < 1L ||
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
  
  if (!is.numeric(num_cores) ||
      length(num_cores) != 1L ||
      !is.finite(num_cores)) {
    stop("num_cores must be a single finite number.")
  }
  
  num_cores <- as.integer(max(1L, num_cores))
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  custom_exports <- c(
    "run_one_replicate_local",
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
  
  # ----------------------------------------------------------
  # Step 1: Optional sanity check before launching parallel jobs
  # ----------------------------------------------------------
  if (isTRUE(run_sanity_check)) {
    cat("Running one sanity-check replicate for ", figure_id, "...\n", sep = "")
    
    sanity_out <- run_one_replicate_local(
      n = n,
      scenario = scenario,
      n_local_groups = local_group_values[1],
      alpha = alpha,
      n_grid = n_grid,
      x2_grid = x2_grid,
      x3_grid = x3_grid,
      x2_breaks = x2_breaks,
      x3_breaks = x3_breaks,
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
    
    print(head(sanity_out))
  }
  
  cat(sprintf(
    "\n===== Running %s: n = %d, local groups = %s =====\n",
    figure_id,
    n,
    paste(local_group_values, collapse = ", ")
  ))
  cat("Using", num_cores, "cores\n")
  flush.console()
  
  # ----------------------------------------------------------
  # Step 2: Start parallel backend
  # ----------------------------------------------------------
  foreach::registerDoSEQ()
  
  cl <- parallel::makeCluster(num_cores, outfile = "")
  doParallel::registerDoParallel(cl)
  parallel::clusterSetRNGStream(cl, iseed = 123)
  
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    try(foreach::registerDoSEQ(), silent = TRUE)
  }, add = TRUE)
  
  all_outputs <- list()
  
  # ----------------------------------------------------------
  # Step 3: Run all replicates for each localization setting
  # ----------------------------------------------------------
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
          run_one_replicate_local(
            n = n,
            scenario = scenario,
            n_local_groups = K_local,
            alpha = alpha,
            n_grid = n_grid,
            x2_grid = x2_grid,
            x3_grid = x3_grid,
            x2_breaks = x2_breaks,
            x3_breaks = x3_breaks,
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
          
          n_cells <- (length(x2_breaks) - 1L) * (length(x3_breaks) - 1L)
          
          x2_mid <- rep(
            (x2_breaks[-length(x2_breaks)] + x2_breaks[-1]) / 2,
            times = length(x3_breaks) - 1L
          )
          x3_mid <- rep(
            (x3_breaks[-length(x3_breaks)] + x3_breaks[-1]) / 2,
            each = length(x2_breaks) - 1L
          )
          
          list(
            data.frame(
              cell_id = seq_len(n_cells),
              x2_mid = x2_mid,
              x3_mid = x3_mid,
              covered = NA_real_,
              length = NA_real_
            )
          )
        }
      )
    }
    
    # --------------------------------------------------------
    # Step 4: Aggregate replicate-level cell summaries
    # --------------------------------------------------------
    cover_mat <- do.call(cbind, lapply(res_list, function(z) z$covered))
    len_mat <- do.call(cbind, lapply(res_list, function(z) z$length))
    
    cell_template <- res_list[[1]][, c("cell_id", "x2_mid", "x3_mid")]
    
    mean_cover <- rowMeans(cover_mat, na.rm = TRUE)
    mean_len <- rowMeans(len_mat, na.rm = TRUE)
    
    sd_cover <- apply(cover_mat, 1, stats::sd, na.rm = TRUE) /
      sqrt(rowSums(!is.na(cover_mat)))
    sd_len <- apply(len_mat, 1, stats::sd, na.rm = TRUE) /
      sqrt(rowSums(!is.na(len_mat)))
    
    out_K <- data.frame(
      figure = figure_id,
      n = n,
      scenario = scenario,
      n_local_groups = K_local,
      cell_template,
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
  
  # ----------------------------------------------------------
  # Step 5: Save figure-level data
  # ----------------------------------------------------------
  final_df <- do.call(rbind, all_outputs)
  rownames(final_df) <- NULL
  
  output_file <- file.path(
    output_dir,
    paste0(figure_id, "_local_data.csv")
  )
  
  write.csv(final_df, file = output_file, row.names = FALSE)
  cat("Saved:", output_file, "\n")
  
  invisible(final_df)
}

# ------------------------------------------------------------
# Run Figures S.4--S.6
# ------------------------------------------------------------

figureS4_data <- run_one_local_figure(
  figure_id = "figureS4",
  n = 300,
  local_group_values = c(1, 10),
  num_simulations_run = 1000,
  num_cores = 7,
  scenario = "Homo",
  train_frac = 0.3,
  S = 1,
  B = 1,
  weight_cap = 30,
  output_dir = results_dir
)

figureS5_data <- run_one_local_figure(
  figure_id = "figureS5",
  n = 100,
  local_group_values = c(1, 5),
  num_simulations_run = 1000,
  num_cores = 7,
  scenario = "Homo",
  train_frac = 0.3,
  S = 1,
  B = 1,
  weight_cap = 30,
  output_dir = results_dir
)

figureS6_data <- run_one_local_figure(
  figure_id = "figureS6",
  n = 500,
  local_group_values = c(1, 20),
  num_simulations_run = 1000,
  num_cores = 7,
  scenario = "Homo",
  train_frac = 0.3,
  S = 1,
  B = 1,
  weight_cap = 30,
  output_dir = results_dir
)