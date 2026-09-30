# ============================================================
# Reproduce Figure 2 Data
# CD4 real-data prediction bands
# HCP-repro
#
# This script generates leave-one-subject-out prediction bands
# for the CD4 real-data example.
#
# Default setting:
#   missing_rate = "20"
#   prediction_type = "pointwise"
#
# Output:
#   results/figure2_cd4_band_data.csv
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
repo_dir <- normalizePath(Sys.getenv("HCP_REPO_DIR", "."), mustWork = TRUE)
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
source(file.path(r_dir, "dwr_region_realdata.R"), chdir = TRUE)
source(file.path(r_dir, "lc_region_realdata.R"), chdir = TRUE)
source(file.path(r_dir, "lmem_region.R"), chdir = TRUE)

# ------------------------------------------------------------
# Main settings
# ------------------------------------------------------------
missing_rate <- "20"
prediction_type <- "pointwise"
alpha <- 0.1
n_grid <- 200
seed <- 123
num_cores <- 7

methods <- c("HCP", "DWR", "LC")

# ------------------------------------------------------------
# Read data and generate artificial missingness
# ------------------------------------------------------------
cd4_dat <- read_cd4_data(file.path(data_dir, "CD4_data.txt"))
cd4_dat <- add_cd4_missingness(
  dat = cd4_dat,
  missing_rate = missing_rate,
  seed = seed
)

test_ids <- sort(unique(cd4_dat$id))
x_cols <- get_cd4_x_cols()


#' Run one CD4 leave-one-subject-out fit and return prediction-band data
#'
#' @description
#' This helper runs one leave-one-subject-out analysis for a selected CD4
#' subject and method. It returns the observed outcomes and prediction-band
#' endpoints used to plot Figure 2.
#'
#' @param dat CD4 data frame with artificial missingness indicator.
#' @param test_id Held-out subject ID.
#' @param method Method used to construct prediction bands.
#' @param prediction_type Prediction type, either `"pointwise"` or
#'   `"simultaneous"`.
#' @param missing_rate Artificial missingness setting, either `"20"` or `"50"`.
#' @param alpha Miscoverage level.
#' @param n_grid Number of response-grid points.
#' @param seed Optional random seed.
#'
#' @return A data frame containing the held-out subject ID, time, observed
#' outcome, lower and upper prediction-band endpoints, method name, missingness
#' setting, and prediction type.
run_one_cd4_band <- function(
    dat,
    test_id,
    method = c("HCP", "DWR", "LC", "LMEM"),
    prediction_type = c("pointwise", "simultaneous"),
    missing_rate = c("20", "50"),
    alpha = 0.1,
    n_grid = 200,
    seed = NULL
) {
  method <- match.arg(method)
  prediction_type <- match.arg(prediction_type)
  missing_rate <- match.arg(as.character(missing_rate), choices = c("20", "50"))
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  x_cols <- get_cd4_x_cols()
  
  dat_test <- dat[dat$id == test_id, , drop = FALSE]
  dat_sample <- dat[dat$id != test_id, , drop = FALSE]
  
  y_true <- as.numeric(dat_test$Y)
  x_test <- as.matrix(dat_test[, x_cols, drop = FALSE])
  storage.mode(x_test) <- "numeric"
  colnames(x_test) <- x_cols
  
  y_grid <- make_cd4_y_grid(
    dat_sample = dat_sample,
    y_col = "Y",
    n_grid = n_grid
  )
  
  alpha_use <- if (prediction_type == "simultaneous") {
    alpha / nrow(dat_test)
  } else {
    alpha
  }
  
  weight_cap_use <- get_cd4_weight_cap(
    method = method,
    missing_rate = missing_rate
  )
  
  res <- switch(
    method,
    
    HCP = hcp_region(
      dat = dat_sample,
      id_col = "id",
      y_col = "Y",
      delta_col = "delta",
      x_cols = x_cols,
      x_test = x_test,
      y_grid = y_grid,
      alpha = alpha_use,
      train_frac = 0.3,
      S = 5,
      B = 5,
      combine_B = "cct",
      combine_S = "cct",
      dens_method = "grf",
      dens_taus = (1:(2^4 - 1)) / (2^4),
      dens_h = NULL,
      dens_h_scenario = NULL,
      enforce_monotone = FALSE,
      tail_decay = TRUE,
      prop_method = "logistic",
      weight_cap = 5,
      seed = seed
    ),
    
    DWR = dwr_region_realdata(
      dat = dat_sample,
      id_col = "id",
      y_col = "Y",
      delta_col = "delta",
      x_cols = x_cols,
      x_test = x_test,
      y_grid = y_grid,
      alpha = alpha_use,
      train_frac = 0.5,
      B = 1,
      seed = seed,
      score_type = "residual",
      reg_method = "grf",
      quant_method = "grf",
      lower_tau = 0.05,
      upper_tau = 0.95
    ),
    
    LC = lc_region_realdata(
      dat = dat_sample,
      id_col = "id",
      y_col = "Y",
      delta_col = "delta",
      x_cols = x_cols,
      x_test = x_test,
      y_grid = y_grid,
      alpha = alpha_use,
      train_frac = 0.5,
      seed = seed,
      quant_method = "grf",
      prop_method = "logistic",
      lower_tau = 0.05,
      upper_tau = 0.95,
      weight_cap = weight_cap_use,
      time_matched = FALSE
    ),
    
    LMEM = lmem_region(
      dat = dat_sample, x_test = dat_test, setting = "cd4",
      alpha = alpha_use, n_sims = 10000L, seed = seed
    )
  )
  
  data.frame(
    id = dat_test$id,
    time = dat_test$time,
    true = y_true,
    lower = res$lo_hi[, "lo"],
    upper = res$lo_hi[, "hi"],
    method = method,
    missing_rate = missing_rate,
    prediction_type = prediction_type,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Parallel computation
# ------------------------------------------------------------
custom_exports <- c(
  "cd4_dat",
  "test_ids",
  "methods",
  "missing_rate",
  "prediction_type",
  "alpha",
  "n_grid",
  "seed",
  "run_one_cd4_band",
  "read_cd4_data",
  "add_cd4_missingness",
  "get_cd4_x_cols",
  "get_cd4_missing_beta",
  "get_cd4_weight_cap",
  "generate_cd4_missing_indicator",
  "make_cd4_y_grid",
  "evaluate_interval_region",
  "hcp_region",
  "dwr_region_realdata",
  "lc_region_realdata",
  "lmem_region", "lmem_random_covariance", "lmem_covariance_root", "lmem_fit_diagnostics",
  "fit_cond_density_qp",
  "fit_propensity_model",
  "quantile_levels"
)

cl <- parallel::makeCluster(num_cores, outfile = "")
doParallel::registerDoParallel(cl)
parallel::clusterSetRNGStream(cl, iseed = seed)

on.exit({
  try(parallel::stopCluster(cl), silent = TRUE)
  try(foreach::registerDoSEQ(), silent = TRUE)
}, add = TRUE)

figure2_data <- foreach(
  ii = seq_along(test_ids),
  .combine = rbind,
  .errorhandling = "pass",
  .export = custom_exports,
  .packages = c("grf", "quantreg", "lme4")
) %dorng% {
  test_id <- test_ids[ii]
  
  do.call(
    rbind,
    lapply(methods, function(method_i) {
      tryCatch(
        run_one_cd4_band(
          dat = cd4_dat,
          test_id = test_id,
          method = method_i,
          prediction_type = prediction_type,
          missing_rate = missing_rate,
          alpha = alpha,
          n_grid = n_grid,
          seed = seed + 1000L * ii + match(method_i, methods)
        ),
        error = function(e) {
          message(sprintf(
            "Failed: test_id=%s, method=%s: %s",
            as.character(test_id), method_i, conditionMessage(e)
          ))
          NULL
        }
      )
    })
  )
}

# ------------------------------------------------------------
# Save output
# ------------------------------------------------------------
output_file <- file.path(results_dir, "figure2_cd4_band_data.csv")

write.csv(
  figure2_data,
  file = output_file,
  row.names = FALSE
)

cat("Saved Figure 2 data to:\n", output_file, "\n")
print(head(figure2_data))