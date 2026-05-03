# ============================================================
# Reproduce Figure S.11 Data
# Gallstones real-data prediction bands
# HCPclust-repro
#
# This script generates the leave-one-subject-out prediction-band
# data for the gallstones real-data analysis. The saved output is
# used by real_data/plot_figureS11_gallstones.R to reproduce
# Figure S.11 in the supplementary materials.
#
# Input:
#   data/gallstones.txt
#
# Output:
#   results/figureS11_gallstones_band_data.csv
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
source(file.path(r_dir, "fit_propensity_model.R"), chdir = TRUE)
source(file.path(r_dir, "fit_cond_density_qp.R"), chdir = TRUE)
source(file.path(r_dir, "tuning_helpers.R"), chdir = TRUE)
source(file.path(r_dir, "evaluation_helpers.R"), chdir = TRUE)

source(file.path(r_dir, "hcp_region.R"), chdir = TRUE)
source(file.path(r_dir, "dwr_region_legacy.R"), chdir = TRUE)
source(file.path(r_dir, "lc_region_legacy.R"), chdir = TRUE)
source(file.path(r_dir, "lmem_region.R"), chdir = TRUE)

source(file.path(r_dir, "gallstones_data_helpers.R"), chdir = TRUE)

# ------------------------------------------------------------
# Main settings
# ------------------------------------------------------------
prediction_type <- "pointwise"
alpha <- 0.1
n_grid <- 200
seed <- 123
impute_seed <- 12345
noise_scale <- 1.2
num_cores <- 7

methods <- c("HCP", "DWR", "LC")
x_cols <- c("T1", "T2", "T3", "Treat")

# ------------------------------------------------------------
# Read data and impute missing outcomes for evaluation
# ------------------------------------------------------------
gall_dat <- load_gallstones_data(file.path(data_dir, "gallstones.txt"))

gall_dat <- impute_gallstones_outcomes(
  gall_dat,
  seed = impute_seed,
  noise_scale = noise_scale
)

cat(sprintf(
  "Original missing rate: %.2f%%\n",
  100 * mean(gall_dat$delta == 0)
))
print(with(gall_dat, table(time, delta)))

test_ids <- sort(unique(gall_dat$id))

# ------------------------------------------------------------
# Run one leave-one-subject-out fit and return prediction-band data
# ------------------------------------------------------------
run_one_gallstones_band <- function(
    dat,
    test_id,
    method = c("HCP", "DWR", "LC", "LMEM"),
    prediction_type = c("pointwise", "simultaneous"),
    alpha = 0.1,
    n_grid = 200,
    x_cols = c("T1", "T2", "T3", "Treat"),
    seed = NULL
) {
  method <- match.arg(method)
  prediction_type <- match.arg(prediction_type)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  dat_test <- dat[dat$id == test_id, , drop = FALSE]
  dat_sample <- dat[dat$id != test_id, , drop = FALSE]
  
  y_true <- as.numeric(dat_test$Y_eval)
  
  x_test <- as.matrix(dat_test[, x_cols, drop = FALSE])
  storage.mode(x_test) <- "numeric"
  colnames(x_test) <- x_cols
  
  y_grid <- make_gallstones_y_grid(
    dat_sample = dat_sample,
    y_col = "Y",
    n_grid = n_grid
  )
  
  alpha_use <- if (prediction_type == "simultaneous") {
    alpha / nrow(dat_test)
  } else {
    alpha
  }
  
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
      S = 12,
      B = 5,
      dens_method = "rq",
      dens_taus = (1:(2^4 - 1)) / (2^4),
      dens_h = NULL,
      dens_h_scenario = NULL,
      prop_method = "logistic",
      weight_cap = 50,
      seed = seed
    ),
    
    DWR = dwr_region_legacy(
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
      score_type = "residual",
      reg_method = "linear",
      seed = seed
    ),
    
    LC = lc_region_legacy(
      dat = dat_sample,
      id_col = "id",
      y_col = "Y",
      delta_col = "delta",
      x_cols = x_cols,
      x_test = x_test,
      y_grid = y_grid,
      alpha = alpha_use,
      train_frac = 0.5,
      quant_method = "linear",
      prop_method = "logistic",
      weight_cap = 50,
      time_matched = FALSE,
      seed = seed
    ),
    
    LMEM = lmem_region(
      dat = dat_sample,
      id_col = "id",
      y_col = "Y",
      delta_col = "delta",
      x_cols = x_cols,
      x_test = x_test,
      y_grid = y_grid,
      alpha = alpha_use,
      fixed_formula = "T1 + T2 + T3 + Treat",
      random_formula = "(1 | id)",
      n_sims = 1000,
      pred_which = "full",
      seed = seed
    )
  )
  
  data.frame(
    id = dat_test$id,
    time = dat_test$time,
    true = y_true,
    lower = res$lo_hi[, "lo"],
    upper = res$lo_hi[, "hi"],
    method = method,
    prediction_type = prediction_type,
    delta = dat_test$delta,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Parallel computation
# ------------------------------------------------------------
custom_exports <- c(
  "gall_dat",
  "test_ids",
  "methods",
  "prediction_type",
  "alpha",
  "n_grid",
  "seed",
  "x_cols",
  "run_one_gallstones_band",
  "make_gallstones_y_grid",
  "hcp_region",
  "dwr_region_legacy",
  "lc_region_legacy",
  "lmem_region",
  "fit_cond_density_qp",
  "fit_propensity_model",
  "quantile_levels",
  "evaluate_interval_region"
)

cl <- parallel::makeCluster(num_cores, outfile = "")
doParallel::registerDoParallel(cl)
parallel::clusterSetRNGStream(cl, iseed = seed)

on.exit({
  try(parallel::stopCluster(cl), silent = TRUE)
  try(foreach::registerDoSEQ(), silent = TRUE)
}, add = TRUE)

figureS11_data <- foreach(
  ii = seq_along(test_ids),
  .combine = rbind,
  .errorhandling = "pass",
  .export = custom_exports,
  .packages = c("quantreg", "lme4", "merTools")
) %dorng% {
  test_id <- test_ids[ii]
  
  do.call(
    rbind,
    lapply(methods, function(method_i) {
      tryCatch(
        run_one_gallstones_band(
          dat = gall_dat,
          test_id = test_id,
          method = method_i,
          prediction_type = prediction_type,
          alpha = alpha,
          n_grid = n_grid,
          x_cols = x_cols,
          seed = seed + 1000L * ii + match(method_i, methods)
        ),
        error = function(e) {
          message(sprintf(
            "Failed: test_id=%s, method=%s: %s",
            as.character(test_id),
            method_i,
            conditionMessage(e)
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
output_file <- file.path(results_dir, "figureS11_gallstones_band_data.csv")

write.csv(
  figureS11_data,
  file = output_file,
  row.names = FALSE
)

cat("Saved Figure S.11 data to:\n", output_file, "\n")
print(head(figureS11_data))