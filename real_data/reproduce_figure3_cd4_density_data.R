# ============================================================
# Reproduce Figure 3 Data
# CD4 conditional density example
# HCP-repro
#
# This script generates leave-one-subject-out conditional density
# estimates for the first 100 CD4 subjects. The saved output is
# used by real_data/plot_figure3_cd4_density.R to reproduce
# Figure 3 in the paper.
#
# Input:
#   data/CD4_data.txt
#
# Output:
#   results/figure3_cd4_density_first100_data.csv
#
# Note:
#   Update repo_dir below to the local path of this repository
#   before running the script.
# ============================================================

suppressPackageStartupMessages({
  library(grf)
  library(quantreg)
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
source(file.path(r_dir, "fit_cond_density_qp.R"), chdir = TRUE)
source(file.path(r_dir, "tuning_helpers.R"), chdir = TRUE)

# ------------------------------------------------------------
# Main settings
# ------------------------------------------------------------
missing_rate <- "50"
seed <- 123

n_subjects_use <- 100
selected_time_index <- c(6, 7)

n_grid <- 200
b_true <- 5
dens_taus <- (1:(2^b_true - 1)) / (2^b_true)

# Density-tail decay factor for the Figure 3 conditional-density estimates.
density_decay_factor <- 0.3

# ------------------------------------------------------------
# Read CD4 data and generate artificial missingness
# ------------------------------------------------------------
cd4_dat <- read_cd4_data(file.path(data_dir, "CD4_data.txt"))
cd4_dat <- add_cd4_missingness(
  dat = cd4_dat,
  missing_rate = missing_rate,
  seed = seed
)

x_cols <- get_cd4_x_cols()
all_ids <- sort(unique(cd4_dat$id))

n_use <- min(n_subjects_use, length(all_ids))
subject_indices <- seq_len(n_use)

# ------------------------------------------------------------
# Run leave-one-subject-out density estimation
# ------------------------------------------------------------
density_list <- vector("list", length(subject_indices))

for (jj in seq_along(subject_indices)) {
  subject_index <- subject_indices[jj]
  test_id <- all_ids[subject_index]
  
  cat(sprintf(
    "Running subject_index = %d, id = %s (%d/%d)\n",
    subject_index, as.character(test_id), jj, length(subject_indices)
  ))
  flush.console()
  
  dat_test <- cd4_dat[cd4_dat$id == test_id, , drop = FALSE]
  dat_sample <- cd4_dat[cd4_dat$id != test_id, , drop = FALSE]
  
  y_grid <- make_cd4_y_grid(
    dat_sample = dat_sample,
    y_col = "Y",
    n_grid = n_grid
  )
  
  set.seed(seed + subject_index)
  
  dens_fit <- fit_cond_density_qp(
    dat = dat_sample,
    y_col = "Y",
    delta_col = "delta",
    x_cols = x_cols,
    taus = dens_taus,
    h = NULL,
    method = "grf",
    enforce_monotone = FALSE,
    tail_decay = TRUE,
    num_extra_points = 20L,
    decay_factor = density_decay_factor,
    seed = seed + subject_index
  )
  
  x_test <- dat_test[, x_cols, drop = FALSE]
  
  dens_mat <- dens_fit$predict_density_grid(
    x_new = x_test,
    y_grid = y_grid
  )
  
  if (!is.matrix(dens_mat) || nrow(dens_mat) != nrow(dat_test)) {
    stop("Unexpected density prediction output dimension for id = ", test_id)
  }
  
  density_list[[jj]] <- data.frame(
    subject_index = rep(subject_index, nrow(dat_test) * length(y_grid)),
    id = rep(test_id, nrow(dat_test) * length(y_grid)),
    time_index = rep(seq_len(nrow(dat_test)), each = length(y_grid)),
    time = rep(dat_test$time, each = length(y_grid)),
    true = rep(dat_test$Y, each = length(y_grid)),
    y_grid = rep(y_grid, times = nrow(dat_test)),
    density = as.vector(t(dens_mat)),
    selected_for_plot = rep(seq_len(nrow(dat_test)) %in% selected_time_index, each = length(y_grid)),
    missing_rate = rep(missing_rate, nrow(dat_test) * length(y_grid)),
    stringsAsFactors = FALSE
  )
}

density_data <- do.call(rbind, density_list)

# ------------------------------------------------------------
# Save output
# ------------------------------------------------------------
output_file <- file.path(results_dir, "figure3_cd4_density_first100_data.csv")

write.csv(
  density_data,
  file = output_file,
  row.names = FALSE
)

cat("Saved Figure 3 density data to:\n", output_file, "\n")
cat("Number of subjects saved:", n_use, "\n")
cat("Missing rate:", missing_rate, "\n")
cat("Selected time indices:", paste(selected_time_index, collapse = ", "), "\n")
print(head(density_data))