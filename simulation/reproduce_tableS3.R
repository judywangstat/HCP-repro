# ============================================================
# Reproduce Table S.3
# 50% missing responses, balanced clusters of size 5
# Input: simulation_dgp.R, fixed replicate seeds 1--1000.
# Output: results/tableS3_final_results_missing50.csv
# Set repo_dir below or HCP_REPO_DIR; run from the repository root.
# Sourcing defines functions. Rscript simulation/reproduce_tableS3.R --run runs the experiment.
# ============================================================
repo_dir <- normalizePath(Sys.getenv('HCP_REPO_DIR', '.'), mustWork = TRUE)
source(file.path(repo_dir, 'R/load_simulation_benchmarks.R'))
load_simulation_benchmarks(repo_dir)
results_dir <- file.path(repo_dir, 'results')

run_tableS3_parallel <- function(
    num_simulations_run = 1000, num_cores = 2L,
    sample_sizes = c(100, 300, 500), scenarios = c('Homo', 'Heter', 'Asym', 'Bimo', 'Bimo-fix'),
    alpha = .1, n_test = 1, n_grid = 200, theta_mean = rep(2, 4), theta_sigma = diag(1, 4),
    x_mean = rep(0, 4), x_sigma = diag(1, 4), beta = c(0, 0, 2, 2, 2),
    weight_cap_rule = function(n) pmin(28, 2.2 * (n / 100)^1.72), output_file = NULL,
    run_sanity_check = FALSE, audit_dir = NULL, return_details = FALSE, workers_log = '') {
  caps <- vapply(sample_sizes, weight_cap_rule, numeric(1))
  if (isTRUE(run_sanity_check)) print(run_one_replicate(sample_sizes[1], scenarios[1],
    alpha = alpha, n_test = n_test, n_grid = n_grid, theta_mean = theta_mean, theta_sigma = theta_sigma,
    x_mean = x_mean, x_sigma = x_sigma, beta = beta, weight_cap = caps[1], seed = 1
  ))
  run <- run_benchmark_settings(repo_dir, num_simulations_run, num_cores, sample_sizes, scenarios,
    alpha, n_test, n_grid, theta_mean, theta_sigma, x_mean, x_sigma, beta, caps,
    rng_kind = "L'Ecuyer-CMRG", return_details = return_details, workers_log = workers_log
  )
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    write.csv(run$summary, output_file, row.names = FALSE)
    if (is.null(audit_dir)) audit_dir <- paste0(tools::file_path_sans_ext(output_file), '_audit')
  }
  if (!is.null(audit_dir)) save_benchmark_audit(run, audit_dir)
  final <- run$summary;attr(final, 'run') <- run
  invisible(final)
}

if (sys.nframe() == 0L && '--run' %in% commandArgs(TRUE)) {
  run_tableS3_parallel(
    num_simulations_run = 1000, num_cores = 7,
    output_file = file.path(results_dir, 'tableS3_final_results_missing50.csv')
  )
}
