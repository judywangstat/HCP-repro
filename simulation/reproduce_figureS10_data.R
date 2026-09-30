# ============================================================
# Reproduce Figure S.10 data
# 20% missing responses, n=500; pointwise data for five methods
# Input: simulation_dgp.R, fixed replicate seeds 1--1000.
# Output: results/figureS10_{HCP,DWR,LC,LMEM,Oracle}_final_mat_n500.csv
# Set repo_dir below or HCP_REPO_DIR; run from the repository root.
# Sourcing defines functions. Rscript simulation/reproduce_figureS10_data.R --run runs the experiment.
# ============================================================
repo_dir <- normalizePath(Sys.getenv('HCP_REPO_DIR', '.'), mustWork = TRUE)
source(file.path(repo_dir, 'R/load_simulation_benchmarks.R'))
load_simulation_benchmarks(repo_dir)
results_dir <- file.path(repo_dir, 'results')

run_figureS10_data <- function(
    num_simulations_run = 1000, num_cores = 7L, n = 500,
    scenarios = c('Homo', 'Heter', 'Asym', 'Bimo', 'Bimo-fix'), alpha = .1, n_test = 1, n_grid = 200,
    theta_mean = rep(2, 4), theta_sigma = diag(1, 4), x_mean = rep(0, 4), x_sigma = diag(1, 4),
    beta = c(3, 0, 2, 2, 2), weight_cap = 30, lmem_n_sims = 500L,
    output_dir = results_dir, audit_dir = file.path(output_dir, 'figureS10_audit'),
    return_details = FALSE, workers_log = '', method_names = simulation_methods()) {
  # This experiment uses Mersenne-Twister for data generation.
  # All five methods use the same generated dataset and held-out subject in
  # each replicate. The default writes all five method files together.
  run <- run_benchmark_settings(repo_dir, num_simulations_run, num_cores, n, scenarios,
    alpha, n_test, n_grid, theta_mean, theta_sigma, x_mean, x_sigma, beta, weight_cap,
    rng_kind = 'Mersenne-Twister', lmem_n_sims = lmem_n_sims,
    return_details = return_details, workers_log = workers_log, method_names = method_names
  )
  by_scenario <- setNames(lapply(scenarios, function(s) {
    r <- Filter(function(z) z$scenario == s, run$records)
    r <- r[order(vapply(r, `[[`, integer(1), 'seed'))]
    do.call(rbind, lapply(r, `[[`, 'pointwise'))
  }), scenarios)
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    for (m in method_names) {
      mat <- do.call(cbind, lapply(scenarios, function(s) {
        a <- by_scenario[[s]][, paste0(m, c('_cov', '_len')), drop = FALSE]
        colnames(a) <- make.names(paste0(s, c('_cov', '_len')));a
      }))
      write.csv(mat, file.path(output_dir, paste0('figureS10_', m, '_final_mat_n', n, '.csv')), row.names = FALSE)
    }
  }
  if (!is.null(audit_dir)) save_benchmark_audit(run, audit_dir)
  attr(by_scenario, 'run') <- run
  invisible(by_scenario)
}

if (sys.nframe() == 0L && '--run' %in% commandArgs(TRUE)) {
  run_figureS10_data(
    num_simulations_run = 1000, num_cores = 7, n = 500, weight_cap = 30,
    output_dir = results_dir
  )
}
