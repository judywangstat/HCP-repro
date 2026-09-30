# All production simulation helpers are sourced identically on master and workers.
load_simulation_benchmarks <- function(repo_dir, envir = .GlobalEnv) {
  files <- c(
    'simulation_dgp.R', 'imbalanced_covariates.R', 'tuning_helpers.R', 'evaluation_helpers.R',
    'fit_propensity_model.R', 'fit_cond_density_qp.R', 'hcp_region.R', 'dwr_region.R',
    'lc_region.R', 'hcp_finite_region.R', 'lmem_region.R', 'oracle_region.R',
    'simulation_method_helpers.R', 'run_one_replicate_pointwise.R', 'run_one_replicate.R',
    'run_benchmark_settings.R'
  )
  for (f in files) sys.source(file.path(repo_dir, 'R', f), envir = envir)
  invisible(files)
}
