# Shared production runner. No statistical model/tuning choices are made here.
# All raw per-method records survive until explicit aggregation and serialization.
run_benchmark_settings <- function(
    repo_dir, num_simulations_run, num_cores, sample_sizes,
    scenarios, alpha, n_test, n_grid, theta_mean, theta_sigma, x_mean, x_sigma, beta, weight_caps,
    rng_kind = "L'Ecuyer-CMRG", lmem_n_sims = 500L, return_details = FALSE, workers_log = '',
    cluster_size_design = 'balanced', method_names = simulation_methods()) {
  if (length(num_simulations_run) != 1L || !is.finite(num_simulations_run) || num_simulations_run < 1 ||
    num_simulations_run != as.integer(num_simulations_run)) stop('Invalid replicate count.')
  if (length(num_cores) != 1L || !is.finite(num_cores) || num_cores < 1) stop('Invalid worker count.')
  if (length(weight_caps) != length(sample_sizes)) stop('One weight cap per sample size required.')
  if (n_test != 1L) stop('The production design holds out exactly one subject.')
  cl <- parallel::makePSOCKcluster(as.integer(num_cores), outfile = workers_log)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, function(repo_dir, rng_kind) {
    source(file.path(repo_dir, 'R/load_simulation_benchmarks.R'), local = .GlobalEnv)
    load_simulation_benchmarks(repo_dir)
    RNGkind(rng_kind, 'Inversion', 'Rejection');options(warn = 1);NULL
  }, repo_dir, rng_kind)
  # Table 1 and Tables S.1/S.3 use L'Ecuyer-CMRG; Figure S.10 uses MT.
  if (rng_kind == "L'Ecuyer-CMRG") parallel::clusterSetRNGStream(cl, iseed = 123)
  records <- list()
  for (ni in seq_along(sample_sizes)) for (s in scenarios) {
    n <- sample_sizes[ni]
    args <- list(
      n = n, scenario = s, alpha = alpha, n_test = n_test, n_grid = n_grid,
      theta_mean = theta_mean, theta_sigma = theta_sigma, x_mean = x_mean, x_sigma = x_sigma,
      beta = beta, weight_cap = weight_caps[ni], lmem_n_sims = lmem_n_sims, return_details = return_details, cluster_size_design = cluster_size_design, method_names = method_names
    )
    cat(sprintf('Running n=%d scenario=%s (%d replicates)\n', n, s, num_simulations_run));flush.console()
    z <- parallel::parLapplyLB(cl, seq_len(num_simulations_run), function(i, args) {
      args$seed <- i
      ans <- tryCatch(do.call(run_one_replicate_pointwise, args), error = function(e) {
        # Only data/setup failures reach this handler. Method failures are already
        # caught inside the replicate, without losing any successful peer values.
        out <- matrix(NA_real_, 5L, 2L * length(args$method_names), dimnames = list(NULL, as.vector(rbind(paste0(args$method_names, '_cov'), paste0(args$method_names, '_len')))))
        attr(out, 'method_diagnostics') <- setNames(lapply(args$method_names, function(m)
          list(
            method = m, success = FALSE, error = conditionMessage(e), error_stage = 'setup',
            fit_failed = FALSE, prediction_failed = FALSE, evaluation_failed = FALSE,
            events = list(list(stage = 'setup', kind = 'error', message = conditionMessage(e))),
            model = NULL, n_coverage = 0L, n_length = 0L, empty_regions = 0L
          )), args$method_names)
        out
      })
      list(n = args$n, scenario = args$scenario, seed = i, pointwise = ans)
    }, args = args)
    records <- c(records, z)
    cat(sprintf(
      'Completed n=%d scenario=%s; method failures=%d\n', n, s,
      sum(vapply(z, function(a) sum(!vapply(attr(a$pointwise, 'method_diagnostics'), `[[`, logical(1), 'success')), integer(1)))
    ))
    flush.console()
  }
  list(
    records = records, summary = summarize_benchmark_records(records),
    config = list(
      num_simulations_run = num_simulations_run, sample_sizes = sample_sizes, scenarios = scenarios,
      alpha = alpha, n_test = n_test, n_grid = n_grid, theta_mean = theta_mean, theta_sigma = theta_sigma,
      x_mean = x_mean, x_sigma = x_sigma, beta = beta, weight_caps = weight_caps, rng_kind = rng_kind,
      lmem_n_sims = lmem_n_sims, cluster_size_design = cluster_size_design, method_names = method_names
    ), session = sessionInfo()
  )
}

summarize_benchmark_records <- function(records) {
  groups <- unique(vapply(records, function(z) paste(z$n, z$scenario, sep = ':'), character(1)))
  final <- list()
  for (key in groups) {
    r <- Filter(function(z) paste(z$n, z$scenario, sep = ':') == key, records)
    mat <- do.call(rbind, lapply(r, function(z) colMeans(z$pointwise, na.rm = TRUE)))
    final[[length(final) + 1L]] <- data.frame(
      sample_size = r[[1]]$n,
      scenario = if (r[[1]]$scenario == 'Bimo-fix') 'Bimo-Fix' else r[[1]]$scenario,
      as.list(colMeans(mat, na.rm = TRUE)), check.names = FALSE
    )
  }
  do.call(rbind, final)
}

benchmark_audit_tables <- function(records) {
  attempts <- events <- oracle <- pointwise <- list()
  nullable <- function(x, default = NA) if (is.null(x) || length(x) == 0L) default else x
  for (z in records) {
    base <- data.frame(n = z$n, scenario = z$scenario, seed = z$seed)
    a <- z$pointwise;dd <- attr(a, 'method_diagnostics')
    pointwise[[length(pointwise) + 1L]] <- cbind(base, test_row = seq_len(nrow(a)), as.data.frame(a))
    od <- attr(a, 'oracle_diagnostics')
    if (!is.null(od)) oracle[[length(oracle) + 1L]] <- cbind(base, od)
    for (m in names(dd)) {
      d <- dd[[m]];model <- d$model
      row <- cbind(base,
        method = m, success = d$success, error_stage = nullable(d$error_stage, ''),
        fit_failure = nullable(d$fit_failed), prediction_failure = nullable(d$prediction_failed),
        evaluation_failure = d$evaluation_failed, singular = nullable(model$singular),
        convergence_warning = nullable(model$convergence_warning),
        convergence_messages = paste(nullable(model$convergence_messages, character()), collapse = ' | '),
        optimizer_code = paste(nullable(model$optimizer_code, character()), collapse = ' | '),
        n_coverage = d$n_coverage, n_length = d$n_length, empty_regions = d$empty_regions,
        warning_count = sum(vapply(d$events, function(e) e$kind == 'warning', logical(1))),
        residual_variance = nullable(model$residual_variance),
        coverage = mean(a[, paste0(m, '_cov')], na.rm = TRUE), length = mean(a[, paste0(m, '_len')], na.rm = TRUE)
      )
      for (v in c('var_intercept', paste0('var_X', 1:4))) row[[v]] <- nullable(model$variances[v])
      attempts[[length(attempts) + 1L]] <- row
      if (length(d$events)) {
        ev <- do.call(rbind, lapply(d$events, as.data.frame))
        count <- aggregate(rep(1L, nrow(ev)), by = ev, FUN = sum);names(count)[4] <- 'count'
        events[[length(events) + 1L]] <- cbind(base, method = m, count)
      }
    }
  }
  a <- do.call(rbind, attempts);summary <- list()
  groups <- unique(a[, c('n', 'scenario', 'method')])
  for (i in seq_len(nrow(groups))) {
    g <- groups[i, ];b <- a[a$n == g$n & a$scenario == g$scenario & a$method == g$method, ]
    ns <- sum(!is.na(b$singular));nc <- sum(!is.na(b$convergence_warning))
    summary[[i]] <- cbind(g,
      attempted = nrow(b), successful = sum(b$success), failed = sum(!b$success),
      failure_proportion = mean(!b$success), fit_failures = sum(b$fit_failure, na.rm = TRUE),
      prediction_failures = sum(b$prediction_failure, na.rm = TRUE),
      singular = sum(b$singular, na.rm = TRUE), fitted_with_singularity_check = ns,
      singular_proportion = if (ns) sum(b$singular, na.rm = TRUE) / ns else NA_real_,
      convergence_warnings = sum(b$convergence_warning, na.rm = TRUE),
      warning_replicates = sum(b$warning_count > 0), valid_coverage_replicates = sum(is.finite(b$coverage)),
      valid_length_replicates = sum(is.finite(b$length)), n_coverage = sum(b$n_coverage), n_length = sum(b$n_length),
      coverage = mean(b$coverage, na.rm = TRUE), length = mean(b$length, na.rm = TRUE),
      coverage_mcse = sd(b$coverage, na.rm = TRUE) / sqrt(sum(is.finite(b$coverage)))
    )
  }
  list(
    attempts = a, events = if (length(events)) do.call(rbind, events) else data.frame(),
    oracle = if (length(oracle)) do.call(rbind, oracle) else data.frame(),
    pointwise = do.call(rbind, pointwise), summary = do.call(rbind, summary)
  )
}

save_benchmark_audit <- function(run, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  tables <- benchmark_audit_tables(run$records)
  for (n in names(tables)) write.csv(tables[[n]], file.path(directory, paste0(n, '.csv')), row.names = FALSE)
  saveRDS(run, file.path(directory, 'run.rds'))
  writeLines(capture.output(str(list(config = run$config, session = run$session))), file.path(directory, 'metadata.txt'))
  invisible(tables)
}
