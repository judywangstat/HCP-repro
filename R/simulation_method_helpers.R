# Non-statistical orchestration: isolate method-level failures and retain evidence.
simulation_methods <- function() c('HCP', 'DWR', 'LC', 'LMEM', 'Oracle')
simulation_result_columns <- function() as.vector(rbind(paste0(simulation_methods(), '_cov'), paste0(simulation_methods(), '_len')))

run_simulation_method <- function(method, build, evaluate, y_true, keep_prediction = FALSE) {
  state <- new.env(parent = emptyenv());state$stage <- 'method';state$fit <- NULL
  events <- list();error <- NULL;error_stage <- NULL;res <- NULL
  record <- function(kind, c) events[[length(events) + 1L]] <<- list(stage = state$stage, kind = kind, message = conditionMessage(c))
  ans <- tryCatch(
    withCallingHandlers(
      {
        res <- build(state);state$stage <- 'evaluation';value <- evaluate(res, y_true)
        if (!identical(dim(value), c(length(y_true), 2L))) stop('Invalid method evaluation shape.')
        value
      },
      warning = function(w) record('warning', w),
      message = function(m) record('message', m)
    ),
    error = function(e) {
      error <<- conditionMessage(e);error_stage <<- state$stage;record('error', e);NULL
    }
  )
  # Warnings/messages are recorded AND emitted. No silent suppression.
  if (is.null(ans)) ans <- matrix(NA_real_, length(y_true), 2, dimnames = list(NULL, c('covered', 'length')))
  model <- lmem_fit_diagnostics(state$fit)
  if (!is.null(model)) {
    fitting_warnings <- vapply(Filter(function(e) e$stage == 'fit' && e$kind == 'warning', events), `[[`, character(1), 'message')
    model$convergence_warning <- model$convergence_warning || any(grepl('converg|gradient|Hessian|eigenvalue|identif', fitting_warnings, ignore.case = TRUE))
  }
  diagnostics <- list(
    method = method, success = is.null(error), error = error, error_stage = error_stage,
    fit_failed = if (method %in% c('LMEM')) identical(error_stage, 'fit') else NA,
    prediction_failed = if (method %in% c('LMEM')) identical(error_stage, 'prediction') else NA,
    evaluation_failed = identical(error_stage, 'evaluation'), events = events, model = model,
    n_coverage = sum(is.finite(ans[, 1])), n_length = sum(is.finite(ans[, 2])),
    empty_regions = if (!is.null(res$region)) sum(lengths(res$region) == 0L) else 0L
  )
  out <- list(evaluation = ans, diagnostics = diagnostics)
  if (method == 'Oracle' && !is.null(res)) out$oracle_diagnostics <- data.frame(
    test_row = seq_along(res$oracles), n_components = res$n_components,
    probability = res$conditional_probability, total_length = res$total_length,
    boundary_abs_error = vapply(res$diagnostics, `[[`, numeric(1), 'boundary_abs_error')
  )
  if (keep_prediction) out$prediction <- res
  out
}
