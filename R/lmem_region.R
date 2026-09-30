# ============================================================
# lmem_region.R
# Linear mixed-effects prediction for a new subject
#
# Simulation: fixed intercept and four fixed slopes, with four independent
# random slopes and no random intercept. CD4 and gallstones: independent
# random intercept and time slope, with the dataset-specific fixed effects.
# Fits use observed responses and REML. Prediction draws include uncertainty
# in beta, new-subject random effects, and independent residual errors.
# ============================================================

#' Collect every random-effect covariance block returned by lmer
lmem_random_covariance <- function(fit, terms) {
  G <- matrix(0, length(terms), length(terms), dimnames = list(terms, terms))
  seen <- character()
  for (block in lme4::VarCorr(fit)) {
    block <- as.matrix(block); nm <- rownames(block)
    if (any(nm %in% seen) || !all(nm %in% terms)) stop("Unexpected random-effect block.")
    G[nm, nm] <- block; seen <- c(seen, nm)
  }
  if (!setequal(seen, terms)) stop("Missing random-effect variance block.")
  if (any(G[row(G) != col(G)] != 0)) stop("The specified random effects must be independent.")
  G
}

#' Square root of a positive semidefinite covariance matrix
lmem_covariance_root <- function(G, cholesky = FALSE) {
  G <- (G + t(G)) / 2
  e <- eigen(G, symmetric = TRUE)
  tol <- 1e-10 * max(1, max(abs(e$values)))
  if (any(e$values < -tol)) stop("Covariance is not positive semidefinite.")
  if (any(e$values < 0)) warning("Tiny negative covariance eigenvalues clipped to zero.")
  if (cholesky && min(e$values) > tol) return(t(chol(G)))
  sweep(e$vectors, 2L, sqrt(pmax(e$values, 0)), "*")
}

#' Fit an LMEM and construct pointwise prediction intervals for a new subject
#'
#' @param dat Training subjects, with id, Y, delta and the required covariates.
#' @param x_test Covariate matrix/data frame for ONE new subject. Gallstones
#'   additionally requires the original time column, in months.
#' @param setting One of simulation, cd4 or gallstones.
#' @param alpha Per-observation miscoverage level. Real-data simultaneous
#'   prediction uses alpha divided by the held-out subject's observation count.
#' @param n_sims Predictive draws: 500 for simulations, 10,000 for real data.
#' @param seed Replicate or held-out-subject seed. Prediction uses the next
#'   L'Ecuyer stream and restores the caller's RNG state.
#' @param return_draws Whether to retain the Monte Carlo components.
#' @param state Optional environment for recording fit/prediction diagnostics.
#' @return lo_hi, predictive_mean, predictive_covariance and fitting diagnostics.
#'   Variance components are plug-in estimates. This is not a full bootstrap
#'   over random-effect or residual variance parameters.
lmem_region <- function(dat, x_test, setting = c("simulation", "cd4", "gallstones"),
                        alpha = .1, n_sims = NULL, seed = NULL,
                        return_draws = FALSE, state = NULL) {
  setting <- match.arg(setting)
  simulation <- setting == "simulation"
  if (is.null(n_sims)) n_sims <- if (simulation) 500L else 10000L
  if (length(n_sims) != 1L || !is.finite(n_sims) || n_sims < 1L || n_sims != as.integer(n_sims)) stop("Invalid n_sims.")
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("Invalid alpha.")
  test <- as.data.frame(x_test)
  if (simulation) {
    if (ncol(test) != 4L) stop("Simulation requires four test covariates.")
    names(test) <- paste0("X", 1:4)
    formula <- Y ~ X1 + X2 + X3 + X4 + (0 + X1 + X2 + X3 + X4 || id)
    terms <- paste0("X", 1:4)
  } else if (setting == "cd4") {
    formula <- Y ~ time + age + smoke + drug + partners + cesd + (1 + time || id)
    terms <- c("(Intercept)", "time")
  } else {
    if (!"time" %in% names(test)) stop("Gallstones prediction requires time in months.")
    dat$t <- (dat$time - 6) / 18; test$t <- (test$time - 6) / 18
    formula <- Y ~ T1 + T2 + T3 + Treat + (1 + t || id)
    terms <- c("(Intercept)", "t")
  }
  observed <- dat[dat$delta == 1 & !is.na(dat$Y), , drop = FALSE]
  if (nrow(observed) < 10L || nrow(test) < 1L) stop("Insufficient observed training or test rows.")
  observed$id <- factor(observed$id)
  if (is.null(state)) state <- new.env(parent = emptyenv())
  state$stage <- "fit"
  fit <- lme4::lmer(formula, data = observed, REML = TRUE)
  state$fit <- fit; state$stage <- "prediction"
  beta <- lme4::fixef(fit)
  Vbeta <- as.matrix(stats::vcov(fit))[names(beta), names(beta), drop = FALSE]
  X <- stats::model.matrix(stats::delete.response(stats::terms(lme4::nobars(formula))), test)
  if (!identical(colnames(X), names(beta))) stop("Fixed-effect design mismatch.")
  G <- lmem_random_covariance(fit, terms)
  Z <- if (simulation) as.matrix(test[, terms, drop = FALSE]) else cbind(1, test[[terms[2]]])
  colnames(Z) <- terms
  sigma2 <- stats::sigma(fit)^2
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  old_kind <- RNGkind(); had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(seed)
  assign(".Random.seed", parallel::nextRNGStream(get(".Random.seed", .GlobalEnv)), .GlobalEnv)
  # Retain each experiment's specified normal-draw order and matrix factorization.
  if (simulation) {
    beta_z <- matrix(rnorm(n_sims * length(beta)), n_sims, length(beta))
    b_z <- matrix(rnorm(n_sims * 4L), n_sims, 4L)
    eps_z <- matrix(rnorm(n_sims * nrow(test)), n_sims, nrow(test))
    beta_draws <- sweep(beta_z %*% t(lmem_covariance_root(Vbeta)), 2L, beta, "+")
    b_draws <- b_z %*% t(lmem_covariance_root(G))
    residual_draws <- sqrt(sigma2) * t(eps_z)
  } else {
    beta_draws <- sweep(matrix(rnorm(n_sims * length(beta)), n_sims) %*% chol(Vbeta), 2L, beta, "+")
    b1 <- rnorm(n_sims)
    residual_draws <- matrix(rnorm(nrow(X) * n_sims, sd = sqrt(sigma2)), nrow(X), n_sims)
    b2 <- rnorm(n_sims)
    b_draws <- cbind(b1, b2) %*% t(lmem_covariance_root(G, cholesky = TRUE))
  }
  # Each draw shares one beta and one random-effect vector across the subject.
  fixed_draws <- X %*% t(beta_draws)
  random_draws <- Z %*% t(b_draws)
  y_draws <- fixed_draws + random_draws + residual_draws
  bounds <- t(apply(y_draws, 1L, stats::quantile, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE, type = 7))
  colnames(bounds) <- c("lo", "hi")
  out <- list(lo_hi = bounds, p_final = NULL,
              predictive_mean = as.vector(X %*% beta),
              predictive_covariance = X %*% Vbeta %*% t(X) + Z %*% G %*% t(Z) + diag(sigma2, nrow(test)),
              diagnostics = lmem_fit_diagnostics(fit))
  if (return_draws) out$draws <- list(y = y_draws, fixed = fixed_draws, random = random_draws,
    residual = residual_draws, b = b_draws, beta = beta_draws, X = X, Z = Z, G = G, Vbeta = Vbeta, sigma2 = sigma2)
  out
}

#' Record singularity, convergence and variance estimates without discarding fits
lmem_fit_diagnostics <- function(fit) {
  if (is.null(fit)) return(NULL)
  V <- as.data.frame(lme4::VarCorr(fit))
  variances <- setNames(rep(NA_real_, 5L), c("var_intercept", paste0("var_X", 1:4)))
  for (term in c("(Intercept)", paste0("X", 1:4))) {
    found <- V$vcov[V$grp != "Residual" & !is.na(V$var1) & V$var1 == term & is.na(V$var2)]
    if (length(found)) variances[if (term == "(Intercept)") "var_intercept" else paste0("var_", term)] <- sum(found)
  }
  messages <- unlist(fit@optinfo$conv$lme4$messages)
  non_singular <- messages[!grepl("boundary.*singular", messages, ignore.case = TRUE)]
  optimizer_code <- as.integer(fit@optinfo$conv$opt)
  if (!length(optimizer_code)) optimizer_code <- 0L
  list(singular = lme4::isSingular(fit, tol = 1e-4), convergence_messages = messages,
       non_singular_convergence_messages = non_singular, optimizer_code = optimizer_code,
       convergence_warning = length(non_singular) > 0L || any(optimizer_code != 0L),
       variances = variances, residual_variance = stats::sigma(fit)^2, n_observed = nobs(fit))
}
