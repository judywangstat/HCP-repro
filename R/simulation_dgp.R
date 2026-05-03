# ============================================================
# simulation_dgp.R
# Core data-generating functions for simulation studies
# ============================================================

#' Generate subject-specific random effects
#'
#' @description
#' Generate subject-level random effects from a multivariate normal distribution.
#'
#' @param n Number of subjects.
#' @param theta_mean Mean vector of the random effects.
#' @param theta_sigma Covariance matrix of the random effects.
#'
#' @return An `n` by `length(theta_mean)` matrix of random effects.
generate_theta <- function(n, theta_mean, theta_sigma) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required.")
  }
  MASS::mvrnorm(n = n, mu = theta_mean, Sigma = theta_sigma)
}

#' Generate covariates for one subject
#'
#' @description
#' Generate the covariate matrix for a single subject. The first covariate is a
#' normalized time variable, and the remaining covariates are generated from a
#' multivariate normal distribution.
#'
#' @param m Number of observations for the subject.
#' @param x_mean Mean vector for the covariates.
#' @param x_sigma Covariance matrix for the covariates.
#'
#' @return An `m` by `length(x_mean)` covariate matrix.
generate_subject_covariates <- function(m, x_mean, x_sigma) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required.")
  }
  
  d <- length(x_mean)
  
  # time variable: permutation of {1,...,m}/m
  x_time <- sample(seq_len(m)) / m
  
  if (d == 1L) {
    return(matrix(x_time, ncol = 1))
  }
  
  X_normal <- MASS::mvrnorm(
    n = m,
    mu = x_mean[-1],
    Sigma = x_sigma[-1, -1, drop = FALSE]
  )
  
  X <- cbind(x_time, X_normal)
  colnames(X) <- paste0("X", seq_len(ncol(X)))
  X
}

#' Generate simulation errors
#'
#' @description
#' Generate outcome errors under one of the simulation scenarios: homogeneous,
#' heterogeneous, asymmetric, bimodal, or fixed-bimodal.
#'
#' @param m Number of observations.
#' @param X Covariate matrix used by scenarios with covariate-dependent errors.
#' @param scenario Simulation scenario. Supported values are `"Homo"`,
#' `"Heter"`, `"Asym"`, `"Bimo"`, and `"Bimo-fix"`.
#'
#' @return A numeric vector of generated errors.
generate_errors <- function(m, X, scenario) {
  X <- as.matrix(X)
  
  if (scenario == "Homo") {
    eps <- rnorm(m, 0, 1)
    
  } else if (scenario == "Heter") {
    eps <- rnorm(m, 0, 1 + 3 * abs(X[, 2]))
    
  } else if (scenario == "Asym") {
    eps <- rgamma(m, shape = 0.1, rate = 0.1)
    
  } else if (scenario %in% c("Bimo", "Bimo-fix")) {
    eps <- ifelse(
      runif(m) < 0.5,
      rnorm(m, -8, 0.5),
      rnorm(m, 8, 0.5)
    )
    
  } else {
    stop("Unknown scenario.")
  }
  
  eps
}

#' Generate outcomes
#'
#' @description
#' Generate outcome values from the linear mixed simulation model
#' `Y = X^T theta + epsilon`.
#'
#' @param X Covariate matrix.
#' @param theta Subject-specific coefficient vector.
#' @param epsilon Error vector.
#'
#' @return A numeric vector of generated outcomes.
generate_outcomes <- function(X, theta, epsilon) {
  rowSums(X * matrix(theta, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)) + epsilon
}

#' Generate missingness indicators
#'
#' @description
#' Generate binary missingness indicators from a logistic missingness model.
#'
#' @param X Covariate matrix.
#' @param beta Coefficient vector for the logistic missingness model, including
#' the intercept.
#'
#' @return An integer vector of missingness indicators, where 1 indicates an
#' observed response and 0 indicates a missing response.
generate_missing_indicator <- function(X, beta) {
  X_tilde <- cbind(1, X)
  log_odds <- X_tilde %*% beta
  prob <- 1 / (1 + exp(-log_odds))
  rbinom(nrow(X), 1, prob)
}

#' Generate full clustered simulation dataset
#'
#' @description
#' Generate a clustered simulation dataset in long format. For each subject, the
#' function generates subject-specific random effects, covariates, errors,
#' outcomes, and missingness indicators under the specified simulation scenario.
#'
#' @param n Number of subjects.
#' @param m Number of observations per subject. Can be a scalar or a vector of
#' length `n`.
#' @param scenario Simulation scenario. Supported values are `"Homo"`,
#' `"Heter"`, `"Asym"`, `"Bimo"`, and `"Bimo-fix"`.
#' @param theta_mean Mean vector of the subject-specific random effects.
#' @param theta_sigma Covariance matrix of the subject-specific random effects.
#' @param x_mean Mean vector used to generate covariates.
#' @param x_sigma Covariance matrix used to generate covariates.
#' @param beta Coefficient vector for the logistic missingness model, including
#' the intercept.
#' @param seed Optional random seed.
#'
#' @return A data frame in long format with columns `id`, `time_index`, `Y`,
#' `delta`, and generated covariates.
generate_simulation_data <- function(
    n,
    m = 5,
    scenario = c("Homo", "Heter", "Asym", "Bimo", "Bimo-fix"),
    theta_mean = rep(2, 4),
    theta_sigma = diag(1, 4),
    x_mean = rep(0, 4),
    x_sigma = diag(1, 4),
    beta = c(3, 0, 2, 2, 2),
    seed = NULL
) {
  
  scenario <- match.arg(scenario)
  
  if (!is.null(seed)) set.seed(seed)
  
  # expand m if scalar
  if (length(m) == 1L) {
    m <- rep(m, n)
  }
  
  d <- length(theta_mean)
  
  # generate random effects
  theta <- generate_theta(n, theta_mean, theta_sigma)
  
  out_list <- vector("list", n)
  
  for (i in seq_len(n)) {
    
    # covariates
    X_i <- generate_subject_covariates(m[i], x_mean, x_sigma)
    
    # handle Bimo-fix
    theta_i <- theta[i, ]
    if (scenario == "Bimo-fix") {
      theta_i <- rep(theta_mean[1], d)
    }
    
    # errors
    eps_i <- generate_errors(m[i], X_i, scenario)
    
    # outcome
    Y_i <- generate_outcomes(X_i, theta_i, eps_i)
    
    # missingness
    delta_i <- generate_missing_indicator(X_i, beta)
    
    # build data.frame
    out_list[[i]] <- data.frame(
      id = i,
      time_index = seq_len(m[i]),
      Y = as.numeric(Y_i),
      delta = as.integer(delta_i),
      X_i,
      check.names = FALSE
    )
  }
  
  dat <- do.call(rbind, out_list)
  rownames(dat) <- NULL
  dat
}

#' Generate outcomes given fixed covariates
#'
#' @description
#' Generate outcomes for a fixed covariate matrix. This helper is used in
#' conditional-coverage experiments where the test covariates are fixed and only
#' the random effects and errors are regenerated.
#'
#' @param X Fixed covariate matrix.
#' @param scenario Simulation scenario. Supported values are `"Homo"`,
#' `"Heter"`, `"Asym"`, `"Bimo"`, and `"Bimo-fix"`.
#' @param theta_mean Mean vector of the subject-specific random effects.
#' @param theta_sigma Covariance matrix of the subject-specific random effects.
#'
#' @return A numeric vector of generated outcomes corresponding to the rows of
#' `X`.
generate_Y_given_X <- function(
    X,
    scenario,
    theta_mean = rep(2, 4),
    theta_sigma = diag(1, 4)
) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  
  theta <- generate_theta(n, theta_mean, theta_sigma)
  
  if (scenario == "Bimo-fix") {
    theta <- matrix(theta_mean[1], nrow = n, ncol = d)
  }
  
  eps <- generate_errors(n, X, scenario)
  as.numeric(generate_outcomes(X, theta, eps))
}