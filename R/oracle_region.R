#' Construct oracle prediction intervals from the true data-generating mechanism
#'
#' @description
#' This function computes oracle prediction intervals for new covariate values
#' under a known simulation data-generating mechanism (DGP). It is intended only
#' for simulation benchmarks, not for real-data analysis.
#'
#' Depending on the scenario, the interval is obtained either from a closed-form
#' distribution or from Monte Carlo approximation to the true conditional
#' distribution of \eqn{Y \mid X = x}.
#'
#' @param x_test New covariate value(s). This can be a numeric vector
#'   (treated as one test point), or a matrix/data.frame with one row per test point
#'   and \code{d} columns.
#' @param scenario Scenario name. Must be one of \code{"Homo"}, \code{"Heter"},
#'   \code{"Asym"}, \code{"Bimo"}, or \code{"Bimo-fix"}.
#' @param alpha Miscoverage level in \eqn{(0,1)}.
#' @param B_mc Number of Monte Carlo samples used when the oracle interval is not
#'   available in closed form.
#' @param theta_mean Mean of the random coefficient distribution.
#' @param theta_sd Standard deviation of the random coefficient distribution.
#' @param seed Optional random seed.
#'
#' @return A numeric matrix with two columns:
#' \describe{
#'   \item{\code{lo}}{Lower oracle prediction bound.}
#'   \item{\code{hi}}{Upper oracle prediction bound.}
#' }
oracle_interval_from_dgp <- function(
    x_test,
    scenario = c("Homo", "Heter", "Asym", "Bimo", "Bimo-fix"),
    alpha = 0.1,
    B_mc = 10000L,
    theta_mean = 2,
    theta_sd = 1,
    seed = NULL
) {
  
  # ---------------------------------------------------------------------------
  # Step 0: Validate and standardize inputs
  # ---------------------------------------------------------------------------
  scenario <- match.arg(scenario)
  
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number in (0,1).")
  }
  
  B_mc <- as.integer(B_mc)
  if (!is.finite(B_mc) || B_mc < 1L) {
    stop("B_mc must be a positive integer.")
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (is.vector(x_test) && !is.list(x_test)) {
    x_test <- matrix(as.numeric(x_test), nrow = 1)
  } else if (is.data.frame(x_test)) {
    x_test <- as.matrix(x_test)
  } else {
    x_test <- as.matrix(x_test)
  }
  
  storage.mode(x_test) <- "numeric"
  
  if (anyNA(x_test) || any(!is.finite(x_test))) {
    stop("x_test must contain only finite numeric values.")
  }
  
  d <- ncol(x_test)
  K <- nrow(x_test)
  
  # ---------------------------------------------------------------------------
  # Step 1: Compute the mean and random-effect variance contribution
  # ---------------------------------------------------------------------------
  mu <- as.numeric(x_test %*% rep(theta_mean, d))
  var_theta <- (theta_sd^2) * rowSums(x_test^2)
  
  # Helper for Gaussian intervals
  normal_interval <- function(mu, var, alpha) {
    z_lo <- stats::qnorm(alpha / 2)
    z_hi <- stats::qnorm(1 - alpha / 2)
    sd <- sqrt(var)
    cbind(lo = mu + sd * z_lo, hi = mu + sd * z_hi)
  }
  
  # ---------------------------------------------------------------------------
  # Step 2: Closed-form oracle intervals when available
  # ---------------------------------------------------------------------------
  if (scenario == "Homo") {
    out <- normal_interval(mu = mu, var = var_theta + 1, alpha = alpha)
    return(out)
  }
  
  if (scenario == "Heter") {
    if (d < 2L) {
      stop("Scenario 'Heter' requires at least two covariates.")
    }
    var_eps <- (1 + 3 * abs(x_test[, 2]))^2
    out <- normal_interval(mu = mu, var = var_theta + var_eps, alpha = alpha)
    return(out)
  }
  
  # ---------------------------------------------------------------------------
  # Step 3: Monte Carlo oracle intervals for non-Gaussian scenarios
  # ---------------------------------------------------------------------------
  out <- matrix(NA_real_, nrow = K, ncol = 2)
  colnames(out) <- c("lo", "hi")
  
  if (scenario == "Bimo-fix") {
    # No random theta in the Bimo-fix DGP
    comp <- sample.int(2, size = B_mc, replace = TRUE)
    means <- ifelse(comp == 1, -8, 8)
    eps <- stats::rnorm(B_mc, mean = means, sd = 0.5)
    
    for (k in seq_len(K)) {
      out[k, ] <- stats::quantile(
        mu[k] + eps,
        probs = c(alpha / 2, 1 - alpha / 2),
        names = FALSE
      )
    }
    return(out)
  }
  
  # For Asym and Bimo, theta remains random
  Theta <- matrix(
    stats::rnorm(B_mc * d, mean = theta_mean, sd = theta_sd),
    nrow = B_mc, ncol = d
  )
  
  MU <- Theta %*% t(x_test)   # B_mc x K
  
  if (scenario == "Asym") {
    EPS <- matrix(
      stats::rgamma(B_mc * K, shape = 0.1, rate = 0.1),
      nrow = B_mc, ncol = K
    )
  } else if (scenario == "Bimo") {
    comp <- matrix(
      sample.int(2, size = B_mc * K, replace = TRUE),
      nrow = B_mc, ncol = K
    )
    means <- ifelse(comp == 1, -8, 8)
    EPS <- matrix(
      stats::rnorm(B_mc * K, mean = as.vector(means), sd = 0.5),
      nrow = B_mc, ncol = K
    )
  } else {
    stop("Unsupported scenario.")
  }
  
  Y_mc <- MU + EPS
  
  for (k in seq_len(K)) {
    out[k, ] <- stats::quantile(
      Y_mc[, k],
      probs = c(alpha / 2, 1 - alpha / 2),
      names = FALSE
    )
  }
  
  out
}

#' Construct an oracle prediction region under a known simulation DGP
#'
#' @description
#' This function constructs oracle prediction regions for new covariate values
#' under a known simulation data-generating mechanism. It is intended only for
#' simulation studies and benchmark comparisons.
#'
#' The oracle interval is computed from the true conditional distribution of
#' \eqn{Y \mid X = x}, either analytically or via Monte Carlo approximation,
#' depending on the scenario.
#'
#' @param x_test New covariate value(s). This can be a numeric vector
#'   (treated as one test point), or a matrix/data.frame with one row per test point.
#' @param y_grid Numeric vector of candidate response values. This is used only
#'   to convert the oracle interval into a grid-based region for consistency with
#'   other methods.
#' @param scenario Scenario name. Must be one of \code{"Homo"}, \code{"Heter"},
#'   \code{"Asym"}, \code{"Bimo"}, or \code{"Bimo-fix"}.
#' @param alpha Miscoverage level in \eqn{(0,1)}.
#' @param B_mc Number of Monte Carlo samples used when the oracle interval is not
#'   available in closed form.
#' @param theta_mean Mean of the random coefficient distribution.
#' @param theta_sd Standard deviation of the random coefficient distribution.
#' @param seed Optional random seed.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{region}}{A list of oracle prediction regions, one for each row of
#'   \code{x_test}, represented as subsets of \code{y_grid}.}
#'   \item{\code{lo_hi}}{A matrix with columns \code{"lo"} and \code{"hi"}
#'   giving the oracle interval endpoints.}
#'   \item{\code{p_final}}{\code{NULL}, since oracle intervals do not arise from
#'   conformal \eqn{p}-values.}
#'   \item{\code{y_grid}}{The candidate grid used to represent the interval as a region.}
#' }
#'
#' @export
oracle_region <- function(
    x_test,
    y_grid,
    scenario = c("Homo", "Heter", "Asym", "Bimo", "Bimo-fix"),
    alpha = 0.1,
    B_mc = 10000L,
    theta_mean = 2,
    theta_sd = 1,
    seed = NULL
) {
  
  # ---------------------------------------------------------------------------
  # Step 0: Validate inputs
  # ---------------------------------------------------------------------------
  scenario <- match.arg(scenario)
  
  y_grid <- as.numeric(y_grid)
  if (anyNA(y_grid) || any(!is.finite(y_grid))) {
    stop("y_grid must contain only finite numeric values.")
  }
  y_grid <- sort(unique(y_grid))
  if (length(y_grid) < 2L) {
    stop("y_grid must contain at least two distinct values.")
  }
  
  # ---------------------------------------------------------------------------
  # Step 1: Compute oracle interval
  # ---------------------------------------------------------------------------
  lo_hi <- oracle_interval_from_dgp(
    x_test = x_test,
    scenario = scenario,
    alpha = alpha,
    B_mc = B_mc,
    theta_mean = theta_mean,
    theta_sd = theta_sd,
    seed = seed
  )
  
  K <- nrow(lo_hi)
  
  # ---------------------------------------------------------------------------
  # Step 2: Convert interval into grid-based region
  # ---------------------------------------------------------------------------
  regions <- vector("list", K)
  for (k in seq_len(K)) {
    regions[[k]] <- y_grid[y_grid >= lo_hi[k, "lo"] & y_grid <= lo_hi[k, "hi"]]
  }
  
  # ---------------------------------------------------------------------------
  # Return final output
  # ---------------------------------------------------------------------------
  list(
    region = regions,
    lo_hi = lo_hi,
    p_final = NULL,
    y_grid = y_grid
  )
}