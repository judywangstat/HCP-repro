#' Run one simulation replicate for the benchmark comparison
#'
#' @description
#' This function runs one full simulation replicate under a specified sample size
#' and scenario. It generates one clustered dataset, holds out one new subject
#' as the test subject, applies all benchmark methods to the remaining sample,
#' evaluates coverage and prediction-region length, and returns one row of
#' replicate-level summary results.
#'
#' The function is intended for simulation studies such as the reproduction of
#' Table 1. It assumes that the required method functions and helper functions
#' have already been sourced into the current R session.
#'
#' @param n Number of subjects in the generated dataset.
#' @param scenario Scenario name. Must be one of \code{"Homo"},
#'   \code{"Heter"}, \code{"Asym"}, \code{"Bimo"}, or \code{"Bimo-fix"}.
#' @param alpha Miscoverage level in \eqn{(0,1)}.
#' @param n_test Number of held-out test subjects. The current implementation
#'   assumes \code{n_test = 1}.
#' @param n_grid Number of grid points used to construct the candidate response grid.
#' @param theta_mean Mean vector of the random-effects distribution.
#' @param theta_sigma Covariance matrix of the random-effects distribution.
#' @param x_mean Mean vector of the covariate distribution.
#' @param x_sigma Covariance matrix of the covariate distribution.
#' @param beta Logistic regression coefficients for the missingness mechanism.
#' @param weight_cap Optional truncation level for inverse-propensity weights.
#'   If \code{NULL}, the default rule \code{get_max_wt(n)} is used.
#' @param seed Optional random seed for reproducibility.
#'
#' @return A named numeric vector containing:
#' \describe{
#'   \item{\code{HCP_cov}}{Coverage of HCP in this replicate.}
#'   \item{\code{HCP_len}}{Prediction-region length of HCP in this replicate.}
#'   \item{\code{DWR_cov}}{Coverage of DWR in this replicate.}
#'   \item{\code{DWR_len}}{Prediction-region length of DWR in this replicate.}
#'   \item{\code{LC_cov}}{Coverage of LC in this replicate.}
#'   \item{\code{LC_len}}{Prediction-region length of LC in this replicate.}
#'   \item{\code{LMEM_cov}}{Coverage of LMEM in this replicate.}
#'   \item{\code{LMEM_len}}{Prediction-region length of LMEM in this replicate.}
#'   \item{\code{Oracle_cov}}{Coverage of Oracle in this replicate.}
#'   \item{\code{Oracle_len}}{Prediction-region length of Oracle in this replicate.}
#' }
#'
#' @export
run_one_replicate <- function(
    n,
    scenario,
    alpha = 0.1,
    n_test = 1,
    n_grid = 200,
    theta_mean = rep(2, 4),
    theta_sigma = diag(1, 4),
    x_mean = rep(0, 4),
    x_sigma = diag(1, 4),
    beta = c(3, 0, 2, 2, 2),
    weight_cap = NULL,
    seed = NULL
) {
  pointwise <- run_one_replicate_pointwise(
    n=n, scenario=scenario, alpha=alpha, n_test=n_test, n_grid=n_grid,
    theta_mean=theta_mean, theta_sigma=theta_sigma, x_mean=x_mean, x_sigma=x_sigma,
    beta=beta, weight_cap=weight_cap, seed=seed)
  out <- colMeans(pointwise, na.rm=TRUE)
  attr(out, "method_diagnostics") <- attr(pointwise, "method_diagnostics")
  attr(out, "oracle_diagnostics") <- attr(pointwise, "oracle_diagnostics")
  attr(out, "input_metadata") <- attr(pointwise, "input_metadata")
  out
}
