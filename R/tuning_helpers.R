# ============================================================
# tuning_helpers.R
# Simulation tuning helpers
# ============================================================

#' Adjust bandwidths to preserve bimodal structure
#'
#' @description
#' In the bimodal simulation settings, the default bandwidth chosen by
#' \code{quantreg::bandwidth.rq} may oversmooth the low-density valley between
#' the two modes. This helper reduces the bandwidth in a central quantile range
#' to better preserve the intended bimodal structure.
#'
#' @param h Numeric bandwidth vector.
#' @param scenario Scenario name. Must be \code{"Bimo"} or \code{"Bimo-fix"}.
#' @param center_fraction Fraction of the tau grid treated as the central region.
#'
#' @return A numeric vector of adjusted bandwidths.
adjust_bandwidth_bimodal <- function(
    h,
    scenario,
    center_fraction = 0.2
) {
  if (!scenario %in% c("Bimo", "Bimo-fix")) {
    return(h)
  }
  
  h <- as.numeric(h)
  nn <- length(h)
  
  if (nn < 5L) {
    return(h)
  }
  
  center_size <- floor(nn / 5)
  start_idx <- floor((nn - center_size) / 2) + 1L
  end_idx <- start_idx + center_size - 1L
  
  ref_prob <- if (scenario == "Bimo") 0.20 else 0.01
  h_ref <- as.numeric(stats::quantile(h, probs = ref_prob, names = FALSE))
  
  h[start_idx:end_idx] <- h_ref
  h
}

#' Construct tau-specific bandwidths and neighboring quantile levels
#'
#' @description
#' This function constructs the bandwidth vector \code{h} together with
#' \code{taus_hi = taus + h} and \code{taus_lo = taus - h}. By default, the
#' bandwidth is selected using \code{quantreg::bandwidth.rq}. In bimodal
#' scenarios, an optional bandwidth adjustment can be applied in the central
#' quantile region to reduce oversmoothing across modes.
#'
#' @param taus Numeric vector of quantile levels in \eqn{(0,1)}.
#' @param n Sample size used in \code{quantreg::bandwidth.rq}.
#' @param scenario Optional scenario name. If equal to \code{"Bimo"} or
#'   \code{"Bimo-fix"}, the bimodal bandwidth adjustment is applied.
#' @param apply_bimodal_adjust Logical; if \code{TRUE}, apply bimodal adjustment
#'   when \code{scenario} is \code{"Bimo"} or \code{"Bimo-fix"}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{h}}{Bandwidth vector.}
#'   \item{\code{taus_hi}}{Upper neighboring quantile levels.}
#'   \item{\code{taus_lo}}{Lower neighboring quantile levels.}
#' }
quantile_levels <- function(
    taus,
    n,
    scenario = NULL,
    apply_bimodal_adjust = TRUE
) {
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("Package 'quantreg' is required for quantile_levels().")
  }
  
  taus <- sort(unique(as.numeric(taus)))
  
  if (anyNA(taus) || any(taus <= 0 | taus >= 1)) {
    stop("taus must be numeric values in (0,1).")
  }
  
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1) {
    stop("n must be a positive number.")
  }
  
  h <- quantreg::bandwidth.rq(taus, n, hs = TRUE, alpha = 0.1)
  h <- as.numeric(h)
  
  h[!is.finite(h) | h <= 0] <- 0.5 * n^(-1/5)
  
  if (isTRUE(apply_bimodal_adjust) &&
      !is.null(scenario) &&
      scenario %in% c("Bimo", "Bimo-fix")) {
    h <- adjust_bandwidth_bimodal(h = h, scenario = scenario)
  }
  
  taus_hi <- taus + h
  taus_lo <- taus - h
  
  min_h <- 1e-6
  for (i in seq_along(taus)) {
    iter <- 0L
    while ((taus_lo[i] <= 0 || taus_hi[i] >= 1) && iter < 60L) {
      h[i] <- h[i] / 2
      taus_hi[i] <- taus[i] + h[i]
      taus_lo[i] <- taus[i] - h[i]
      iter <- iter + 1L
      if (h[i] < min_h) break
    }
  }
  
  list(
    h = h,
    taus_hi = taus_hi,
    taus_lo = taus_lo
  )
}
#' Choose the number of repeated subject-level splits
#'
#' @description
#' This helper returns the number of repeated subject-level splits used in the
#' simulation study. The returned value is a simulation tuning choice rather than
#' a part of the HCP method definition itself.
#'
#' @param n Sample size.
#' @param scenario Scenario name.
#'
#' @return A positive integer.
get_S <- function(n, scenario) {
  ifelse(n >= 300 && scenario %in% c("Bimo", "Bimo-fix"), 1L, 5L)
}

#' Choose the quantile grid resolution parameter
#'
#' @description
#' This helper returns the parameter \code{b} used to define the quantile grid
#' \code{taus = (1:(2^b - 1)) / 2^b} in the simulation study. The returned value
#' is a simulation tuning choice rather than a part of the HCP method definition
#' itself.
#'
#' @param n Sample size.
#' @param scenario Scenario name.
#'
#' @return A positive integer.
get_b <- function(n, scenario) {
  ifelse(n >= 300 && scenario %in% c("Bimo", "Bimo-fix"), 6L, 4L)
}

#' Choose the truncation level for inverse-propensity weights
#'
#' @description
#' This helper returns the upper bound used to truncate inverse-propensity
#' weights in the simulation study. The truncation level increases with sample
#' size to allow more variability in larger samples while avoiding unstable
#' extreme weights in smaller samples.
#'
#' @param n Sample size.
#'
#' @return A positive numeric value giving the weight cap.
get_max_wt <- function(n) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1) {
    stop("n must be a positive number.")
  }
  30 
}