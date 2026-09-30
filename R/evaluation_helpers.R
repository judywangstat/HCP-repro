# ============================================================
# evaluation_helpers.R
# Helpers for simulation evaluation
# ============================================================

#' Evaluate a single-interval prediction region
#'
#' @description
#' This function evaluates coverage and interval length for methods that return
#' a single interval per test point through \code{lo_hi}.
#'
#' It is intended for interval-based methods such as DWR, LC, LMEM, and Oracle.
#' It can also be used for HCP in scenarios where the retained region is treated
#' as a single connected interval.
#'
#' @param res A method output list containing \code{lo_hi}.
#' @param y_true Numeric vector of true outcomes, one for each test point.
#'
#' @return A numeric matrix with two columns:
#' \describe{
#'   \item{\code{covered}}{Indicator of whether the true outcome is covered.}
#'   \item{\code{length}}{Interval length \code{hi - lo}.}
#' }
#'
#' @export
evaluate_interval_region <- function(res, y_true) {
  if (is.null(res$lo_hi)) {
    stop("res must contain 'lo_hi'.")
  }
  
  lo_hi <- as.matrix(res$lo_hi)
  y_true <- as.numeric(y_true)
  
  if (ncol(lo_hi) != 2L) {
    stop("res$lo_hi must have exactly two columns.")
  }
  
  if (length(y_true) != nrow(lo_hi)) {
    stop("length(y_true) must equal nrow(res$lo_hi).")
  }
  
  if (!is.null(colnames(lo_hi)) && all(c("lo", "hi") %in% colnames(lo_hi))) {
    lo <- lo_hi[, "lo"]
    hi <- lo_hi[, "hi"]
  } else {
    lo <- lo_hi[, 1]
    hi <- lo_hi[, 2]
  }
  
  covered <- logical(length(y_true))
  length_out <- numeric(length(y_true))
  
  for (k in seq_along(y_true)) {
    if (is.na(lo[k]) || is.na(hi[k])) {
      covered[k] <- FALSE
      length_out[k] <- NA_real_
    } else {
      covered[k] <- (y_true[k] >= lo[k] && y_true[k] <= hi[k])
      length_out[k] <- hi[k] - lo[k]
    }
  }
  
  cbind(
    covered = as.numeric(covered),
    length = length_out
  )
}

#' Evaluate a bimodal grid-based prediction region
#'
#' @description
#' This function evaluates coverage and region length for grid-based prediction
#' regions that may contain two disconnected components, as in the bimodal
#' simulation scenarios.
#'
#' If the retained region contains more than two apparent components, the largest
#' gap is used as the main split, so that the region is treated as two parts.
#' This matches the simulation setting in which the target distribution is at
#' most bimodal.
#'
#' @param res A method output list containing \code{region} and \code{y_grid}.
#' @param y_true Numeric vector of true outcomes, one for each test point.
#'
#' @return A numeric matrix with two columns:
#' \describe{
#'   \item{\code{covered}}{Indicator of whether the true outcome is covered.}
#'   \item{\code{length}}{Total length of the retained bimodal region.}
#' }
#'
#' @export
evaluate_bimodal_region <- function(res, y_true) {
  if (is.null(res$region) || is.null(res$y_grid)) {
    stop("res must contain 'region' and 'y_grid'.")
  }
  
  y_true <- as.numeric(y_true)
  K <- length(res$region)
  
  if (length(y_true) != K) {
    stop("length(y_true) must equal length(res$region).")
  }
  
  y_grid <- sort(unique(as.numeric(res$y_grid)))
  if (length(y_grid) < 2L) {
    stop("res$y_grid must contain at least two distinct values.")
  }
  
  grid_step <- min(diff(y_grid))
  tol <- grid_step / 2
  
  covered <- logical(K)
  length_out <- numeric(K)
  
  for (k in seq_len(K)) {
    reg <- as.numeric(res$region[[k]])
    
    if (length(reg) == 0L || all(!is.finite(reg))) {
      covered[k] <- FALSE
      length_out[k] <- NA_real_
      next
    }
    
    reg <- sort(unique(reg))
    
    # one-point region
    if (length(reg) == 1L) {
      covered[k] <- abs(y_true[k] - reg[1]) <= tol
      length_out[k] <- 0
      next
    }
    
    gaps <- diff(reg)
    candidate_idx <- which(gaps > 1.5 * grid_step)
    
    # no major gap: treat as one connected interval
    if (length(candidate_idx) == 0L) {
      lo <- min(reg)
      hi <- max(reg)
      covered[k] <- (y_true[k] >= lo && y_true[k] <= hi)
      length_out[k] <- hi - lo
      
    } else {
      # use the largest gap as the main split
      split_idx <- candidate_idx[which.max(gaps[candidate_idx])]
      
      reg1 <- reg[1:split_idx]
      reg2 <- reg[(split_idx + 1):length(reg)]
      
      lo1 <- min(reg1)
      hi1 <- max(reg1)
      lo2 <- min(reg2)
      hi2 <- max(reg2)
      
      covered[k] <-
        (y_true[k] >= lo1 && y_true[k] <= hi1) ||
        (y_true[k] >= lo2 && y_true[k] <= hi2)
      
      length_out[k] <- (hi1 - lo1) + (hi2 - lo2)
    }
  }
  
  cbind(
    covered = as.numeric(covered),
    length = length_out
  )
}
# Continuous union evaluator for the HPD Oracle. Coverage uses component
# membership and length sums component measures, without filling the gaps.
evaluate_oracle_union <- function(res,y_true) {
  if(is.null(res$oracles)||length(res$oracles)!=length(y_true))
    stop('One HPD Oracle object is required for each test outcome.')
  out<-t(vapply(seq_along(y_true),function(k)
    as.numeric(evaluate_hpd_union(res$oracles[[k]],y_true[k])[1,]),numeric(2)))
  colnames(out)<-c('covered','length');out
}
