# ============================================================
# lmem_region.R
# LMEM-based prediction region
# ============================================================

#' Construct an LMEM-based prediction region using observed clustered data
#'
#' @description
#' This function fits a linear mixed-effects model to the observed clustered data
#' and constructs prediction intervals for new covariate values using
#' \code{merTools::predictInterval()}.
#'
#' The function is intended to implement the LMEM benchmark under the correctly
#' specified random-effects model. Only observed outcomes are used for model
#' fitting.
#'
#' @param dat A data frame containing clustered observations.
#' @param id_col Name of the subject or cluster identifier column.
#' @param y_col Name of the outcome column.
#' @param delta_col Name of the missingness indicator column, where 1 indicates
#'   an observed outcome and 0 indicates a missing outcome.
#' @param x_cols Character vector giving the names of the covariate columns.
#' @param x_test New covariate value(s). This can be a numeric vector
#'   (treated as one test point), or a matrix/data.frame with one row per test point
#'   and \code{length(x_cols)} columns.
#' @param y_grid Numeric vector of candidate response values. This is used only
#'   to convert the model-based prediction interval into a grid-based region for
#'   consistency with other methods.
#' @param alpha Miscoverage level in \eqn{(0,1)}.
#' @param fixed_formula Right-hand side of the fixed-effects formula.
#' @param random_formula Right-hand side of the random-effects formula.
#' @param level Optional confidence level for prediction intervals. If \code{NULL},
#'   it is set to \code{1 - alpha}.
#' @param n_sims Number of simulations used in \code{merTools::predictInterval()}.
#' @param pred_which Which uncertainty component to include in
#'   \code{merTools::predictInterval()}. Default is \code{"random"} to preserve
#'   the current behavior. Use \code{"full"} for legacy gallstones reproduction.
#' @param seed Optional random seed.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{region}}{A list of prediction regions, one for each row of
#'   \code{x_test}, represented as subsets of \code{y_grid}.}
#'   \item{\code{lo_hi}}{A matrix with columns \code{"lo"} and \code{"hi"}
#'   giving the model-based prediction interval endpoints.}
#'   \item{\code{p_final}}{\code{NULL}, since LMEM does not produce conformal
#'   \eqn{p}-values.}
#'   \item{\code{y_grid}}{The candidate grid used to represent the interval as a region.}
#' }
#'
#' @export
lmem_region <- function(
    dat,
    id_col,
    y_col = "Y",
    delta_col = "delta",
    x_cols,
    x_test,
    y_grid,
    alpha = 0.1,
    fixed_formula = NULL,
    random_formula = NULL,
    level = NULL,
    n_sims = 1000L,
    pred_which = "random",
    seed = NULL
) {
  
  # ---------------------------------------------------------------------------
  # Step 0: Validate inputs
  # ---------------------------------------------------------------------------
  if (!is.data.frame(dat)) {
    stop("dat must be a data.frame.")
  }
  
  if (!all(c(id_col, y_col, delta_col) %in% names(dat))) {
    stop("id_col, y_col, and delta_col must all be columns in dat.")
  }
  
  if (!is.character(x_cols) || length(x_cols) < 1L) {
    stop("x_cols must be a non-empty character vector.")
  }
  
  if (!all(x_cols %in% names(dat))) {
    stop("Some columns specified in x_cols are missing from dat.")
  }
  
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number in (0,1).")
  }
  
  if (is.null(level)) {
    level <- 1 - alpha
  }
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("level must be a single number in (0,1).")
  }
  
  n_sims <- as.integer(n_sims)
  if (!is.finite(n_sims) || n_sims < 1L) {
    stop("n_sims must be a positive integer.")
  }
  
  if (!is.character(pred_which) || length(pred_which) != 1L ||
      !pred_which %in% c("full", "fixed", "random", "all")) {
    stop("pred_which must be one of 'full', 'fixed', 'random', or 'all'.")
  }
  
  y_grid <- as.numeric(y_grid)
  if (anyNA(y_grid) || any(!is.finite(y_grid))) {
    stop("y_grid must contain only finite numeric values.")
  }
  y_grid <- sort(unique(y_grid))
  if (length(y_grid) < 2L) {
    stop("y_grid must contain at least two distinct values.")
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # ---------------------------------------------------------------------------
  # Step 1: Standardize x_test input
  # ---------------------------------------------------------------------------
  if (is.vector(x_test) && !is.list(x_test)) {
    x_test <- matrix(as.numeric(x_test), nrow = 1)
  } else if (is.data.frame(x_test)) {
    x_test <- as.matrix(x_test)
  } else {
    x_test <- as.matrix(x_test)
  }
  
  if (ncol(x_test) != length(x_cols)) {
    stop("x_test must have length(x_cols) columns.")
  }
  
  storage.mode(x_test) <- "numeric"
  colnames(x_test) <- x_cols
  
  if (anyNA(x_test) || any(!is.finite(x_test))) {
    stop("x_test contains NA or non-finite values after coercion.")
  }
  
  K <- nrow(x_test)
  
  # ---------------------------------------------------------------------------
  # Helper: Robust binary check for observed outcomes
  # ---------------------------------------------------------------------------
  is_delta1 <- function(d_raw) {
    d_chr <- toupper(trimws(as.character(d_raw)))
    d_num <- suppressWarnings(as.numeric(d_chr))
    (!is.na(d_num) & d_num == 1) | (d_chr %in% c("1", "TRUE", "T"))
  }
  
  # ---------------------------------------------------------------------------
  # Step 2: Restrict to observed data
  # ---------------------------------------------------------------------------
  obs_idx <- is_delta1(dat[[delta_col]]) & !is.na(dat[[y_col]])
  dat_obs <- dat[obs_idx, c(id_col, y_col, x_cols), drop = FALSE]
  
  if (nrow(dat_obs) < 10L) {
    stop("Too few observed outcomes are available for LMEM fitting.")
  }
  
  dat_obs[[id_col]] <- as.factor(dat_obs[[id_col]])
  
  # ---------------------------------------------------------------------------
  # Step 3: Build model formula
  # ---------------------------------------------------------------------------
  if (is.null(fixed_formula)) {
    fixed_formula <- paste(x_cols, collapse = " + ")
  }
  
  if (is.null(random_formula)) {
    random_formula <- paste0("(1 | ", id_col, ")")
  }
  
  full_formula <- stats::as.formula(
    paste(y_col, "~", fixed_formula, "+", random_formula)
  )
  
  # ---------------------------------------------------------------------------
  # Step 4: Fit linear mixed-effects model
  # ---------------------------------------------------------------------------
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required for lmem_region().")
  }
  if (!requireNamespace("merTools", quietly = TRUE)) {
    stop("Package 'merTools' is required for lmem_region().")
  }
  
  fit <- lme4::lmer(
    formula = full_formula,
    data = dat_obs,
    control = lme4::lmerControl(
      check.conv.singular = lme4::.makeCC(action = "ignore", tol = 1e-4)
    )
  )
  
  # ---------------------------------------------------------------------------
  # Step 5: Build test data for a new subject
  # ---------------------------------------------------------------------------
  new_group_label <- "__new_subject__"
  while (new_group_label %in% levels(dat_obs[[id_col]])) {
    new_group_label <- paste0(new_group_label, "_x")
  }
  
  newdata <- data.frame(rep(new_group_label, K))
  names(newdata) <- id_col
  newdata[[id_col]] <- factor(
    newdata[[id_col]],
    levels = c(levels(dat_obs[[id_col]]), new_group_label)
  )
  
  x_test_df <- as.data.frame(x_test)
  colnames(x_test_df) <- x_cols
  newdata[x_cols] <- x_test_df
  
  # ---------------------------------------------------------------------------
  # Step 6: Construct prediction intervals
  # ---------------------------------------------------------------------------
  pred_intervals <- merTools::predictInterval(
    merMod = fit,
    newdata = newdata,
    level = level,
    n.sims = n_sims,
    which = pred_which
  )
  
  lo_hi <- cbind(
    lo = pred_intervals$lwr,
    hi = pred_intervals$upr
  )
  
  # ---------------------------------------------------------------------------
  # Step 7: Convert model-based intervals into grid-based regions
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