# ============================================================
# cd4_data_helpers.R
# Helper functions for the CD4 real-data analysis
# ============================================================

#' Read and standardize the CD4 longitudinal dataset
#'
#' @description
#' This function reads the CD4 longitudinal dataset, checks that the required
#' variables are present, renames the outcome variable from `CD4` to `Y`,
#' converts variables to the expected numeric types, and orders the data by
#' subject ID and time.
#'
#' @param file_path Character string. Path to the CD4 data file.
#'
#' @return A data frame containing the standardized CD4 data with columns
#' `id`, `Y`, and the CD4 covariates used in the analysis.
read_cd4_data <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("CD4 data file does not exist: ", file_path)
  }
  
  dat <- read.table(file_path, header = TRUE, quote = "\"")
  
  required_cols <- c(
    "id", "CD4",
    "time", "age", "smoke", "drug", "partners", "cesd"
  )
  
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0L) {
    stop(
      "The CD4 data file is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  dat <- dat[, required_cols, drop = FALSE]
  names(dat)[names(dat) == "CD4"] <- "Y"
  
  dat$id <- as.integer(dat$id)
  dat$Y <- as.numeric(dat$Y)
  
  x_cols <- get_cd4_x_cols()
  for (nm in x_cols) {
    dat[[nm]] <- as.numeric(dat[[nm]])
  }
  
  dat <- dat[order(dat$id, dat$time), , drop = FALSE]
  rownames(dat) <- NULL
  
  dat
}


#' Return the covariate names used in the CD4 analysis
#'
#' @description
#' This function returns the covariate names used as predictors in the CD4
#' real-data analysis.
#'
#' @return A character vector of CD4 covariate names.
get_cd4_x_cols <- function() {
  c("time", "age", "smoke", "drug", "partners", "cesd")
}


#' Return missingness coefficients for the CD4 analysis
#'
#' @description
#' This function returns the logistic missingness-model coefficients used to
#' generate artificial missingness in the CD4 analysis.
#'
#' @param missing_rate Character or numeric value indicating the missingness
#' setting. Supported values are `"20"` and `"50"`.
#'
#' @return A numeric vector of coefficients for the logistic missingness model,
#' including the intercept.
get_cd4_missing_beta <- function(missing_rate = c("20", "50")) {
  missing_rate <- match.arg(as.character(missing_rate), choices = c("20", "50"))
  
  switch(
    missing_rate,
    "20" = c(7, 0, 0, -3, 0, 0, 0),
    "50" = c(1, 0, 0, -3, 0, 0, 0)
  )
}


#' Generate a missingness indicator from covariates
#'
#' @description
#' This function generates binary missingness indicators from a logistic model
#' using the supplied covariate matrix and coefficient vector.
#'
#' @param X Numeric matrix or data frame of covariates.
#' @param beta Numeric coefficient vector. Its length must be equal to
#' `ncol(X) + 1`, including the intercept.
#'
#' @return An integer vector of missingness indicators generated from a
#' Bernoulli distribution.
generate_cd4_missing_indicator <- function(X, beta) {
  X <- as.matrix(X)
  storage.mode(X) <- "numeric"
  
  if (anyNA(X) || any(!is.finite(X))) {
    stop("X contains NA or non-finite values.")
  }
  
  beta <- as.numeric(beta)
  if (length(beta) != ncol(X) + 1L) {
    stop("length(beta) must be equal to ncol(X) + 1.")
  }
  
  X_tilde <- cbind(1, X)
  log_odds <- as.numeric(X_tilde %*% beta)
  prob <- 1 / (1 + exp(-log_odds))
  
  stats::rbinom(nrow(X), size = 1, prob = prob)
}


#' Add artificial missingness to the CD4 dataset
#'
#' @description
#' This function adds a binary missingness indicator `delta` to the CD4 dataset
#' using the CD4-specific logistic missingness model.
#'
#' @param dat Data frame. Standardized CD4 dataset containing `id`, `Y`, and
#' the CD4 covariates.
#' @param missing_rate Character or numeric value indicating the missingness
#' setting. Supported values are `"20"` and `"50"`.
#' @param seed Integer or `NULL`. Random seed used to generate missingness.
#'
#' @return The input data frame with an additional column `delta`.
add_cd4_missingness <- function(
    dat,
    missing_rate = c("20", "50"),
    seed = 123
) {
  missing_rate <- match.arg(as.character(missing_rate), choices = c("20", "50"))
  
  if (!is.data.frame(dat)) {
    stop("dat must be a data.frame.")
  }
  
  required_cols <- c("id", "Y", get_cd4_x_cols())
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0L) {
    stop(
      "dat is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  x_cols <- get_cd4_x_cols()
  beta <- get_cd4_missing_beta(missing_rate)
  
  dat$delta <- generate_cd4_missing_indicator(
    X = dat[, x_cols, drop = FALSE],
    beta = beta
  )
  
  dat
}


#' Return CD4-specific inverse-propensity weight truncation level
#'
#' @description
#' This function returns the inverse-propensity weight truncation level used in
#' the CD4 real-data analysis. The truncation level depends on the method and
#' the missingness setting.
#'
#' @param method Character string. Method name. Supported values are `"HCP"`,
#' `"LC"`, `"DWR"`, and `"LMEM"`.
#' @param missing_rate Character or numeric value indicating the missingness
#' setting. Supported values are `"20"` and `"50"`.
#'
#' @return A numeric weight cap. Returns `Inf` for methods without truncation.
get_cd4_weight_cap <- function(
    method = c("HCP", "LC", "DWR", "LMEM"),
    missing_rate = c("20", "50")
) {
  method <- match.arg(method)
  missing_rate <- match.arg(as.character(missing_rate), choices = c("20", "50"))
  
  if (method == "HCP") {
    return(if (missing_rate == "20") 8 else 13)
  }
  
  if (method == "LC") {
    return(if (missing_rate == "20") 8 else 13)
  }
  
  Inf
}


#' Build a response grid for CD4 leave-one-subject-out analysis
#'
#' @description
#' This function constructs an equally spaced response grid over the observed
#' finite outcome range in a CD4 sample. The grid is used to evaluate prediction
#' regions.
#'
#' @param dat_sample Data frame containing the sample used to define the grid.
#' @param y_col Character string. Name of the outcome column.
#' @param n_grid Integer. Number of grid points.
#'
#' @return A numeric vector of response-grid values.
make_cd4_y_grid <- function(
    dat_sample,
    y_col = "Y",
    n_grid = 200
) {
  if (!is.data.frame(dat_sample)) {
    stop("dat_sample must be a data.frame.")
  }
  if (!y_col %in% names(dat_sample)) {
    stop("y_col is not a column in dat_sample.")
  }
  
  n_grid <- as.integer(n_grid)
  if (!is.finite(n_grid) || n_grid < 2L) {
    stop("n_grid must be an integer >= 2.")
  }
  
  y <- as.numeric(dat_sample[[y_col]])
  y <- y[is.finite(y)]
  
  if (length(y) < 2L) {
    stop("dat_sample must contain at least two finite outcome values.")
  }
  
  seq(
    from = min(y),
    to = max(y),
    length.out = n_grid
  )
}


#' Summarize CD4 prediction results for one held-out subject
#'
#' @description
#' This function summarizes evaluation results for one held-out CD4 subject by
#' computing coverage and average prediction-region length. For pointwise
#' prediction, coverage is averaged over observations. For simultaneous
#' prediction, coverage is one only if all observations are covered.
#'
#' @param eval_mat Matrix or matrix-like object containing columns `covered`
#' and `length`.
#' @param prediction_type Character string. Type of prediction summary.
#' Supported values are `"pointwise"` and `"simultaneous"`.
#'
#' @return A numeric vector with entries `coverage` and `length`.
summarize_cd4_leaveout <- function(
    eval_mat,
    prediction_type = c("pointwise", "simultaneous")
) {
  prediction_type <- match.arg(prediction_type)
  
  if (!is.matrix(eval_mat)) {
    eval_mat <- as.matrix(eval_mat)
  }
  
  if (!all(c("covered", "length") %in% colnames(eval_mat))) {
    stop("eval_mat must contain columns named 'covered' and 'length'.")
  }
  
  covered <- as.numeric(eval_mat[, "covered"])
  len <- as.numeric(eval_mat[, "length"])
  
  if (prediction_type == "pointwise") {
    cov_out <- mean(covered, na.rm = TRUE)
  } else {
    cov_out <- as.numeric(length(covered) > 0L && all(covered == 1))
  }
  
  c(
    coverage = cov_out,
    length = mean(len, na.rm = TRUE)
  )
}