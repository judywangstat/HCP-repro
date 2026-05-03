# ============================================================
# gallstones_data_helpers.R
# Helpers for loading and preparing the gallstones data
# ============================================================

#' Load and preprocess the gallstones data
#'
#' @description
#' Reads the gallstones data from a tab-delimited text file and converts it into
#' the long-format data structure used by the real-data reproduction scripts.
#' The covariates are three visit-time indicators and the treatment indicator.
#'
#' @param data_file Path to the gallstones data file.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{\code{id}}{Subject identifier.}
#'   \item{\code{time}}{Original visit time.}
#'   \item{\code{Y}}{Outcome value, with missing values retained as NA.}
#'   \item{\code{delta}}{Missingness indicator; 1 if Y is observed and 0 otherwise.}
#'   \item{\code{T1}, \code{T2}, \code{T3}}{Visit-time indicators.}
#'   \item{\code{Treat}}{Treatment indicator.}
#' }
#'
#' @export
load_gallstones_data <- function(data_file) {
  if (!file.exists(data_file)) {
    stop("data_file does not exist: ", data_file)
  }
  
  raw_dat <- utils::read.delim(
    data_file,
    na.strings = c("#NULL!", "NA", "")
  )
  
  required_cols <- c("id", "time", "trt", "y")
  if (!all(required_cols %in% names(raw_dat))) {
    stop(
      "The gallstones data must contain columns: ",
      paste(required_cols, collapse = ", ")
    )
  }
  
  Y_raw <- as.numeric(raw_dat$y)
  dat <- data.frame(
    id = raw_dat$id,
    time = raw_dat$time,
    Y = Y_raw,
    delta = as.integer(!is.na(Y_raw)),
    T1 = as.integer(raw_dat$time == 12),
    T2 = as.integer(raw_dat$time == 20),
    T3 = as.integer(raw_dat$time == 24),
    Treat = as.numeric(raw_dat$trt)
  )
  
  dat <- dat[order(dat$id, dat$time), , drop = FALSE]
  rownames(dat) <- NULL
  
  dat
}


#' Impute missing gallstones outcomes for evaluation
#'
#' @description
#' Imputes missing outcomes in the gallstones data using the same strategy as the
#' original reproduction code. A random forest estimates the conditional mean,
#' and a quantile forest estimates an interquartile-range scale. Missing outcomes
#' are then filled by adding Gaussian noise to the random-forest prediction.
#'
#' This imputed outcome is used only as a pseudo-truth for evaluating coverage in
#' the leave-one-subject-out analysis. The original observed/missing indicator is
#' kept unchanged.
#'
#' @param dat A data frame returned by \code{load_gallstones_data()}.
#' @param y_col Outcome column name.
#' @param x_cols Covariate column names used for imputation.
#' @param seed Random seed used for reproducible imputation.
#' @param noise_scale Multiplier applied to the estimated interquartile range.
#'
#' @return The input data frame with an additional column \code{Y_eval}, where
#' missing outcomes have been imputed.
#'
#' @export
impute_gallstones_outcomes <- function(
    dat,
    y_col = "Y",
    x_cols = c("T1", "T2", "T3", "Treat"),
    seed = 123,
    noise_scale = 1.2
) {
  if (!is.data.frame(dat)) {
    stop("dat must be a data.frame.")
  }
  
  if (!all(c(y_col, x_cols) %in% names(dat))) {
    stop("Some columns specified by y_col or x_cols are missing from dat.")
  }
  
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("Package 'randomForest' is required for gallstones imputation.")
  }
  
  if (!requireNamespace("grf", quietly = TRUE)) {
    stop("Package 'grf' is required for gallstones imputation.")
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  X <- as.matrix(dat[, x_cols, drop = FALSE])
  storage.mode(X) <- "numeric"
  
  Y <- as.numeric(dat[[y_col]])
  obs_idx <- !is.na(Y)
  miss_idx <- is.na(Y)
  
  if (!any(miss_idx)) {
    dat$Y_eval <- Y
    return(dat)
  }
  
  if (sum(obs_idx) < 10L) {
    stop("Too few observed outcomes for imputation.")
  }
  
  mean_fit <- randomForest::randomForest(
    x = X[obs_idx, , drop = FALSE],
    y = Y[obs_idx]
  )
  
  mean_pred <- stats::predict(
    mean_fit,
    X[miss_idx, , drop = FALSE]
  )
  
  q_fit <- grf::quantile_forest(
    X = X[obs_idx, , drop = FALSE],
    Y = Y[obs_idx],
    quantiles = c(0.25, 0.75)
  )
  
  q_pred <- stats::predict(
    q_fit,
    newdata = X[miss_idx, , drop = FALSE]
  )$predictions
  
  scale_pred <- q_pred[, 2] - q_pred[, 1]
  noise <- stats::rnorm(length(mean_pred))
  
  Y[miss_idx] <- mean_pred + noise_scale * scale_pred * noise
  
  dat$Y_eval <- as.numeric(Y)
  dat
}


#' Create a response grid for gallstones prediction regions
#'
#' @description
#' Constructs the candidate response grid from the training sample, matching the
#' original gallstones reproduction code.
#'
#' @param dat_sample Training/calibration data after excluding the held-out
#'   subject.
#' @param y_col Outcome column name.
#' @param n_grid Number of grid points.
#'
#' @return A numeric vector of grid values.
#'
#' @export
make_gallstones_y_grid <- function(
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
  
  seq(
    from = min(y, na.rm = TRUE),
    to = max(y, na.rm = TRUE),
    length.out = n_grid
  )
}