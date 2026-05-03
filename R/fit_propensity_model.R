#' Fit a missingness propensity model from pooled observations
#'
#' @description
#' This function fits the missingness propensity
#' \eqn{\pi(x) = \mathbb{P}(\delta = 1 \mid x)}
#' using pooled observations under a marginal missingness model.
#' The missingness mechanism can be estimated using logistic regression,
#' generalized random forests, or gradient boosting.
#'
#' Both continuous and categorical covariates are supported. Character variables
#' are converted to factors, and categorical covariates are expanded into dummy
#' variables using \code{model.matrix()}.
#'
#' @param dat A data frame containing the missingness indicator and covariates.
#' @param delta_col Name of the missingness indicator column. The indicator should
#'   represent whether the outcome is observed (\code{1}) or missing (\code{0}).
#' @param x_cols Character vector giving the names of the covariate columns used
#'   to model missingness.
#' @param method Estimation method. Must be one of \code{"logistic"},
#'   \code{"grf"}, or \code{"boosting"}.
#' @param eps Small symmetric clipping constant used to truncate the estimated
#'   missingness propensity away from both 0 and 1, that is, to
#'   \eqn{[\epsilon, 1-\epsilon]}.
#' @param ... Additional arguments passed to the selected learner.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{predict}}{A function that evaluates the estimated missingness
#'   propensity \eqn{\hat\pi(x)} at new covariate values.}
#'   \item{\code{fit_object}}{The fitted missingness model object.}
#'   \item{\code{method}}{The estimation method used.}
#' }
#'
#' @export
fit_propensity_model <- function(
    dat,
    delta_col = "delta",
    x_cols,
    method = c("logistic", "grf", "boosting"),
    eps = 1e-6,
    ...
) {
  
  # ---------------------------------------------------------------------------
  # Step 0: Match arguments and perform basic input checks
  # ---------------------------------------------------------------------------
  method <- match.arg(method)
  
  if (!is.data.frame(dat)) {
    stop("dat must be a data.frame.")
  }
  
  if (!is.character(delta_col) || length(delta_col) != 1L) {
    stop("delta_col must be a single character string.")
  }
  
  if (!delta_col %in% names(dat)) {
    stop("delta_col is not a column in dat.")
  }
  
  if (!is.character(x_cols) || length(x_cols) < 1L) {
    stop("x_cols must be a non-empty character vector.")
  }
  
  if (anyDuplicated(x_cols)) {
    stop("x_cols contains duplicated column names.")
  }
  
  if (!all(x_cols %in% names(dat))) {
    stop("Some columns specified in x_cols are missing from dat.")
  }
  
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0 || eps >= 0.5) {
    stop("eps must be a single finite number in (0, 0.5).")
  }
  
  dots <- list(...)
  
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  # ---------------------------------------------------------------------------
  # Helper: Parse a binary missingness indicator into {0, 1}
  # ---------------------------------------------------------------------------
  parse_binary_indicator <- function(z, name = "delta") {
    z_chr <- toupper(trimws(as.character(z)))
    out <- rep(NA_integer_, length(z_chr))
    
    out[z_chr %in% c("1", "TRUE", "T")] <- 1L
    out[z_chr %in% c("0", "FALSE", "F")] <- 0L
    
    suppressWarnings(z_num <- as.numeric(z_chr))
    idx_num <- is.na(out) & !is.na(z_num)
    out[idx_num] <- as.integer(z_num[idx_num])
    
    if (anyNA(out)) {
      stop(name, " contains missing or unparseable values after binary parsing.")
    }
    if (any(!out %in% c(0L, 1L))) {
      stop(name, " must be binary after parsing.")
    }
    
    out
  }
  
  # ---------------------------------------------------------------------------
  # Step 1: Parse the missingness indicator
  # ---------------------------------------------------------------------------
  delta <- parse_binary_indicator(dat[[delta_col]], name = delta_col)
  
  # ---------------------------------------------------------------------------
  # Step 2: Build the pooled design matrix for the missingness model
  # ---------------------------------------------------------------------------
  x_df <- dat[, x_cols, drop = FALSE]
  
  # Convert character variables to factors so that model.matrix handles them properly.
  for (nm in names(x_df)) {
    if (is.character(x_df[[nm]])) {
      x_df[[nm]] <- as.factor(x_df[[nm]])
    }
  }
  
  # Store factor levels from the training data for future prediction.
  x_levels <- lapply(x_df, function(v) if (is.factor(v)) levels(v) else NULL)
  
  # Build a design matrix with an intercept and dummy variables as needed.
  terms_obj <- stats::terms(stats::as.formula("~ ."), data = x_df)
  X <- stats::model.matrix(terms_obj, data = x_df)
  
  if (nrow(X) != length(delta)) {
    stop("Internal error: design matrix and delta have incompatible dimensions.")
  }
  
  # Make column names safe for downstream learners.
  colnames(X) <- make.names(colnames(X), unique = TRUE)
  
  # ---------------------------------------------------------------------------
  # Step 3A: Fit logistic regression
  # ---------------------------------------------------------------------------
  fit_logistic <- function() {
    df_glm <- as.data.frame(X)
    df_glm$.delta <- delta
    
    dots_local <- dots
    if (is.null(dots_local$family)) {
      dots_local$family <- stats::binomial()
    }
    
    do.call(
      stats::glm,
      c(list(formula = .delta ~ ., data = df_glm), dots_local)
    )
  }
  
  # ---------------------------------------------------------------------------
  # Step 3B: Fit generalized random forest
  # ---------------------------------------------------------------------------
  fit_grf <- function() {
    if (!requireNamespace("grf", quietly = TRUE)) {
      stop("Package 'grf' is required for method = 'grf'.")
    }
    
    y_fac <- factor(delta, levels = c(0, 1))
    
    do.call(
      grf::probability_forest,
      c(list(X = X, Y = y_fac), dots)
    )
  }
  
  # ---------------------------------------------------------------------------
  # Step 3C: Fit gradient boosting model
  # ---------------------------------------------------------------------------
  fit_boosting <- function() {
    if (!requireNamespace("xgboost", quietly = TRUE)) {
      stop("Package 'xgboost' is required for method = 'boosting'.")
    }
    
    dtrain <- xgboost::xgb.DMatrix(data = X, label = delta)
    
    dots_local <- dots
    params <- dots_local$params
    if (is.null(params)) {
      params <- list()
    }
    params$objective <- params$objective %||% "binary:logistic"
    params$eval_metric <- params$eval_metric %||% "logloss"
    
    nrounds <- dots_local$nrounds %||% 400L
    dots_local$params <- NULL
    dots_local$nrounds <- NULL
    
    do.call(
      xgboost::xgb.train,
      c(list(params = params, data = dtrain, nrounds = nrounds, verbose = 0), dots_local)
    )
  }
  
  # ---------------------------------------------------------------------------
  # Step 4: Fit the selected missingness model
  # ---------------------------------------------------------------------------
  fit_object <- switch(
    method,
    logistic = fit_logistic(),
    grf = fit_grf(),
    boosting = fit_boosting()
  )
  
  # ---------------------------------------------------------------------------
  # Helper: Build a prediction design matrix for new covariate values
  # ---------------------------------------------------------------------------
  build_X_new <- function(x_new) {
    if (is.matrix(x_new)) {
      x_new <- as.data.frame(x_new)
    }
    
    if (!is.data.frame(x_new)) {
      stop("x_new must be a data.frame or matrix.")
    }
    
    if (!all(x_cols %in% names(x_new))) {
      stop("x_new must contain columns: ", paste(x_cols, collapse = ", "))
    }
    
    x_new <- x_new[, x_cols, drop = FALSE]
    
    for (nm in names(x_new)) {
      if (is.character(x_new[[nm]])) {
        x_new[[nm]] <- as.factor(x_new[[nm]])
      }
      
      if (is.factor(x_new[[nm]]) && !is.null(x_levels[[nm]])) {
        x_new[[nm]] <- factor(as.character(x_new[[nm]]), levels = x_levels[[nm]])
      }
    }
    
    X_new <- stats::model.matrix(terms_obj, data = x_new)
    colnames(X_new) <- make.names(colnames(X_new), unique = TRUE)
    
    # Add any missing training columns as zeros.
    miss <- setdiff(colnames(X), colnames(X_new))
    if (length(miss) > 0) {
      add <- matrix(0, nrow = nrow(X_new), ncol = length(miss))
      colnames(add) <- miss
      X_new <- cbind(X_new, add)
    }
    
    X_new <- X_new[, colnames(X), drop = FALSE]
    X_new
  }
  
  # ---------------------------------------------------------------------------
  # Output function: Predict missingness propensity at new covariate values
  # ---------------------------------------------------------------------------
  predict_fun <- function(x_new) {
    X_new <- build_X_new(x_new)
    
    p_hat <- switch(
      method,
      logistic = {
        as.numeric(
          stats::predict(fit_object, newdata = as.data.frame(X_new), type = "response")
        )
      },
      grf = {
        pred_obj <- stats::predict(fit_object, newdata = X_new)
        pr <- pred_obj$predictions
        
        if (is.matrix(pr)) {
          if (ncol(pr) == 1) {
            as.numeric(pr[, 1])
          } else {
            as.numeric(pr[, 2])
          }
        } else {
          as.numeric(pr)
        }
      },
      boosting = {
        dtest <- xgboost::xgb.DMatrix(data = X_new)
        as.numeric(stats::predict(fit_object, newdata = dtest))
      }
    )
    
    pmax(pmin(p_hat, 1 - eps), eps)
  }
  
  # ---------------------------------------------------------------------------
  # Return fitted model object and user-facing prediction function
  # ---------------------------------------------------------------------------
  list(
    method = method,
    fit_object = fit_object,
    predict = predict_fun
  )
}