#' Fit a conditional density estimator via a quantile-process approach
#'
#' @description
#' Fits a conditional quantile model using pooled observed outcomes and converts
#' the fitted quantile curve into a conditional density estimator through a
#' quotient-based approximation.
#'
#' @param dat A data frame in long format.
#' @param y_col Outcome column name.
#' @param delta_col Missingness indicator column name; 1 means observed.
#' @param x_cols Character vector of covariate column names.
#' @param taus Quantile levels in (0,1).
#' @param h Optional bandwidth; either NULL, a scalar, or a vector of length
#'   \code{length(taus)}.
#' @param method Quantile engine: \code{"rq"}, \code{"qrf"}, or \code{"grf"}.
#'   \code{"qrf"} uses \code{quantregForest}; \code{"grf"} uses
#'   \code{grf::quantile_forest()}.
#' @param enforce_monotone Logical; if TRUE, apply isotonic adjustment instead
#'   of the crossing-removal rule.
#' @param tail_decay Logical; if TRUE, add decaying tail points for interpolation.
#' @param num_extra_points Number of tail points on each side.
#' @param decay_factor Tail decay factor.
#' @param seed Optional random seed.
#' @param ... Additional arguments passed to the selected quantile engine.
#'
#' @return A list with elements \code{taus}, \code{fit_object},
#'   \code{predict_Q}, \code{predict_density}, and \code{predict_density_grid}.
fit_cond_density_qp <- function(
    dat,
    y_col = "Y",
    delta_col = "delta",
    x_cols,
    taus = seq(0.05, 0.95, by = 0.01),
    h = NULL,
    method = c("rq", "qrf", "grf"),
    enforce_monotone = FALSE,
    tail_decay = TRUE,
    num_extra_points = 20L,
    decay_factor = 0.7,
    seed = NULL,
    ...
) {
  method <- match.arg(method)
  
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (!is.data.frame(dat)) {
    stop("dat must be a data.frame.")
  }
  if (!all(c(y_col, delta_col, x_cols) %in% names(dat))) {
    stop("Some columns specified by y_col, delta_col, or x_cols are missing from dat.")
  }
  
  taus <- sort(unique(as.numeric(taus)))
  if (anyNA(taus) || any(taus <= 0 | taus >= 1)) {
    stop("taus must lie strictly inside (0,1).")
  }
  if (length(taus) < 5L) {
    stop("taus must contain at least 5 distinct values.")
  }
  
  if (!is.logical(enforce_monotone) || length(enforce_monotone) != 1L) {
    stop("enforce_monotone must be TRUE or FALSE.")
  }
  if (!is.logical(tail_decay) || length(tail_decay) != 1L) {
    stop("tail_decay must be TRUE or FALSE.")
  }
  
  num_extra_points <- as.integer(num_extra_points)
  if (num_extra_points < 0L) {
    stop("num_extra_points must be nonnegative.")
  }
  if (!is.numeric(decay_factor) || length(decay_factor) != 1L ||
      !is.finite(decay_factor) || decay_factor <= 0 || decay_factor >= 1) {
    stop("decay_factor must be a finite number in (0,1).")
  }
  
  parse_delta1 <- function(z) {
    z_chr <- toupper(trimws(as.character(z)))
    z_num <- suppressWarnings(as.numeric(z_chr))
    (!is.na(z_num) & z_num == 1) | (z_chr %in% c("1", "TRUE", "T"))
  }
  
  assert_numeric_x <- function(x, x_cols, name = "x_new") {
    if (is.vector(x) && !is.list(x)) {
      xm <- matrix(as.numeric(x), nrow = 1L)
      if (length(xm) != length(x_cols) || anyNA(xm)) {
        stop(name, " must be a numeric vector of length length(x_cols).")
      }
    } else if (is.data.frame(x)) {
      if (ncol(x) != length(x_cols)) {
        stop(name, " must have length(x_cols) columns.")
      }
      bad <- names(x)[!vapply(x, is.numeric, logical(1))]
      if (length(bad) > 0L) {
        stop(name, " has non-numeric columns: ", paste(bad, collapse = ", "), ".")
      }
      xm <- as.matrix(x)
    } else {
      xm <- as.matrix(x)
      if (!is.numeric(xm)) {
        suppressWarnings(storage.mode(xm) <- "numeric")
      }
      if (anyNA(xm) || ncol(xm) != length(x_cols)) {
        stop(name, " must be coercible to a numeric matrix with length(x_cols) columns.")
      }
    }
    colnames(xm) <- x_cols
    xm
  }
  
  # ---------------------------------------------------------------------------
  # Step 1: observed outcomes only
  # ---------------------------------------------------------------------------
  dobs <- dat[parse_delta1(dat[[delta_col]]) & !is.na(dat[[y_col]]), , drop = FALSE]
  n_obs <- nrow(dobs)
  if (n_obs < 20L) {
    stop("Too few observed outcomes after filtering delta == 1.")
  }
  
  x_df <- dobs[, x_cols, drop = FALSE]
  bad_x <- names(x_df)[!vapply(x_df, is.numeric, logical(1))]
  if (length(bad_x) > 0L) {
    stop("Non-numeric covariates found in x_cols: ", paste(bad_x, collapse = ", "), ".")
  }
  
  X_obs <- as.matrix(x_df)
  y_obs <- as.numeric(dobs[[y_col]])
  
  # ---------------------------------------------------------------------------
  # Step 2: bandwidth construction
  # ---------------------------------------------------------------------------
  if (is.null(h)) {
    if (!requireNamespace("quantreg", quietly = TRUE)) {
      stop("Package 'quantreg' is required when h = NULL.")
    }
    h_vec <- tryCatch(
      quantreg::bandwidth.rq(taus, n_obs, hs = TRUE, alpha = 0.1),
      error = function(e) rep(0.5 * n_obs^(-1 / 5), length(taus))
    )
  } else if (length(h) == 1L) {
    h_vec <- rep(as.numeric(h), length(taus))
  } else {
    h_vec <- as.numeric(h)
    if (length(h_vec) != length(taus)) {
      stop("If h is a vector, it must have length equal to length(taus).")
    }
  }
  
  h_vec[!is.finite(h_vec) | h_vec <= 0] <- 0.5 * n_obs^(-1 / 5)
  
  taus_hi <- taus + h_vec
  taus_lo <- taus - h_vec
  min_h <- 1e-6
  
  for (i in seq_along(taus)) {
    iter <- 0L
    while ((taus_lo[i] <= 0 || taus_hi[i] >= 1) && iter < 60L) {
      h_vec[i] <- h_vec[i] / 2
      taus_hi[i] <- taus[i] + h_vec[i]
      taus_lo[i] <- taus[i] - h_vec[i]
      iter <- iter + 1L
      if (h_vec[i] < min_h) break
    }
  }
  
  ok_tau <- is.finite(h_vec) & (h_vec > min_h) & (taus_lo > 0) & (taus_hi < 1)
  if (sum(ok_tau) < 5L) {
    stop("Too few tau values remain after enforcing tau +/- h in (0,1).")
  }
  
  taus0 <- taus[ok_tau]
  h0 <- h_vec[ok_tau]
  taus_minus <- taus_lo[ok_tau]
  taus_plus <- taus_hi[ok_tau]
  taus_pm <- sort(unique(c(taus_minus, taus0, taus_plus)))
  
  idx_m <- match(taus_minus, taus_pm)
  idx_0 <- match(taus0, taus_pm)
  idx_p <- match(taus_plus, taus_pm)
  
  if (anyNA(idx_m) || anyNA(idx_0) || anyNA(idx_p)) {
    stop("Internal tau-grid mismatch.")
  }
  
  # ---------------------------------------------------------------------------
  # Step 3: fit quantile engine
  # ---------------------------------------------------------------------------
  dots <- list(...)
  rq_method <- dots$rq_method %||% "br"
  
  fit_rq_engine <- function() {
    if (!requireNamespace("quantreg", quietly = TRUE)) {
      stop("Package 'quantreg' is required for method = 'rq'.")
    }
    
    rq_dat <- as.data.frame(X_obs)
    rq_dat[[y_col]] <- y_obs
    
    fit_obj <- quantreg::rq(
      formula = stats::as.formula(paste(y_col, "~ .")),
      data = rq_dat,
      tau = taus_pm,
      method = rq_method
    )
    
    predict_Q <- function(x_new, taus_use = taus_pm) {
      x_new <- assert_numeric_x(x_new, x_cols, name = "x_new")
      pred <- as.matrix(stats::predict(fit_obj, newdata = as.data.frame(x_new)))
      
      idx <- match(taus_use, taus_pm)
      if (anyNA(idx)) {
        stop("taus_use must be a subset of the fitted tau grid.")
      }
      pred[, idx, drop = FALSE]
    }
    
    list(object = fit_obj, predict_Q = predict_Q)
  }
  
  fit_qrf_engine <- function() {
    if (!requireNamespace("quantregForest", quietly = TRUE)) {
      stop("Package 'quantregForest' is required for method = 'qrf'.")
    }
    
    dots_qrf <- dots
    dots_qrf$rq_method <- NULL
    
    qrf_fit <- do.call(
      quantregForest::quantregForest,
      c(list(x = as.data.frame(X_obs), y = y_obs), dots_qrf)
    )
    
    predict_Q <- function(x_new, taus_use = taus_pm) {
      x_new <- assert_numeric_x(x_new, x_cols, name = "x_new")
      Q <- as.matrix(stats::predict(
        qrf_fit,
        newdata = as.data.frame(x_new),
        what = taus_use
      ))
      if (ncol(Q) != length(taus_use)) {
        stop("Unexpected output dimension from quantile random forest prediction.")
      }
      Q
    }
    
    list(object = qrf_fit, predict_Q = predict_Q)
  }
  
  fit_grf_engine <- function() {
    if (!requireNamespace("grf", quietly = TRUE)) {
      stop("Package 'grf' is required for method = 'grf'.")
    }
    
    dots_grf <- dots
    dots_grf$rq_method <- NULL
    
    grf_fit <- do.call(
      grf::quantile_forest,
      c(list(X = X_obs, Y = y_obs), dots_grf)
    )
    
    predict_Q <- function(x_new, taus_use = taus_pm) {
      x_new <- assert_numeric_x(x_new, x_cols, name = "x_new")
      
      Q <- stats::predict(
        grf_fit,
        newdata = x_new,
        quantiles = taus_use
      )$predictions
      
      Q <- as.matrix(Q)
      
      if (ncol(Q) != length(taus_use)) {
        stop("Unexpected output dimension from grf quantile_forest prediction.")
      }
      
      Q
    }
    
    list(object = grf_fit, predict_Q = predict_Q)
  }
  
  engine <- switch(
    method,
    rq = fit_rq_engine(),
    qrf = fit_qrf_engine(),
    grf = fit_grf_engine()
  )
  
  # ---------------------------------------------------------------------------
  # Step 4: matrix density helper
  # ---------------------------------------------------------------------------
  dens_y_given_x_matrix <- function(
    y_val,
    y_given_x_hat,
    hi_y_given_x,
    lo_y_given_x,
    h_use,
    enforce_monotone = FALSE,
    tail_decay = TRUE,
    num_extra_points = 20L,
    decay_factor = 0.7
  ) {
    if (!identical(dim(y_given_x_hat), dim(hi_y_given_x)) ||
        !identical(dim(hi_y_given_x), dim(lo_y_given_x))) {
      stop("The quantile matrices do not have the same dimensions.")
    }
    
    if (length(h_use) != ncol(y_given_x_hat)) {
      stop("The length of h_use must match the number of columns in y_given_x_hat.")
    }
    
    n <- nrow(y_given_x_hat)
    f_val <- matrix(0, nrow = n, ncol = length(y_val))
    
    cross <- y_given_x_hat[, -1, drop = FALSE] -
      y_given_x_hat[, -ncol(y_given_x_hat), drop = FALSE] <= 0
    
    for (ii in seq_len(n)) {
      q0 <- y_given_x_hat[ii, ]
      qh <- hi_y_given_x[ii, ]
      ql <- lo_y_given_x[ii, ]
      h_ii <- h_use
      
      if (isTRUE(enforce_monotone)) {
        q0 <- as.numeric(stats::isoreg(q0)$yf)
        qh <- as.numeric(stats::isoreg(qh)$yf)
        ql <- as.numeric(stats::isoreg(ql)$yf)
      } else {
        if (any(cross[ii, ])) {
          ind <- which(cross[ii, ]) + 1L
          keep <- setdiff(seq_along(q0), ind)
          q0 <- q0[keep]
          qh <- qh[keep]
          ql <- ql[keep]
          h_ii <- h_ii[keep]
        }
      }
      
      if (length(q0) == 0L) {
        f_val[ii, ] <- NA_real_
        next
      }
      
      pdf_y_given_x <- (2 * h_ii) / (qh - ql - 1e-8)
      pdf_y_given_x[is.infinite(pdf_y_given_x)] <- NA_real_
      pdf_y_given_x[pdf_y_given_x < 1e-10] <- 1e-10
      pdf_y_given_x[pdf_y_given_x > 1] <- 1
      
      valid_idx <- !is.na(pdf_y_given_x)
      x <- q0[valid_idx]
      f_x <- pdf_y_given_x[valid_idx]
      
      if (length(x) == 0L) {
        f_val[ii, ] <- NA_real_
        next
      }
      
      ord <- order(x)
      x <- x[ord]
      f_x <- f_x[ord]
      
      keep <- !duplicated(x)
      x <- x[keep]
      f_x <- f_x[keep]
      
      if (length(x) == 1L) {
        x <- c(x, x + 0.1)
        f_x <- c(f_x, f_x)
      }
      
      x_all <- x
      f_all <- f_x
      
      if (isTRUE(tail_decay) && num_extra_points > 0L) {
        gap <- stats::quantile(diff(sort(x)), 0.5, names = FALSE, na.rm = TRUE)
        if (is.na(gap) || !is.finite(gap) || gap <= 0) {
          gap <- 2
        }
        
        extra_points_left <- seq(min(x) - num_extra_points * gap, min(x) - gap, by = gap)
        extra_points_right <- seq(max(x) + gap, max(x) + num_extra_points * gap, by = gap)
        extra_density_left <- f_x[1] * decay_factor ^ seq(num_extra_points, 1)
        extra_density_right <- f_x[length(f_x)] * decay_factor ^ seq(1, num_extra_points)
        
        x_all <- c(extra_points_left, x, extra_points_right)
        f_all <- c(extra_density_left, f_x, extra_density_right)
      }
      
      f_all[f_all < 1e-10] <- 1e-10
      f_all[f_all > 1] <- 1
      
      ord_all <- order(x_all)
      x_all <- x_all[ord_all]
      f_all <- f_all[ord_all]
      
      keep_all <- !duplicated(x_all)
      x_all <- x_all[keep_all]
      f_all <- f_all[keep_all]
      
      if (length(x_all) == 1L) {
        x_all <- c(x_all, x_all + 0.1)
        f_all <- c(f_all, f_all)
      }
      
      fy <- stats::approx(
        x = x_all,
        y = f_all,
        xout = y_val,
        method = "linear",
        rule = 2
      )$y
      
      fy <- pmin(pmax(fy, 1e-10), 1)
      
      if (any(y_val < min(x_all))) {
        fy[y_val < min(x_all)] <- f_all[1]
      }
      if (any(y_val > max(x_all))) {
        fy[y_val > max(x_all)] <- f_all[length(f_all)]
      }
      
      f_val[ii, ] <- fy
    }
    
    f_val
  }
  
  # ---------------------------------------------------------------------------
  # Step 5: prepare quantile matrices at new x
  # ---------------------------------------------------------------------------
  build_quantile_mats <- function(x_new) {
    x_new <- assert_numeric_x(x_new, x_cols, name = "x_new")
    Qall <- engine$predict_Q(x_new, taus_pm)
    
    list(
      Q0 = Qall[, idx_0, drop = FALSE],
      Qhi = Qall[, idx_p, drop = FALSE],
      Qlo = Qall[, idx_m, drop = FALSE]
    )
  }
  
  # ---------------------------------------------------------------------------
  # Step 6: matrix prediction on a full y-grid
  # ---------------------------------------------------------------------------
  predict_density_grid <- function(x_new, y_grid) {
    x_new <- assert_numeric_x(x_new, x_cols, name = "x_new")
    y_grid <- as.numeric(y_grid)
    
    if (anyNA(y_grid) || any(!is.finite(y_grid))) {
      stop("y_grid must contain only finite numeric values.")
    }
    
    qm <- build_quantile_mats(x_new)
    
    dens_y_given_x_matrix(
      y_val = y_grid,
      y_given_x_hat = qm$Q0,
      hi_y_given_x = qm$Qhi,
      lo_y_given_x = qm$Qlo,
      h_use = h0,
      enforce_monotone = enforce_monotone,
      tail_decay = tail_decay,
      num_extra_points = num_extra_points,
      decay_factor = decay_factor
    )
  }
  
  # ---------------------------------------------------------------------------
  # Step 7: paired prediction
  # ---------------------------------------------------------------------------
  predict_density <- function(x_new, y_new, rule = 2) {
    x_new <- assert_numeric_x(x_new, x_cols, name = "x_new")
    y_new <- as.numeric(y_new)
    
    if (any(!is.finite(y_new))) {
      stop("y_new must contain only finite numeric values.")
    }
    if (nrow(x_new) != length(y_new)) {
      stop("x_new and y_new must have the same number of rows/elements.")
    }
    
    row_keys <- apply(signif(x_new, 14), 1, paste, collapse = "\r")
    key_unique <- unique(row_keys)
    
    out <- numeric(length(y_new))
    
    for (key in key_unique) {
      idx <- which(row_keys == key)
      x_one <- x_new[idx[1], , drop = FALSE]
      y_sub <- y_new[idx]
      dens_sub <- predict_density_grid(x_one, y_sub)
      out[idx] <- as.numeric(dens_sub[1, ])
    }
    
    out
  }
  
  # ---------------------------------------------------------------------------
  # Step 8: user-facing quantile prediction
  # ---------------------------------------------------------------------------
  predict_Q_user <- function(x_new, taus_use = taus0) {
    x_new <- assert_numeric_x(x_new, x_cols, name = "x_new")
    engine$predict_Q(x_new, taus_use)
  }
  
  list(
    taus = taus0,
    fit_object = engine$object,
    predict_Q = predict_Q_user,
    predict_density = predict_density,
    predict_density_grid = predict_density_grid
  )
}