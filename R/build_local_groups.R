#' Build localized groups for calibration and test covariates
#'
#' @description
#' This function constructs localized groups for calibration and test covariate
#' values using a profile-based representation of the fitted conditional density.
#'
#' The partition is learned from the training covariates. For each covariate
#' vector \eqn{x}, the fitted conditional density \eqn{\widehat\pi(y \mid x)} is
#' evaluated on a grid of \eqn{y}-values, and a profile representation is built
#' from the conditional distribution of \eqn{\widehat\pi(Y \mid x)}. K-means is
#' then applied to the training profiles, and the calibration and test profiles
#' are assigned to the nearest cluster centers.
#'
#' @param X_train Training covariate matrix or data.frame.
#' @param X_cal Calibration covariate matrix or data.frame.
#' @param X_test Test covariate matrix or data.frame.
#' @param dens_fit A fitted conditional density object returned by
#'   \code{\link{fit_cond_density_qp}}.
#' @param y_grid Numeric grid of response values used to evaluate the fitted
#'   conditional density.
#' @param n_local_groups Number of localized groups (clusters).
#' @param t_num Number of grid points used in the profile representation.
#' @param seed Optional random seed for k-means.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{group_cal}}{Cluster labels for the calibration covariates.}
#'   \item{\code{group_test}}{Cluster labels for the test covariates.}
#'   \item{\code{kmeans_fit}}{The fitted k-means object on the training profiles.}
#' }
#'
#' @export
build_local_groups <- function(
    X_train,
    X_cal,
    X_test,
    dens_fit,
    y_grid,
    n_local_groups = 5L,
    t_num = 200L,
    seed = NULL
) {
  
  # ---------------------------------------------------------------------------
  # Step 0: Validate and standardize inputs
  # ---------------------------------------------------------------------------
  to_numeric_matrix <- function(x, name) {
    if (is.data.frame(x)) {
      x <- as.matrix(x)
    } else {
      x <- as.matrix(x)
    }
    
    storage.mode(x) <- "numeric"
    
    if (anyNA(x) || any(!is.finite(x))) {
      stop(name, " must contain only finite numeric values.")
    }
    
    x
  }
  
  X_train <- to_numeric_matrix(X_train, "X_train")
  X_cal   <- to_numeric_matrix(X_cal,   "X_cal")
  X_test  <- to_numeric_matrix(X_test,  "X_test")
  
  if (ncol(X_train) != ncol(X_cal) || ncol(X_train) != ncol(X_test)) {
    stop("X_train, X_cal, and X_test must have the same number of columns.")
  }
  
  y_grid <- sort(unique(as.numeric(y_grid)))
  if (length(y_grid) < 2L || anyNA(y_grid) || any(!is.finite(y_grid))) {
    stop("y_grid must contain at least two distinct finite numeric values.")
  }
  
  n_local_groups <- as.integer(n_local_groups)
  if (!is.finite(n_local_groups) || n_local_groups < 1L) {
    stop("n_local_groups must be a positive integer.")
  }
  
  t_num <- as.integer(t_num)
  if (!is.finite(t_num) || t_num < 10L) {
    stop("t_num must be an integer greater than or equal to 10.")
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n_train <- nrow(X_train)
  n_cal   <- nrow(X_cal)
  n_test  <- nrow(X_test)
  
  if (n_local_groups > n_train) {
    stop("n_local_groups cannot exceed the number of training points.")
  }
  
  # If only one group is requested, return a trivial partition.
  if (n_local_groups == 1L) {
    return(list(
      group_cal = rep(1L, n_cal),
      group_test = rep(1L, n_test),
      kmeans_fit = NULL
    ))
  }
  
  # ---------------------------------------------------------------------------
  # Helper: Evaluate fitted densities on a common y-grid
  # ---------------------------------------------------------------------------
  density_on_grid <- function(X_mat, y_grid, dens_fit) {
    nx <- nrow(X_mat)
    ny <- length(y_grid)
    
    X_big <- X_mat[rep(seq_len(nx), each = ny), , drop = FALSE]
    y_big <- rep(y_grid, times = nx)
    
    X_big_df <- as.data.frame(X_big)
    
    f_big <- as.numeric(dens_fit$predict_density(X_big_df, y_big))
    if (anyNA(f_big) || any(!is.finite(f_big))) {
      stop("dens_fit$predict_density returned NA or non-finite values.")
    }
    
    matrix(f_big, nrow = nx, ncol = ny, byrow = TRUE)
  }
  
  # ---------------------------------------------------------------------------
  # Helper: Construct profile representation from a fitted density matrix
  # ---------------------------------------------------------------------------
  # For each x, we approximate
  #   G_{hat pi}(t | x) = P{ hat pi(Y|x) <= t | x }
  # using the density evaluated on y_grid.
  profile_from_density <- function(dens_mat, y_grid, t_grid) {
    nx <- nrow(dens_mat)
    ny <- ncol(dens_mat)
    nt <- length(t_grid)
    
    if (ny != length(y_grid)) {
      stop("dens_mat must have length(y_grid) columns.")
    }
    
    # Build approximate cell widths on the y-grid.
    y_mid_left  <- c(y_grid[1], (y_grid[-length(y_grid)] + y_grid[-1]) / 2)
    y_mid_right <- c((y_grid[-length(y_grid)] + y_grid[-1]) / 2, y_grid[length(y_grid)])
    cell_width  <- y_mid_right - y_mid_left
    
    prof <- matrix(0, nrow = nx, ncol = nt)
    
    for (i in seq_len(nx)) {
      f_row <- dens_mat[i, ]
      f_row[!is.finite(f_row)] <- 0
      f_row[f_row < 0] <- 0
      
      # Approximate probability mass associated with each y-grid point.
      mass_row <- f_row * cell_width
      total_mass <- sum(mass_row, na.rm = TRUE)
      
      if (!is.finite(total_mass) || total_mass <= 0) {
        stop("Degenerate density profile encountered during profile construction.")
      }
      
      mass_row <- mass_row / total_mass
      
      ord <- order(f_row)
      f_sorted <- f_row[ord]
      mass_sorted <- mass_row[ord]
      cdf_sorted <- cumsum(mass_sorted)
      
      idx <- findInterval(t_grid, f_sorted, left.open = FALSE)
      g <- numeric(nt)
      pos <- idx > 0
      g[pos] <- cdf_sorted[idx[pos]]
      g[!pos] <- 0
      
      prof[i, ] <- g
    }
    
    prof
  }
  
  # ---------------------------------------------------------------------------
  # Helper: Assign new profiles to nearest k-means centers
  # ---------------------------------------------------------------------------
  assign_to_centers <- function(profile_new, centers) {
    n_new <- nrow(profile_new)
    n_centers <- nrow(centers)
    out <- integer(n_new)
    
    for (i in seq_len(n_new)) {
      d2 <- rowSums(
        (centers - matrix(profile_new[i, ], nrow = n_centers,
                          ncol = ncol(centers), byrow = TRUE))^2
      )
      out[i] <- which.min(d2)
    }
    
    out
  }
  
  # ---------------------------------------------------------------------------
  # Step 1: Evaluate fitted densities for training, calibration, and test points
  # ---------------------------------------------------------------------------
  dens_train <- density_on_grid(X_train, y_grid, dens_fit)
  dens_cal   <- density_on_grid(X_cal,   y_grid, dens_fit)
  dens_test  <- density_on_grid(X_test,  y_grid, dens_fit)
  
  # ---------------------------------------------------------------------------
  # Step 2: Construct a common t-grid for the profile representation
  # ---------------------------------------------------------------------------
  max_dens <- max(c(dens_train, dens_cal, dens_test), na.rm = TRUE)
  if (!is.finite(max_dens) || max_dens <= 0) {
    stop("Unable to construct profile grid because the fitted densities are degenerate.")
  }
  
  t_grid <- seq(0, max_dens, length.out = t_num)
  
  # ---------------------------------------------------------------------------
  # Step 3: Build profile representations
  # ---------------------------------------------------------------------------
  profile_train <- profile_from_density(dens_train, y_grid, t_grid)
  profile_cal   <- profile_from_density(dens_cal,   y_grid, t_grid)
  profile_test  <- profile_from_density(dens_test,  y_grid, t_grid)
  
  # ---------------------------------------------------------------------------
  # Step 4: Fit k-means on the training profiles
  # ---------------------------------------------------------------------------
  kmeans_fit <- stats::kmeans(profile_train, centers = n_local_groups)
  
  # ---------------------------------------------------------------------------
  # Step 5: Assign calibration and test profiles to training-derived clusters
  # ---------------------------------------------------------------------------
  group_cal  <- assign_to_centers(profile_cal,  kmeans_fit$centers)
  group_test <- assign_to_centers(profile_test, kmeans_fit$centers)
  
  # ---------------------------------------------------------------------------
  # Return localized grouping information
  # ---------------------------------------------------------------------------
  list(
    group_cal = group_cal,
    group_test = group_test,
    kmeans_fit = kmeans_fit
  )
}