# ============================================================
# imbalanced_covariates.R
# Table S.1: cluster sizes 5 or 50 with identical covariate support
# X1 is {0.2,0.4,0.6,0.8,1}, once or ten times per subject.
# X2--X4 and every outcome/missingness draw use simulation_dgp.R.
# ============================================================

generate_imbalanced_covariates <- function(m, x_mean, x_sigma) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required.")
  }

  d <- length(x_mean)

  # Five time values, each repeated m/5 times, in random order.
  if (!m %in% c(5L, 50L)) stop("This experiment permits only M=5 or M=50.")
  x_time <- if (m == 5L) sample(seq_len(m)) / m else
    rep(seq_len(5L) / 5L, each=10L)[sample.int(50L)]

  if (d == 1L) {
    return(matrix(x_time, ncol = 1))
  }

  X_normal <- MASS::mvrnorm(
    n = m,
    mu = x_mean[-1],
    Sigma = x_sigma[-1, -1, drop = FALSE]
  )

  X <- cbind(x_time, X_normal)
  colnames(X) <- paste0("X", seq_len(ncol(X)))
  X
}
