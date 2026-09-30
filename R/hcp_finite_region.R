# ============================================================
# hcp_finite_region.R
# Event-based representation of HCP sets on a fixed finite domain
#
# Used for Table 1, Table S.1, Table S.3 and Figure S.10.
# Dependencies: hcp_region.R and its density/propensity helpers.
# Uses the supplied domain endpoints, fitted scores, weights, ranks and CCT rule.
# ============================================================

hcp_combine_scores <- function(p_arr, method) {
  p_arr <- pmin(pmax(p_arr, 1e-15), 1 - 1e-15)
  t_arr <- tan((0.5 - p_arr) * pi)
  tbar <- apply(t_arr, c(2, 3), mean, na.rm = TRUE)
  bad <- !is.finite(tbar)
  if (any(bad)) {
    tbar[is.nan(tbar)] <- 0
    inf <- is.infinite(tbar)
    tbar[inf] <- sign(tbar[inf]) * .Machine$double.xmax
  }
  p <- 0.5 - atan(tbar) / pi
  pmin(pmax(p, 0), 1)
}

hcp_capture_interpolants <- function(fit, X, grid) {
  saved <- vector('list', nrow(X));count <- 0L
  pe <- environment(fit$predict_density_grid)
  helper <- get('dens_y_given_x_matrix', pe)
  en <- new.env(parent = environment(helper))
  en$.save_knots <- function(i, x, f) saved[[i]] <<- list(x = x, f = f)
  rewrite <- function(x) {
    if (!is.call(x)) return(x)
    if (identical(x[[1]], as.name('{'))) {
      nodes <- list(as.name('{'))
      for (a in as.list(x)[-1L]) {
        nodes[[length(nodes) + 1L]] <- rewrite(a)
        if (is.call(a) && identical(a[[1]], as.name('<-')) && identical(a[[2]], quote(f_val[ii, ])) && identical(a[[3]], as.name('fy'))) {
          nodes[[length(nodes) + 1L]] <- quote(.save_knots(ii, x_all, f_all));count <<- count + 1L
        }
      }
      return(as.call(nodes))
    }
    for (i in seq_along(x)[-1L]) if (!identical(x[[i]], quote(expr = ))) x[i] <- list(rewrite(x[[i]]))
    x
  }
  body(helper) <- rewrite(body(helper));environment(helper) <- en;stopifnot(count == 1L)
  pred <- fit$predict_density_grid;pen <- new.env(parent = pe);pen$dens_y_given_x_matrix <- helper;environment(pred) <- pen
  check <- pred(X, grid)
  stopifnot(all(vapply(saved, function(z) !is.null(z) && length(z$x) >= 2L && all(diff(z$x) > 0), logical(1))))
  rebuilt <- t(vapply(saved, function(z) pmin(pmax(approx(z$x, z$f, xout = grid, rule = 2)$y, 1e-10), 1), numeric(length(grid))))
  stopifnot(identical(unname(check), unname(rebuilt)))
  list(knots = saved, prediction = check)
}

hcp_p_from_scores <- function(fits, scores) {
  S <- length(fits);K <- nrow(scores[[1]]);Ny <- ncol(scores[[1]]);B <- length(fits[[1]]$cal)
  ps <- array(NA_real_, c(S, K, Ny))
  for (s in seq_len(S)) {
    pb <- array(NA_real_, c(B, K, Ny));ft <- fits[[s]]
    for (b in seq_len(B)) {
      cb <- ft$cal[[b]]
      if (isTRUE(cb$empty)) next
      for (k in seq_len(K)) {
        j <- findInterval(scores[[s]][k, ], cb$R_sorted, left.open = TRUE) + 1L
        ge <- ifelse(j <= length(cb$W_ge), cb$W_ge[j], 0)
        pb[b, k, ] <- (ge + ft$weight[k]) / (cb$sum_w + ft$weight[k])
      }
    }
    ps[s, , ] <- hcp_combine_scores(pb, 'cct')
  }
  hcp_combine_scores(ps, 'cct')
}
hcp_cache_p_grid <- function(fits, y) {
  scores <- lapply(fits, function(ft) -t(vapply(ft$knots, function(z)
    pmin(pmax(approx(z$x, z$f, xout = y, rule = 2)$y, 1e-10), 1), numeric(length(y)))))
  hcp_p_from_scores(fits, scores)
}
hcp_cache_p_point <- function(fits, k, y) {
  fs <- lapply(fits, function(ft) {
    ft$weight <- ft$weight[k];ft$knots <- ft$knots[k];ft
  })
  as.numeric(hcp_cache_p_grid(fs, y))
}

hcp_union_length <- function(x) if (nrow(x)) sum(x[, 2] - x[, 1]) else 0
hcp_union_member <- function(x, y) vapply(y, function(v) any(v >= x[, 1] & v <= x[, 2]), logical(1))
hcp_numerical_control <- list(
  alpha = .1, y_relative = 1e-8, total_band_relative = 1e-6,
  certificate_margin = 1e-12, max_refine = 64L, max_expand = 64L, max_events = 100000L, max_p_evaluations = 1000000L
)

hcp_merge_bands <- function(a) {
  if (!length(a) || !nrow(a)) return(matrix(numeric(), ncol = 2))
  a <- a[order(a[, 1], a[, 2]), , drop = FALSE];out <- vector('list', nrow(a));j <- 1L;out[[j]] <- a[1, ]
  if (nrow(a) > 1L) for (i in 2:nrow(a)) if (a[i, 1] <= out[[j]][2]) out[[j]][2] <- max(out[[j]][2], a[i, 2]) else {
    j <- j + 1L;out[[j]] <- a[i, ]
  }
  do.call(rbind, out[seq_len(j)])
}
hcp_point_fits <- function(fits, k) lapply(fits, function(ft) {
  ft$knots <- ft$knots[k];ft$weight <- ft$weight[k];ft
})
hcp_density_at <- function(ft, y) pmin(pmax(approx(ft$knots[[1]]$x, ft$knots[[1]]$f, xout = y, rule = 2)$y, 1e-10), 1)
hcp_set_membership <- function(set, y) {
  a <- set$intervals;hit <- rep(FALSE, length(y))
  if (nrow(a)) for (i in seq_len(nrow(a))) hit <- hit | ((y > a$lower[i] | (a$lower_closed[i] & y == a$lower[i])) &
    (y < a$upper[i] | (a$upper_closed[i] & y == a$upper[i])))
  uncertain <- hcp_union_member(set$bands, y)
  list(member = hit, uncertain = uncertain)
}
hcp_extract_atoms <- function(events, cells, points) {
  m <- length(events);on <- logical(2 * m + 1L);on[seq(1, 2 * m + 1L, 2)] <- cells;on[seq(2, 2 * m, 2)] <- points
  lower <- upper <- rep(NA_real_, length(on));lc <- uc <- rep(FALSE, length(on))
  lower[seq(1, 2 * m + 1L, 2)] <- c(-Inf, events);upper[seq(1, 2 * m + 1L, 2)] <- c(events, Inf)
  lower[seq(2, 2 * m, 2)] <- upper[seq(2, 2 * m, 2)] <- events;lc[seq(2, 2 * m, 2)] <- uc[seq(2, 2 * m, 2)] <- TRUE
  st <- which(on & !c(FALSE, head(on, -1)));en <- which(on & !c(tail(on, -1), FALSE))
  data.frame(lower = lower[st], upper = upper[en], lower_closed = lc[st], upper_closed = uc[en])
}

# Find a verified density-threshold bracket, including any floating-point
# equality plateau. This never interpolates or smooths the conformal p-value.
hcp_bracket_event <- function(meta, ft, y, tol, cfg) {
  x0 <- meta$x0;x1 <- meta$x1;q <- meta$q;sgn <- sign(meta$f1 - meta$f0)
  fun <- function(x) sgn * (hcp_density_at(ft, x) - q)
  radius <- 64 * .Machine$double.eps * max(1, abs(y), abs(x0), abs(x1));calls <- 0L
  for (it in 0:cfg$max_refine) {
    lo <- max(x0, y - radius);hi <- min(x1, y + radius);v <- fun(c(lo, hi));calls <- calls + 2L
    valid <- (v[1] < 0 || (lo == x0 && v[1] == 0)) && (v[2] > 0 || (hi == x1 && v[2] == 0))
    if (valid) break
    radius <- radius * 2
  }
  if (!valid) return(list(lo = x0, hi = x1, ok = FALSE, calls = calls, iterations = it))
  iterations <- it
  if (hi - lo > tol) {
    # Lower and upper limits of the density==q plateau, each to tol/4.
    al <- bl <- lo;ah <- bh <- hi
    for (j in seq_len(cfg$max_refine)) {
      if (ah - al <= tol / 4) break
      mid <- al + (ah - al) / 2;if (mid <= al || mid >= ah) break
      v <- fun(mid);calls <- calls + 1L;if (v >= 0) ah <- mid else al <- mid
    }
    for (jj in seq_len(cfg$max_refine)) {
      if (bh - bl <= tol / 4) break
      mid <- bl + (bh - bl) / 2;if (mid <= bl || mid >= bh) break
      v <- fun(mid);calls <- calls + 1L;if (v > 0) bh <- mid else bl <- mid
    }
    lo <- al;hi <- bh;iterations <- max(iterations, j, jj)
  }
  lo <- min(lo, y);hi <- max(hi, y)
  list(lo = lo, hi = hi, ok = hi - lo <= tol, calls = calls, iterations = iterations)
}

hcp_construct_events <- function(fits, k, original_grid, cfg = hcp_numerical_control) {
  started <- proc.time()[3];ft <- hcp_point_fits(fits, k);L <- min(original_grid);U <- max(original_grid);W <- U - L
  scale <- max(1, W);ytol <- cfg$y_relative * scale;btol <- cfg$total_band_relative * scale
  count <- 0L;density_queries <- 0L
  pfun <- function(y) {
    stopifnot(all(is.finite(y)), all(y >= L), all(y <= U))
    count <<- count + length(y)
    if (count > cfg$max_p_evaluations) stop('P_EVALUATION_BUDGET')
    as.numeric(hcp_cache_p_grid(ft, y))
  }
  minscore <- lapply(ft, function(z) matrix(-min(z$knots[[1]]$f), 1, 1))
  lower_bound <- as.numeric(hcp_p_from_scores(ft, minscore))
  proof <- lower_bound > cfg$alpha + cfg$certificate_margin
  empty_intervals <- data.frame(lower = numeric(), upper = numeric(), lower_closed = logical(), upper_closed = logical())
  if (proof) {
    ints <- data.frame(lower = -Inf, upper = Inf, lower_closed = FALSE, upper_closed = FALSE)
    return(list(
      intervals = ints, bands = matrix(numeric(), ncol = 2), boundaries = data.frame(), length = Inf, components = 1L,
      status = 'resolved', precision_pass = TRUE, component_certified = TRUE, whole_line = TRUE, whole_line_certified = TRUE,
      lower_bound = lower_bound, left_unbounded = TRUE, right_unbounded = TRUE, domain_c = NA_real_,
      events = 0L, p_evaluations = 1L, density_queries = 0L, uncertainty_length = 0, max_boundary_width = 0,
      elapsed = unname(proc.time()[3] - started), method = 'global_lower_bound', issues = character()
    ))
  }
  knotmin <- min(vapply(ft, function(z) min(z$knots[[1]]$x), numeric(1)));knotmax <- max(vapply(ft, function(z) max(z$knots[[1]]$x), numeric(1)))
  stopifnot(all(is.finite(c(knotmin, knotmax, L, U))), W > 0)
  cc <- 0
  dl <- L;du <- U
  roots <- list();flatbands <- list();ri <- fi <- 0L
  for (s in seq_along(ft)) {
    z <- ft[[s]]$knots[[1]];threshold <- sort(unique(-unlist(lapply(ft[[s]]$cal, `[[`, 'R_sorted'))))
    for (j in seq_len(length(z$x) - 1L)) {
      f0 <- z$f[j];f1 <- z$f[j + 1L];if (f0 == f1) next
      lo <- min(f0, f1);hi <- max(f0, f1)
      a <- findInterval(lo, threshold, left.open = TRUE) + 1L;b <- findInterval(hi, threshold)
      if (a > b) next
      qq <- threshold[a:b];x0 <- z$x[j];x1 <- z$x[j + 1L]
      if (abs(f1 - f0) <= 1024 * .Machine$double.eps * max(abs(f0), abs(f1), 1e-10)) {
        fi <- fi + 1L;flatbands[[fi]] <- c(x0, x1)
      }
      yy <- x0 + (qq - f0) / (f1 - f0) * (x1 - x0);yy[qq == f0] <- x0;yy[qq == f1] <- x1
      yy <- pmax(x0, pmin(x1, yy))
      ri <- ri + 1L;roots[[ri]] <- data.frame(y = yy, s = s, x0 = x0, x1 = x1, f0 = f0, f1 = f1, q = qq)
    }
  }
  roots <- if (length(roots)) do.call(rbind, roots) else data.frame(y = numeric(), s = integer(), x0 = numeric(), x1 = numeric(), f0 = numeric(), f1 = numeric(), q = numeric())
  ev <- sort(unique(c(dl, du, unlist(lapply(ft, function(z) z$knots[[1]]$x)), roots$y)))
  ev <- ev[ev >= L & ev <= U]
  flatbands <- lapply(Filter(function(a) a[2] >= L && a[1] <= U, flatbands), function(a) c(max(a[1], L), min(a[2], U)))
  if (length(ev) > cfg$max_events) stop('EVENT_BUDGET')
  mid <- head(ev, -1) + diff(ev) / 2;unrepresentable <- which(mid <= head(ev, -1) | mid >= tail(ev, -1))
  pp <- pfun(ev);pm <- pfun(mid);ptail <- pfun(c(dl, du))
  points <- pp > cfg$alpha;cells <- c(ptail[1] > .1, pm > cfg$alpha, ptail[2] > .1)
  active <- which(points != head(cells, -1) | points != tail(cells, -1) | head(cells, -1) != tail(cells, -1))
  boundary_tol <- min(ytol, btol / (2 * max(1, length(active))))
  root_index <- match(roots$y, ev);groups <- split(seq_len(nrow(roots)), root_index)
  bands <- boundaries <- list();bi <- 0L;ok <- TRUE;issues <- character()
  for (j in active) {
    y <- ev[j];rr <- groups[[as.character(j)]];lo <- hi <- y;iterations <- 0L;certified <- TRUE
    if (length(rr)) for (ir in rr) {
      br <- hcp_bracket_event(roots[ir, ], ft[[roots$s[ir]]], y, boundary_tol, cfg)
      lo <- min(lo, br$lo);hi <- max(hi, br$hi);certified <- certified && br$ok
      density_queries <- density_queries + br$calls;iterations <- max(iterations, br$iterations)
    }
    rad <- 64 * .Machine$double.eps * max(1, abs(y), scale)
    lo <- min(lo, y - rad);hi <- max(hi, y + rad)
    bi <- bi + 1L;bands[[bi]] <- c(lo, hi)
    boundaries[[bi]] <- data.frame(
      event_index = j, y = y, lo = lo, hi = hi, width = hi - lo, point_accepted = points[j],
      left_accepted = cells[j], right_accepted = cells[j + 1L], density_bracket_verified = certified, iterations = iterations
    )
    ok <- ok && certified && hi - lo <= ytol
  }
  # Close but distinct events cannot be silently merged. Keep a bounded
  # uncertainty region around cells without a representable midpoint.
  if (length(unrepresentable)) {
    for (j in unrepresentable) {
      bi <- bi + 1L;bands[[bi]] <- c(ev[j], ev[j + 1L])
    }
    issues <- c(issues, 'unrepresentable_event_cell')
  }
  if (length(flatbands)) {
    for (v in flatbands) {
      bi <- bi + 1L;bands[[bi]] <- v
    }
    issues <- c(issues, 'numerically_flat_threshold_segment');ok <- FALSE
  }
  bm <- if (length(bands)) hcp_merge_bands(do.call(rbind, bands)) else matrix(numeric(), ncol = 2)
  bound <- if (length(boundaries)) do.call(rbind, boundaries) else data.frame()
  overlapping <- if (nrow(bound) > 1L) any(head(bound$hi, -1) >= tail(bound$lo, -1)) else FALSE
  if (overlapping) issues <- c(issues, 'overlapping_boundary_brackets')
  uncertainty <- hcp_union_length(bm);ok <- ok && uncertainty <= btol
  if (!ok) issues <- c(issues, 'accuracy_budget_not_met')
  intervals <- hcp_extract_atoms(ev, cells, points)
  component_certified <- !length(unrepresentable) && !length(flatbands) && !overlapping
  whole <- length(cells) > 0 && all(cells) && all(points)
  list(
    intervals = intervals, bands = bm, boundaries = bound, length = if (nrow(intervals)) sum(intervals$upper - intervals$lower) else 0,
    components = nrow(intervals), status = if (ok && component_certified) 'resolved' else if (ok) 'accuracy_resolved_topology_uncertain' else 'unresolved',
    precision_pass = ok, component_certified = component_certified, whole_line = whole, whole_line_certified = whole && ok && component_certified,
    lower_bound = lower_bound, left_unbounded = cells[1], right_unbounded = tail(cells, 1), domain_c = cc, domain = c(dl, du),
    events = length(ev), p_evaluations = count, density_queries = density_queries, uncertainty_length = uncertainty,
    max_boundary_width = if (nrow(bound)) max(bound$width) else 0, elapsed = unname(proc.time()[3] - started), method = 'score_events', issues = issues
  )
}

hcp_clip_set <- function(s, L, U) {
  s$raw_status <- s$status;s$raw_components <- s$components;s$raw_unrestricted_length <- s$length
  a <- s$intervals
  if (nrow(a)) {
    a$lower_closed <- ifelse(a$lower < L, TRUE, a$lower_closed);a$upper_closed <- ifelse(a$upper > U, TRUE, a$upper_closed)
    a$lower <- pmax(a$lower, L);a$upper <- pmin(a$upper, U)
    a <- a[a$lower < a$upper | (a$lower == a$upper & a$lower_closed & a$upper_closed), , drop = FALSE];rownames(a) <- NULL
  }
  s$intervals <- a;s$components <- nrow(a);s$length <- if (nrow(a)) sum(a$upper - a$lower) else 0
  if (nrow(s$bands)) {
    b <- s$bands;b[, 1] <- pmax(b[, 1], L);b[, 2] <- pmin(b[, 2], U);s$bands <- hcp_merge_bands(b[b[, 1] <= b[, 2], , drop = FALSE])
  }
  if (nrow(s$boundaries)) s$boundaries <- s$boundaries[s$boundaries$y >= L & s$boundaries$y <= U, , drop = FALSE]
  s$uncertainty_length <- hcp_union_length(s$bands);s$domain <- c(L, U);s$domain_c <- 0;s$left_unbounded <- s$right_unbounded <- FALSE
  s$whole_line <- s$whole_line_certified <- FALSE;s$entire_finite_domain <- nrow(a) == 1L && a$lower == L && a$upper == U && a$lower_closed && a$upper_closed
  stopifnot(is.finite(s$length), s$length <= U - L + 1e-10, all(a$lower >= L), all(a$upper <= U));s
}
hcp_construct_finite <- function(fits, k, grid) hcp_clip_set(hcp_construct_events(fits, k, grid), min(grid), max(grid))
#' Represent the fitted HCP acceptance rule on the supplied finite domain
#'
#' Calls hcp_region once. Observation callbacks retain fitted density knots,
#' calibration scores and weights; they neither refit nor consume randomness.
#' The event constructor keeps every connected component with strict p > alpha.
#' @param dat Training/calibration subjects in long format.
#' @param x_test Covariates of one held-out subject.
#' @param y_grid Candidate grid whose endpoints define the finite domain.
#' @param alpha Miscoverage level (the manuscript uses 0.1).
#' @param ... Statistical arguments passed unchanged to hcp_region.
#' @return HCP output plus finite_sets, total_length, and score_cache.
hcp_finite_region <- function(dat, x_test, y_grid, alpha = .1, ...) {
  if (!identical(as.numeric(alpha), .1)) stop("This numerical implementation is validated for alpha = 0.1.")
  X <- as.matrix(x_test); fits <- list()
  observer <- function(stage, e) {
    s <- get("s", e)
    if (stage == "test") {
      captured <- hcp_capture_interpolants(get("dens_fit", e), X, y_grid)
      stopifnot(identical(unname(captured$prediction), unname(-get("R_test_mat", e))))
      fits[[s]] <<- list(knots = captured$knots, weight = get("w_test_vec", e))
    } else if (stage == "split") {
      fits[[s]]$cal <<- get("sub_cache", e)
    }
  }
  ans <- hcp_region(
    dat = dat, x_test = X, y_grid = y_grid, alpha = alpha,
    score_observer = observer, ...
  )
  stopifnot(identical(unname(hcp_cache_p_grid(fits, y_grid)), unname(ans$p_final)))
  ans$finite_sets <- lapply(seq_len(nrow(X)), function(k) hcp_construct_finite(fits, k, y_grid))
  if (any(!vapply(ans$finite_sets, `[[`, logical(1), "precision_pass")))
    stop("Finite-domain set construction exceeded the numerical accuracy budget.")
  ans$total_length <- vapply(ans$finite_sets, `[[`, numeric(1), "length")
  ans$score_cache <- fits
  ans
}

#' Evaluate membership and total Lebesgue measure of finite-domain HCP sets
evaluate_hcp_finite <- function(res, y_true) {
  stopifnot(length(res$finite_sets) == length(y_true))
  cbind(
    covered = vapply(seq_along(y_true), function(k)
      as.numeric(hcp_set_membership(res$finite_sets[[k]], y_true[k])$member), numeric(1)),
    length = res$total_length
  )
}
