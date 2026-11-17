#' Maximum rank correlation fit with a smoothed objective
#'
#' Maximizes the smoothed area \code{\link{smooth_auc}} of the linear predictor
#' \code{X[,1] + X[,-1] \%*\% eta + Z \%*\% gamma} over \code{(eta, gamma)}.
#' The coefficient of the first established marker is fixed at one for
#' identifiability (Han 1987; Sherman 1993). Pass \code{Z = NULL} for the
#' reduced model.
#'
#' @param y binary outcome.
#' @param X numeric matrix (n by p) of established markers, p >= 1. The first
#'   column should be oriented so that larger values indicate a case.
#' @param Z numeric matrix (n by q) of new markers, or \code{NULL}.
#' @param h bandwidth.
#' @param start optional starting values of length p - 1 + q; by default the
#'   logistic estimates scaled to the normalization are supplied by
#'   \code{\link{aucnull_test}}, and a grid search is used when none are given.
#' @param reltol convergence tolerance passed to \code{optim}. With one free
#'   coefficient the maximization is a bounded scalar search (Brent) on an
#'   interval of half width two (on the scale of the anchor) around the start,
#'   widened while the maximum lies on the boundary; with several free
#'   coefficients a quasi Newton (BFGS) search is used.
#' @return a list with \code{theta} (named coefficient vector, excluding the
#'   fixed first coefficient), \code{auc} (maximized smoothed area), the
#'   objective function \code{f} (returns the smoothed area at a parameter
#'   vector) and \code{convergence}.
#' @export
mrc_fit <- function(y, X, Z = NULL, h, start = NULL, reltol = 1e-10) {
  X <- as.matrix(X); n <- nrow(X); p <- ncol(X)
  if (!is.null(Z)) Z <- as.matrix(Z)
  q <- if (is.null(Z)) 0L else ncol(Z)
  W <- cbind(if (p > 1) X[, -1, drop = FALSE], Z)   # free coefficients
  k <- ncol(W); if (is.null(k)) k <- 0L
  x1 <- X[, 1]
  f <- function(theta) {
    s <- if (k > 0) x1 + as.vector(W %*% theta) else x1
    smooth_auc(s, y, h)
  }
  if (k == 0L) {
    return(list(theta = numeric(0), auc = f(numeric(0)), f = f, convergence = 0L))
  }
  if (is.null(start)) {
    start <- rep(0, k)
    if (k == 1L) {
      grid <- seq(-2, 2, length.out = 41) * sd(x1) / (sd(W[, 1]) + 1e-12)
      vals <- vapply(grid, f, numeric(1))
      start <- grid[which.max(vals)]
    }
  }
  if (k == 1L) {
    # bounded scalar search on an interval of half width two around the start, widened while the maximum lies on the boundary
    half <- 2 * sd(x1) / (sd(W[, 1]) + 1e-12)
    repeat {
      o <- optim(start, function(th) -f(th), method = "Brent", lower = start - half, upper = start + half,
                 control = list(reltol = reltol, maxit = 2000))
      if (abs(o$par - start) < 0.999 * half || half > 1e3) break
      start <- o$par; half <- 2 * half
    }
  } else {
    o <- optim(start, function(th) -f(th), method = "BFGS", control = list(reltol = reltol, maxit = 2000))
  }
  theta <- o$par
  nm <- c(if (p > 1) paste0("eta", seq_len(p - 1)), if (q > 0) paste0("gamma", seq_len(q)))
  names(theta) <- nm
  list(theta = theta, auc = -o$value, f = f, convergence = o$convergence)
}

#' Numerical Hessian of the smoothed area at a parameter vector
#' @keywords internal
num_hessian <- function(f, theta, eps = 0.05) {
  k <- length(theta); H <- matrix(0, k, k)
  for (i in seq_len(k)) for (j in seq_len(k)) {
    ei <- ej <- rep(0, k); ei[i] <- eps; ej[j] <- eps
    H[i, j] <- (f(theta + ei + ej) - f(theta + ei - ej) - f(theta - ei + ej) + f(theta - ei - ej)) / (4 * eps * eps)
  }
  (H + t(H)) / 2
}
