#' Heller test of the null hypothesis of equal AUCs for new markers in nested models
#'
#' Tests whether adding the new markers \code{Z} to the established markers
#' \code{X} improves the area under the ROC curve of the best linear predictor.
#' The statistic is \code{2 n (A_full - A_reduced)} with both areas obtained
#' by maximum rank correlation on a kernel smoothed objective. Under the null
#' hypothesis that the best linear ordering does not use \code{Z}, the
#' statistic converges to a weighted sum of independent chi-square variables
#' with one degree of freedom (Heller et al. 2017, Theorem 1); the weights are
#' the eigenvalues of \code{-V_gamma \%*\% solve(D^{gamma gamma})}, where
#' \code{V_gamma} is the asymptotic covariance of the rank correlation estimate
#' of the coefficients of \code{Z} and \code{D} is the Hessian of the smoothed
#' area. Because the weights do not involve the link, the test keeps its asymptotic level for
#' any strictly increasing link, which is not true of the likelihood ratio test
#' when the logistic link is misspecified.
#'
#' Two estimators of \code{V_gamma} are offered. \code{"boot"} (recommended)
#' uses \code{n} times the bootstrap covariance of the rank correlation
#' estimate. \code{"plugin"} uses the sandwich form \code{solve(D) Sigma_g solve(D)}
#' with \code{Sigma_g} estimated from the per observation gradients of the
#' Hajek projection of the smoothed area, computed at the bandwidth
#' \code{sd(X[,1]) n^(-1/3)} regardless of \code{expo}.
#'
#' Clustered data (for example several lesions per patient) are handled through
#' \code{cluster}: the bootstrap then resamples clusters, the plugin covariance
#' sums the projections within clusters, and the comparison tests use a cluster
#' robust covariance. The limiting distribution is unchanged; the asymptotics are
#' in the number of clusters, so the weights are poorly estimated when there are
#' few clusters. Rows with missing values in \code{y}, \code{X}, \code{Z} or
#' \code{cluster} are dropped with a message.
#'
#' @param y binary outcome (0/1 or logical).
#' @param X established markers: a numeric vector or matrix (n by p). The first
#'   column is the anchor whose coefficient is fixed at one; it must be
#'   continuous. Its orientation is chosen by fitting the reduced model with the
#'   anchor and with its negative and keeping the orientation with the larger
#'   smoothed area (for p = 1 this amounts to an empirical AUC of at least 0.5).
#'   Factors must be converted to numeric columns beforehand.
#' @param Z new markers: a numeric vector or matrix (n by q).
#' @param method \code{"boot"} or \code{"plugin"} for the covariance of the
#'   rank correlation estimate.
#' @param B number of bootstrap replicates for \code{method = "boot"}.
#' @param expo bandwidth exponent for the statistic and the curvature; the
#'   bandwidth is \code{sd(X[,1]) * n^(-expo)}.
#' @param n_sim number of Monte Carlo draws used to evaluate the weighted
#'   chi-square distribution when \code{q > 1}.
#' @param seed optional seed for the bootstrap; the random number generator state
#'   is restored on exit.
#' @param compare logical; if \code{TRUE}, the likelihood ratio test of the
#'   logistic model, the Wald test with a robust (cluster robust when
#'   \code{cluster} is given) covariance, the Wald test of \code{Z} after a
#'   restricted cubic spline in each established marker with the same robust
#'   covariance, and the DeLong test on the fitted linear predictors are also
#'   returned for comparison. With \code{cluster} the likelihood ratio and
#'   DeLong tests, which assume independent observations, are not returned.
#' @param cluster optional vector of cluster identifiers of length n.
#' @param hessian_threshold identifiability threshold on the standardized scale:
#'   if the smallest eigenvalue of the curvature of the smoothed area in the
#'   coefficients of \code{Z} after profiling the other coefficients,
#'   \code{-solve(solve(D)[gamma, gamma])}, is below \code{-hessian_threshold}
#'   (that is, below 0.01 by default; for one established and one new marker
#'   this is the condition that the Hessian exceeds -0.01), the AUC criterion is
#'   treated as unidentified and the P value is \code{NA}.
#' @return an object of class \code{"aucnull_test"}: a list with the statistic,
#'   the weights \code{lambda}, the p value, the rank correlation estimates,
#'   the two smoothed areas, the number of clusters, and (if \code{compare = TRUE})
#'   the comparison tests.
#' @references
#' Heller G, Seshan VE, Moskowitz CS, Gonen M (2017). Inference for the
#' difference in the area under the ROC curve derived from nested binary
#' regression models. Biostatistics 18, 260-274.
#'
#' Han AK (1987). Non-parametric analysis of a generalized regression model.
#' Journal of Econometrics 35, 303-316.
#'
#' Sherman RP (1993). The limiting distribution of the maximum rank correlation
#' estimator. Econometrica 61, 123-137.
#'
#' Demler OV, Pencina MJ, D'Agostino RB (2012). Misuse of DeLong test to compare
#' AUCs for nested models. Statistics in Medicine 31, 2577-2587.
#' @examples
#' set.seed(1)
#' n <- 200
#' x <- rnorm(n); z <- rnorm(n)
#' y <- rbinom(n, 1, plogis(1.2 * x))           # z carries no information
#' aucnull_test(y, x, z, B = 20)                 # use the default B = 200 in practice
#' \donttest{
#' # clustered data: two observations per cluster
#' id <- rep(1:100, each = 2)
#' aucnull_test(y, x, z, B = 20, cluster = id)
#' }
#' @export
aucnull_test <- function(y, X, Z, method = c("boot", "plugin"), B = 200, expo = 1 / 3,
                         n_sim = 1e5, seed = NULL, compare = TRUE, cluster = NULL, hessian_threshold = -0.01) {
  method <- match.arg(method)
  y <- as.integer(as.logical(y))
  X <- as.matrix(X); Z <- as.matrix(Z)
  n <- length(y); p <- ncol(X); q <- ncol(Z)
  if (nrow(X) != n || nrow(Z) != n) stop("y, X and Z must have the same number of rows")
  if (!is.null(cluster) && length(cluster) != n) stop("cluster must have length n")
  ok <- !is.na(y) & stats::complete.cases(X) & stats::complete.cases(Z) & (if (is.null(cluster)) TRUE else !is.na(cluster))
  if (any(!ok)) {
    message(sum(!ok), " rows with missing values were dropped")
    y <- y[ok]; X <- X[ok, , drop = FALSE]; Z <- Z[ok, , drop = FALSE]
    if (!is.null(cluster)) cluster <- cluster[ok]
    n <- length(y)
  }
  if (sum(y) < 2 || sum(1 - y) < 2) stop("both outcome classes need at least two observations")
  cl <- if (is.null(cluster)) seq_len(n) else as.integer(factor(cluster))
  K <- max(cl); members <- split(seq_len(n), cl)
  sc <- function(M) sweep(sweep(M, 2, colMeans(M)), 2, apply(M, 2, sd), "/")
  Xs <- sc(X); Zs <- sc(Z)
  h <- default_bandwidth(Xs[, 1], n, expo)
  # starting values: logistic estimates scaled to the normalization (coefficient of X1 equal to one)
  xn <- paste0("x", seq_len(p)); zn <- paste0("z", seq_len(q))
  d0 <- data.frame(y = y, Xs, Zs); names(d0) <- c("y", xn, zn)
  fr <- glm(y ~ ., data = d0[, c("y", xn)], family = binomial())
  ff <- glm(y ~ ., data = d0, family = binomial())
  scaled_start <- function(cf, orient) {
    b1 <- orient * cf[xn[1]]
    if (!is.finite(b1) || abs(b1) < 1e-8) return(NULL)
    st <- c(if (p > 1) cf[xn[-1]], cf[zn[zn %in% names(cf)]]) / b1
    if (all(is.finite(st))) unname(st) else NULL
  }
  # orientation of the anchor: the reduced model is fitted with X1 and with -X1, and the orientation with the larger
  # smoothed area is kept (the chart that contains the true index attains the larger population area)
  fit_reduced <- function(orient) {
    Xo <- Xs; Xo[, 1] <- orient * Xo[, 1]
    list(orient = orient, Xs = Xo, red = mrc_fit(y, Xo, NULL, h, start = scaled_start(coef(fr), orient)))
  }
  cand <- list(fit_reduced(1), fit_reduced(-1))
  best <- if (cand[[2]]$red$auc > cand[[1]]$red$auc) cand[[2]] else cand[[1]]
  orient <- best$orient; Xs <- best$Xs; red <- best$red
  full <- mrc_fit(y, Xs, Zs, h, start = scaled_start(coef(ff), orient))
  stat <- 2 * n * (full$auc - red$auc)
  k <- length(full$theta); gi <- (k - q + 1):k      # indices of gamma
  D <- num_hessian(full$f, full$theta)
  Dinv <- tryCatch(solve(D), error = function(e) NULL)
  # identifiability: curvature of the smoothed criterion in gamma after profiling the other coefficients,
  # C = -solve(Dinv[gamma, gamma]); for p = q = 1 this is -D itself
  curv_gamma <- if (is.null(Dinv)) NULL else tryCatch(-solve(Dinv[gi, gi, drop = FALSE]), error = function(e) NULL)
  unidentified <- is.null(curv_gamma) || any(!is.finite(curv_gamma)) ||
    min(eigen((curv_gamma + t(curv_gamma)) / 2, symmetric = TRUE, only.values = TRUE)$values) <= -hessian_threshold
  if (unidentified) {
    warning("the largest eigenvalue of the Hessian of the smoothed area exceeds the threshold; the new markers may be nearly monotone in the established index, the AUC criterion is treated as unidentified, and no P value is reported")
  }
  if (method == "boot") {
    if (!is.null(seed)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        old_seed <- get(".Random.seed", envir = globalenv())
        on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
      }
      set.seed(seed)
    }
    G <- matrix(NA_real_, B, q)
    for (b in seq_len(B)) {
      idx <- if (is.null(cluster)) sample.int(n, n, replace = TRUE) else unlist(members[sample.int(K, K, replace = TRUE)], use.names = FALSE)
      fb <- mrc_fit(y[idx], Xs[idx, , drop = FALSE], Zs[idx, , drop = FALSE], h, start = full$theta)
      G[b, ] <- fb$theta[gi]
    }
    Vg <- n * cov(G)
  } else if (!unidentified) {
    h2 <- default_bandwidth(Xs[, 1], n, 1 / 3)
    Sg <- sigma_g(y, Xs, Zs, full$theta, h2, cl)
    V <- Dinv %*% Sg %*% Dinv
    Vg <- V[gi, gi, drop = FALSE]
  }
  if (unidentified) {
    lambda <- rep(NA_real_, q); pv <- NA_real_
  } else {
    Dgg <- Dinv[gi, gi, drop = FALSE]
    lambda <- Re(eigen(-Vg %*% solve(Dgg), only.values = TRUE)$values)
    pv <- pweightedchisq(stat, lambda, n_sim = n_sim)
  }
  out <- list(statistic = stat, lambda = lambda, p.value = pv, q = q, p = p, n = n, n_clusters = K,
              clustered = !is.null(cluster), orientation = orient, gamma = full$theta[gi], theta_full = full$theta,
              theta_reduced = red$theta, auc_full = full$auc, auc_reduced = red$auc, bandwidth = h,
              method = method, B = if (method == "boot") B else NA, scaled = TRUE)
  if (compare) {
    out$robust_wald <- robust_wald(ff, zn, cl)
    # restricted cubic spline (four knots) in each established marker, Wald test of Z with the robust covariance
    Bs <- do.call(cbind, lapply(seq_len(p), function(j) rcs_basis(X[, j])))
    if (is.null(Bs)) {
      out$spline_wald <- out$robust_wald
    } else {
      ds <- data.frame(y = y, X, Bs, Z); names(ds) <- c("y", xn, paste0("s", seq_len(ncol(Bs))), zn)
      fs <- glm(y ~ ., data = ds, family = binomial())
      out$spline_wald <- robust_wald(fs, zn, cl)
    }
    if (is.null(cluster)) {
      lrt <- as.numeric(2 * (logLik(ff) - logLik(fr)))
      out$lrt <- list(statistic = lrt, df = q, p.value = pchisq(lrt, q, lower.tail = FALSE),
                      gamma = coef(ff)[zn])
      out$delong <- delong_test(predict(ff), predict(fr), y)
    }
  }
  class(out) <- "aucnull_test"
  out
}

#' Wald test of a set of coefficients with a (cluster) robust covariance
#' @keywords internal
robust_wald <- function(fit, coefs, cl) {
  keep <- !is.na(coef(fit))
  Xm <- stats::model.matrix(fit)[, keep, drop = FALSE]; mu <- stats::fitted(fit); yv <- fit$y
  sc <- Xm * (yv - mu)
  G <- rowsum(sc, cl); K <- nrow(G)
  bread <- stats::vcov(fit, complete = FALSE)
  V <- bread %*% crossprod(G) %*% bread * K / max(K - 1, 1)
  b <- coef(fit)[coefs]; Vb <- V[coefs, coefs, drop = FALSE]
  w <- as.numeric(t(b) %*% solve(Vb) %*% b)
  list(statistic = w, df = length(coefs), p.value = pchisq(w, length(coefs), lower.tail = FALSE), gamma = b,
       n_clusters = K)
}

#' Restricted cubic spline basis (Harrell) without the linear term
#' @keywords internal
rcs_basis <- function(x, probs = c(0.05, 0.35, 0.65, 0.95)) {
  k <- unique(stats::quantile(x, probs, names = FALSE)); K <- length(k)
  if (K < 3) return(NULL)
  out <- sapply(seq_len(K - 2), function(j) {
    t <- function(a) pmax(x - a, 0)^3
    (t(k[j]) - t(k[K - 1]) * (k[K] - k[j]) / (k[K] - k[K - 1]) + t(k[K]) * (k[K - 1] - k[j]) / (k[K] - k[K - 1])) / (k[K] - k[1])^2
  })
  matrix(out, ncol = K - 2)
}

#' Gradient covariance of the smoothed area from the Hoeffding projections (plugin method)
#'
#' Per observation contributions are summed within clusters; with one observation
#' per cluster this is the usual two sample U statistic variance.
#' @keywords internal
sigma_g <- function(y, Xs, Zs, theta, h, cl = seq_along(y)) {
  x1 <- Xs[, 1]; p <- ncol(Xs)
  W <- cbind(if (p > 1) Xs[, -1, drop = FALSE], Zs)
  s <- x1 + as.vector(W %*% theta)
  i1 <- y == 1L; i0 <- y == 0L
  d <- outer(s[i1], s[i0], "-") / h
  phi <- dnorm(d) / h
  n1 <- sum(i1); n0 <- sum(i0); n <- length(y); k <- ncol(W)
  U <- matrix(0, n1, k); V <- matrix(0, n0, k)
  for (m in seq_len(k)) {
    dw <- outer(W[i1, m], W[i0, m], "-")
    U[, m] <- rowMeans(phi * dw); V[, m] <- colMeans(phi * dw)
  }
  E <- matrix(0, n, k)
  E[i1, ] <- sweep(U, 2, colMeans(U)) / n1
  E[i0, ] <- sweep(V, 2, colMeans(V)) / n0
  G <- rowsum(E, cl); K <- nrow(G)
  n * crossprod(G) * K / max(K - 1, 1)
}

#' Distribution function of a weighted sum of chi-square(1) variables
#'
#' Returns \code{P(sum lambda_j chi2_1j > x)}. Exact for a single weight;
#' Monte Carlo otherwise.
#' @param x value of the statistic.
#' @param lambda vector of weights.
#' @param n_sim number of Monte Carlo draws when \code{length(lambda) > 1}.
#' @return upper tail probability.
#' @export
pweightedchisq <- function(x, lambda, n_sim = 1e5) {
  lambda <- lambda[is.finite(lambda)]
  if (length(lambda) == 0L) return(NA_real_)
  if (length(lambda) == 1L) {
    if (lambda <= 0) return(NA_real_)
    return(pchisq(x / lambda, 1, lower.tail = FALSE))
  }
  Zm <- matrix(rnorm(n_sim * length(lambda))^2, n_sim)
  mean(as.vector(Zm %*% lambda) > x)
}

#' DeLong test for two correlated AUCs (paired design)
#'
#' Provided for comparison only. It is valid for two fixed scores on the same
#' subjects and is not valid for the linear predictors of nested models fitted
#' on the same data (Demler et al. 2012).
#' @param s1,s2 scores.
#' @param y binary outcome.
#' @return list with the two empirical areas, their difference, the standard
#'   error, the z statistic and the two sided p value.
#' @export
delong_test <- function(s1, s2, y) {
  y <- as.integer(as.logical(y))
  comp <- function(s) {
    a <- s[y == 1L]; b <- s[y == 0L]; m <- length(a); k <- length(b)
    r <- rank(c(a, b)); ra <- rank(a); rb <- rank(b)
    list(v10 = (r[1:m] - ra) / k, v01 = 1 - (r[(m + 1):(m + k)] - rb) / m)
  }
  c1 <- comp(s1); c2 <- comp(s2)
  m <- length(c1$v10); k <- length(c1$v01)
  S10 <- cov(cbind(c1$v10, c2$v10)); S01 <- cov(cbind(c1$v01, c2$v01))
  S <- S10 / m + S01 / k
  d <- mean(c1$v10) - mean(c2$v10)
  v <- S[1, 1] + S[2, 2] - 2 * S[1, 2]
  z <- if (v > 0) d / sqrt(v) else 0
  list(auc1 = mean(c1$v10), auc2 = mean(c2$v10), diff = d, se = sqrt(max(v, 0)), z = z,
       p.value = 2 * pnorm(-abs(z)))
}

#' Print method for \code{aucnull_test} objects
#'
#' @param x an object returned by \code{\link{aucnull_test}}.
#' @param digits number of significant digits (currently unused; four are shown).
#' @param ... further arguments, ignored.
#' @return \code{x}, invisibly.
#' @export
print.aucnull_test <- function(x, digits = 4, ...) {
  cat("Heller test of equal AUCs for new markers in nested models (Heller et al. 2017)\n")
  cat(sprintf("  n = %d%s, q = %d, bandwidth = %.4g, weights by %s%s\n", x$n,
              if (x$clustered) sprintf(" in %d clusters", x$n_clusters) else "", x$q, x$bandwidth, x$method,
              if (x$method == "boot") sprintf(" (B = %d%s)", x$B, if (x$clustered) ", cluster bootstrap" else "") else ""))
  cat(sprintf("  smoothed AUC: reduced %.4f, full %.4f\n", x$auc_reduced, x$auc_full))
  if (is.na(x$p.value)) cat(sprintf("  statistic 2n(A_full - A_red) = %.4g; the AUC criterion is unidentified (Hessian threshold), no P value\n", x$statistic)) else
  cat(sprintf("  statistic 2n(A_full - A_red) = %.4g, lambda = %s, p = %.4g\n",
              x$statistic, paste(signif(x$lambda, 4), collapse = ", "), x$p.value))
  cat("  rank correlation estimate of gamma (standardized scale): ", paste(signif(x$gamma, 4), collapse = ", "), "\n")
  if (!is.null(x$robust_wald)) {
    rb <- if (x$clustered) "cluster robust" else "robust"
    cat("  for comparison:\n")
    if (!is.null(x$lrt)) cat(sprintf("    likelihood ratio test of Z in the logistic model: chi2(%d) = %.4g, p = %.4g\n", x$lrt$df, x$lrt$statistic, x$lrt$p.value))
    cat(sprintf("    %s Wald test of Z in the logistic model: chi2(%d) = %.4g, p = %.4g\n", rb, x$robust_wald$df, x$robust_wald$statistic, x$robust_wald$p.value))
    cat(sprintf("    %s Wald test of Z after a restricted cubic spline in X: chi2(%d) = %.4g, p = %.4g\n", rb, x$spline_wald$df, x$spline_wald$statistic, x$spline_wald$p.value))
    if (!is.null(x$delong)) cat(sprintf("    DeLong test on the fitted linear predictors: diff = %.4g, p = %.4g (not valid for nested models)\n", x$delong$diff, x$delong$p.value))
  }
  invisible(x)
}
