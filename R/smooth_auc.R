#' Kernel smoothed area under the ROC curve
#'
#' Computes the kernel smoothed Mann and Whitney statistic
#' \eqn{A_h(s) = (n_0 n_1)^{-1} \sum_{i:y_i=1} \sum_{j:y_j=0} \Phi((s_i - s_j)/h)},
#' as used by Heller et al. (2017) to obtain a differentiable objective.
#'
#' @param score numeric vector of scores (higher means more likely to be a case).
#' @param y binary outcome (0/1, or logical).
#' @param h bandwidth. With \code{h = 0} the empirical (unsmoothed) area is returned.
#' @return a single number in (0, 1).
#' @export
smooth_auc <- function(score, y, h) {
  y <- as.integer(y)
  s1 <- score[y == 1L]; s0 <- score[y == 0L]
  if (h <= 0) {
    d <- outer(s1, s0, "-")
    return(mean((d > 0) + 0.5 * (d == 0)))
  }
  mean(pnorm(outer(s1, s0, "-") / h))
}

#' Default bandwidth
#'
#' \code{sd(x1) * n^(-expo)}, where \code{x1} is the first column of the
#' established marker matrix (whose coefficient is fixed at one).
#' @param x1 numeric vector.
#' @param n sample size.
#' @param expo exponent; the default 1/3 satisfies the rate conditions of the paper (0.2 is the choice of Heller et al. 2017).
#' @return a bandwidth.
#' @keywords internal
default_bandwidth <- function(x1, n, expo = 1 / 3) {
  sd(x1) * n^(-expo)
}
