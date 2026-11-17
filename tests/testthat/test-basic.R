test_that("smooth_auc reduces to the empirical area at h = 0", {
  set.seed(1); y <- rbinom(100, 1, 0.4); s <- rnorm(100) + y
  a0 <- smooth_auc(s, y, 0)
  w <- wilcox.test(s[y == 1], s[y == 0])$statistic / (sum(y) * sum(1 - y))
  expect_equal(as.numeric(a0), as.numeric(w), tolerance = 1e-10)
})

test_that("null case: p value is not tiny and statistic is nonnegative", {
  set.seed(2); n <- 300; x <- rnorm(n); z <- rnorm(n); y <- rbinom(n, 1, plogis(1.2 * x))
  r <- aucnull_test(y, x, z, B = 30, seed = 3)
  expect_true(r$statistic >= -1e-8)
  expect_true(r$p.value > 0.001)
  expect_true(all(r$lambda > 0))
})

test_that("alternative case: p value is small", {
  set.seed(4); n <- 400; x <- rnorm(n); z <- rnorm(n); y <- rbinom(n, 1, plogis(1.2 * x + 1.0 * z))
  r <- aucnull_test(y, x, z, B = 30, seed = 5)
  expect_true(r$p.value < 0.01)
  expect_true(r$lrt$p.value < 0.01)
})

test_that("plugin and boot give weights of the same order", {
  set.seed(6); n <- 300; x <- rnorm(n); z <- rnorm(n); y <- rbinom(n, 1, plogis(1.2 * x))
  rb <- aucnull_test(y, x, z, method = "boot", B = 40, seed = 7, compare = FALSE)
  rp <- aucnull_test(y, x, z, method = "plugin", compare = FALSE)
  expect_true(rp$lambda / rb$lambda > 0.3 && rp$lambda / rb$lambda < 3)
})

test_that("q = 2 runs and returns two weights", {
  set.seed(8); n <- 300; x <- rnorm(n); z <- matrix(rnorm(2 * n), n); y <- rbinom(n, 1, plogis(1.2 * x))
  r <- aucnull_test(y, x, z, B = 20, seed = 9, n_sim = 2e4, compare = FALSE)
  expect_length(r$lambda, 2)
  expect_true(is.finite(r$p.value))
})

test_that("clustered data: cluster bootstrap and cluster robust Wald run, iid comparators are dropped", {
  set.seed(10); K <- 150; id <- rep(seq_len(K), each = 2); n <- 2 * K
  a <- rnorm(K)[id]; b <- 0.8 * rnorm(K)[id]
  x <- (a + rnorm(n)) / sqrt(2); z <- rnorm(n); y <- rbinom(n, 1, plogis(1.2 * x + b))
  r <- aucnull_test(y, x, z, B = 20, seed = 11, cluster = id)
  expect_equal(r$n_clusters, K)
  expect_true(r$clustered)
  expect_null(r$lrt); expect_null(r$delong)
  expect_true(is.finite(r$robust_wald$p.value) && is.finite(r$spline_wald$p.value))
  rp <- aucnull_test(y, x, z, method = "plugin", cluster = id, compare = FALSE)
  expect_true(is.finite(rp$p.value))
  expect_error(aucnull_test(y, x, z, cluster = id[-1]))
})

test_that("rows with missing values are dropped with a message", {
  set.seed(12); n <- 200; x <- rnorm(n); z <- rnorm(n); y <- rbinom(n, 1, plogis(x)); z[1:3] <- NA
  expect_message(r <- aucnull_test(y, x, z, B = 10, compare = FALSE), "3 rows")
  expect_equal(r$n, n - 3)
})

test_that("anchor orientation with two established markers: the chart with the correct sign is chosen", {
  set.seed(13); n <- 600
  x1 <- rnorm(n); x2 <- 0.9 * x1 + sqrt(1 - 0.81) * rnorm(n)   # correlation 0.9
  z <- rnorm(n)
  y <- rbinom(n, 1, plogis(x1 - 1.5 * x2 + 1.0 * z))            # conditional coefficient of x1 positive, marginal AUC of x1 below 0.5
  r <- aucnull_test(y, cbind(x1, x2), z, B = 20, seed = 14, compare = FALSE)
  expect_equal(r$orientation, 1)
  expect_true(r$theta_full["eta1"] < 0)                          # coefficient of x2 relative to x1 is about -1.5
  expect_true(r$p.value < 0.01)
})

test_that("an unidentified criterion gives a warning and an NA P value", {
  set.seed(15); n <- 300; x <- rnorm(n); y <- rbinom(n, 1, plogis(1.5 * x)); z <- exp(x)   # z is a monotone transform of x
  expect_warning(r <- aucnull_test(y, x, z, B = 10, compare = FALSE), "unidentified")
  expect_true(is.na(r$p.value))
})
