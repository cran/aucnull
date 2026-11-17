# Clustered data example: K patients with two lesions each, as in Design V of the paper (Section 6) except that the
# new marker carries information here. A patient random intercept enters the risk and a shared patient component
# enters both markers. The first call ignores the clustering (ordinary bootstrap, model based covariances); the second
# resamples patients in the bootstrap and uses cluster robust covariances in the Wald comparisons.
library(aucnull)
set.seed(1)
K <- 150; id <- rep(seq_len(K), each = 2); n <- 2 * K
a <- rnorm(K)[id]; b <- 0.8 * rnorm(K)[id]; cc <- rnorm(K)[id]
x <- (a + rnorm(n)) / sqrt(2)
z <- (cc + rnorm(n)) / sqrt(2)
y <- rbinom(n, 1, plogis(1.2 * x + 0.5 * z + b))
print(aucnull_test(y, x, z, B = 200, seed = 2))                 # clustering ignored
print(aucnull_test(y, x, z, B = 200, seed = 2, cluster = id))   # patients resampled, cluster robust Wald
