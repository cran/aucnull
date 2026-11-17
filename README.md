# aucnull

Testing whether a new marker improves the area under the ROC curve beyond established markers in nested models.

## Why

Applying the DeLong test to the linear predictors of two nested logistic models is not valid: under the null the difference in areas is not asymptotically normal and the test is extremely conservative (Demler, Pencina and D'Agostino 2012). When the logistic link is correct, the likelihood ratio test answers the same question with more power (Pepe et al. 2013). When the link is misspecified, the null hypothesis of no added association and the null hypothesis of equal AUCs differ: the likelihood test can reject with probability tending to one although the new marker does not improve the ordering of risks. This package implements the test of the null hypothesis of equal AUCs of Heller, Seshan, Moskowitz and Gonen (2017), the Heller test, with the weights estimated by the bootstrap of the smoothed maximum rank correlation estimator (default) or by the sandwich formula; the test keeps its asymptotic level under any strictly increasing link.

## Install

From CRAN:

```r
install.packages("aucnull")
```

From the source archive:

```r
install.packages("remotes")
remotes::install_local("path/to/aucnull")
```

## Use

```r
library(aucnull)
set.seed(1)
n <- 300
x <- rnorm(n); z <- rnorm(n)
y <- rbinom(n, 1, plogis(1.2 * x))     # z carries no information
aucnull_test(y, x, z, B = 200)
```

The print method shows the statistic, the weights, the p value, and, for comparison, the likelihood ratio test, the Wald test of `Z` with a robust covariance, the same Wald test after a restricted cubic spline in each established marker, and the (invalid) DeLong test on the fitted linear predictors.

Arguments: `y` binary outcome; `X` established markers (vector or matrix; the first column is the anchor whose coefficient is fixed at one, and its orientation is chosen by fitting the reduced model with the anchor and with its negative and keeping the larger smoothed area); `Z` new markers (vector or matrix); `method = "boot"` (recommended) or `"plugin"` (sandwich formula); `B` bootstrap resamples (default 200); `expo` bandwidth exponent (default 1/3, bandwidth sd(X[,1]) n^(-expo)); `cluster` optional cluster identifiers.

### Clustered data

When several observations come from one unit, for example two lesions per patient, pass the unit identifier:

```r
id <- rep(1:150, each = 2)
aucnull_test(y, x, z, B = 200, cluster = id)
```

The bootstrap then resamples units, the plugin covariance sums the projections within units, and the Wald comparisons use a cluster robust covariance. The limiting distribution of the statistic is unchanged. The likelihood ratio test and the DeLong test, which assume independent observations, are not reported. The asymptotics are in the number of units, so the weights are poorly estimated when there are few units. Rows with missing values are dropped with a message.

## Scope and warnings

- The test addresses the hypothesis that the linear predictor with the largest population AUC does not use `Z`. It requires the identifiability conditions of Han (1987). If `Z` is a monotone transform of the anchor, the area is flat in the coefficient of `Z`, the Hessian is singular, and the test is undefined; when the curvature of the smoothed area in the coefficients of `Z` after profiling the other coefficients falls below 0.01 on the standardized scale (`hessian_threshold = -0.01`; for one established and one new marker this is the condition that the Hessian exceeds -0.01) a warning is issued and the P value is `NA`.
- The weights are estimated from `B` bootstrap resamples; the asymptotic validity of the test refers to an exact bootstrap covariance, and with `B = 200` the Monte Carlo error of the covariance alone raises the level of a nominal 5 percent test to about 0.051. Starting values are the logistic estimates scaled to the normalization; with one free coefficient the smoothed area is maximized by a bounded scalar search, otherwise by a quasi Newton search.
- The comparison tests reported with `compare = TRUE` are the likelihood ratio test, a Wald test with a robust covariance, a Wald test of `Z` with a robust covariance after a restricted cubic spline in each established marker, and the DeLong test; the robust Wald tests address the logistic null hypothesis, not the hypothesis of equal AUCs.
- All markers are standardized internally; the reported `gamma` is on the standardized scale.
- Run time is dominated by the bootstrap (`B` refits of the rank correlation estimate). For n up to a few thousand and `B = 200` this takes seconds to minutes.

## Reproducing the real data example

`inst/examples/wdbc_example.R` reproduces the two marker pairs of Section 7 of the paper (Wisconsin diagnostic breast cancer data, downloaded from the UCI repository). `inst/examples/screen_r_datasets.R` screens marker pairs in other public data sets shipped with CRAN packages (pROC, MASS, survival).

## References

Heller G, Seshan VE, Moskowitz CS, Gonen M (2017). Biostatistics 18, 260-274.
Han AK (1987). Journal of Econometrics 35, 303-316.
Sherman RP (1993). Econometrica 61, 123-137.
Demler OV, Pencina MJ, D'Agostino RB (2012). Statistics in Medicine 31, 2577-2587.
Pepe MS, Kerr KF, Longton G, Wang Z (2013). Statistics in Medicine 32, 1467-1482.
