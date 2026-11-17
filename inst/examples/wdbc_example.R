# Real data example: Wisconsin diagnostic breast cancer data (public; UCI Machine Learning Repository,
# Wolberg, Street and Mangasarian). n = 569, outcome malignant (212) vs benign (357).
# Main example of the paper (Section 7.3). Established marker X: standard error of concavity (empirical AUC 0.781).
# New marker Z: standard error of compactness (correlation 0.80 with X).
# Values reported in Section 7.3 of the paper (bandwidth sd(X) n^(-1/3), B = 500): likelihood ratio test p = 0.012;
# Heller test statistic 0.25, weight about 1.0, p about 0.6; MRC coefficient -0.05 (bootstrap SE 0.10);
# empirical AUC of the fitted full logistic model 0.748, below the 0.781 of X alone; DeLong on the fitted linear
# predictors p < 0.0001 in favour of the reduced model; likelihood ratio test for Z after a restricted cubic spline
# in X p = 0.68. R results will differ slightly because of the bootstrap and the optimizer.
library(aucnull)
url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
d <- read.csv(url, header = FALSE)
feat <- c("radius", "texture", "perimeter", "area", "smoothness", "compactness", "concavity",
          "concave_points", "symmetry", "fractal_dimension")
names(d) <- c("id", "diagnosis", paste0("mean_", feat), paste0("se_", feat), paste0("worst_", feat))
y <- as.integer(d$diagnosis == "M")
X <- d$se_concavity
Z <- d$se_compactness
r <- aucnull_test(y, X, Z, B = 500, seed = 1)
print(r)

# Mechanism checks (Theorem 6): nonlinearity of E[Z | X] and misspecification of the logistic link in X.
dec <- cut(X, quantile(X, 0:10 / 10), include.lowest = TRUE)
print(round(tapply(Z, dec, mean), 3))                       # E[Z | X] by decile of X
fit1 <- glm(y ~ X, family = binomial())
print(round(cbind(observed = tapply(y, dec, mean), logistic = tapply(fitted(fit1), dec, mean)), 2))
anova(lm(Z ~ X), lm(Z ~ poly(X, 3)))                        # curvature in E[Z | X]

# Agreeing example of the paper (Section 7.2): X = standard error of concave points (AUC 0.792), Z = standard error of texture.
# Expected: likelihood ratio p = 0.001; Heller test statistic 9.2, weight about 0.8, p about 0.0006; DeLong p = 0.39.
X2 <- d$se_concave_points
Z2 <- d$se_texture
r2 <- aucnull_test(y, X2, Z2, B = 500, seed = 1)
print(r2)
