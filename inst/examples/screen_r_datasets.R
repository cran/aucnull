# Screen public R datasets for marker pairs where the likelihood ratio test rejects
# but the Heller test of the AUC difference does not. Run on a machine with R and the suggested packages:
#   install.packages(c("pROC", "MASS", "survival", "mlbench"))
# Outputs a CSV in the working directory.

library(aucnull)

screen_pairs <- function(y, M, label, min_abs_cor = 0.9, B = 100, seed = 1) {
  M <- as.matrix(M); nm <- colnames(M); out <- list()
  keep <- complete.cases(M) & !is.na(y); y <- y[keep]; M <- M[keep, , drop = FALSE]
  C <- cor(M)
  for (i in seq_len(ncol(M))) for (j in seq_len(ncol(M))) {
    if (i == j || abs(C[i, j]) > min_abs_cor) next
    r <- tryCatch(aucnull_test(y, M[, i], M[, j], B = B, seed = seed), error = function(e) NULL)
    if (is.null(r)) next
    out[[length(out) + 1]] <- data.frame(data = label, X = nm[i], Z = nm[j], n = r$n,
      auc_reduced = r$auc_reduced, auc_full = r$auc_full, stat = r$statistic, lambda = r$lambda[1],
      p_aucnull = r$p.value, p_lrt = r$lrt$p.value, p_delong = r$delong$p.value, gamma_mrc = unname(r$gamma))
    cat(sprintf("%-10s %-14s + %-14s AUC %.3f -> %.3f  AUCnull p=%.3f  LRT p=%.2e  DeLong p=%.3f\n",
                label, nm[i], nm[j], r$auc_reduced, r$auc_full, r$p.value, r$lrt$p.value, r$delong$p.value))
  }
  do.call(rbind, out)
}

res <- list()

if (requireNamespace("pROC", quietly = TRUE)) {
  data(aSAH, package = "pROC")
  y <- as.integer(aSAH$outcome == "Poor")
  M <- data.frame(wfns = as.integer(aSAH$wfns), s100b = aSAH$s100b, ndka = aSAH$ndka, age = aSAH$age)
  res$aSAH <- screen_pairs(y, M, "aSAH")
}

if (requireNamespace("MASS", quietly = TRUE)) {
  data(birthwt, package = "MASS")
  M <- data.frame(age = birthwt$age, lwt = birthwt$lwt, ptl = birthwt$ptl, ftv = birthwt$ftv)
  res$birthwt <- screen_pairs(birthwt$low, M, "birthwt")
  data(Pima.tr, package = "MASS"); data(Pima.te, package = "MASS")
  P <- rbind(Pima.tr, Pima.te)
  M <- P[, c("npreg", "glu", "bp", "skin", "bmi", "ped", "age")]
  res$Pima <- screen_pairs(as.integer(P$type == "Yes"), M, "Pima")
}

if (requireNamespace("survival", quietly = TRUE)) {
  data(pbc, package = "survival")
  pbc5 <- pbc[!is.na(pbc$time) & (pbc$status == 2 | pbc$time >= 5 * 365.25), ]
  y <- as.integer(pbc5$status == 2 & pbc5$time < 5 * 365.25)
  M <- pbc5[, c("age", "bili", "albumin", "protime", "copper", "platelet", "alk.phos", "ast", "chol")]
  res$pbc <- screen_pairs(y, M, "pbc5y")
  data(flchain, package = "survival")
  fl <- flchain[!is.na(flchain$creatinine), ]
  fl <- fl[fl$death == 1 | fl$futime >= 5 * 365.25, ]
  y <- as.integer(fl$death == 1 & fl$futime < 5 * 365.25)
  M <- fl[, c("age", "kappa", "lambda", "creatinine")]
  set.seed(2); idx <- sample(nrow(fl), min(3000, nrow(fl)))   # keep the run time moderate
  res$flchain <- screen_pairs(y[idx], M[idx, ], "flchain5y", B = 60)
}

all <- do.call(rbind, res)
write.csv(all, "aucnull_screen_results.csv", row.names = FALSE)
cat("\nPairs with LRT p < 0.01 and Heller test p > 0.10:\n")
print(all[all$p_lrt < 0.01 & all$p_aucnull > 0.10, ])
