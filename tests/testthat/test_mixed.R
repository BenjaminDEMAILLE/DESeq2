context("mixed")
test_that("nbinomMixedWaldTest fits random effects and works with results()", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("MASS")
  set.seed(1)
  # 8 individuals in two groups, two conditions within individual,
  # plus four individuals contributing a single sample
  n <- 40
  coldata <- DataFrame(grp=factor(c(rep(c("X","Y"), each=8), rep(c("X","Y"), each=2))),
                       ind=factor(c(rep(1:8, each=2), 9:12)),
                       cnd=factor(c(rep(c("A","B"), 8), rep("A", 4))))
  m <- nrow(coldata)
  indEffect <- matrix(rnorm(n * 12, 0, 0.5), n, 12)[, as.integer(coldata$ind)]
  grpEffect <- outer(c(rep(1, 10), rep(0, n - 10)), as.numeric(coldata$grp == "Y"))
  mu <- 200 * exp(indEffect + grpEffect)
  counts <- matrix(rnbinom(n * m, mu=mu, size=1/0.05), n, m)
  dds <- DESeqDataSetFromMatrix(counts, coldata, ~ cnd + grp)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds, quiet=TRUE)

  expect_error(nbinomMixedWaldTest(dds, random = ~ ind))
  expect_error(nbinomMixedWaldTest(dds, random = ~ (1 | foo)))

  dds <- nbinomMixedWaldTest(dds, random = ~ (1 | ind), quiet=TRUE)
  expect_equal(resultsNames(dds), c("Intercept","cnd_B_vs_A","grp_Y_vs_X"))
  res <- results(dds, name="grp_Y_vs_X")
  expect_true(all(res$pvalue >= 0 & res$pvalue <= 1, na.rm=TRUE))
  expect_true("randomSD_ind_.Intercept." %in% names(mcols(dds)))
  # the simulated group effect (log2 = 1.44) is recovered for the first 10 genes
  expect_equal(mean(res$log2FoldChange[1:10]), 1/log(2), tolerance=0.3)
  expect_true(mean(res$padj[1:10] < 0.1) >= 0.7)
  expect_true(mean(res$padj[11:40] < 0.1) <= 0.2)

  # character contrast with the reference level pulls the stored coefficient
  res2 <- results(dds, contrast=c("grp","Y","X"))
  expect_equal(res2$log2FoldChange, res$log2FoldChange)
  res3 <- results(dds, contrast=c("grp","X","Y"))
  expect_equal(res3$log2FoldChange, -res$log2FoldChange)
  # numeric contrast uses the stored covariance matrix
  res4 <- results(dds, contrast=c(0,1,-1))
  b <- as.matrix(mcols(dds)[, resultsNames(dds)])
  expect_equal(res4$log2FoldChange, as.vector(b %*% c(0,1,-1)))
  v1 <- mcols(dds)$mixedCov_2_2; v2 <- mcols(dds)$mixedCov_3_3; v12 <- mcols(dds)$mixedCov_2_3
  expect_equal(res4$lfcSE, sqrt(v1 + v2 - 2 * v12))

  # shrinkage: only ashr, which works from the results table
  expect_error(lfcShrink(dds, coef="grp_Y_vs_X", type="normal"))
  expect_error(lfcShrink(dds, coef="grp_Y_vs_X", type="apeglm"))

  # refitting without random effects clears the mixed-model state
  dds <- nbinomWaldTest(dds, quiet=TRUE)
  expect_null(attr(dds, "mixedRandom"))
  expect_false("mixedCov_1_1" %in% names(mcols(dds)))
})

test_that("nbinomMixedWaldTest agrees with nbinomWaldTest without individual variation", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("MASS")
  set.seed(2)
  dds <- makeExampleDESeqDataSet(n=30, m=12, interceptMean=6)
  dds$ind <- factor(rep(1:6, each=2))
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds, quiet=TRUE)
  fixed <- results(nbinomWaldTest(dds, quiet=TRUE))
  mixed <- results(nbinomMixedWaldTest(dds, random = ~ (1 | ind), quiet=TRUE))
  # with a near-zero random effect variance the two fits are close
  expect_equal(cor(fixed$log2FoldChange, mixed$log2FoldChange), 1, tolerance=0.02)
  expect_equal(cor(fixed$stat, mixed$stat), 1, tolerance=0.05)
})
