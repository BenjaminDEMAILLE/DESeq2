#' Wald test for a negative binomial GLMM with random effects
#'
#' Fits, for each gene, a negative binomial generalized linear mixed model
#' (GLMM) with the fixed effects of \code{design(object)} and the random
#' effects given in \code{random}, for example \code{~ (1 | individual)}.
#' This is meant for comparisons between individuals when some individuals
#' contribute several samples (several tissues or time points), possibly
#' alongside individuals with a single sample. In such designs a fixed effect
#' for the individual makes the comparison inestimable, see the vignette
#' section on individuals nested within groups.
#'
#' The size factors (or normalization factors) and the dispersion estimates
#' of \code{object} are used as fixed inputs: the gene-wise dispersion is
#' not re-estimated, and enters the model through the negative binomial
#' family \code{MASS::negative.binomial(theta=1/dispersion)}. Each gene is fit
#' by \code{lme4::glmer} (Laplace approximation). As the dispersion is
#' estimated by \code{\link{estimateDispersions}} with the fixed effects only,
#' it may also absorb part of the variation between individuals, which makes
#' the test somewhat conservative.
#'
#' The coefficients, standard errors, Wald statistics and p-values are stored
#' in \code{mcols(object)} as with \code{\link{nbinomWaldTest}}, on the log2
#' scale, so that \code{\link{results}}, \code{\link{resultsNames}} and
#' \code{summary} work as usual, including with \code{contrast}: the
#' covariance matrix of the fixed effects is stored for each gene for this
#' purpose. Cook's distances are not defined for the mixed model, so no
#' outlier filtering is performed. For shrinkage of log2 fold changes, use
#' \code{lfcShrink} with \code{type="ashr"}: the \code{"apeglm"} and
#' \code{"normal"} estimators refit a model without the random effects.
#'
#' The standard deviation of each random effect is stored in
#' \code{mcols(object)} in columns starting with \code{randomSD_}, and
#' \code{mcols(object)$betaConv} is \code{FALSE} for genes where
#' \code{glmer} failed or reported a convergence problem. As with very
#' low counts in one group \code{glmer} can return extreme coefficients
#' with small standard errors, the Wald p-value is set to \code{NA} for
#' coefficients of 30 or more in absolute value on the log2 scale.
#'
#' @param object a DESeqDataSet with size factors and dispersion estimates,
#' and a design containing the fixed effects only
#' @param random a one-sided formula with the random-effect terms, in
#' \code{lme4} syntax, e.g. \code{~ (1 | individual)}. The variables must be
#' columns of \code{colData(object)}
#' @param control an optional \code{lme4::glmerControl()} object
#' @param quiet whether to print messages
#' @param parallel if FALSE, no parallelization. if TRUE, parallel
#' execution using \code{BiocParallel}, see \code{\link{DESeq}}
#' @param BPPARAM an optional parameter object passed internally
#' to \code{\link{bplapply}} when \code{parallel=TRUE}
#'
#' @return a DESeqDataSet with results columns accessible
#' with the \code{\link{results}} function. The coefficients and standard
#' errors are reported on a log2 scale.
#'
#' @seealso \code{\link{nbinomWaldTest}}, \code{\link{results}}
#'
#' @examples
#'
#' # two groups of three individuals, two conditions per individual
#' dds <- makeExampleDESeqDataSet(n=50, m=12)
#' dds$grp <- factor(rep(c("X","Y"), each=6))
#' dds$ind <- factor(rep(1:6, each=2))
#' dds$cnd <- factor(rep(c("A","B"), 6))
#' design(dds) <- ~ cnd + grp
#' dds <- estimateSizeFactors(dds)
#' dds <- estimateDispersions(dds)
#' if (requireNamespace("lme4", quietly=TRUE)) {
#'   dds <- nbinomMixedWaldTest(dds, random = ~ (1 | ind))
#'   res <- results(dds, name="grp_Y_vs_X")
#' }
#'
#' @export
nbinomMixedWaldTest <- function(object, random, control=NULL, quiet=FALSE,
                                parallel=FALSE, BPPARAM=bpparam()) {
  if (!requireNamespace("lme4", quietly=TRUE) |
      !requireNamespace("MASS", quietly=TRUE)) {
    stop("nbinomMixedWaldTest() requires installing the packages 'lme4' and 'MASS'")
  }
  if (is.null(dispersions(object))) {
    stop("testing requires dispersion estimates, first call estimateDispersions()")
  }
  # findbars() moved from lme4 to reformulas in recent versions
  findbars <- if (requireNamespace("reformulas", quietly=TRUE)) {
    reformulas::findbars
  } else {
    lme4::findbars
  }
  if (missing(random) || !is(random, "formula") ||
      length(findbars(random)) == 0) {
    stop("'random' should be a one-sided formula with random-effect terms, e.g. ~ (1 | individual)")
  }
  if (length(random) != 2) {
    stop("'random' should be a one-sided formula, e.g. ~ (1 | individual)")
  }
  randomVars <- all.vars(random)
  if (!all(randomVars %in% names(colData(object)))) {
    stop("the variables in 'random' should be columns of colData(object)")
  }
  if ("weights" %in% assayNames(object)) {
    stop("observation weights are not supported by nbinomMixedWaldTest()")
  }
  if (is.null(control)) {
    control <- lme4::glmerControl()
  }
  object <- sanitizeRowRanges(object)

  if ("results" %in% mcols(mcols(object))$type) {
    if (!quiet) message("found results columns, replacing these")
    object <- removeResults(object)
  }
  if (is.null(mcols(object)$allZero)) {
    object <- getBaseMeansAndVariances(object)
  }

  # only continue on the rows with non-zero row mean
  objectNZ <- object[!mcols(object)$allZero,,drop=FALSE]

  # fixed effects: the model matrix of the design, with the column
  # names used by nbinomWaldTest()
  modelAsFormula <- is(design(object), "formula")
  modelMatrix <- getModelMatrix(object)
  stopifnot(all(MatrixGenerics::colSums(abs(modelMatrix)) > 0))
  if (qr(modelMatrix)$rank < ncol(modelMatrix)) {
    stop("the model matrix of the fixed effects is not full rank")
  }
  modelMatrixNames <- colnames(modelMatrix)
  modelMatrixNames[modelMatrixNames == "(Intercept)"] <- "Intercept"
  modelMatrixNames <- make.names(modelMatrixNames)
  if (modelAsFormula && attr(terms(design(object)),"intercept") == 1) {
    convertNames <- renameModelMatrixColumns(colData(object), design(object))
    convertNames <- convertNames[convertNames$from %in% modelMatrixNames,,drop=FALSE]
    modelMatrixNames[match(convertNames$from, modelMatrixNames)] <- convertNames$to
  }
  colnames(modelMatrix) <- modelMatrixNames
  p <- ncol(modelMatrix)

  # the fixed effects enter glmer as one matrix term, so that the
  # coefficients come back in the order of the model matrix columns
  data <- as.data.frame(colData(object))[, randomVars, drop=FALSE]
  data$.X <- unname(modelMatrix)
  bars <- vapply(findbars(random), function(b) paste0("(", deparse(b), ")"), character(1))
  formula <- as.formula(paste(".y ~ 0 + .X + offset(.offset) +", paste(bars, collapse=" + ")))

  counts <- counts(objectNZ)
  logNF <- log(getSizeOrNormFactors(objectNZ))
  alpha <- dispersions(objectNZ)

  fitRows <- function(rows) {
    lapply(rows, function(i) {
      fitMixedRow(counts[i,], logNF[i,], alpha[i], data, formula, control, p)
    })
  }
  if (!quiet) message("fitting negative binomial GLMMs with lme4::glmer, ",
                      nrow(objectNZ), " genes")
  if (!parallel) {
    fits <- fitRows(seq_len(nrow(objectNZ)))
  } else {
    nworkers <- getNworkers(BPPARAM)
    idx <- factor(sort(rep(seq_len(nworkers), length.out=nrow(objectNZ))))
    fits <- do.call(c, bplapply(levels(idx), function(l) {
      fitRows(which(idx == l))
    }, BPPARAM=BPPARAM))
  }

  # random effect standard deviations: named after the first successful fit
  ok <- vapply(fits, function(f) !is.null(f$beta), logical(1))
  if (!any(ok)) {
    stop("lme4::glmer failed for all genes")
  }
  sdNames <- names(fits[[which(ok)[1]]]$randomSD)
  nsd <- length(sdNames)

  betaMatrix <- t(vapply(fits, function(f) if (is.null(f$beta)) rep(NA_real_, p) else f$beta, numeric(p)))
  covUpper <- upper.tri(diag(p), diag=TRUE)
  ncov <- sum(covUpper)
  covMatrix <- t(vapply(fits, function(f) {
    if (is.null(f$beta)) rep(NA_real_, ncov) else f$vcov[covUpper]
  }, numeric(ncov)))
  if (ncov == 1) covMatrix <- t(covMatrix)
  if (p == 1) betaMatrix <- t(betaMatrix)
  randomSD <- t(vapply(fits, function(f) {
    if (is.null(f$beta) || length(f$randomSD) != nsd) rep(NA_real_, nsd) else unname(f$randomSD)
  }, numeric(nsd)))
  if (nsd == 1) randomSD <- t(randomSD)
  betaConv <- vapply(fits, function(f) f$conv, logical(1))
  deviance <- vapply(fits, function(f) f$deviance, numeric(1))

  # convert from natural log to log2 scale
  betaMatrix <- log2(exp(1)) * betaMatrix
  covMatrix <- log2(exp(1))^2 * covMatrix
  betaSE <- sqrt(covMatrix[, which(diag(p)[covUpper] == 1), drop=FALSE])
  colnames(betaMatrix) <- modelMatrixNames
  colnames(betaSE) <- paste0("SE_", modelMatrixNames)
  WaldStatistic <- betaMatrix/betaSE
  colnames(WaldStatistic) <- paste0("WaldStatistic_", modelMatrixNames)
  WaldPvalue <- 2*pnorm(abs(WaldStatistic), lower.tail=FALSE)
  # with very low counts in one group, glmer can return coefficients
  # beyond the range fitNbinomGLMsOptim() treats as unstable (30 on the
  # log2 scale) with a deceptively small standard error: no p-value there
  WaldPvalue[abs(betaMatrix) >= 30] <- NA
  colnames(WaldPvalue) <- paste0("WaldPvalue_", modelMatrixNames)
  covIdx <- which(covUpper, arr.ind=TRUE)
  colnames(covMatrix) <- paste0("mixedCov_", covIdx[,1], "_", covIdx[,2])
  colnames(randomSD) <- make.names(paste0("randomSD_", sdNames))

  if (any(!betaConv)) {
    if (!quiet) message(paste(sum(!betaConv), "rows did not converge or failed in lme4::glmer, labelled in mcols(object)$betaConv"))
  }

  attr(object,"betaPrior") <- FALSE
  attr(object,"betaPriorVar") <- rep(1e6, p)
  attr(object,"modelMatrix") <- modelMatrix
  attr(object,"modelMatrixType") <- if (modelAsFormula) "standard" else "user-supplied"
  attr(object,"dispModelMatrix") <- modelMatrix
  attr(object,"test") <- "Wald"
  attr(object,"mixedRandom") <- random

  resultsList <- c(matrixToList(betaMatrix),
                   matrixToList(betaSE),
                   matrixToList(WaldStatistic),
                   matrixToList(WaldPvalue),
                   matrixToList(randomSD),
                   matrixToList(covMatrix),
                   list(betaConv = betaConv,
                        deviance = deviance,
                        maxCooks = rep(NA_real_, nrow(objectNZ))))
  WaldResults <- buildDataFrameWithNARows(resultsList, mcols(object)$allZero)

  modelMatrixNamesSpaces <- gsub("_"," ",modelMatrixNames)
  covNames <- paste0(modelMatrixNamesSpaces[covIdx[,1]], ", ", modelMatrixNamesSpaces[covIdx[,2]])
  mcols(WaldResults) <- DataFrame(
    type = rep("results", ncol(WaldResults)),
    description = c(paste("log2 fold change (MLE):", modelMatrixNamesSpaces),
                    paste("standard error:", modelMatrixNamesSpaces),
                    paste("Wald statistic:", modelMatrixNamesSpaces),
                    paste("Wald test p-value:", modelMatrixNamesSpaces),
                    paste("random effect standard deviation (natural log scale):", sdNames),
                    paste("covariance of coefficients (log2 scale):", covNames),
                    "convergence of lme4::glmer fit",
                    "deviance for the fitted mixed model",
                    "maximum Cook's distance for row (not defined for mixed models)"))
  mcols(object) <- cbind(mcols(object), WaldResults)
  object
}

# fit one row with lme4::glmer, returning coefficients and their
# covariance on the natural log scale, or NULL coefficients on failure
fitMixedRow <- function(y, offset, alpha, data, formula, control, p) {
  data$.y <- y
  data$.offset <- offset
  failed <- list(beta=NULL, conv=FALSE, deviance=NA_real_)
  fit <- tryCatch(suppressMessages(suppressWarnings(
    lme4::glmer(formula, data=data, control=control,
                family=MASS::negative.binomial(theta=1/alpha)))),
    error=function(e) NULL)
  if (is.null(fit)) return(failed)
  beta <- unname(lme4::fixef(fit))
  vcov <- unname(as.matrix(stats::vcov(fit)))
  if (length(beta) != p || any(!is.finite(vcov))) return(failed)
  vc <- as.data.frame(lme4::VarCorr(fit))
  vc <- vc[is.na(vc$var2),,drop=FALSE]
  randomSD <- stats::setNames(vc$sdcor, paste(vc$grp, vc$var1, sep="_"))
  conv <- is.null(fit@optinfo$conv$lme4$messages) ||
    all(grepl("singular", fit@optinfo$conv$lme4$messages))
  list(beta=beta, vcov=vcov, randomSD=randomSD, conv=conv,
       deviance=-2 * as.numeric(stats::logLik(fit)))
}

# contrasts for the results of nbinomMixedWaldTest(), using the stored
# covariance matrix of the coefficients (log2 scale)
getContrastMixed <- function(object, contrast) {
  coefColumns <- resultsNames(object)
  p <- length(coefColumns)
  stopifnot(length(contrast) == p)
  beta <- as.matrix(mcols(object)[, coefColumns, drop=FALSE])
  covIdx <- which(upper.tri(diag(p), diag=TRUE), arr.ind=TRUE)
  covColumns <- paste0("mixedCov_", covIdx[,1], "_", covIdx[,2])
  covMatrix <- as.matrix(mcols(object)[, covColumns, drop=FALSE])
  # c' V c, summing the upper triangle with off-diagonal terms counted twice
  w <- ifelse(covIdx[,1] == covIdx[,2], 1, 2) * contrast[covIdx[,1]] * contrast[covIdx[,2]]
  contrastEstimate <- as.vector(beta %*% contrast)
  contrastSE <- sqrt(as.vector(covMatrix %*% w))
  contrastStatistic <- contrastEstimate / contrastSE
  contrastPvalue <- 2*pnorm(abs(contrastStatistic), lower.tail=FALSE)
  DataFrame(log2FoldChange=contrastEstimate,
            lfcSE=contrastSE,
            stat=contrastStatistic,
            pvalue=contrastPvalue)
}
