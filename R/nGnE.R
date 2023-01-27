#' @title
#' here comes the title
#'
#' @description
#' \loadmathjax
#' here comes the description
#'
#' @param GnEOut ...
#' @param K ...
#' @param dNew ...
#' @param pMethod ...
#' @param corType type of correlation: Pearson (default) or spearman rank sum
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{parGenoAll}{...}
#'   \item{predNew}{...},
#'   \item{testAccuracyEnv}{...}
#' }
#'
#' @export
nGnE <- function(GnEOut,
                 K,
                 dNew,
                 pMethod = c("parametric", "nonparametric"),
                 corType = c("pearson", "spearman")) {
  pMethod <- match.arg(pMethod)
  corType <- match.arg(corType)
  # check:
  # target trait in d and dNew must be the same
  # K should have row and colnames, and of class matrix
  # change names dNew....
  # to do : add env main
  E <- GnEOut$E
  G <- GnEOut$G
  Y <- GnEOut$Y
  indices <- GnEOut$indices
  if (GnEOut$quadratic == TRUE) {
    originalIndices <- indices[setdiff((1:length(indices)),
                                       grep(x = indices, pattern = "_quad"))]
    dNewQuad <- (dNew[, originalIndices]) ^ 2
    colnames(dNewQuad) <- paste0(originalIndices, "_quad")
    dNew <- cbind(dNew, dNewQuad)
  }
  dNew <- droplevels(dNew)
  nIndices <- ncol(GnEOut$parGeno) - 1
  stopifnot(all(GnEOut$genotypes %in% colnames(K)))
  stopifnot(all(levels(dNew$G) %in% colnames(K)))
  trainGeno <- rownames(GnEOut$parGeno) #GnEOut$genotypes
  testGeno <- setdiff(levels(dNew[, G]), trainGeno) #GnEOut$genotypes)
  stopifnot(length(intersect(trainGeno, testGeno)) == 0)
  stopifnot(length(testGeno) > 0)
  genoAll <- c(trainGeno, testGeno)
  K <- K[genoAll, genoAll]
  nAll <- ncol(K)
  parGenoAll <- matrix(NA, nAll, nIndices + 1)
  colnames(parGenoAll) <- colnames(GnEOut$parGeno)
  rownames(parGenoAll) <- genoAll
  parGeno <- data.frame(GnEOut$parGeno, G = rownames(GnEOut$parGeno))
  if (pMethod == "parametric") {
    for (j in 1:(nIndices + 1)) {
      b <- rrBLUP::kin.blup(data = parGeno, geno = "G",
                            pheno = names(parGeno)[j], K = K)
      parGenoAll[names((b$pred)), j] <- as.numeric(b$pred)
    }
  }
  predNew <- data.frame(dNew, pred = as.numeric(parGenoAll[dNew[, G], 1]),
                        predMain = as.numeric(parGenoAll[dNew[, G], 1]))
  for (j in 2:(nIndices + 1)) {
    predNew$pred <- predNew$pred +
      as.numeric(parGenoAll[dNew[, G], j]) * dNew[, indices[j - 1]]
  }
  ## Compute statistics for test data.
  testAccuracyEnv <- getAccuracyEnv(datNew = dNew[, Y],
                                    datPred = predNew[["pred"]],
                                    datPredMain = predNew[["predMain"]],
                                    datE = dNew[, E],
                                    corType = corType)
  return(list(parGenoAll = parGenoAll,
              predNew = predNew,
              testAccuracyEnv = testAccuracyEnv))
}
