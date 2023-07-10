#' @title
#' Penalized factorial regression on the GBLUPs
#'
#' @description
#' Based on multi-environment field trials, fits the factorial regression model
#' \deqn{G_{ij} = \mu + e_j + g_i + \sum_{k=1}^s \beta_{ik} x_{ij} + \epsilon_{ij},}
#' where .... with environmental main effects $e_j$, genotypic main effects \eqn{g_{i}} and
#' genotype-specific environmental sensitivities \eqn{\beta_{ik}}. See e.g. Millet
#' et al 2019 and van Eeuwijk / Denis 1996 (?). There are \eqn{s} environmental
#' indices with values \eqn{x_{ij}}. Optionally, predictions can be made for a set
#' of test environments, for which environmental indices are available. The new
#' environments must contain the same set of genotypes, or a subset
#' (For genomic prediction of new genotypes, the function ...).
#' Constraints: the first level of the factor genotype is put to 0.
#' Penalization: model (...) is fitted using glmnet, simultaneously penalizing
#' \eqn{e_j + g_i} and \eqn{\beta_{ik}}. If penG = FALSE (or penE = FALSE) the
#' main effects \eqn{g_i} (respectively \eqn{e_j}) are not been penalised.
#' Predictions for the test environments are first constructed using the
#' estimated main effects and sensitivities; next, predicted environmental
#' main effects are added. The latter are obtained by regressing the estimated
#' environmental main effects for the training environment on the average values
#' of the indices in these environments, as in Millet et al 2019.
#'
#' @inheritParams GnE
#'
#' @param K kinship matrix, or a list of such matrices. (TO DO: checks)
#' @param method (TO DO) "factReg" (default) or "perGeno"
#' @param type Character vector of length nrow(dat), which whould indicate the
#' type of observation, which should be one of "GE", "GnE", "nGE", "nGnE"
#' @param target character; should be either "GnE" or "nGnE" (default)
#' @param gpWeight vector of length \code{nrow(dat)}; default \code{NULL}
#' @param useRes If \code{method} is "perGeno": Indicates whether the
#' genotype-specific regressions are to be fitted on the residuals of a model
#' with main effects. 0 means no (i.e. fit the regressions on the original
#' data). 1: residuals of a model with environmental main effects.
#' 2 (default): genotypic and environmental main effects.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{predTrain}{Vector with predictions for the training set
#'   (to do: Add the factors genotype and environment; make a dataframe)}
#'   \item{predTest}{Vector with predictions for the test set (to do: Add the
#'   factors genotype and environment; make a dataframe). To do: add estimated
#'   environmental main effects, not only predicted environmental main effects}
#' }
#'
#' @export
frBLUP <- function(dat,
                   Y,
                   G,
                   E,
                   indices = NULL,
                   indicesData = NULL,
                   K,
                   method = c("factReg", "perGeno"),
                   testEnv = NULL,
                   type = NULL,
                   target = c("nGnE", "GnE"),
                   gpWeight = NULL,
                   useRes = 2,
                   corType = c("pearson", "spearman"),
                   partition = NULL,
                   alpha = 1,
                   lambda = NULL,
                   penE = 0,
                   scaling = c("train", "all", "no"),
                   quadratic = FALSE) {
  ## Input checks.
  if (!inherits(dat, "data.frame")) {
    stop("dat should be a data.frame.\n")
  }
  ## Get column names.
  traitName <- if (is.numeric(Y)) names(dat)[Y] else Y
  genoName <- if (is.numeric(G)) names(dat)[G] else G
  envName <- if (is.numeric(E)) names(dat)[E] else E
  ## Rename data columns for Y, G and E. - this includes checks.
  dat <- renameYGE(dat = dat, Y = Y, G = G, E = E)
  stopifnot(penE >= 0)
  scaling <- match.arg(scaling)
  corType <- match.arg(corType)
  ## Either indices or indicesData should be provided.
  if ((is.null(indices) && is.null(indicesData)) ||
      (!is.null(indices) && !is.null(indicesData))) {
    stop("Either indices or indicesData should be provided.\n")
  }
  if (!is.null(indices) && (!is.character(indices) ||
                            length(indices) <= 1 ||
                            !all(hasName(x = dat, name = indices)))) {
    stop("indices should be a vector of length > 1 of columns in dat.\n")
  }
  if (!is.null(indicesData)) {
    if (!inherits(indicesData, "data.frame") ||
        !all(levels(dat$E) %in% rownames(indicesData))) {
      stop("indicesData should be a data.frame with all environments in its ",
           "rownames.\n")
    }
    if (!all(rownames(indicesData) %in% levels(dat$E))) {
      stop("All environments in indicesData should be in dat.\n")
    }
    presCols <- colnames(indicesData)[colnames(indicesData) %in%
                                        colnames(dat)]
    if (length(presCols) > 0) {
      warning("The following columns in indicesDat are already in dat. Values ",
              "in dat will be overwritten:\n",
              paste(presCols, collapse = ", "), ".\n")
    }
  }
  ## Check testEnv.
  if (!is.null(testEnv) && (!is.character(testEnv) || length(testEnv) < 1 ||
                            !all(testEnv %in% levels(dat$E)))) {
    stop("testEnv should be a vector of environments present in dat.\n")
  }
  ## Check kinship.
  if (!is.null(K)) {
    if (!is.matrix(K) || !setequal(levels(dat$G), colnames(K)) ||
        !setequal(levels(dat$G), rownames(K))) {
      stop("K should be a matrix with all genotypes in dat in its row and ",
           "column names.\n")
    }
  }
  # either : type = NULL + testEnv specified
  # or     : type and target specified
  rownames(dat) <- paste0(dat$E, "_X_", dat$G)
  dat <- droplevels(dat)
  if (is.null(testEnv)) {
    dat$type <- type
    tre <- as.character(unique(dat$E[dat$type == "GE"]))
    tst <- setdiff(levels(dat$E), tre)
  } else {
    tst <- testEnv
    tre <- setdiff(levels(dat$E), tst)
  }
  genoAcc <- NULL
  if (target == "nGnE") {
    treG <- as.character(unique(dat$G[dat$type == "GE"]))
    tstG <- setdiff(levels(dat$G), treG)
    datG <- droplevels(dat[dat$type %in% c("GE", "nGE"), ])
    datG$Y[datG$G %in% tstG] <- NA
    gblups <- nGE(dat = datG, Y = "Y", G = "G", E = "E", K = K)
    gblups$pred <- gblups$pred[rownames(gblups$pred) %in% rownames(dat), ]
    # for the new environments, having yield itself is useful for evaluating nGnE
    dat$gblup <- dat$Y
    dat$weight <- NA
    dat[rownames(gblups$pred), "gblup"] <- gblups$pred$value
    # or without sqrt ???
    dat[rownames(gblups$pred), "weight"] <- 1 / sqrt(gblups$pred$PEV)
  } else if (target == "GnE") {
    datG <- droplevels(dat[dat$type == "GE", ])
    gblups <- nGE(dat = datG, Y = "Y", G = "G", E = "E", K = K)
    gblups$pred <- gblups$pred[rownames(gblups$pred) %in% rownames(dat), ]
    # for the new environments, having yield itself is useful for evaluating nGnE
    dat$gblup <- dat$Y
    dat$weight <- NA
    dat[rownames(gblups$pred), "gblup"]  <- gblups$pred$value
    # or without sqrt ???
    dat[rownames(gblups$pred), "weight"] <- 1 / sqrt(gblups$pred$PEV)
  }
  names(dat)[names(dat) == "Y"] <- "Yobs"
  if (sum(is.na(dat[dat$type %in% c("GE", "nGE"), "weight"])) > 0) {
    ww <- NULL
  } else {
    ww <- dat$weight
  }
  if (method == "factReg") {
    res <- GnE(dat = dat,
               Y = "gblup",
               G = "G",
               E = "E",
               scaling = scaling,
               indices = indices,
               indicesData = indicesData,
               weight = ww,
               partition = partition,
               testEnv = tst,
               quadratic = quadratic,
               penG = 0,
               penE = penE,
               alpha = alpha,
               lambda = lambda,
               corType = corType)
  } else if (method == "perGeno") {
    res <- perGeno(dat = dat,
                   Y = "gblup",
                   G = "G",
                   E = "E",
                   indices = indices,
                   indicesData = indicesData,
                   weight = ww,
                   partition = partition,
                   testEnv = tst,
                   corType = corType,
                   useRes = useRes)
  }
  genPred <- dat[c("E", "G", "gblup")]
  gpWeight <- dat[c("E", "G", "weight")]
  rownames(genPred) <- rownames(gpWeight) <- NULL
  res$genPred <- genPred
  res$gpWeight <- gpWeight
  res$Y <- traitName
  res$trainAccuracyEnv <- NULL
  res$RMSEtrain <- NULL
  return(res)
}
