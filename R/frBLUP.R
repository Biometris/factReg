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
#' @param genPred vector of length \code{nrow(dat)}; default \code{NULL}
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
                   mainCovariates = NULL,
                   mainCovPen = 1,
                   K,
                   method = c("factReg", "perGeno"),
                   testEnv = NULL,
                   type = NULL,
                   target = c("generic", "nGnE", "GnE"),
                   genPred = NULL,
                   gpWeight = NULL,
                   useRes = 2,
                   outputFile = "results_glmnet",
                   outputDir = getwd(),
                   corType = c("pearson", "spearman"),
                   partition = NULL,
                   alpha = 1,
                   lambda = NULL,
                   penE = 0,
                   scaling = c("no", "train", "all"),
                   postLasso = FALSE,
                   quadratic = FALSE,
                   verbose = FALSE) {
  ###############
  scaling <- match.arg(scaling)
  corType <- match.arg(corType)
  method  <- match.arg(method)
  target  <- match.arg(target)
  if (target == "generic") {
    stopifnot(!is.null(indicesData))
    # if target = 'generic', indices can only be at the level of env's
  }
  traitName <- ifelse(is.numeric(Y), names(dat)[Y], Y)
  ## Rename data columns for Y, G and E.
  dat <- renameYGE(dat = dat, Y = Y, G = G, E = E)
  # to do: check if K has row and column names.
  #
  # should be contained in levels(dat$G)
  # checks on target and dat$type
  #
  # either : type = NULL + testEnv specified
  # or     : type and target specified
  ##############
  rownames(dat) <- paste0(dat$E, "_X_", dat$G)
  dat <- droplevels(dat)
  if (is.null(testEnv)) {
    dat$type <- type
    stopifnot(sum(is.na(dat$G)) == 0)
    #dat <- dat[!is.na(dat$G),]

    #d1 <- droplevels(dat[dat$type %in% c('GE','GnE'),])
    #tre <- as.character(unique(d1$E[d1$type=='GE']))
    tre <- as.character(unique(dat$E[dat$type == "GE"]))
    tst <- setdiff(levels(dat$E), tre)
  } else {
    tst <- testEnv
    tre <- setdiff(levels(dat$E), tst)
  }
  ###################
  genoAcc <- NULL
  if (!is.null(genPred)) {
    dat$gblup  <- genPred
    dat$weight <- gpWeight
    if (target == "nGnE") {
      treG <- as.character(unique(dat$G[dat$type == "GE"]))
      tstG <- setdiff(levels(dat$G), treG)
      genoAcc <- tstG
    }
  } else if (target == "nGnE") {
    treG <- as.character(unique(dat$G[dat$type == "GE"]))
    tstG <- setdiff(levels(dat$G), treG)
    d3 <- droplevels(dat[dat$type %in% c("GE", "nGE"), ])
    d3$Y[d3$G %in% tstG] <- NA
    gblups <- nGE(dat = d3,
                  Y = "Y",
                  G = "G",
                  E = "E",
                  K = K)
    gblups$pred <- gblups$pred[rownames(gblups$pred) %in% rownames(dat), ]
    # for the new environments, having yield itself is useful for evaluating nGnE
    dat$gblup <- dat$Y #rep(NA, nrow(dat))
    dat$weight <- NA
    dat[rownames(gblups$pred), "gblup"]  <- gblups$pred$value
    # or without sqrt ???
    dat[rownames(gblups$pred), "weight"] <- 1 / sqrt(gblups$pred$PEV) # 1 / (gblups$pred$PEV)
    genoAcc <- tstG
  } else if (target == "GnE") {
    d3 <- droplevels(dat[dat$type == "GE", ])
    gblups <- nGE(dat = d3,
                  Y = "Y",
                  G = "G",
                  E = "E",
                  K = K)
    gblups$pred <- gblups$pred[rownames(gblups$pred) %in% rownames(dat), ]
    # for the new environments, having yield itself is useful for evaluating nGnE
    dat$gblup  <- dat$Y #rep(NA, nrow(dat))
    dat$weight <- NA
    dat[rownames(gblups$pred), "gblup"]  <- gblups$pred$value
    # or without sqrt ???
    dat[rownames(gblups$pred), "weight"] <- 1 / sqrt(gblups$pred$PEV) # 1 / (gblups$pred$PEV)
  } else if (target == "generic") {
    stopifnot(!is.null(indicesData))
    d3 <- droplevels(dat[!dat$E %in% tst, ])
    gblups <- nGE(dat = d3,
                  Y = "Y",
                  G = "G",
                  E = "E",
                  K = K)
    # now combine gblups (all genotypes, but only training env's) and
    # dat (not necessarily all genotypes, but also test env's)
    # we will add the unnique 'test' part from dat to gblups
    temp2 <- matrix(nrow(gblups$pred), ncol(dat))
    colnames(temp2) <- colnames(dat)
    rownames(temp2) <- rownames(gblups$pred)
    temp2 <- as.data.frame(temp2)
    temp2$G <- gblups$pred$G
    temp2$E <- gblups$pred$E
    temp2$G <- factor(temp2$G, levels = levels(dat$G))
    temp2$E <- factor(temp2$E, levels = levels(dat$E))
    temp2[rownames(dat) %in% rownames(temp2),
          setdiff(colnames(temp2), c("G", "E"))] <-
      dat[rownames(dat) %in% rownames(temp2),
          setdiff(colnames(dat), c("G", "E"))]
    dat <- rbind(dat[which(!(rownames(dat) %in% rownames(temp2))), ], temp2)
    # for the new environments, having yield itself is useful for evaluating nGnE
    dat$gblup  <- dat$Y #rep(NA, nrow(dat))
    dat$weight <- NA
    dat[rownames(gblups$pred), "gblup"]  <- gblups$pred$value
    # or without sqrt ???
    dat[rownames(gblups$pred), "weight"] <- 1 / sqrt(gblups$pred$PEV) # 1 / (gblups$pred$PEV)
  }
  names(dat)[names(dat) == "Y"] <- "Yobs"
  ################
  # sum(is.na(dat$gblup)); sum(is.na(dat[dat$E %in% tst, 'gblup']))
  # as.numeric(table(dat[which(!(dat$E %in% tst) & !is.na(dat$gblup)), 'G']))
  if (sum(is.na(dat$weight)) > 0) {
    ww <- NULL
  } else {
    ww <- dat$weight
  }
  if (method == "factReg") {
    out1 <- GnE(dat = dat,
                Y = "gblup",
                G = "G",
                E = "E",
                scaling = scaling,
                indices = indices,
                indicesData = indicesData,
                mainCovariates = mainCovariates,
                weight = ww,
                partition = partition,
                testEnv = tst,
                quadratic = quadratic,
                penG = 0,
                penE = penE,
                genoAcc = genoAcc,
                postLasso = postLasso,
                alpha = alpha,
                lambda = lambda,
                corType = corType)
  } else if (method == "perGeno") {
    out1 <- perGeno(dat = dat,
                    Y = "gblup",
                    G = "G",
                    E = "E",
                    indices = indices,
                    indicesData = indicesData,
                    weight = ww, #dat$weight,
                    partition = partition,
                    testEnv = tst,
                    genoAcc = genoAcc,
                    corType = corType,
                    useRes = useRes)
  }
  if (is.null(genPred)) {
    out1$genPred <- dat$gblup
    out1$gpWeight <- dat$weight
  } else {
    out1$genPred <- genPred
    out1$gpWeight <- gpWeight
  }
  out1$Y <- traitName
  out1$trainAccuracyEnv <- NULL
  out1$RMSEtrain <- NULL
  return(out1)
}
