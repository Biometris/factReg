#' @title
#' Penalized factorial regression using glmnet
#'
#' @description
#' \loadmathjax
#' Based on multi-environment field trials, fits the factorial regression model
#' \mjeqn{Y_{ij} = \mu + e_j + g_i + \sum_{k=1}^s \beta_{ik} x_{ij} + \epsilon_{ij},}{ascii}
#' with environmental main effects \mjeqn{e_j}{ascii}, genotypic main effects \mjeqn{g_{i}}{ascii} and
#' genotype-specific environmental sensitivities \mjeqn{\beta_{ik}}{ascii}. See e.g. Millet
#' et al 2019 and Bustos-Korts et al 2019. There are \mjeqn{s}{ascii}
#' environmental indices with values \mjeqn{x_{ij}}{ascii}. Optionally, predictions can be made for a set
#' of test environments, for which environmental indices are available. The new
#' environments must contain the same set of genotypes, or a subset
#' (For genomic prediction of new genotypes, see the function nGnE).
#' Constraints: the genotypic main effect for the first genotype is put to 0.
#'
#' Penalization: model (refer to equation above) is fitted using glmnet, simultaneously penalizing
#' \mjeqn{e_j}{ascii}, \mjeqn{g_i}{ascii} and \mjeqn{\beta_{ik}}{ascii}.
#' If penG = 0 and penE = 0, the main effects \mjeqn{g_i}{ascii} and
#' \mjeqn{e_j}{ascii} are not penalized. If these parameters are 1, the
#' the main effects are penalized to the same degree as the sensitivities. Any
#' nonnegative values are allowed.
#'
#' Another extension is the option to inlcude common sensitivities that are
#' the same for each genotype, replacing model (...) by
#' \mjeqn{Y_{ij} = \mu + e_j + g_i + \sum_{k=1}^s (\beta_{k} + \beta_{ik}) x_{ij} + \epsilon_{ij},}{ascii}
#'  where, for the the kth environmental covariate,
#'  \mjeqn{\beta_{k}}{ascii} is the general and
#'  \mjeqn{\beta_{ik}}{ascii} the genotype-specific sensitivity.
#'
#' Predictions for the test environments are first constructed using the
#' estimated genotypic main effects and sensitivities; next, predicted environmental
#' main effects are added. The latter are obtained by regressing the estimated
#' environmental main effects for the training environment on the average values
#' of the indices in these environments, as in Millet et al 2019.
#'
#' @param dat A \code{data.frame} with data from multi-environment trials.
#' Each row corresponds to a particular genotype in a particular environment.
#' The data do not need to be balanced, i.e. an environment does not need to
#' contain all genotypes. \code{dat} should contain the training as well as the
#' test environments (see testEnv)
#' @param Y The trait to be analyzed: either of type character, in which case
#' it should be one of the column names in \code{dat}, or numeric, in which
#' case the Yth column of \code{dat} will be analyzed.
#' @param G The column in \code{dat} containing the factor genotype (either
#' character or numeric).
#' @param E The column in \code{dat} containing the factor environment
#' (either character or numeric).
#' @param indices The columns in \code{dat} containing the environmental
#' indices (vector of type character). Alternatively, if the indices are always
#' constant within environments (i.e. not genotype dependent), the
#' environmental data can also be provided using the argument \code{indicesData}
#' (see below).
#' @param indicesData An optional \code{data.frame} containing environmental
#' indices (covariates); one value for each environment and index. It should
#' have the environment names as row names (corresponding to the names
#' contained in \code{dat$E}); the column names are the indices. If
#' \code{indices} (see before) is also provided, the latter will be ignored.
#' @param mainCovariates Optional. The columns in \code{dat} containing main
#' covariates, i.e. with effects that are constant across genotypes. May
#' overlap with indices.
#' @param mainCovPen Numeric vector corresponding to \code{mainCovariates},
#' indicating the relative penalties used in glmnet. Default: all 1.
#' @param testEnv vector (character). Data from these environments are not used
#' for fitting the model. Accuracy is evaluated for training and test
#' environments separately. The default is \code{NULL}, i.e. no test
#' environments, in which case the whole dataset is training. It is also
#' possible that there are test environments, but without any data; in this
#' case, no accuracy is reported for test environments (CHECK correctness).
#' @param weight Numeric vector of length \code{nrow(dat)}, specifying the
#' weight (inverse variance) of each observation, used in glmnet. Default
#' \code{NULL}, giving constant weights.
#' @param outputFile The file name of the output files, without .csv extension,
#' which is added by the function
#' @param outputDir The directory to which output-files are to be
#' written.
#' @param corType type of correlation: Pearson (default) or spearman rank sum.
#' @param partition \code{data.frame} with columns E and partition. The column
#' E should contain the training environments (type character); partition
#' should be of type integer. Environments in the same fold should have
#' the same integer value. Default is \code{data.frame()}, in which case the
#' function uses a leave-one-environment out cross-validation. If \code{NULL},
#' the (inner) training sets used for cross-validation will be drawn randomly
#' from all observations, ignoring the environment structure. In the latter
#' case, the number of folds (nfolds) can be specified.
#' @param nfolds Default \code{NULL}. If \code{partition == NULL}, this can be
#' used to specify the number of folds to be used in glmnet.
#' @param alpha Type of penalty, as in glmnet (1 = LASSO, 0 = ridge; in between
#'  = elastic net). Default is 1.
#' @param lambda Numeric vector; defines the grid over which the penalty lambda
#' is optimised in cross validation.
#' Default: NULL (defined by glmnet).
#' Important special case: lambda = 0 (no penalty).
#' May be better handled with lm rather than glmnet (ask Vahe)
#' @param penG numeric; default 0. If 1, genotypic main effects are
#' penalized. If 0, they are not. Any non negative real number is allowed.
#' @param penE numeric; default 0. If 1, environmental main effects are
#' penalized. If 0, they are not. Any non negative real number is allowed.
#' @param scaling determines how the environmental variables are scaled.
#' "train" : all data (test and training environments) are scaled
#' using the mean and and standard deviation in the training environments.
#' "all" : using the mean and standard deviation of all environments.
#' "no" : No scaling.
#' @param postLasso boolean; default \code{FALSE}. Indicates whether to refit
#' the model without penalty, using only the variables selected by
#' lasso. To do: Require that alpha = 1.
#' @param quadratic boolean; default \code{FALSE}. If \code{TRUE}, quadratic
#' terms (i.e., squared indices) are added to the model. Only for those indices
#' that were selected (CLARIFY).
#' @param genoAcc character. The accuracies per environment are evaluated
#' using these genotypes. It must be a subset of the genotypes in \code{dat}
#' @param verbose boolean; default \code{FALSE}. If \code{TRUE}, the accuracies
#' per environment are printed on screen.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{predTrain}{Vector with predictions for the training set (to do: Add
#'   the factors genotype and environment; make a dataframe)}
#'   \item{predTest}{Vector with predictions for the test set (to do: Add the
#'   factors genotype and environment; make a dataframe). To do: add estimated
#'   environmental main effects, not only predicted environmental main effects}
#'   \item{resTrain}{Vector with residuals for the training set}
#'   \item{resTest}{Vector with residuals for the test set}
#'   \item{mu}{the estimated overall (grand) mean}
#'   \item{envInfoTrain}{The estimated environmental main effects, and the
#'   predicted effects, obtained when the former are regressed on the averaged
#'   indices, using penalized regression.}
#'   \item{envInfoTest}{The predicted environmental main effects for the test
#'   environments, obtained from penalized regression using the estimated
#'   main effects for the training environments and the averaged indices.}
#'   \item{parGeno}{dataframe containing the estimated genotypic main effects
#'   (first column) and sensitivities (subsequent columns)}
#'   \item{main.par}{If the option mainCovariates is used, this vector contains
#'   estimates of the main effects of the indices specified in mainCovariates.}
#'   \item{trainAccuracyEnv}{a data-frame with the accuracy (r) for each
#'   training environment, as well as the root mean squre error (RMSE), mean
#'   abosolute deviation (MAD) and rank (the latter is a proportion: how many
#'   of the best 5 genotypes are in the top 10). To be removed or further
#'   developed. All these quantities are also evaluated for a model with only
#'   genotypic and environmental main effects (columns rMain, RMSEmain and
#'   rankMain).}
#'   \item{testAccuracyEnv}{A data-frame with the accuracy for each test
#'   environment, with the same columns as trainAccuracyEnv.}
#'   \item{trainAccuracyGeno}{a data-frame with the accuracy (r) for each
#'   genotype, averaged over the training environments}
#'   \item{testAccuracyGeno}{a data-frame with the accuracy (r) for each
#'   genotype, averaged over the test environments}
#'   \item{lambda}{The value of lambda selected using cross validation.}
#'   \item{lambdaSequence}{...}
#'   \item{RMSEtrain}{The root mean squared error on the training environments}
#'   \item{RMSEtest}{The root mean squared error on the test environments}
#'   \item{Y}{The name of the trait that was predicted, i.e. the column name
#'   in \code{dat} that was used}
#'   \item{G}{The genotype label that was used, i.e. the argument G that was
#'   used}
#'   \item{E}{The environment label that was used, i.e. the argument E that
#'   was used}
#'   \item{indices}{The indices that were used, i.e. the argument indices that
#'   was used}
#'   \item{postLasso}{The postLasso option that was used}
#'   \item{quadratic}{The quadratic option that was used}
#' }
#'
#' @export
GnE <- function(dat,
                Y,
                G,
                E,
                indices = NULL,
                indicesData = NULL,
                mainCovariates = NULL,
                mainCovPen = 1,
                testEnv = NULL,
                weight = NULL,
                outputFile = "results_glmnet",
                outputDir = getwd(),
                corType = c("pearson", "spearman"),
                partition = data.frame(),
                nfolds = NULL,
                alpha = 1,
                lambda = NULL,
                penG = 0,
                penE = 0,
                scaling = c("train", "all", "no"),
                postLasso = FALSE,
                quadratic = FALSE,
                genoAcc = NULL,
                verbose = FALSE) {

  stopifnot(penG >= 0) # also > 1 is possible!
  stopifnot(penE >= 0)

  scaling <- match.arg(scaling)
  corType <- match.arg(corType)

  if (is.null(lambda)) {
    lambdaProvided <- FALSE
  } else {
    lambdaProvided <- TRUE
    lambda <- sort(lambda, decreasing = TRUE)
  }
  ## Get traitName.
  traitName <- ifelse(is.numeric(Y), names(dat)[Y], Y)
  ## Rename data columns for Y, G and E.
  dat <- renameYGE(dat = dat, Y = Y, G = G, E = E)
  ## Get training envs.
  trainEnv <- setdiff(levels(dat$E), testEnv)
  ## Remove missing values from training envs.
  dat <- dat[!(is.na(dat$Y) & dat$E %in% trainEnv),]
  dat <- droplevels(dat)

  ##################

  if (!is.null(indicesData)) {

    indices <- colnames(indicesData)
    nIndices <- ncol(indicesData)

    stopifnot(all(levels(dat$E) %in% rownames(indicesData)))
    dat <- dat[, setdiff(names(dat), indices)]

    qw <- matrix(NA, nrow(dat), nIndices)
    colnames(qw) <- indices
    dat <- data.frame(dat, qw)

    for (ind in indices) {
      for (env in levels(dat$E)) {
        dat[which(dat$E == env), ind] <- indicesData[[env, ind]]
      }
    }
  }
  ###################

  ## Scale environmental variables.
  if (scaling == "train") {
    muTr <- colMeans(dat[dat$E %in% trainEnv, indices])
    sdTr <- sapply(X = dat[dat$E %in% trainEnv, indices], sd)
    dat[, indices] <- scale(dat[, indices], center = muTr, scale = sdTr)
  } else if (scaling == "all") {
    dat[, indices] <- scale(dat[, indices])
  }

  if (is.null(weight)) {
    dat$W <- rep(1, nrow(dat))
  } else {
    stopifnot(length(weight)==nrow(dat))
    dat$W <- weight
  }

  ## Split dat into training and test set.
  dTrain <- dat[dat$E %in% trainEnv, ]

  redundantGeno <- names(which(table(dTrain$G[!is.na(dTrain$Y)]) < 10))
  if (length(redundantGeno) > 0) {
    warning("the following genotypes have < 10 observations, and are removed:",
            "\n\n", paste(redundantGeno, collapse = ", "), "\n\n",
            "See also the documentation of the function nGnE.\n")
    dTrain <- dTrain[!(dTrain$G %in% redundantGeno), ]
  }
  dTrain <- droplevels(dTrain)

  if (!is.null(testEnv)) {
    dTest <- dat[dat$E %in% testEnv, ]
    dTest <- droplevels(dTest)
    nGenoTest <- nlevels(dTest$G)
  } else{
    dTest <- NULL
  }

  nEnvTest <- length(testEnv)
  nEnvTrain <- nlevels(dat$E) - nEnvTest
  nGenoTrain <- nlevels(dTrain$G)

  ## When partition == data.frame() (the default),
  ## do leave one environment out cross-validation.
  if (!is.null(partition)) {
    if (ncol(partition) == 0) {
      partition <- unique(data.frame(E = dTrain$E,
                                     partition = as.numeric(dTrain$E)))
    } else {
      stopifnot(all(c("E", "partition") %in% colnames(partition)))
      stopifnot(all(dTrain$E %in% partition$E))
      partition <- partition[partition$E %in% trainEnv, ]
    }
    ## Construct foldid from partition
    foldDat <- merge(dTrain, partition)
    foldid <- foldDat$partition
    nfolds <- length(unique(foldid))
  } else {
    foldid <- NULL
    if (is.null(nfolds)) {nfolds <- 10}
  }

  # just to het the names ...
  w2 <- aggregate(x = dTrain$Y, by = list(dTrain$G),
                  FUN = mean, na.rm = TRUE)
  e3 <- aggregate(x = dTrain$Y, by = list(dTrain$E),
                  FUN = mean, na.rm = TRUE)

  mm <- Matrix::sparse.model.matrix(Y ~ E + G, data = dTrain)
  modelMain <- glmnet::glmnet(y = dTrain$Y, x = mm, thresh = 1e-18, lambda = 0)
  cf <- c(modelMain$a0, modelMain$beta[-1])
  names(cf) <- c("(Intercept)", paste0('E', e3$Group.1[-1]),
                 paste0('G', w2$Group.1[-1]))

  mainOnly <- rep(0, nGenoTrain)
  names(mainOnly) <- levels(dTrain$G)
  mainTemp <- cf[(nEnvTrain + 1):(nEnvTrain + nGenoTrain - 1)]
  names(mainTemp) <- substring(names(mainTemp), first = 2)
  mainOnly[names(mainTemp)] <- mainTemp

  ## Define the design matrix for the factorial regression model.
  geFormula <- as.formula(paste0("Y ~ -1 + E + G +",
                                 paste(paste0(indices, ":G"),
                                       collapse = " + ")))
  ## Construct design matrix for training + test set.
  # important when there are NA's in the test set.
  opts <- options(na.action = "na.pass")
  on.exit(options(opts), add = TRUE)
  ma <- Matrix::sparse.model.matrix(geFormula, data = rbind(dTrain, dTest))
  colnames(ma)[1:(nEnvTrain + nGenoTrain + nEnvTest - 1)] <-
    substring(colnames(ma)[1:(nEnvTrain + nGenoTrain + nEnvTest - 1)], first = 2)

  if (!is.null(mainCovariates)) {
    ncolMaOld <- ncol(ma)
    ma <- cbind(ma,
                Matrix::Matrix(as.matrix(rbind(dTrain, dTest)[, mainCovariates])))
  }

  if (quadratic) {
    ## Add quadratic columns to design matrix.
    mQuad <- ma[, (nEnvTrain + nEnvTest + nGenoTrain):ncol(ma)] ^ 2
    colnames(mQuad) <- paste0(colnames(mQuad), "_quad")
    ma <- cbind(ma, mQuad)
    ## Add the quadratic columns to the indices.
    # from this point: only used to estimate env.main effects for test env.
    # If there are mainCovariates, these are not included here
    indices <- c(indices, paste0(indices, "_quad"))
  }

  # define the vector indicating to which extent parameters are to be penalized.
  penaltyFactorA <- rep(1, ncol(ma))
  penaltyFactorA[1:(nEnvTrain + nEnvTest)] <- penE
  penaltyFactorA[(nEnvTrain + nEnvTest + 1):
                   (nEnvTrain + nEnvTest + nGenoTrain - 1)] <- penG

  if (!is.null(mainCovariates)) {
    if (quadratic) {
      tempInd <- c((ncolMaOld + 1):(ncolMaOld + length(mainCovariates)),
                   (ncol(ma) - length(mainCovariates) + 1):ncol(ma))
      #(ncol(ma) - 2 * length(mainCovariates) + 1):ncol(ma)
    } else {
      tempInd <- (ncol(ma) - length(mainCovariates) + 1):ncol(ma)
    }
    penaltyFactorA[tempInd] <- mainCovPen
  } else {
    tempInd <- integer()
  }

  # note: even if unpenalized, the estimated main effects change,
  # depending on lambda!

  # run glmnet, either (if provided) for a single value of lambda
  # (using glmnet), or using cv.glmnet

  thr <- 1e-07

  if (lambdaProvided && length(lambda) == 1) {
    if (lambda == 0) {
      thr <- 1e-11
    }
    glmnetOutA <- glmnet::glmnet(x = ma[1:nrow(dTrain),],
                                 y = dTrain$Y,
                                 lambda = lambda,
                                 weights = dTrain$W,
                                 alpha = alpha, standardize = TRUE,
                                 penalty.factor = penaltyFactorA,
                                 intercept = TRUE,
                                 thres = thr)
    lambdaIndex <- 1
    lambdaSequence <- lambda
    cfe <- glmnetOutA$beta
    mu <- as.numeric(glmnetOutA$a0)

    if (!is.null(mainCovariates)) {
      main.par <- glmnetOutA$beta[tempInd]
      names(main.par) <- colnames(ma)[tempInd]
    } else {
      main.par <- NULL
    }
  } else {
    glmnetOutA <- glmnet::cv.glmnet(x = ma[1:nrow(dTrain),], y = dTrain$Y,
                                    lambda = lambda,
                                    weights = dTrain$W,
                                    foldid = foldid, nfolds = nfolds,
                                    alpha = alpha, standardize = TRUE,
                                    penalty.factor = penaltyFactorA,
                                    intercept = TRUE)
    lambda <- glmnetOutA$lambda
    lambdaSequence <- lambda
    lambdaIndex <- which(lambda == glmnetOutA$lambda.min)
    cfe <- glmnetOutA$glmnet.fit$beta
    mu <- glmnetOutA$glmnet.fit$a0[lambdaIndex]

    if (postLasso) {
      lassoModel <- glmnet::glmnet(x = ma , y = dTrain$Y, alpha = 1,
                                   weights = dTrain$W,
                                   family = "gaussian",
                                   lambda = glmnetOutA$lambda.min)

      ## Post-Selection, with selected significant variables: Refitting.
      ## Select the non-zero coefficients of env. variables.
      coeffic <- coef(lassoModel)[-1]
      ## Keep the subsample of the significant variables.
      positions <- which(coeffic != 0)
      # before :
      # positions <- union(1:(nEnvTrain + nGenoTrain - 1), positions)
      positions <- union(1:(nEnvTrain + nEnvTest + nGenoTrain - 1), positions)
      ## Refitting with the selected variables, without penalty.
      ma <- ma[, positions]
      glmnetOutA <- glmnet::glmnet(x = ma, y = dTrain$Y, lambda = 0,
                                   family = "gaussian",
                                   weights = dTrain$W,
                                   thresh = 1e-11)
      lambdaIndex <- 1
      cfe <- glmnetOutA$beta
      mu <- as.numeric(glmnetOutA$a0)
    }
    if (!is.null(mainCovariates)) {
      main.par <- glmnetOutA$glmnet.fit$beta[tempInd, lambdaIndex]
    } else {
      main.par <- NULL
    }
  }
  ## Extract from the output the estimated environmental main effects
  envMain <- as.matrix(cfe[1:nEnvTrain, lambdaIndex, drop = FALSE])
  envMain2 <- envMain
  envMain2[,1] <- cf[1:nEnvTrain] #c(0,cf)
  ## Create a data.frame that will contain the estimated genotypic
  ## main effects (first column), and the estimated environmental
  ## sensitivities (subsequent) columns.
  tempInd2 <- setdiff((nEnvTrain + nEnvTest + nGenoTrain):ncol(ma), tempInd)
  ## assign the genotypic main effects
  parGeno <-
    data.frame(main = c(0, cfe[(nEnvTrain + nEnvTest + 1):
                                 (nGenoTrain + nEnvTrain + nEnvTest - 1),
                               lambdaIndex]), row.names = levels(dTrain$G))
  ## assign the genotype specific sensitivities (subsequent columns)
  #cfeIndRows <- rownames(cfe[(nGenoTrain + nEnvTrain + nEnvTest):(nrow(cfe) - length(mainCovariates)), , drop = FALSE])
  cfeIndRows <- rownames(cfe[tempInd2, , drop = FALSE])
  cfeGenotypes <- sapply(X = cfeIndRows, FUN = function(cfeIndRow) {
    substring(strsplit(x = cfeIndRow, split = ":")[[1]][1], first = 2)
  })
  cfeIndices <- sapply(X = cfeIndRows, FUN = function(cfeIndRow) {
    strsplit(x = cfeIndRow, split = ":")[[1]][2]
  })
  # cfeDf <- data.frame(geno = cfeGenotypes, index = cfeIndices,
  #                     val = cfe[(nGenoTrain + nEnvTrain + nEnvTest):(nrow(cfe) - length(mainCovariates)), lambdaIndex],
  #                     stringsAsFactors = FALSE)
  cfeDf <- data.frame(geno = cfeGenotypes, index = cfeIndices,
                      val = cfe[tempInd2, lambdaIndex],
                      stringsAsFactors = FALSE)
  cfeDf$geno <- factor(cfeDf$geno, levels = rownames(parGeno))
  cfeDf$index <- factor(cfeDf$index, levels = indices)
  parGenoIndices <- reshape2::dcast(cfeDf, geno ~ index, fill = 0, drop = FALSE,
                                    value.var = "val")
  parGeno <- merge(parGeno, parGenoIndices, by.x = "row.names", by.y = "geno")
  rownames(parGeno) <- parGeno[["Row.names"]]
  parGeno <- parGeno[, c("main", indices)]
  ## Make predictions for training set.
  predTrain <- as.numeric(predict(object = glmnetOutA, newx = ma[1:nrow(dTrain), ],
                                  s = "lambda.min"))
  resTrain <- dTrain$Y - predTrain
  if (!is.null(testEnv)) {
    predTest <- as.numeric(
      predict(object = glmnetOutA,
              newx = ma[(nrow(dTrain)+ 1):(nrow(dTrain) + nrow(dTest)), ],
              s = "lambda.min"))
    resTest <- dTest$Y - predTest
  } else {
    predTest <- NULL
    resTest <- NULL
  }
  ## Compute the mean of each environmental index, in each environment
  if (!quadratic) {
    indFrame <- aggregate(dat[, indices], by = list(E = dat$E), FUN = mean)
  } else {
    indFrame <- merge(aggregate(dat[, indices[1:(length(indices) / 2)]],
                                by = list(E = dat$E), FUN = mean),
                      aggregate(dat[, indices[1:(length(indices) / 2)]] ^ 2,
                                by = list(E = dat$E), FUN = mean), by = "E")
    colnames(indFrame)[-1] <- indices
  }
  rownames(indFrame) <- indFrame$E
  indFrameTrain <- indFrame[indFrame$E %in% trainEnv, ]
  if (!is.null(testEnv)) {
    indFrameTest  <- indFrame[indFrame$E %in% testEnv, ]
  }
  indFrameTrain2 <- merge(indFrameTrain, envMain2,
                          by.x = "E", by.y = "row.names")
  indFrameTrain <- merge(indFrameTrain, envMain, by.x = "E", by.y = "row.names")
  colnames(indFrameTrain)[ncol(indFrameTrain)] <- "envMainFitted"
  colnames(indFrameTrain2)[ncol(indFrameTrain2)] <- "envMainFitted"
  if (!is.null(partition)) {
    indFrameTrain <- merge(indFrameTrain, partition)
    indFrameTrain2 <- merge(indFrameTrain2, partition)
  }
  rownames(indFrameTrain) <- indFrameTrain$E
  rownames(indFrameTrain2) <- indFrameTrain2$E
  ## predict env. main effects.
  ## (to do : quadratic terms. NO post-LASSO ?)
  ## note : parEnvTrain and parEnvTest will now be matrices; not vectors.
  if (is.null(partition)) {
    glmnetOut <- glmnet::cv.glmnet(x = as.matrix(indFrameTrain[, indices]),
                                   y = indFrameTrain$envMainFitted,
                                   alpha = alpha, nfolds = nfolds)
    glmnetOut2 <- glmnet::cv.glmnet(x = as.matrix(indFrameTrain2[, indices]),
                                    y = indFrameTrain2$envMainFitted,
                                    alpha = alpha, nfolds = nfolds)
  } else {
    glmnetOut <- glmnet::cv.glmnet(x = as.matrix(indFrameTrain[, indices]),
                                   y = indFrameTrain$envMainFitted,
                                   alpha = alpha,
                                   foldid = indFrameTrain$partition)
    glmnetOut2 <- glmnet::cv.glmnet(x = as.matrix(indFrameTrain2[, indices]),
                                    y = indFrameTrain2$envMainFitted,
                                    alpha = alpha,
                                    foldid = indFrameTrain2$partition)
  }
  # names(which(glmnetOut$glmnet.fit$beta[,which(glmnetOut$lambda == glmnetOut$lambda.min)]>0))
  parEnvTrain <- predict(object = glmnetOut,
                         newx = as.matrix(indFrameTrain[, indices]),
                         s = "lambda.min")
  parEnvTrain2 <- predict(object = glmnetOut2,
                          newx = as.matrix(indFrameTrain2[, indices]),
                          s = "lambda.min")
  if (!is.null(testEnv)) {
    parEnvTest  <- predict(object = glmnetOut,
                           newx = as.matrix(indFrameTest[, indices]),
                           s = "lambda.min")
    parEnvTest2 <- predict(object = glmnetOut2,
                           newx = as.matrix(indFrameTest[, indices]),
                           s = "lambda.min")
  }
  indicesTest <- NULL
  if (!is.null(testEnv)) {
    ## add the predicted environmental main effects:
    predTest <- predTest + as.numeric(parEnvTest[as.character(dTest$E), ])
    indicesTest <- data.frame(indFrameTest,
                              envMainPred = as.numeric(parEnvTest))
  }
  indicesTrain <- data.frame(indFrameTrain,
                             envMainPred = as.numeric(parEnvTrain))
  if (!is.null(genoAcc)) {
    predTrain <- predTrain[which(dTrain$G %in% genoAcc)]
    dTrain <- dTrain[dTrain$G %in% genoAcc,]
    dTrain <- droplevels(dTrain)
    if (!is.null(testEnv)) {
      predTest <- predTest[which(dTest$G %in% genoAcc)]
      dTest <- dTest[dTest$G %in% genoAcc,]
      dTest <- droplevels(dTest)
    }
  }
  ## Split raw data and predictions by environment for computing statistics.
  s1 <- split(dTrain$Y, dTrain$E)
  s2 <- split(predTrain, dTrain$E)
  s1G <- split(dTrain$Y, dTrain$G)
  s2G <- split(predTrain, dTrain$G)
  fm <- indFrameTrain$envMainFitted
  names(fm) <- as.character(indFrameTrain$E)
  #s3 <- split(as.numeric(mainOnly[as.character(dTrain$G)]) + as.numeric(fm[as.character(dTrain$E)]), dTrain$E)
  s3 <- split(as.numeric(predict(modelMain, newx = mm)), dTrain$E)

  #s3 <- split(as.numeric(mainOnly[as.character(dTrain$G)]) + as.numeric(envMain[dTrain$E,1]), dTrain$E)
  #  + fm[as.character(dTrain$E)] ??
  # plot(cfe[1:25], cf[1:25])

  ## Compute statistics for training data.
  trainAccuracyEnv <-
    data.frame(Env = levels(dTrain$E),
               r = mapply(FUN = cor, s1, s2,
                          MoreArgs = list(use = "na.or.complete",
                                          method = corType)),
               rMain = mapply(FUN = cor, s1, s3,
                              MoreArgs = list(use = "na.or.complete",
                                              method = corType)),
               RMSE = mapply(FUN = function(x, y) {
                 sqrt(mean((x - y) ^ 2, na.rm = TRUE))
                 }, s1, s2),
               RMSEmain = mapply(FUN = function(x, y) {
                 sqrt(mean((x - y) ^ 2, na.rm = TRUE))
                 }, s1, s3),
               MAD = mapply(FUN = function(x, y) {
                 mean(abs(x - y), na.rm = TRUE)
                 }, s1, s2),
               rank = mapply(FUN = exRank, s1, s2),
               rankMain = mapply(FUN = exRank, s1, s3),
               row.names = NULL)
  trainAccuracyGeno <-
    data.frame(Geno = levels(dTrain$G),
               r = mapply(FUN = cor, s1G, s2G,
                          MoreArgs = list(use = "na.or.complete",
                                          method = corType)))
  if (!is.null(testEnv)) {
    ## Split raw data and predictions by environment for computing statistics.
    s1t <- split(dTest$Y, dTest$E)
    s2t <- split(predTest, dTest$E)
    s1tG <- split(dTest$Y, dTest$G)
    s2tG <- split(predTest, dTest$G)
    s3t <- split(mainOnly[as.character(dTest$G)] +
                   as.numeric(parEnvTest2[as.character(dTest$E), ]), dTest$E)
    ## Compute statistics for test data.
    testAccuracyEnv <-
      data.frame(Env = levels(dTest$E),
                 r = mapply(FUN = cor, s1t, s2t,
                            MoreArgs = list(use = "na.or.complete",
                                            method = corType)),
                 rMain = mapply(FUN = cor, s1t, s3t,
                                MoreArgs = list(use = "na.or.complete",
                                                method = corType)),
                 RMSE = mapply(FUN = function(x, y) {
                   sqrt(mean((x - y) ^ 2, na.rm = TRUE))
                   }, s1t, s2t),
                 RMSEmain = mapply(FUN = function(x, y) {
                   sqrt(mean((x - y) ^ 2, na.rm = TRUE))
                   }, s1t, s3t),
                 MAD = mapply(FUN = function(x, y) {
                   mean(abs(x - y), na.rm = TRUE)
                   }, s1t, s2t),
                 rank = mapply(FUN = exRank, s1t, s2t),
                 rankMain = mapply(FUN = exRank, s1t, s3t),
                 row.names = NULL)
    ##################
    if (is.null(partition)) {
      glmnetOut <- glmnet::cv.glmnet(x = as.matrix(indFrameTrain[, indices]),
                                     y = trainAccuracyEnv$r, alpha = alpha,
                                     foldid = nfolds)
    } else {
      glmnetOut <- glmnet::cv.glmnet(x = as.matrix(indFrameTrain[, indices]),
                                     y = trainAccuracyEnv$r, alpha = alpha,
                                     foldid = indFrameTrain$partition)
    }
    rTest  <- predict(object = glmnetOut,
                      newx = as.matrix(indFrameTest[, indices]),
                      s = "lambda.min")
    testAccuracyEnv$rEst <- as.numeric(rTest)
    ##################
    testAccuracyGeno <-
      data.frame(Geno = levels(dTest$G),
                 r = mapply(FUN = cor, s1tG, s2tG,
                            MoreArgs = list(use = "na.or.complete",
                                            method = corType)))
  } else {
    testAccuracyEnv  <- NULL
    testAccuracyGeno <- NULL
  }
  ## Create RMSE.
  RMSEtrain <- sqrt(mean((dTrain$Y - predTrain) ^ 2, na.rm = TRUE))
  if (!is.null(testEnv)) {
    RMSEtest  <- sqrt(mean((dTest$Y - predTest) ^ 2, na.rm = TRUE))
  } else {
    RMSEtest <- NULL
  }
  ## Write output to csv.
  resultFilePerEnvTrain <- file.path(outputDir,
                                     paste0(outputFile, "_perEnv_Train.csv"))
  write.csv(trainAccuracyEnv, file = resultFilePerEnvTrain,
            row.names = FALSE, quote = FALSE)
  if (!is.null(testEnv)) {
    resultFilePerEnvTest <- file.path(outputDir,
                                      paste0(outputFile, "_perEnv_Test.csv"))
    write.csv(testAccuracyEnv, file = resultFilePerEnvTest,
              row.names = FALSE, quote = FALSE)
  }
  if (verbose == TRUE) {
    ## Print output to console.
    if (!is.null(testEnv)) {
      cat("\n\n", "Test environments (", traitName, ")", "\n\n")
      print(format(testAccuracyEnv, digits = 2, nsmall = 2))
    }
    cat("\n\n", "Training environments (", traitName,")", "\n\n")
    print(format(trainAccuracyEnv, digits = 2, nsmall = 2))
  }
  ## Create output object.
  out <- list(predTrain = predTrain,
              predTest = predTest,
              resTrain = resTrain,
              resTest = resTest,
              mu = mu,
              envInfoTrain = indicesTrain,
              envInfoTest  = indicesTest,
              parGeno = parGeno,
              main.par = main.par,
              trainAccuracyEnv = trainAccuracyEnv,
              testAccuracyEnv = testAccuracyEnv,
              trainAccuracyGeno = trainAccuracyGeno,
              testAccuracyGeno = testAccuracyGeno,
              lambda = lambda[lambdaIndex],
              lambdaSequence = lambdaSequence,
              RMSEtrain = RMSEtrain,
              RMSEtest = RMSEtest,
              Y = Y,
              G = G,
              E = E,
              indices = indices,
              postLasso = postLasso,
              quadratic = quadratic)
  #dTrain = dTrain, dTest = dTest)
  return(out)
}
