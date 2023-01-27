#' @title
#' Genomic prediction using glmnet, with a genotype-specific penalized
#' regression model.
#'
#' @description
#' \loadmathjax
#' .... These models can be fitted either for the original
#' data, or on the residuals of a model with only main effects.
#'
#' @inheritParams GnE
#'
#' @param useRes Indicates whether the genotype-specific regressions are to be
#' fitted on the residuals of a model with main effects. 0 means no
#' (i.e. fit the regressions on the original data). 1:
#' residuals of a model with environmental main effects. 2 (default): genotypic and
#' environmental main effects.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{predTrain}{Vector with predictions for the training set (to do: Add
#'   the factors genotype and environment; make a dataframe)}
#'   \item{predTest}{Vector with predictions for the test set (to do: Add the
#'   factors genotype and environment; make a dataframe). To do: add estimated
#'   environmental main effects, not only predicted environmental main effects}
#'   \item{mu}{the estimated overall (grand) mean}
#'   \item{envInfoTrain}{The estimated environmental main effects, and the
#'   predicted effects, obtained when the former are regressed on the averaged
#'   indices, using penalized regression.}
#'   \item{envInfoTest}{The predicted environmental main effects for the test
#'   environments, obtained from penalized regression using the estimated
#'   main effects for the training environments and the averaged indices.}
#'   \item{parGeno}{dataframe containing the estimated genotypic
#'   main effects (first column) and sensitivities (subsequent columns)}
#'   \item{testAccuracyEnv}{a \code{data.frame} with the accuracy (r) for each
#'   test environment}
#'   \item{trainAccuracyEnv}{a \code{data.frame} with the accuracy (r) for each
#'   training environment}
#'   \item{trainAccuracyGeno}{a \code{data.frame} with the accuracy (r) for
#'   each genotype, averaged over the training environments}
#'   \item{testAccuracyGeno}{a \code{data.frame} with the accuracy (r) for
#'   each genotype, averaged over the test environments}
#'   \item{RMSEtrain}{The root mean squared error on the training environments}
#'   \item{RMSEtest}{The root mean squared error on the test environments}
#'   \item{Y}{The name of the trait that was predicted, i.e. the column name
#'   in dat that was used}
#'   \item{G}{The genotype label that was used, i.e. the argument G that was
#'   used}
#'   \item{E}{The environment label that was used, i.e. the argument E that was
#'   used}
#'   \item{indices}{The indices that were used, i.e. the argument indices that
#'   was used}
#'   \item{lambdaOpt}{}
#'   \item{pargeno}{}
#'   \item{quadratic}{}
#' }
#'
#' @export
perGeno <- function(dat,
                    Y,
                    G,
                    E,
                    indices = NULL,
                    indicesData = NULL,
                    testEnv = NULL,
                    weight = NULL,
                    useRes = 2,
                    outputFile = "results_glmnet",
                    outputDir = getwd(),
                    corType = c("pearson", "spearman"),
                    partition = data.frame(),
                    nfolds = NULL,
                    scaling = c( "no", "train", "all"),
                    genoAcc = NULL,
                    verbose = FALSE,
                    alpha = 1) {
  scaling <- match.arg(scaling)
  corType <- match.arg(corType)
  ## Get traitName.
  traitName <- ifelse(is.numeric(Y), names(dat)[Y], Y)
  ## Rename data columns for Y, G and E.
  dat <- renameYGE(dat = dat, Y = Y, G = G, E = E)
  ## Get training envs.
  trainEnv <- setdiff(levels(dat$E), testEnv)
  #remove missing values from training env's
  dat <- dat[!(is.na(dat$Y) & dat$E %in% trainEnv),]
  dat <- droplevels(dat)
  if (!is.null(indicesData)) {
    indices <- colnames(indicesData)
    nIndices <- ncol(indicesData)
    stopifnot(all(levels(dat$E) %in% rownames(indicesData)))
    dat <- dat[, setdiff(names(dat), indices)]
    qw <- matrix(nrow = nrow(dat), ncol = nIndices,
                 dimnames = list(NULL, indices))
    dat <- data.frame(dat, qw)
    for (ind in indices) {
      for (env in levels(dat$E)) {
        dat[which(dat$E == env), ind] <- indicesData[env, ind]
      }
    }
  }
  ## Scale environmental variables.
  if (scaling == "train") {
    muTr <- colMeans(dat[dat$E %in% trainEnv, indices])
    sdTr <- sapply(X = dat[dat$E %in% trainEnv, indices], sd)
    dat[, indices] <- scale(dat[, indices], center = muTr, scale = sdTr)
    # for (gn in levels(dTrain$G)) {
    #   muTr <- colMeans(dat[(dat$G == gn) & (dat$E %in% trainEnv), indices])
    #   sdTr <- sapply(X = dat[(dat$G == gn) & (dat$E %in% trainEnv), indices], sd)
    #   dat[dat$G == gn, indices] <- scale(dat[dat$G == gn, indices], center = muTr, scale = sdTr)
    # }
  } else if (scaling == "all") {
    dat[, indices] <- scale(dat[, indices])
  }
  if (is.null(weight)) {
    dat$W <- 1
  } else {
    stopifnot(length(weight) == nrow(dat))
    dat$W <- weight
  }
  trainEnv <- setdiff(levels(dat$E), testEnv)
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
    if (is.null(nfolds)) {
      nfolds <- 10
    }
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
  ## For comparison, fit the unpenalized linear model with main effects only.
  #modelMain <- lm(Y ~ E + G, data = dTrain)
  #cf <- coef(modelMain)
  # modelMain.lme4 <- lme4::lmer(Y ~ E + (1|G), data = dTrain)
  # dim(coef(modelMain.lme4)$G)
  # plot(coef(modelMain.lme4)$G[,1], c(0,cf[-(1:25)]))
  if (useRes == 1) {
    modelMain2 <- lm(Y ~ E, data = dTrain)
    cf2 <- coef(modelMain2)
  }
  # Even if useRes is 0 or 1, mainOnly will still be used for comparison with a
  # main effects only model
  mainOnly <- rep(0, nGenoTrain)
  names(mainOnly) <- levels(dTrain$G)
  mainTemp <- cf[(nEnvTrain + 1):(nEnvTrain + nGenoTrain - 1)]
  names(mainTemp) <- substring(names(mainTemp), first = 2)
  mainOnly[names(mainTemp)] <- mainTemp
  ## Extract from the output the estimated environmental main effects
  if (useRes == 1) {
    envMain <- as.matrix(c(0, cf2[2:nEnvTrain]))
    rownames(envMain) <- levels(dTrain$E)
  } else {
    envMain <- as.matrix(c(0, cf[2:nEnvTrain]))
    rownames(envMain) <- levels(dTrain$E)
  }
  ##########################
  if (useRes == 1) {
    dTrain$yRes <- residuals(modelMain2)
  } else {
    #dTrain$yRes <- residuals(modelMain)
    #
    dTrain$yRes <- dTrain$Y -
      as.numeric(predict(modelMain, newx = mm, thresh = 1e-18, lambda = 0))
  }
  #dTrain$yMainFitted <- fitted(modelMain)
  #dTest$yMainFitted <- predict.lm(object = modelMain, newdata = dTest)

  # qwe = predict.lm(modelMain, newdata = dTrain); sum(abs(dTrain$yMainFitted - qwe))

  ## Define the design matrix for the factorial regression model.
  geFormula <- as.formula(paste0("yRes ~ -1 +",
                                 paste(paste0(indices, ":G"),
                                       collapse = " + ")))
  ## Construct design matrix for training set.
  m <- Matrix::sparse.model.matrix(geFormula, data = dTrain)
  #m[which(is.nan(as.matrix(m)), arr.ind = T)] <- 0
  if (!is.null(testEnv)) {
    ## Construct design matrix for test set.
    geFormula <- as.formula(paste0(" ~ -1 + ",
                                   paste(paste0(indices, ":G"),
                                         collapse = " + ")))
    #mTest[which(is.nan(as.matrix(mTest)), arr.ind = T)] <- 0
    mTest <- Matrix::sparse.model.matrix(geFormula, data = dTest)
  }
  #####################################
  predTrain <- rep(NA, nrow(dTrain))
  names(predTrain) <- rownames(dTrain)
  if (!is.null(testEnv)) {
    predTest <- rep(NA, nrow(dTest))
    names(predTest) <- rownames(dTest)
  } else {
    predTest <- NULL
  }
  sdTr <- muTr <- rep(NA, ncol(m))
  names(sdTr) <- names(muTr) <- colnames(m)
  nIndices <- length(indices)

  ## Create a data.frame that will contain the estimated genotypic
  ## main effects (first column), and the estimated environmental
  ## sensitivities (subsequent) columns.
  parGeno <- matrix(nrow = nGenoTrain, ncol = nIndices + 1,
                    dimnames = list(levels(dTrain$G), c("main", indices)))
  parGeno <- as.data.frame(parGeno)
  lambdaOpt <- rep(NA, nlevels(dTrain$G))
  names(lambdaOpt) <- levels(dTrain$G)
  for (gg in levels(dTrain$G)) {
    # gg <- levels(dTrain$G)[2]
    ggInd <- which(levels(dTrain$G) == gg)
    cat(ggInd / length(levels(dTrain$G)), "\n")
    # wrong!
    #cn <- grep(pattern = gg, x = colnames(m))
    cn <- ggInd + (0:(nIndices - 1)) * nlevels(dTrain$G)
    scn <- (colnames(m))[cn]
    gnInd <- which(dTrain$G == gg)
    #mgg <- as.matrix(m[gnInd, scn, drop = FALSE])
    mgg <- m[gnInd, scn, drop = FALSE]
    mugg <- colMeans(mgg)
    sdgg <- apply(X = mgg, MARGIN = 2, FUN = sd)
    muTr[scn] <- mugg
    sdTr[scn] <- sdgg
    m[gnInd, scn] <- scale(mgg, center = mugg, scale = sdgg)
    #rfOut <- randomForest::randomForest(x = as.matrix(m[gnInd, scn]),
    #                                    y = dTrain$yRes[gnInd])
    if (useRes > 0) {
      yTemp <- dTrain$yRes[gnInd]
    } else {
      yTemp <- dTrain$Y[gnInd]
    }
    #qwe=(as.matrix(m[gnInd, scn, drop = FALSE])); image(cor(qwe))
    glmnetOut <- glmnet::cv.glmnet(x = as.matrix(m[gnInd, scn, drop = FALSE]),
                                   y = yTemp,
                                   weights = dTrain$W[gnInd],
                                   foldid = foldid[gnInd],
                                   nfolds = nfolds,
                                   standardize = TRUE,
                                   intercept = TRUE,
                                   alpha = alpha)
    lambda <- glmnetOut$lambda
    lambdaSequence <- lambda
    lambdaIndex <- which(lambda == glmnetOut$lambda.min)
    cfe <- as.numeric(glmnetOut$glmnet.fit$beta[, lambdaIndex])
    mu <- as.numeric(glmnetOut$glmnet.fit$a0[lambdaIndex])
    parGeno[ggInd, ] <- c(mu, cfe)
    lambdaOpt[ggInd] <- glmnetOut$lambda.min
    #rfOut <- ranger::ranger(x = as.matrix(m[gnInd, scn]),
    #                        y = dTrain$yRes[gnInd])

    # for rf: as.numeric(rfOut$predicted)
    predTrain[gnInd] <- as.numeric(predict(object = glmnetOut,
                                           newx = as.matrix(m[gnInd, scn]),
                                           s = "lambda.min"))
    if (useRes == 2) {
      predTrain[gnInd] <- predTrain[gnInd]  + mainOnly[gg]
    }
    #######  !!!!!!!!! #######

    # labels predTest/Tr correct ?
    if (!is.null(testEnv) & (gg %in% as.character(dTest$G))) {
      mggTest <- as.matrix(mTest[which(dTest$G == gg), scn, drop = FALSE])
      mTest[which(dTest$G == gg), scn] <-
        scale(mggTest, center = mugg, scale = sdgg)
      #rfPred  <- predict(object = rfOut,
      #                   newdata = mTest[which(dTest$G == gg), scn])
      ##rfPred  <- predict(object = rfOut, data = mTest[which(dTest$G == gg), scn])

      #predTest[which(dTest$G == gg)] <- rfPred + mainOnly[gg]
      glmnetPred  <- as.numeric(predict(object = glmnetOut,
                                        newx = mTest[dTest$G == gg,
                                                     scn, drop = FALSE],
                                        s = "lambda.min"))
      predTest[which(dTest$G == gg)] <- glmnetPred
      if (useRes == 2) {
        predTest[which(dTest$G == gg)] <-
          predTest[which(dTest$G == gg)] + mainOnly[gg]
      }
    }
  }
  ## Compute the mean of each environmental index, in each environment
  indFrame <- aggregate(dat[, indices], by = list(E = dat$E), FUN = mean)
  rownames(indFrame) <- indFrame$E
  indFrameTrain <- indFrame[indFrame$E %in% trainEnv, ]
  if (!is.null(testEnv)) {
    indFrameTest  <- indFrame[indFrame$E %in% testEnv, ]
  }
  indFrameTrain <- merge(indFrameTrain, envMain, by.x = "E", by.y = "row.names")
  colnames(indFrameTrain)[ncol(indFrameTrain)] <- "envMainFitted"
  indFrameTrain <- merge(indFrameTrain, partition)
  rownames(indFrameTrain) <- indFrameTrain$E
  glmnetOut <- glmnet::cv.glmnet(x = as.matrix(indFrameTrain[, indices]),
                                 y = indFrameTrain$envMainFitted,
                                 alpha = alpha,
                                 foldid = indFrameTrain$partition)
  parEnvTrain <- predict(object = glmnetOut,
                         newx = as.matrix(indFrameTrain[, indices]),
                         s = "lambda.min")
  if (!is.null(testEnv)) {
    parEnvTest  <- predict(object = glmnetOut,
                           newx = as.matrix(indFrameTest[, indices]),
                           s = "lambda.min")
  }
  indicesTest <- NULL
  #if (useRes == TRUE) {
  ## add the predicted environmental main effects:
  predTrain <- predTrain + as.numeric(parEnvTrain[as.character(dTrain$E), ])
  indicesTrain <- data.frame(indFrameTrain,
                             envMainPred = as.numeric(parEnvTrain))
  if (!is.null(testEnv)) {
    ## add the predicted environmental main effects:
    predTest <- predTest + as.numeric(parEnvTest[as.character(dTest$E), ])
    indicesTest <- data.frame(indFrameTest,
                              envMainPred = as.numeric(parEnvTest))
  }
  #}
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
  ## Compute statistics for training data.
  predMain <- mainOnly[as.character(dTrain$G)]
  trainAccuracyEnv <- getAccuracyEnv(datNew = dTrain[, "Y"],
                                     datPred = predTrain,
                                     datPredMain = predMain,
                                     datE = dTrain[, "E"],
                                     corType = corType,
                                     rank = TRUE)
  ## Compute accuracies for genotypes.
  trainAccuracyGeno <- getAccuracyGeno(datNew = dTrain[, "Y"],
                                       datPred = predTrain,
                                       datG = dTrain[, "G"],
                                       corType = corType)
  if (!is.null(testEnv)) {
    ## Compute statistics for test data.
    predMainTest <- mainOnly[as.character(dTest$G)]
    testAccuracyEnv <- getAccuracyEnv(datNew = dTest[, "Y"],
                                      datPred = predTest,
                                      datPredMain = predMainTest,
                                      datE = dTest[, "E"],
                                      corType = corType,
                                      rank = TRUE)
    ## Compute accuracies for genotypes.
    testAccuracyGeno <- getAccuracyGeno(datNew = dTest[, "Y"],
                                        datPred = predTest,
                                        datG = dTest[, "G"],
                                        corType = corType)
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
  if (verbose) {
    ## Print output to console.
    if (!is.null(testEnv)) {
      cat("\n\n", "Test environments (", traitName, ")", "\n\n")
      print(format(testAccuracyEnv, digits = 2, nsmall = 2))
    }
    cat("\n\n", "Training environments (", traitName,")", "\n\n")
    print(format(trainAccuracyEnv, digits = 2, nsmall = 2))
  }
  # needed for use in nGnE
  quadratic <- FALSE
  # to do: quadratic terms in perGeno
  ## Create output object.
  out <- list(predTrain = predTrain,
              predTest = predTest,
              envInfoTrain = indicesTrain,
              envInfoTest  = indicesTest,
              testAccuracyEnv = testAccuracyEnv,
              trainAccuracyEnv = trainAccuracyEnv,
              trainAccuracyGeno = trainAccuracyGeno,
              testAccuracyGeno = testAccuracyGeno,
              RMSEtrain = RMSEtrain,
              RMSEtest = RMSEtest,
              Y = Y,
              G = G,
              E = E,
              indices = indices,
              genotypes = levels(dTrain$G),
              lambdaOpt = lambdaOpt,
              parGeno = parGeno,
              quadratic = quadratic)
  return(out)
}
