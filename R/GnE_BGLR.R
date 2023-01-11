#' @title
#' here comes the title
#'
#' @description
#' \loadmathjax
#' here comes the description
#'
#' @inheritParams GnE_glmnet
#'
#' @param K ...
#' @param nIter ...
#' @param burnIn ...
#' @param thin ...
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{yHat}{...}
#' }
#'
#' @export
GnE_BGLR <- function(dat,
                     Y,
                     G,
                     E,
                     K,
                     indices = NULL,
                     indicesData = NULL,
                     testEnv = NULL,
                     corType = c("pearson", "spearman"),
                     scaling = c("train", "all", "no"),
                     nIter = 50,
                     burnIn = 10,
                     thin = 1) {

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

  ##################

  if (!is.null(indicesData)) {

    indices <- colnames(indicesData)
    nIndices <- ncol(indicesData)

    stopifnot(all(levels(dat$E) %in% rownames(indicesData)))
    dat <- dat[, setdiff(names(dat), indices)]

    for (ind in indices) {                             # scale the index
      indicesData[,ind] <- as.numeric(scale(indicesData[, ind]))
    }

    qw <- matrix(NA, nrow(dat), nIndices)
    colnames(qw) <- indices
    dat <- data.frame(dat, qw)

    for (ind in indices) {
      for (env in levels(dat$E)) {
        dat[which(dat$E == env), ind] <- indicesData[env, ind]
      }
    }

    indicesData <- as.matrix(indicesData)
    covECJarquin <- indicesData %*% t(indicesData)
    covECJarquin2 <- covECJarquin[dat$E, dat$E]
  } else {

    covECJarquin <- as.matrix(dat[, indices])
    covECJarquin2 <- tcrossprod(covECJarquin)
  }
  ###################

  ## Scale environmental variables.
  if (scaling == "train") {
    muTr <- colMeans(dat[dat$E %in% trainEnv, indices])
    sdTr <- sapply(X = dat[dat$E %in% trainEnv, indices], FUN = sd)
    dat[, indices] <- scale(dat[, indices], center = muTr, scale = sdTr)
  } else if (scaling == "all") {
    dat[, indices] <- scale(dat[, indices])
  }

  nEnv <- nlevels(dat$E)
  nGeno <- nlevels(dat$G)

  # to do: checks
  K <- K[levels(dat$G), levels(dat$G)]

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

  ########################
  # 1/ Design matrices
  ########################

  # ### Create the design matrix for genotypes
  # Z <- matrix(0, length(dat$Y), nGeno)
  # for (i in 1:nrow(Z)) {
  #   Z[i, dat$G[i]] <- 1
  # }

  ### Create the design matrix for environments
  Ze = matrix(0, length(dat$Y), nEnv)
  for (i in 1:nrow(Ze)) {
    Ze[i, dat$E[i]] <- 1
  }

  ########################
  # 2/ Creating the GW and GE matrices
  ########################

  ## creating the matrix for G
  #K2 <- kronecker(X = matrix(1, nEnv, nEnv), Y = K)
  K2 <- K[dat$G, dat$G]

  ## creating the matrix for GW
  #cov_GWJarquin2 <- Z %*% K %*% t(Z) * covECJarquin2
  cov_GWJarquin2 <- K2 * covECJarquin2

  # creating the matrix for GE
  #cov_GE <- Z %*% K %*% t(Z) * Ze %*% t(Ze)
  cov_GE <- K2 * Ze %*% t(Ze)

  ########################
  # 3/ Decompose into eigenvalues (long!)
  ########################

  EVD_cov_GWJarquin2 <- eigen(cov_GWJarquin2)

  EVD_cov_GE <- eigen(cov_GE)

  EVD_K2 <- eigen(K2)

  ########################
  # 4/ Writing the model with the BGLR format
  ########################
  #
  # Model : EG_GxW_GxE
  ETA_Jarquin <- list(list(~factor(dat$E),
                           model = "BRR"),
                      list(V = EVD_K2$vectors,
                           d = EVD_K2$values,
                           model = "RKHS"),
                      list(V = EVD_cov_GWJarquin2$vectors,
                           d = EVD_cov_GWJarquin2$values,
                           model = "RKHS"),
                      list(V = EVD_cov_GE$vectors,
                           dat = EVD_cov_GE$values,
                           model = "RKHS")
                      )

  ### Run BGLR
  set.seed(1)

  dat$Y[dat$E %in% testEnv] <- NA

  res <- BGLR::BGLR(y = dat$Y,
                    response_type = "gaussian",
                    ETA = ETA_Jarquin,
                    nIter = nIter,  # low value for test, actual value was 50000
                    burnIn = burnIn, # low value for test, actual value was 10000
                    thin = thin,
                    saveAt = "essai",
                    S0 = NULL,
                    df0 =5,
                    R2 = 0.5,
                    weights = NULL,
                    verbose = FALSE,
                    rmExistingFiles = TRUE,
                    groups = NULL)

  ### Get the predictions
  yHat <- res$yHat

  predTrain <- yHat[!(dat$E %in% testEnv)]
  if (!is.null(testEnv)) {
    predTest <- yHat[dat$E %in% testEnv]
  }

  ## Split raw data and predictions by environment for computing statistics.
  s1 <- split(dTrain$Y, dTrain$E)
  s2 <- split(predTrain, dTrain$E)

  s1G <- split(dTrain$Y, dTrain$G)
  s2G <- split(predTrain, dTrain$G)

  ## Compute statistics for training data.
  trainAccuracyEnv <-
    data.frame(Env = levels(dTrain$E),
               r = mapply(FUN = cor, s1, s2,
                          MoreArgs = list(use = "na.or.complete",
                                          method = corType)),
               RMSE = mapply(FUN = function(x, y) {sqrt(mean((x - y) ^ 2, na.rm = TRUE))},
                             s1, s2),
               MAD = mapply(FUN = function(x, y) {mean(abs(x - y), na.rm = TRUE)},
                            s1, s2),
               rank = mapply(FUN = exRank,
                             s1, s2),
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

    ## Compute statistics for test data.
    testAccuracyEnv <-
      data.frame(Env = levels(dTest$E),
                 r = mapply(FUN = cor, s1t, s2t,
                            MoreArgs = list(use = "na.or.complete",
                                            method = corType)),
                 RMSE = mapply(FUN = function(x, y) {sqrt(mean((x - y) ^ 2, na.rm = TRUE))},
                               s1t, s2t),
                 MAD = mapply(FUN = function(x, y) {mean(abs(x - y), na.rm = TRUE)},
                              s1t, s2t),
                 rank = mapply(FUN = exRank,
                               s1t, s2t),
                 row.names = NULL)
  }
  ## Create RMSE.
  RMSEtrain <- sqrt(mean((dTrain$Y - predTrain) ^ 2, na.rm = TRUE))
  if (!is.null(testEnv)) {
    RMSEtest  <- sqrt(mean((dTest$Y - predTest) ^ 2, na.rm = TRUE))
  } else {
    RMSEtest <- NULL
  }

  ## Create output object.
  out <- list(predTrain = predTrain, predTest = predTest,
              trainAccuracyEnv = trainAccuracyEnv,
              testAccuracyEnv = testAccuracyEnv,
              #trainAccuracyGeno = trainAccuracyGeno,
              #testAccuracyGeno = testAccuracyGeno,
              RMSEtrain = RMSEtrain,
              RMSEtest = RMSEtest,
              Y = Y, G = G, E = E,
              indices = indices)
  return(out)

}
