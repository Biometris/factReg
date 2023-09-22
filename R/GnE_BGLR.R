#' @title
#' here comes the title
#'
#' @description
#' \loadmathjax
#' here comes the description
#'
#' @inheritParams GnE
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
  # to do: checks
  redundantGeno <- names(which(table(dat$G[!is.na(dat$Y) &
                                             dat$E %in% trainEnv]) < 10))
  if (length(redundantGeno) > 0) {
    warning("the following genotypes have < 10 observations, and are removed:",
            "\n\n", paste(redundantGeno, collapse = ", "), "\n")
    dat <- dat[!(dat$G %in% redundantGeno), ]
    if (nrow(dat) == 0) {
      stop("No data left after removing genotypes with < 10 observations.\n")
    }
  }
  if (!is.null(indicesData)) {
    indices <- colnames(indicesData)
    ## Remove columns from dat that are also in indices and then merge indices.
    dat <- dat[, setdiff(names(dat), indices)]
    dat <- merge(dat, indicesData, by.x = "E", by.y = "row.names")
    covECJarquin <- tcrossprod(as.matrix(indicesData))
    covECJarquin2 <- covECJarquin[dat$E, dat$E]
  } else {
    covECJarquin <- as.matrix(dat[, indices])
    covECJarquin2 <- tcrossprod(covECJarquin)
  }
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
  K <- K[levels(dat$G), levels(dat$G)]
  ## Split dat into training and test set.
  dTrain <- dat[dat$E %in% trainEnv, ]
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
  ### Create the design matrix for environments
  Ze = matrix(0, length(dat$Y), nEnv)
  for (i in 1:nrow(Ze)) {
    Ze[i, dat$E[i]] <- 1
  }
  ########################
  # 2/ Creating the GW and GE matrices
  ########################
  ## creating the matrix for G
  K2 <- K[dat$G, dat$G]
  ## creating the matrix for GW
  cov_GWJarquin2 <- K2 * covECJarquin2
  # creating the matrix for GE
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
  # Model : EG_GxW_GxE
  ETA_Jarquin <- list(list(~factor(dat$E), model = "BRR"),
                      list(V = EVD_K2$vectors, d = EVD_K2$values,
                           model = "RKHS"),
                      list(V = EVD_cov_GWJarquin2$vectors,
                           d = EVD_cov_GWJarquin2$values, model = "RKHS"),
                      list(V = EVD_cov_GE$vectors, dat = EVD_cov_GE$values,
                           model = "RKHS"))

  ### Run BGLR
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
  ## Compute statistics for training data.
  trainAccuracyEnv <-
    data.frame(Env = levels(dTrain$E),
               r = mapply(FUN = corFun, s1, s2,
                          MoreArgs = list(corType = corType)),
               RMSE = mapply(FUN = RMSEFun, s1, s2),
               MAD = mapply(FUN = MADFun, s1, s2),
               row.names = NULL)
  if (!is.null(testEnv)) {
    ## Split raw data and predictions by environment for computing statistics.
    s1t <- split(dTest$Y, dTest$E)
    s2t <- split(predTest, dTest$E)
    ## Compute statistics for test data.
    testAccuracyEnv <-
      data.frame(Env = levels(dTest$E),
                 r = mapply(FUN = corFun, s1t, s2t,
                            MoreArgs = list(corType = corType)),
                 RMSE = mapply(FUN = RMSEFun, s1t, s2t),
                 MAD = mapply(FUN = MADFun, s1t, s2t),
                 row.names = NULL)
  }
  ## Create RMSE.
  RMSEtrain <- sqrt(mean((dTrain$Y - predTrain) ^ 2, na.rm = TRUE))
  if (!is.null(testEnv)) {
    RMSEtest <- sqrt(mean((dTest$Y - predTest) ^ 2, na.rm = TRUE))
  } else {
    RMSEtest <- NULL
  }
  ## Create output object.
  out <- list(predTrain = predTrain,
              predTest = predTest,
              trainAccuracyEnv = trainAccuracyEnv,
              testAccuracyEnv = testAccuracyEnv,
              RMSEtrain = RMSEtrain,
              RMSEtest = RMSEtest,
              Y = Y,
              G = G,
              E = E,
              indices = indices)
  return(out)

}
