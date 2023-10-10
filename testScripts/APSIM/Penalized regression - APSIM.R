### Required packages and functions ###
library(factReg)
library(BGLR)
library(tictoc)


# Help function to calculate correlations
corType <- "pearson"
corFun <- function(x, y, corType) {
  cor(x, y, use = "na.or.complete", method = corType)
}

# Function to calculate accuracies per environment
calcAcc <- function(yobs, ypred, env) {
  sObs <- split(yobs, env)
  sPred <- split(ypred, env)
  accuracy <- data.frame(
    Env = levels(env),
    r = mapply(FUN = corFun, sObs, sPred, MoreArgs = list(corType = corType))
  )
  accuracy
}

# Help function to calculate RMSE
rmseFun <- function(x, y, na.rm = TRUE) {
  sqrt(mean((x - y)^2, na.rm = na.rm))
}

# Function to calculate RMSE per environment
calcRMSE <- function(yobs, ypred, env) {
  sObs <- split(yobs, env)
  sPred <- split(ypred, env)
  RMSE <- data.frame(
    Env = levels(env),
    RMSE = mapply(FUN = rmseFun, sObs, sPred)
  )
  RMSE
}

### Load data and subset to type "GE" and "GnE" ####
load("../factorialregression/APSIM data/apsim.RData")

dat.apsim <- droplevels(all_apsim[all_apsim$type %in% c("GE", "GnE"), ])

# Create training and test set ####
train <- droplevels(all_apsim[all_apsim$type %in% "GE", ])
test <- droplevels(all_apsim[all_apsim$type %in% "GnE", ])

# Define testing and training environments.
testenv <- levels(test$Env)
trainenv <- levels(train$Env)

# Define testing and training genotypes.
testgen <- levels(test$geno)
traingen <- levels(train$geno)

# Define indices
ind <- colnames(dat.apsim)[c(5, 7:62)]

### Load the apsim data from the simulation (without noise) ####
# apsim_noiseless <- read.table("R/unused/data/Cov_windows_wide_248envs_199geno.txt",
#                               header = TRUE)
# apsim_noiseless$Env <- as.factor(apsim_noiseless$Env)
# apsim_noiseless$geno <- as.factor(apsim_noiseless$geno)

# Create another training set and GnE (without noise). Only used for model validation! ####
train_2 <- droplevels(GE_apsim[GE_apsim$Env %in% train$Env, ])
train_2 <- droplevels(GE_apsim[GE_apsim$geno %in% train$geno, ])
GnE_apsim_noiseless <- droplevels(GE_apsim[GE_apsim$geno %in% train$geno, ])


### Penalized factorial regression ####
# Additive model, only main effects (set the penalty parameter to a large value)
set.seed(1234)
Additive_model <- GnE(dat.apsim, Y = "yld_noise", lambda = 100000,
                      G = "geno", E = "Env", testEnv = testenv,
                      indices = ind, penG = FALSE, penE = FALSE,
                      alpha = 0.5, scaling = "train")

# Full model, no penalization (set the penalty parameter to zero).
set.seed(1234)
Full_model <- GnE(dat.apsim, Y = "yld_noise", lambda = 0,
                  G = "geno", E = "Env", testEnv = testenv,
                  indices = ind, penG = FALSE, penE = FALSE,
                  alpha = 0.5, scaling = "train")

# Elastic Net model, set alpha parameter to 0.5.
set.seed(1234)
tic()
Elnet_model <- GnE(dat.apsim, Y = "yld_noise", lambda = NULL,
                   G = "geno", E = "Env", testEnv = testenv,
                   indices = ind, penG = FALSE, penE = FALSE,
                   alpha = 0.5, scaling = "train")
toc()

# Lasso model, set alpha parameter to 1.
set.seed(1234)
tic()
Lasso_model <- GnE(dat.apsim, Y = "yld_noise", lambda = NULL,
                   G = "geno", E = "Env", testEnv = testenv,
                   indices = ind, penG = FALSE, penE = FALSE,
                   alpha = 1, scaling = "train")
toc()

# Ridge model, set alpha parameter to 0.
set.seed(1234)
tic()
Ridge_model <- GnE(dat.apsim, Y = "yld_noise", lambda = NULL,
                   G = "geno", E = "Env", testEnv = testenv,
                   indices = ind, penG = FALSE, penE = FALSE,
                   alpha = 0, scaling = "train")

toc()

### Get accuracies for training set manually ####
train.obs <- train_2$yld
cor.additive <- calcAcc(train.obs, Additive_model$predTrain$pred, train_2$Env)
RMSE.additive <- calcRMSE(train.obs, Additive_model$predTrain$pred, train_2$Env)
Ov.additive <- rmseFun(train.obs,Additive_model$predTrain$pred)

cor.full <- calcAcc(train.obs, Full_model$predTrain$pred, train_2$Env)
RMSE.full <- calcRMSE(train.obs, Full_model$predTrain$pred, train_2$Env)
Ov.full <- rmseFun(train.obs,Full_model$predTrain$pred)

cor.elnet <- calcAcc(train.obs, Elnet_model$predTrain$pred, train_2$Env)
RMSE.elnet <- calcRMSE(train.obs, Elnet_model$predTrain$pred, train_2$Env)
Ov.elnet <- rmseFun(train.obs, Elnet_model$predTrain$pred)

cor.lasso <- calcAcc(train.obs, Lasso_model$predTrain$pred, train_2$Env)
RMSE.lasso <- calcRMSE(train.obs, Lasso_model$predTrain$pred, train_2$Env)
Ov.lasso <- rmseFun(train.obs, Lasso_model$predTrain$pred)

cor.ridge <- calcAcc(train.obs, Ridge_model$predTrain$pred, train_2$Env)
RMSE.ridge <- calcRMSE(train.obs, Ridge_model$predTrain$pred, train_2$Env)
Ov.ridge <- rmseFun(train.obs, Ridge_model$predTrain$pred)

### Results from factreg methods ####
APCOR_test_results <- c(mean(Additive_model$testAccuracyEnv$r),
                        mean(Full_model$testAccuracyEnv$r),
                        mean(Ridge_model$testAccuracyEnv$r),
                        mean(Lasso_model$testAccuracyEnv$r),
                        mean(Elnet_model$testAccuracyEnv$r))

APCOR_train_results <- c(mean(cor.additive$r),
                         mean(cor.full$r),
                         mean(cor.ridge$r),
                         mean(cor.lasso$r),
                         mean(cor.elnet$r))

ARMSE_test_results <- c(mean(Additive_model$testAccuracyEnv$RMSE),
                        mean(Full_model$testAccuracyEnv$RMSE),
                        mean(Ridge_model$testAccuracyEnv$RMSE),
                        mean(Lasso_model$testAccuracyEnv$RMSE),
                        mean(Elnet_model$testAccuracyEnv$RMSE))

ARMSE_train_results <- c(mean(RMSE.additive$RMSE),
                         mean(RMSE.full$RMSE),
                         mean(RMSE.ridge$RMSE),
                         mean(RMSE.lasso$RMSE),
                         mean(RMSE.elnet$RMSE))

results <- cbind(APCOR_test_results,
                 APCOR_train_results,
                 ARMSE_test_results,
                 ARMSE_train_results)

rownames(results) <- c("Additive",
                       "Full",
                       "Ridge",
                       "Lasso",
                       "El-Net")

colnames(results) <- c("APCOR_Env Test",
                       "APCOR_Env Train",
                       "ARMSE_Env Test",
                       "ARMSE_Env Train")

View(round(results, 3))

### BRR model ####
# Boolean indicator:
IsTrainSet <- dat.apsim$Env %in% trainenv
IsTestSet <- !IsTrainSet

# First standardize and centre the covariates:
indices <- ind
muTr <- colMeans(dat.apsim[IsTrainSet, indices])
sdTr <- pmax(sapply(X = dat.apsim[IsTrainSet, indices], sd), 1e-8)
dat.apsim[, ind] <- scale(dat.apsim[, ind], center = muTr, scale = sdTr)

# Make the design matrix for cov*G interactions:
fText <- paste0(ind, ":geno", collapse = "+")
f <- as.formula(paste0("~", fText))
Z <- model.matrix(f, data = dat.apsim)
Z <- Z[,-1]
q <- ncol(Z)
dim(Z)
dat <- data.frame(Z, dat.apsim)
dat1 <- dat[IsTrainSet, ]
dat1 <- droplevels(dat1)

# Set the argument for BGLR function
nIter <- 10000
burnIn <- 0
thin <- 5
saveAt <- ""
S0 <- NULL;
weights <- NULL;
R2 <- 0.5;

# Setting the linear predictor
ETA <- list( list(~factor(Env),
                  data = dat.apsim, model = "FIXED"),
             list(~factor(geno),
                  data = dat.apsim, model = "FIXED"),
             list(X = Z, model = "BRR")
)

title <- paste("Computation time BRR", nIter) ##Computation time: 2330.52 sec elapsed

y <- dat$yld_noise
y[!IsTrainSet] <- NA

tic(title)
BRR_model <- BGLR(y = y, ETA = ETA, nIter = nIter, burnIn = burnIn, thin = thin,
                  saveAt = saveAt, df0 = 5, S0 = S0, weights = weights, R2 = R2,
                  verbose = FALSE)
toc()


# this is a trick to remove environmental effects from predictions
# first level is equal to zero:
env_eff_corr <- c(0, BRR_model$ETA[[1]]$b)
env_eff_corr[(nlevels(dat$Env) %in% trainenv)] <- 0
env_eff_corr <- as.vector(env_eff_corr)

Xenv <- model.matrix(~Env - 1, data = dat)

ypred_BGLR <- BRR_model$yHat - Xenv %*% env_eff_corr

yobs <- dat$yld_noise
env <- dat$Env

# Get accuracies for test and training set BRR:
acc.bglr <- calcAcc(yobs, ypred_BGLR, env)
acc.bglr.train <- acc.bglr[acc.bglr$Env %in% trainenv,]
acc.bglr.test <- acc.bglr[acc.bglr$Env %in% testenv,]

### Jarquin model ####
# apsim_K_split <- apsim_K[rownames(apsim_K) %in% traingen,
#                          rownames(apsim_K) %in% traingen]
# Jarquin_model <- GnE_BGLR(dat = dat.apsim,
#                           Y = "yld",
#                           G = "geno",
#                           E = "Env",
#                           K = apsim_K_split,
#                           indices = ind,
#                           indicesData = NULL,
#                           testEnv = testenv,
#                           corType = "pearson",
#                           scaling = "train",
#                           nIter = 10000,
#                           burnIn = 0,
#                           thin = 5)

