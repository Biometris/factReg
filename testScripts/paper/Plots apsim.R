#### Figures

library(factReg)
library(ggplot2)

corType <- "pearson"
# Help function to calculate correlations
corFun <- function(x, y, corType) {
  cor(x, y, use = "na.or.complete", method = corType)
}

# Calculate accuracies per environment
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

## Read Aalt-Jan results on deep learning, testing results
## GnE2, GnE3, GnE4 and GnE6 contain results for ML models W1 to W4.
yhats1 <- read.csv("testScripts/paper/apsim_Oct2023/GnE2_test_exp_pred2")
yhats2 <- read.csv("testScripts/paper/apsim_Oct2023/GnE3_test_exp_pred2")
yhats3 <- read.csv("testScripts/paper/apsim_Oct2023/GnE5_test_exp_pred2")
yhats4 <- read.csv("testScripts/paper/apsim_Oct2023/GnE6_test_exp_pred2")

## Load APSIM data.
load("data-raw/apsim.RData")
GEGnE_apsim <- droplevels(all_apsim[all_apsim$type %in% c("GE", "GnE"), ])

## Column X0 in results contains inputs, column X0.1 model predictions.
## Combine everything in one data.frame.
DLRes <- data.frame(Env = GnE_apsim$Env, yld = GnE_apsim$yld,
                    yld_noise = yhats1$X0, W1 = yhats1$X0.1,
                    W2 = yhats2$X0.1, W3 = yhats3$X0.1, W4 = yhats4$X0.1)

## Compute APCOR per environment for each of the DL models.
DLAPCOR <- sapply(X = DLRes[c("W1", "W2", "W3", "W4")], FUN = function(W) {
  calcAcc(W, DLRes$yld, DLRes$Env)$r
})
## Compute RMSE per environment for each of the DL models.
DLRMSE <- sapply(X = DLRes[c("W1", "W2", "W3", "W4")], FUN = function(W) {
  calcRMSE(W, DLRes$yld, DLRes$Env)$RMSE
})
## Compute overall RMSE for each of the DL models.
DLRMSE <- sapply(X = DLRes[c("W1", "W2", "W3", "W4")], FUN = function(W) {
  rmseFun(W, DLRes$yld)
})

## Compute average APCOR per environment for each of the models.
DL_APCORTest <- colMeans(DLAPCOR)
## Compute average RMSE per environment for each of the models.
DL_RMSETest <- colMeans(DLRMSE)

## Add the Deep learning to factreg results for plotting.
testenv <- levels(GnE_apsim$Env)
trainenv <- levels(GE_apsim$Env)
ind <- colnames(all_apsim)[c(5, 7:62)]

set.seed(1234)
Additive_model <- GnE(GEGnE_apsim, Y = "yld_noise", lambda = 100000,
                      G = "geno", E = "Env", testEnv = testenv,
                      indices = ind, penG = FALSE, penE = FALSE,
                      alpha = 0.5, scaling = "train")

set.seed(1234)
Full_model <- GnE(GEGnE_apsim, Y = "yld_noise", lambda = 0,
                  G = "geno", E = "Env", testEnv = testenv,
                  indices = ind, penG = FALSE, penE = FALSE,
                  alpha = 0.5, scaling = "train")

set.seed(1234)
Lasso_model <- GnE(GEGnE_apsim, Y = "yld_noise", lambda = NULL,
                   G = "geno", E = "Env", testEnv = testenv,
                   indices = ind, penG = FALSE, penE = FALSE,
                   partition = NULL, alpha = 1, scaling = "train")

## Combine all results.
APCORTestTot <- data.frame(Env = Additive_model$testAccuracyEnv$Env,
                           add = Additive_model$testAccuracyEnv$r,
                           full = Full_model$testAccuracyEnv$r,
                           lasso = Lasso_model$testAccuracyEnv$r,
                           dl = DLAPCOR[, "W3"])
## Add year and environment type.
APCORTestTot <- merge(APCORTestTot,
                      unique(GEGnE_apsim[c("Env", "year", "groups")]))

## Order withing groups by increasing mean accuracy over all models.
APCORTestTot$meanAcc <- rowMeans(APCORTestTot[c("add", "full", "lasso", "dl")])

## First split by ET, order and then combine again.
APCORTestTotLst <- split(APCORTestTot, APCORTestTot$groups)
APCORTestTotLst <- lapply(X = APCORTestTotLst, FUN = function(dat) {
  dat[order(dat$meanAcc), ]
})
APCORTestTotSort <- do.call(rbind, APCORTestTotLst)

colnames(APCORTestTotSort) <- c("Env", "Additive", "Full","Lasso",
                                "Deep learning", "Year", "Group", "MeanAcc")


## The actual plot.
## Accuracies per environment are plotted for the four models.
par(mar = c(6, 6, 1, 6))

plot(APCORTestTotSort$Additive, type = "b", xaxt = "n", xlab = "",
     ylim = c(0, 1.05), ylab = "Accuracy", col = "blue", pch = 19, lwd = 2,
     main = "")
lines(APCORTestTotSort$Full, type = "b", pch = 19, lwd = 2, col = "green")
lines(APCORTestTotSort$Lasso, type = "b", pch = 19, lwd = 2, col = "red")
lines(APCORTestTotSort$`Deep learning`, type = "b", pch = 19, lwd = 2, col = "black")
abline(v = 1:40, lty = 3, col = "grey")
axis(1, at = seq_along(APCORTestTotSort$Env), APCORTestTotSort$Env, las = 3,
     cex.axis = 0.7, labels = FALSE)
text(x = seq_along(APCORTestTotSort$Env), labels = APCORTestTotSort$Env,
     xpd = NA, srt = 90, y = -0.19, cex = 0.55)

text(x = 4.5, labels = "ET1", y = 1, cex = 0.95)
abline(v = 9.5)
text(x = 13.5, labels = "ET2", y = 1, cex = 0.9)
abline(v = 17.5)
text(x = 23.5, labels = "ET3", y = 1, cex =0.9)
abline(v = 29.5)
text(x = 35.5, labels= "ET4", y = 1, cex = 0.9)

## Box plot with the accuracies per ET for each of the methods.
my_colors <- c("blue", "black", "green", "red")

## ggplot requires data in a long format, so pivot first.
APCORTestTotSortLong <- tidyr::pivot_longer(APCORTestTotSort,
                                            cols = -c(Env, Year, Group, MeanAcc),
                                            names_to = "Model", values_to = "r")

ggplot(APCORTestTotSortLong, aes(x = Model , y = r, color = Model)) +
  geom_boxplot() +
  scale_color_manual(values = my_colors) +
  theme() +
  ylim(0, 1) +
  theme_bw() +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank())+
  facet_wrap(~Group)
ggsave("testScripts/paper/plots/TestingBoxplotApsim (update 8 Dec).jpeg",
       width = 9, height = 6)


## Similar for training data.

##Aalt-Jan results on deep learning, training results

## Read Aalt-Jan results on deep learning, testing results
## GnE2, GnE3, GnE4 and GnE6 contain results for ML models W1 to W4.
yhats_train1 <- read.csv("testScripts/paper/apsim_Oct2023/GnE2_train_exp_pred2")
yhats_train2 <- read.csv("testScripts/paper/apsim_Oct2023/GnE3_train_exp_pred2")
yhats_train3 <- read.csv("testScripts/paper/apsim_Oct2023/GnE5_train_exp_pred2")
yhats_train4 <- read.csv("testScripts/paper/apsim_Oct2023/GnE6_train_exp_pred2")

## Column X0 in results contains inputs, column X0.1 model predictions.
## Combine everything in one data.frame.
DLRes_train <- data.frame(Env = GE_apsim$Env, yld = GE_apsim$yld,
                          yld_noise = yhats_train1$X0, W1 = yhats_train1$X0.1,
                          W2 = yhats_train2$X0.1, W3 = yhats_train3$X0.1,
                          W4 = yhats_train4$X0.1)

## Compute APCOR per environment for each of the DL models.
DLAPCOR_train <- sapply(X = DLRes_train[c("W1", "W2", "W3", "W4")],
                        FUN = function(W) {
                          calcAcc(W, DLRes_train$yld, DLRes_train$Env)$r
                        })

## Compute average APCOR per environment for each of the models.
DL_APCORTrain <- colMeans(DLAPCOR_train)

## Combine all results.
train.obs <- GE_apsim$yld
train.env <- GE_apsim$Env

## The accuracies have to be computed based on the noiseless yield
## This means they cannot be taken straight from the output.
cor.additive <- calcAcc(train.obs, Additive_model$predTrain$pred, train.env)
cor.full <- calcAcc(train.obs, Full_model$predTrain$pred, train.env)
cor.lasso <- calcAcc(train.obs, Lasso_model$predTrain$pred, train.env)

APCORTrainTot <- data.frame(Env = Additive_model$trainAccuracyEnv$Env,
                            add = cor.additive$r,
                            full = cor.full$r,
                            lasso = cor.lasso$r,
                            dl = DLAPCOR_train[, "W3"])
## Add year and environment type.
APCORTrainTot <- merge(APCORTrainTot,
                       unique(GEGnE_apsim[c("Env", "year", "groups")]))

## Order withing groups by increasing mean accuracy over all models.
APCORTrainTot$meanAcc <- rowMeans(APCORTrainTot[c("add", "full", "lasso", "dl")])

## First split by ET, order and then combine again.
APCORTrainTotLst <- split(APCORTrainTot, APCORTrainTot$groups)
APCORTrainTotLst <- lapply(X = APCORTrainTotLst, FUN = function(dat) {
  dat[order(dat$meanAcc), ]
})
APCORTrainTotSort <- do.call(rbind, APCORTrainTotLst)

colnames(APCORTrainTotSort) <- c("Env", "Additive", "Full","Lasso",
                                 "Deep learning", "Year", "Group", "MeanAcc")

## The actual plot.
## Accuracies per environment are plotted for the four models.
par(mar = c(6, 6, 1, 6))

plot(APCORTrainTotSort$Additive, type = "b", xaxt = "n", xlab = "",
     ylim = c(0, 1.05), ylab = "Accuracy", col = "blue", pch = 19, lwd = 2,
     main = "")
lines(APCORTrainTotSort$Full, type = "b", pch = 19, lwd = 2, col = "green")
lines(APCORTrainTotSort$Lasso, type = "b", pch = 19, lwd = 2, col = "red")
lines(APCORTrainTotSort$`Deep learning`, type = "b", pch = 19, lwd = 2, col = "black")
abline(v = 1:208, lty = 3, col = "grey")
axis(1, at = seq_along(APCORTrainTotSort$Env), APCORTrainTotSort$Env, las = 3,
     cex.axis = 0.7, labels = FALSE)
text(x = seq_along(APCORTrainTotSort$Env), labels = APCORTrainTotSort$Env,
     xpd = NA, srt = 90, y = -0.19, cex = 0.55)

text(x = 33.5, labels = "ET1", y = 1, cex = 0.95)
abline(v = 66.5)
text(x = 82.5, labels = "ET2", y = 1, cex = 0.9)
abline(v = 115.5)
text(x = 137.5, labels = "ET3", y = 1, cex = 0.9)
abline(v = 159)
text(x = 184.5, labels = "ET4", y = 1, cex = 0.9)


## Box plot with the accuracies per ET for each of the methods.
my_colors <- c("blue", "black", "green", "red")

## ggplot requires data in a long format, so pivot first.
APCORTrainTotSortLong <- tidyr::pivot_longer(APCORTrainTotSort,
                                             cols = -c(Env, Year, Group, MeanAcc),
                                             names_to = "Model", values_to = "r")

ggplot(APCORTrainTotSortLong, aes(x = Model, y = r, color = Model)) +
  geom_boxplot() +
  scale_color_manual(values = my_colors) +
  theme() +
  ylim(0, 1) +
  theme_bw() +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank())+
  facet_wrap(~Group)
ggsave("testScripts/paper/plots/TrainingBoxplotApsim (update 8 Dec).jpeg",
       width = 9, height = 6)
