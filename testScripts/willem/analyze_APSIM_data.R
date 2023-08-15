
# below, /mnt/data/M/research/STATISTICAL_GENETICS/corteva/
# is the location of our corteva git repository

# we use 'indices' and 'environmental covariates' interchangeably
# (to be updated)

#options(warn=2, error=recover)
#options(warn=0)

library(factReg)
library(ggplot2)

set.seed(123)


# load deep learning results provided by aalt-jan, for nGnE
# Vahe: please look in the files if also the GnE results are present; if not
#         please ask aalt-jan
# rkw <- read.csv("/mnt/data/M/research/STATISTICAL_GENETICS/corteva/users/willem/results_aaltJan/predicties_KW_yield.csv")[1:40,]
# rownames(rkw) <- rkw$Env

# load the data, which are contained in the package

load("./R/unused/apsim_K.rda")
# contains the kinship matrix
#data("apsim_info")  # contains the names of the indices
load("./R/unused/all_apsim.rda")

# contains training and test data; note the column 'type'



# re-name, just to have a shorter name
d           <- all_apsim
# define meaningful row-names
rownames(d) <- paste0(d$Env, '_X_', d$geno)
d           <- droplevels(d)

##################

d0 <- droplevels(d[d$type %in% c('GE','nGE'),])
d1 <- droplevels(d[d$type %in% c('GE','GnE'),])
d2 <- droplevels(d[d$type == 'nGnE',])

# define vectors with the  names of the training and test environments
tre <- as.character(unique(d1$Env[d1$type=='GE']))
tst <- setdiff(levels(d1$Env), tre)

# define vectors with the  names of the training and test genotypes
treG <- as.character(unique(d$geno[d$type=='GE']))
tstG <- setdiff(levels(d$geno), treG)

# define which training environments will be together in the same fold,
# when doing inner cross validation
partitionDF <- data.frame(E = tre,
                          partition = sample(cut(1:length(tre),
                                                 breaks = 12, labels = FALSE)))

indices <- names(all_apsim)[7:62]

#######################################

# unpenalized factorial regression (GE --> GnE)
out0 <- GnE(dat = d1, lambda = 0,
            Y = 'yld', G = "geno", E = "Env",
            indices = indices,
            partition = partitionDF,
            testEnv = tst)

# Penalized factorial regression (GE --> GnE)
out1 <- GnE(dat = d1,
            Y = 'yld', G = "geno", E = "Env",
            indices = indices,
            partition = partitionDF,
            testEnv = tst,
            penG = 0,
            alpha = 1)


# per genotype regression, using residuals (GE --> GnE)
out2 <- perGeno(dat = d1,
                Y = 'yld', G = "geno", E = "Env",
                indices = indices,
                partition = partitionDF,
                testEnv = tst,
                useRes = 2,
                alpha = 1)

# per genotype regression, NOT using residuals (GE --> GnE)
out3 <- perGeno(dat = d1,
                Y = 'yld', G = "geno", E = "Env",
                indices = indices,
                partition = partitionDF,
                testEnv = tst,
                useRes = 0,
                alpha = 1)

# # Penalized factorial regression (GE --> nGnE)
out1nG <- nGnE(GnEOut = out1, K = apsim_K, dNew = d2)
#
# #    GE --> nGnE with per-genotype regression
# out2nG <- nGnE(GnEOut = out2, K = apsim_K, dNew = d2)
#
# summary(out1$testAccuracyEnv$r - out1$testAccuracyEnv$rMain)
# summary(out1b$testAccuracyEnv$r - out1$testAccuracyEnv$rMain)
# summary(out1b$testAccuracyEnv$r - out1$testAccuracyEnv$r)
# summary(out1$testAccuracyEnv$r - out0$testAccuracyEnv$r)
# summary(out2$testAccuracyEnv$r - out1$testAccuracyEnv$r)
# summary(out2$testAccuracyEnv$r - out0$testAccuracyEnv$r)
# summary(out3$testAccuracyEnv$r - out0$testAccuracyEnv$r)
# summary(out3$testAccuracyEnv$r - out2$testAccuracyEnv$r)
# ###############################################################33
# # Example plots:
#
# # example plot 1: penalized vs unpenalized
# qw <- data.frame(FR_unpenalized =  out0$testAccuracyEnv$r,
#                  FR_penalized = out1$testAccuracyEnv$r)
# pdf(file = 'apsim_GnE_penalized_vs_unpenalized.pdf')
# ggplot(qw, aes(x=FR_unpenalized, y=FR_penalized)) + geom_point(size = 3)
# + geom_segment(aes(x = -0.2, y = -0.2, xend = 0.6, yend = 0.6)) + theme(axis.text=element_text(size=18), axis.title=element_text(size=25,face="bold"))
# dev.off()
#
# # example plot 2: penalized vs main effects only
# qw <- data.frame(mainOnly =  out0$testAccuracyEnv$rMain,
#                  FR_penalized = out1$testAccuracyEnv$r)
# pdf(file = 'apsim_GnE_penalized_vs_mainOnly.pdf')
# ggplot(qw, aes(x=mainOnly, y=FR_penalized)) + geom_point(size = 3) +
#   geom_segment(aes(x = -0.2, y = -0.2, xend = 1, yend = 1)) + theme(axis.text=element_text(size=18), axis.title=element_text(size=25,face="bold"))
# dev.off()
#
# # example plot 3 (without ggplot2): penalized vs genotype-specific
# plot(x = out1$testAccuracyEnv$r, y = out2$testAccuracyEnv$r,
#      main = 'apsim GnE', xlab = 'factReg', ylab = 'perGeno',
#      pch = 20); abline(a=0, b=1)
#
#
# # example plot 4: fact. regr. (penalized) for nGnE (rather than GnE),
# #                 vs deep learning (DL)
# qw <- data.frame(DL = rkw[out1$testAccuracyEnv$Env, "X.2"],
#                  factReg = out1nG$testAccuracyEnv$r)
# dev.off()
# pdf(file = 'apsim_nGnE_factReg_vs_DL.pdf')
# ggplot(qw, aes(x=DL, y=factReg)) + geom_point(size = 3) + geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1)) + theme(axis.text=element_text(size=18), axis.title=element_text(size=25,face="bold"))
# dev.off()


#############################################################################
# The following two analyses are to be used in a later paper:
#   fact. regr. on predicted genetic effects:

out4 <- frBLUP(dat = d, type = d$type, target = 'nGnE',
               method = 'perGeno',
               Y = 'yld', G = "geno", E = "Env",
               indices = indices,
               partition = partitionDF,
               alpha = 1,
               K = apsim_K,
               useRes = 1)

out5 <- frBLUP(dat = d, type = d$type, target = 'nGnE',
               method = 'factReg',
               Y = 'yld', G = "geno", E = "Env",
               indices = indices,
               partition = partitionDF,
               alpha = 1,
               K = apsim_K)

# dat1 <- data.frame(G = out1$dTrain$G, E = out1$dTrain$E, rr = out1$resTrain)
# aa <- reshape(dat1, idvar = "G", timevar = "E", direction = "wide", v.names = 'rr')
# for (j in 2:ncol(aa)) {aa[,j] <- as.numeric(aa[,j])}
# heatmap(as.matrix(aa[,-1]))

out5$testAccuracyEnv$r
range(out5$testAccuracyEnv$r)
mean(out5$testAccuracyEnv$r)

pred5 <- out5$genPred
comp5 <- merge(out5$data, pred5, by.x = c("G", "E"), by.y = c("G", "E"))

tst <- comp5[comp5$type == "nGnE", ]
cor(tst$Yobs, tst$Y)
tst$Yobs - tst$Y



range(tst$Yobs - tst$gblup)
