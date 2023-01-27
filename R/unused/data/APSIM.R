set.seed(123)

## Read genotypic data.
apsim_geno <- read.table(file = "./data-raw/genoAus199_maf02_code012.txt")
## Convert to matrix.
apsim_geno <- t(as.matrix(apsim_geno))

## Compute kinship matrix.
M <- apsim_geno
m <- colMeans(M)
for (i in 1:ncol(M)) {
  M[, i] <- (M[,i] - m[i]) / sqrt(m[i] * (2 - m[i]))
}
apsim_K <- tcrossprod(M) / ncol(M)

## Read phenotypic data and indices.
apsimRaw <- read.table(file = "./data-raw/Cov_windows_wide_248envs_199geno.txt",
                       header = TRUE)
## Convert some columns to factor.
apsimRaw$groups <- as.factor(apsimRaw$groups)
apsimRaw$Env <- as.factor(apsimRaw$Env)
apsimRaw$trial <- as.factor(apsimRaw$trial)
apsimRaw$geno <- as.factor(apsimRaw$geno)

## Add year column.
apsimRaw$year <- as.numeric(sapply(X = strsplit(as.character(apsimRaw$trial), split = '_'),
                                   FUN = `[`, 2))
## Convert sowing to 0-1 (June 15 = 1).
apsimRaw$sowing <- ifelse(apsimRaw$sowing == "15june", 1, 0)


## Select genotypic covariates.
geno.cov <- c("p1_biomassmean", "p2_biomassmean", "p3_biomassmean",
              "p4_biomassmean", "p5_biomassmean","p5_flowering")

## Construct vector of indices.
indices  <- setdiff(colnames(apsimRaw)[c(4, 7:62)], geno.cov)

## Construct test and training sets for genotypes and environments.
testGeno <- sample(unique(apsimRaw$geno), 39)
trainGeno <- setdiff(unique(apsimRaw$geno), testGeno)
testEnv <- unique(apsimRaw$Env[apsimRaw$year > 2008])
trainEnv <- setdiff(unique(apsimRaw$Env), testEnv)

## Compute the GBLUP for each env.
## GBLUPs themselves are not needed; only the estimated variance components.
## Depending on these, we will next add noise such that the simulated
## data will have heritability of approximately 0.5, in each environment.
gBlups <- GnEglmnet::nGE(dat = apsimRaw, Y = "yld", G = "geno", E = "Env",
                         keepObserved = FALSE, K = apsim_K)

## Add noise.
for (env in trainEnv) {
  j <- unique(apsimRaw$Env) == env
  n <- sum(apsimRaw$Env == env)
  noise <- rnorm(n = n, sd = 0.5 * sqrt(gBlups$Vg[j]))
  apsimRaw[apsimRaw$Env == env, "yld"] <-
    apsimRaw[apsimRaw$Env == env, "yld"] + noise
}

## Construct subsets.
all_apsim <- apsimRaw

GE_apsim <- all_apsim[(all_apsim$geno %in% trainGeno) & (all_apsim$Env %in% trainEnv), ]
nGE_apsim <- all_apsim[(all_apsim$geno %in% testGeno) & (all_apsim$Env %in% trainEnv), ]
GnE_apsim <- all_apsim[(all_apsim$geno %in% trainGeno) & (all_apsim$Env %in% testEnv), ]
nGnE_apsim <- all_apsim[(all_apsim$geno %in% testGeno) & (all_apsim$Env %in% testEnv), ]

GE_apsim <- droplevels(GE_apsim)
nGE_apsim <- droplevels(nGE_apsim)
GnE_apsim <- droplevels(GnE_apsim)
nGnE_apsim <- droplevels(nGnE_apsim)

## Add type to all_apsim.
all_apsim$type <- "GE"
all_apsim$type[(all_apsim$geno %in% trainGeno) &
                 (all_apsim$Env %in% testEnv)] <- "GnE"
all_apsim$type[(all_apsim$geno %in% testGeno) &
                 (all_apsim$Env %in% trainEnv)] <- "nGE"
all_apsim$type[(all_apsim$geno %in% testGeno) &
                 (all_apsim$Env %in% testEnv)] <- "nGnE"
all_apsim$type <- as.factor(all_apsim$type)

usethis::use_data(indices, geno.cov, apsim_K, apsim_geno, all_apsim,
                  GE_apsim, nGE_apsim, GnE_apsim, nGnE_apsim, overwrite = TRUE)

