data(drops_GE)
data(drops_GnE)
data(drops_nGnE)
data(drops_K)

## Raw data files extracted from
## https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/IASSTN

## Read BLUEs.
BLUEs <- read.csv(file = "data-raw/2b-GrainYield_components_BLUEs_level-1.csv",
                  stringsAsFactors = TRUE)

## One record in BLUEs is duplicated.
BLUEs <- BLUEs[!duplicated(BLUEs), ]

BLUEs <- BLUEs[, -which(colnames(BLUEs) %in% c("geno.panel", "type"))]
rownames(BLUEs) <- paste0(BLUEs$Experiment, '_X_', BLUEs$Variety_ID)

## Remove experiments with less than 100 observations.
BLUEsObs <- table(BLUEs$Experiment)
BLUEsGE <- BLUEs[BLUEs$Experiment %in% names(BLUEsObs[BLUEsObs > 100]), ]

## Read indices.
indexDat <- read.csv(file = "data-raw/3-Indices_GenoEnv_level-1.csv",
                     stringsAsFactors = TRUE)

## Merge indices to BLUEs.
drops_GEnw <- merge(BLUEsGE, indexDat)
drops_GEnw <- droplevels(drops_GEnw)
drops_GEnw$type <- "GE"
rownames(drops_GEnw) <- paste0(drops_GEnw$Experiment, '_X_', drops_GEnw$Variety_ID)


tst <- drops_GE[rownames(drops_GEnw), colnames(drops_GEnw)]
all.equal(tst, drops_GEnw)

## Create GnE data.
BLUEsGnE <- BLUEs[!rownames(BLUEs) %in% rownames(drops_GE), ]

## Merge indices to BLUEs.
drops_GnEnw <- merge(BLUEsGnE, indexDat)
drops_GnEnw <- drops_GnEnw[drops_GnEnw$type == "GnE_ext", ]
drops_GnEnw <- droplevels(drops_GnEnw)
drops_GnEnw$type <- "GnE"
rownames(drops_GnEnw) <- paste0(drops_GnEnw$Experiment, '_X_', drops_GnEnw$Variety_ID)

tst2 <- drops_GnE[rownames(drops_GnEnw), colnames(drops_GnEnw)]
all.equal(tst2, drops_GnEnw)




## Create nGnE data.
BLUEsnGnE <- BLUEs[!rownames(BLUEs) %in% rownames(drops_GE), ]

## Merge indices to BLUEs.
drops_nGnEnw <- merge(BLUEsnGnE, indexDat)
drops_nGnEnw <- drops_nGnEnw[drops_nGnEnw$type == "nGnE_ext", ]
drops_nGnEnw <- droplevels(drops_nGnEnw)
drops_nGnEnw$type <- "nGnE"
rownames(drops_nGnEnw) <- paste0(drops_nGnEnw$Experiment, '_X_', drops_nGnEnw$Variety_ID)

tst3 <- drops_nGnE[rownames(drops_nGnEnw),
                   colnames(drops_nGnEnw)[colnames(drops_nGnEnw) %in% colnames(drops_nGnE)]]
all.equal(tst3, drops_nGnEnw[, colnames(drops_nGnEnw) %in% colnames(drops_nGnE)])


comp <- read.csv("c:/Projects/R_packages/corteva/users/willem/drops_data/Final_files_Amaizing/7a-Genotyping_50K_41722_Amaizing_56Geno.csv",
                 row.names = 1)

load("c:/Projects/R_packages/corteva/users/willem/drops_data/inra_website/7a-Genotyping_50K_41722.RData")

rownames(x) <- x[, 1]
x <- x[, -1]
x <- as.matrix(x)

all.equal(x, as.matrix(dropsGeno))

## Read genotypic information.
dropsGeno <- read.csv("data-raw/7a-Genotyping_50K_41722.csv",
                      row.names = 1)

dropsGenoNG <- read.csv("data-raw/7c-Genotyping_50K_41722_nG.csv",
                        row.names = 1)

dropsGenoTot <- as.matrix(rbind(dropsGeno, dropsGenoNG))


m <- colMeans(dropsGeno)

Z <- dropsGenoTot

for (j in 1:ncol(dropsGenoTot)) {
  Z[, j] <- (Z[, j] - m[j]) / sqrt(m[j] * (2 - m[j]))
}

drops_Knw <- Z %*% t(Z) / ncol(M)

drops_Knw <- statgenGWAS::kinship(dropsGenoTot) / 2

all.equal(drops_K, drops_Knw)

round(drops_K[1:10,1:10] - drops_Knw[1:10,1:10], 2)


