## Raw data files extracted from
## https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/IASSTN

## Read BLUEs.
BLUEs <- read.csv(file = "data-raw/2b-GrainYield_components_BLUEs_level-1.csv",
                  stringsAsFactors = TRUE)

## Remove parent1 column to match data in other statgen packages.
BLUEs <- BLUEs[, colnames(BLUEs) != "parent1"]

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
## Remove parent1 column to match data in other statgen packages.
indexDat <- indexDat[, colnames(indexDat) != "parent1"]

## Merge indices to BLUEs.
drops_GE <- merge(BLUEsGE, indexDat)
drops_GE <- droplevels(drops_GE)
drops_GE$type <- "GE"
rownames(drops_GE) <-
  paste0(drops_GE$Experiment, '_X_', drops_GE$Variety_ID)

## Create GnE data.
BLUEsGnE <- BLUEs[!rownames(BLUEs) %in% rownames(drops_GE), ]

## Merge indices to BLUEs.
drops_GnE <- merge(BLUEsGnE, indexDat)
drops_GnE <- drops_GnE[drops_GnE$type == "GnE_ext", ]
drops_GnE <- droplevels(drops_GnE)
drops_GnE$type <- "GnE"
rownames(drops_GnE) <-
  paste0(drops_GnE$Experiment, '_X_', drops_GnE$Variety_ID)

## Create nGnE data.
BLUEsnGnE <- BLUEs[!rownames(BLUEs) %in% rownames(drops_GE), ]

## Merge indices to BLUEs.
drops_nGnE <- merge(BLUEsnGnE, indexDat)
drops_nGnE <- drops_nGnE[drops_nGnE$type == "nGnE_ext", ]
drops_nGnE <- droplevels(drops_nGnE)
drops_nGnE$type <- "nGnE"
rownames(drops_nGnE) <-
  paste0(drops_nGnE$Experiment, '_X_', drops_nGnE$Variety_ID)

## Read genotypic information.
dropsGeno <- read.csv("data-raw/7a-Genotyping_50K_41722.csv",
                      row.names = 1)

dropsGenoNG <- read.csv("data-raw/7c-Genotyping_50K_41722_nG.csv",
                        row.names = 1)

dropsGenoTot <- as.matrix(rbind(dropsGeno, dropsGenoNG))

drops_K <- statgenGWAS::kinship(dropsGenoTot)

usethis::use_data(drops_GE, drops_GnE, drops_nGnE, drops_K, overwrite = TRUE)
