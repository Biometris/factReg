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
rownames(drops_GEnw) <-
  paste0(drops_GEnw$Experiment, '_X_', drops_GEnw$Variety_ID)

## Create GnE data.
BLUEsGnE <- BLUEs[!rownames(BLUEs) %in% rownames(drops_GE), ]

## Merge indices to BLUEs.
drops_GnEnw <- merge(BLUEsGnE, indexDat)
drops_GnEnw <- drops_GnEnw[drops_GnEnw$type == "GnE_ext", ]
drops_GnEnw <- droplevels(drops_GnEnw)
drops_GnEnw$type <- "GnE"
rownames(drops_GnEnw) <-
  paste0(drops_GnEnw$Experiment, '_X_', drops_GnEnw$Variety_ID)

## Create nGnE data.
BLUEsnGnE <- BLUEs[!rownames(BLUEs) %in% rownames(drops_GE), ]

## Merge indices to BLUEs.
drops_nGnEnw <- merge(BLUEsnGnE, indexDat)
drops_nGnEnw <- drops_nGnEnw[drops_nGnEnw$type == "nGnE_ext", ]
drops_nGnEnw <- droplevels(drops_nGnEnw)
drops_nGnEnw$type <- "nGnE"
rownames(drops_nGnEnw) <-
  paste0(drops_nGnEnw$Experiment, '_X_', drops_nGnEnw$Variety_ID)

## Read genotypic information.
dropsGeno <- read.csv("data-raw/7a-Genotyping_50K_41722.csv",
                      row.names = 1)

dropsGenoNG <- read.csv("data-raw/7c-Genotyping_50K_41722_nG.csv",
                        row.names = 1)

dropsGenoTot <- as.matrix(rbind(dropsGeno, dropsGenoNG))

drops_K <- statgenGWAS::kinship(dropsGenoTot)

usethis::use_data(drops_GE, drops_GnE, drops_nGnE, drops_K, overwrite = TRUE)
