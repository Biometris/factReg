orig_GE <- drops_GE


## Raw data files extracted from
## https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/IASSTN

## Read BLUEs.
BLUEs <- read.csv(file = "data-raw/2b-GrainYield_components_BLUEs_level-1.csv")

## One record in BLUEs is duplicated.
BLUEs <- BLUEs[!duplicated(BLUEs), ]

rownames(BLUEs) <- paste0(BLUEs$Experiment, '_X_', BLUEs$Variety_ID)
BLUEs <- BLUEs[, -which(colnames(BLUEs) %in% c("geno.panel", "type"))]

## Remove experiments with less than 100 observations.
BLUEsObs <- table(BLUEs$Experiment)
BLUEs <- BLUEs[BLUEs$Experiment %in% names(BLUEsObs[BLUEsObs > 100]), ]

## Read indices.
indexDat <- read.csv(file = "data-raw/3-Indices_GenoEnv_level-1.csv")
rownames(indexDat) <- paste0(indexDat$Experiment, '_X_', indexDat$Variety_ID)

## Merge indices to BLUEs.
drops_GEnw <- merge(BLUEs, indexDat)




d2 <- data.frame(d2, d1[rownames(d2), -(1:6)])



d3 <- read.csv(file = "inra_website/8-Info_Maize_variety.csv")

de <- read.csv(file = "files_emilie_6jan_2021/DROPS-data_V.3.0/3-Indices_GenoEnv_level.csv")
#d3$parent1 <- factor(unlist(strsplit(as.character(d3$parent1), split = '_uh')))

rownames(de) <- paste0(de$Experiment, '_X_', de$Variety_ID)

dd <- read.csv(file = "files_emilie_6jan_2021/DROPS-data_V.3.0/2b-GrainYield_components_BLUEs_level.csv")
#dd <- dd[dd$type=='GnE_ext',]
rownames(dd) <- paste0(dd$Experiment, '_X_', dd$Variety_ID)

#d3$parent1 <- factor(unlist(strsplit(as.character(d3$parent1), split = '_uh')))

rownames(de) <- paste0(de$Experiment, '_X_', de$Variety_ID)



drops_GnE  <- de[de$type == 'GnE_ext', ]


drops_GnE  <- cbind(drops_GnE[,1:16], dd[rownames(drops_GnE), 7:15])

drops_GnE$type  <- 'GnE'

drops_GnE  <- droplevels(drops_GnE)

rownames(d3) <- d3$Variety_ID

# x, d1-d3 : variety_ID
# a, d4-d5 : also variety_ID
#intersect(rownames(a), d4$Variety_ID)

##########################################################3

a <- read.csv(file = "Final_files_Amaizing/7a-Genotyping_50K_41722_Amaizing_56Geno.csv", row.names = 1)

load(file = "inra_website/7a-Genotyping_50K_41722.RData")

rownames(x) <- d3[as.character(x$Ind),'Variety_ID']
x <- x[,-1]

#xx <- t(as.matrix(GWAS.obj$markers))
#colnames(xx) <- gsub(x = colnames(xx), pattern = "-", repl = ".")

#length(intersect(colnames(xx), colnames(a)))

#x1 <- xx[, intersect(colnames(xx), colnames(a))]
#a1 <- a[, intersect(colnames(xx), colnames(a))]

x1 <- x[, intersect(colnames(x), colnames(a))]
a1 <- a[, intersect(colnames(x), colnames(a))]


length(intersect(rownames(x1), rownames(a1)))

all(rownames(x1) %in% as.character(d3$parent1))

#rownames(a1) <- d3[rownames(a1), "parent1"]

x <- as.matrix(rbind(x1, a1))

n <- nrow(x)
p <- ncol(x)

m <- apply(x1,2,mean)

M <- x

for (j in 1:p) {
  M[,j] <- (M[,j] - m[j]) / sqrt(m[j] * (2 - m[j]))
}

# for (j in 1:p) {
#   M[,j] <- (M[,j] - m[j]) / sd(M[,j])
# }

M <- as.matrix(M)

drops_KK2 <- M %*% t(M) / ncol(M)

save(drops_K, file = 'drops_amaizing_kinship_50k.RData')

###############################


d4 <- read.csv(file = "Final_files_Amaizing/2b-GrainYield_components_BLUEs_level_Amaizing_32+56Geno.csv")
d4 <- d4[-5,]

rownames(d4) <- paste0(d4$Experiment, '_X_', d4$Variety_ID)

#d5 <- read.csv(file = "Final_files_Amaizing/3-Indices_GenoEnv_level_Amaizing_32Geno.csv")
d5 <- read.csv(file = "Final_files_Amaizing/3-Indices_GenoEnv_level_Amaizing_56Geno.csv")

rownames(d5) <- paste0(d5$Experiment, '_X_', d5$Variety_ID)

d5 <- data.frame(d5, d4[rownames(d5), -(1:6)])



drops_GE <- d2
drops_nGnE  <- d5

drops_GE$type <- 'GE'
drops_nGnE$type <- 'nGnE'

#drops_GnE  <- d1[d1$Experiment %in% setdiff(levels(d1$Experiment), levels(d2$Experiment)) & d1$Variety_ID %in% as.character(d1$Variety_ID[d1$Experiment=='Cam11W']), ]
#drops_GnE  <- cbind(drops_GnE, de[rownames(drops_GnE), ])

#levels(dTrain$parent1)

#levels(d4$parent1) %in% levels(dTrain$parent1)

#save(drops_GE, drops_nGnE, drops_GnE, file = 'drops_amaizing_pheno_and_indices.RData')

setwd(data.folder)
save(drops_GE, file = 'drops_GE.RData')
save(drops_GnE, file = 'drops_GnE.RData')
save(drops_nGnE, file = 'drops_nGnE.RData')
save(drops_K, file = 'drops_K.RData')

