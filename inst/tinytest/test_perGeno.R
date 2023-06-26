## load the data.
data(drops_GE)
data(drops_GnE)
data(drops_K)

## Remove identifiers that we don't need.
drops_GE_GnE <- rbind(drops_GE[, -c(2, 3, 5)], drops_GnE[, -c(2, 3, 5)])

## Restrict to 10 genotypes.
testDat <- drops_GE_GnE[drops_GE_GnE$Variety_ID %in%
                          levels(drops_GE_GnE$Variety_ID)[1:10], ]
testDat <- droplevels(testDat)

testK <- drops_K[levels(testDat$Variety_ID), levels(testDat$Variety_ID)]

## Restrict to 11 environments.
# At least 10 are needed for training to prevent removal of all genotypes.
# 1 used for testing.
testDat <- testDat[testDat$Experiment %in% levels(testDat$Experiment)[1:11], ]
testDat <- droplevels(testDat)

indices <- c("Tnight.Early", "Tnight.Flo")

indicesDat <- testDat[testDat$Variety_ID == "A3",
                      c("Tnight.Early", "Tnight.Flo")]
rownames(indicesDat) <- substring(rownames(indicesDat), first = 1, last = 6)

### Check that input checks work as expected.
expect_error(perGeno(dat = 1, Y = 1, G = 1, E = 1),
             "dat should be a data.frame")
expect_error(perGeno(dat = testDat, Y = "a", G = 1, E = 1),
             "a should be a column in dat")
expect_error(perGeno(dat = testDat, Y = "grain.yield", G = "a", E = 1),
             "a should be a column in dat")
expect_error(perGeno(dat = testDat, Y = "grain.yield", G = "Variety_ID", E = "a"),
             "a should be a column in dat")
expect_error(perGeno(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                     E = "Experiment"),
             "Either indices or indicesData should be provided")
expect_error(perGeno(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                     E = "Experiment"),
             "Either indices or indicesData should be provided")
expect_error(perGeno(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                     E = "Experiment", indices = 1, indicesData = 1),
             "Either indices or indicesData should be provided")
expect_error(perGeno(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                     E = "Experiment", indices = 1),
             "indices should be a vector of length > 1 of columns in dat")
expect_error(perGeno(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                     E = "Experiment", indices = indices[1]),
             "indices should be a vector of length > 1 of columns in dat")


## Check that columns can be refered to by names and by numbers.
modBase <- perGeno(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                   E = "Experiment", indices = indices)
mod1a <- perGeno(dat = testDat, Y = 14, G = 2, E = 1, indices = indices)

## Check output structure.
expect_inherits(modBase, "list")
expect_equal(names(modBase),
             c("predTrain", "predTest", "envInfoTrain", "envInfoTest",
               "testAccuracyEnv", "trainAccuracyEnv", "trainAccuracyGeno",
               "testAccuracyGeno", "RMSEtrain", "RMSEtest", "Y", "G", "E",
               "indices", "genotypes", "lambdaOpt", "parGeno", "quadratic"))

## Check full output object.
expect_equal_to_reference(modBase, "modBasePerGeno")

## Check that test environments can be added.

expect_error(perGeno(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                     E = "Experiment", indices = indices, testEnv = 1),
             "testEnv should be a vector of environments present in dat")
expect_error(perGeno(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                     E = "Experiment", indices = indices, testEnv = "bla"),
             "testEnv should be a vector of environments present in dat")


modTest <- perGeno(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                   E = "Experiment", indices = indices, testEnv = "Cam12R")

## Check full output object.
expect_equal_to_reference(modTest, "modTestPerGeno")


## Check that indicesData can be used.
expect_error(perGeno(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                     E = "Experiment", indicesData = 1),
             "indicesData should be a data.frame with all environments")
expect_error(perGeno(dat = testDat[testDat$Experiment != "Cam12R", ],
                     Y = "grain.yield", G = "Variety_ID", E = "Experiment",
                     indicesData = indicesDat),
             "All environments in indicesData should be in dat")

expect_warning(modIndDat <- perGeno(dat = testDat, Y = "grain.yield",
                                    G = "Variety_ID", E = "Experiment",
                                    indicesData = indicesDat),
               "The following columns in indicesDat are already in dat")

expect_equal(mean(modIndDat$trainAccuracyEnv$r), 0.723234644385349)


## Check that redundant genotypes are removed.
expect_warning(perGeno(dat = testDat[!testDat$Experiment %in% c("Gai12W", "Gai13R"), ],
                       Y = "grain.yield", G = "Variety_ID",
                       E = "Experiment", indices = indices),
               "the following genotypes have < 10 observations, and are removed")
expect_error(perGeno(dat = testDat[!testDat$Experiment %in% c("Gai12W", "Gai13R"), ],
                     Y = "grain.yield", G = "Variety_ID",
                     E = "Experiment", indices = indices),
             "No data left in training set")


## Check that option quadratic works correctly.
modQuad <- perGeno(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                   E = "Experiment", indices = indices, quadratic = TRUE)
expect_equal(modQuad$indices, c("Tnight.Early", "Tnight.Flo",
                                "Tnight.Early_quad", "Tnight.Flo_quad"))
expect_equal(mean(modQuad$trainAccuracyEnv$r), 0.883214087481033)

