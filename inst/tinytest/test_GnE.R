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

### Check that input checks work as expected.
expect_error(GnE(dat = 1, Y = 1, G = 1, E = 1),
             "dat should be a data.frame")
expect_error(GnE(dat = testDat, Y = "a", G = 1, E = 1),
             "a should be a column in dat")
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "a", E = 1),
             "a should be a column in dat")
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID", E = "a"),
             "a should be a column in dat")
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment"),
             "Either indices or indicesData should be provided")
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment"),
             "Either indices or indicesData should be provided")
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment", indices = 1, indicesData = 1),
             "Either indices or indicesData should be provided")
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment", indices = 1),
             "indices should be a vector of columns in dat")
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment", indicesData = 1),
             "indicesData should be a data.frame with all environments")


## Check that columns can be refered to by names and by numbers.
modBase <- GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
               E = "Experiment", indices = indices)
mod1a <- GnE(dat = testDat, Y = 14, G = 2, E = 1, indices = indices)

## Check output structure.
expect_inherits(modBase, "list")
expect_equal(names(modBase),
             c("predTrain", "predTest", "resTrain", "resTest", "mu",
               "envInfoTrain", "envInfoTest", "parGeno", "trainAccuracyEnv",
               "testAccuracyEnv", "trainAccuracyGeno", "testAccuracyGeno",
               "lambda", "lambdaSequence", "RMSEtrain", "RMSEtest", "Y", "G",
               "E", "indices", "quadratic"))

## Check full output object.
expect_equal_to_reference(modBase, "modBase")

## Check that test environments can be added.
modTest <- GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
               E = "Experiment", indices = indices, testEnv = "Cam12R")

## Check full output object.
expect_equal_to_reference(modBase, "modTest")




