load("./data-raw/drops_GE.RData")
load("./data-raw/drops_GnE.RData")
load("./data-raw/drops_nGnE.RData")
load("./data-raw/drops_K.RData")

usethis::use_data(drops_GE, drops_GnE, drops_nGnE, drops_K, overwrite = TRUE)

