#' @title
#' Compute means of environmental covariates or indices
#'
#' @description
#' Given a data-frame in long format with columns genotype, environment, and a 
#' number of environmental variables, compute the mean for each of these variables 
#' in each environment. This is useful when environmental variables have 
#' genotype-specific values within environments. This occurs especially if 
#' environmental variables are corrected for developmental timing.
#'
#' @inheritParams GnE
#'
#' @param useY If set to TRUE, all rows for which the trait of interest has a 
#' missing value are removed
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{indexMeans}{A data-frame with means per environment 
#'   (environments in the rows; environmental variables in the columns)}
#' }
#'
#' @export
indicesMeans <- function(dat,
                         Y,
                         G,
                         E,
                         indices,
                         useY = FALSE) {
  ## Get traitName.
  traitName <- ifelse(is.numeric(Y), names(dat)[Y], Y)
  ## Rename data columns for Y, G and E.
  dat <- renameYGE(dat = dat, Y = Y, G = G, E = E)
  dat <- droplevels(dat)
  ## Remove missing values
  if (useY) {
    #stopifnot(!is.null(Y))
    dat <- dat[!(is.na(dat$Y)),]
  }
  dat <- droplevels(dat)
  indexMeans <- matrix(nrow = nlevels(dat$E), ncol = length(indices),
                       dimnames = list(levels(dat$E), indices))
  indexMeans <- as.data.frame(indexMeans)
  for (ind in indices) {
    for (er in levels(dat$E)) {
      rs <- which(dat$E == er)
      mu <- mean(dat[rs, ind], na.rm = TRUE)
      indexMeans[er, ind] <- mu
    }
  }
  return(indexMeans)
}
