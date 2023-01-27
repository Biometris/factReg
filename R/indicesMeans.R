#' @title
#' compute means of indices
#'
#' @description
#' compute means of indices
#'
#' @inheritParams GnE
#'
#' @param useY ...
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{indexMeans}{...}
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
