#' @title
#' compute means of indices
#'
#' @description
#' compute means of indices
#'
#' @inheritParams GnE
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
                use.Y = FALSE) {

  ## Get traitName.
  traitName <- ifelse(is.numeric(Y), names(dat)[Y], Y)
  ## Rename data columns for Y, G and E.
  dat <- renameYGE(dat = dat, Y = Y, G = G, E = E)

  dat <- droplevels(dat)

  ## Remove missing values
  if (use.Y) {
    #stopifnot(!is.null(Y))
    dat <- dat[!(is.na(dat$Y)),]
  }

  dat <- droplevels(dat)

  indexMeans <- matrix(NA, nlevels(dat$E), length(indices))
  rownames(indexMeans) <- levels(dat$E)
  colnames(indexMeans) <- indices
  indexMeans <- as.data.frame(indexMeans)

  for (ind in indices) {
    for (er in levels(dat$E)) {
      rs <- which(dat$E == er)
      mu <- mean(dat[rs, ind], na.rm = TRUE)
      indexMeans[er, ind] <- mu
    }
  }

  return(list(indexMeans = indexMeans))
}
