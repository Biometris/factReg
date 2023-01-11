#' Helper function for renaming Y, G, and E
#'
#' Helper function for renaming columns Y, G and Y in dat
#'
#' @noRd
#' @keywords internal
renameYGE <- function(dat,
                      Y,
                      G,
                      E) {
  traitName <- ifelse(is.numeric(Y), names(dat)[Y], Y)
  for (var in c("Y", "G", "E")) {
    varVal <- get(var)
    if (is.numeric(varVal)) {
      names(dat)[varVal] <- var
    } else if (is.character(varVal)) {
      stopifnot(varVal %in% names(dat))
      names(dat)[names(dat) == varVal] <- var
    }
  }
  ## Convert G and E to factor.
  if (!is.factor(dat[["G"]])) {
    dat[["G"]] <- as.factor(dat[["G"]])
  }
  if (!is.factor(dat[["E"]])) {
    dat[["E"]] <- as.factor(dat[["E"]])
  }
  dat <- droplevels(dat)
  return(dat)
}

#' Helper function
#'
#' Helper function for ....
#'
#' @noRd
#' @keywords internal
exRank <- function(y0,
                   yPred,
                   topOnly = TRUE,
                   k = 5,
                   m = 10) {
  b <- sort(y0, decreasing = TRUE)[k]
  topK <- which(y0 >= b)

  a <- sort(yPred, decreasing = TRUE)[m]
  topM <- which(yPred >= a)

  return(length(intersect(topK, topM)) / k)
}

