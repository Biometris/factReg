#' @title
#' this is the title
#'
#' @description
#' \loadmathjax
#' here comes the description
#'
#' @inheritParams GnE
#'
#' @param K kinship matrix, or a list of such matrices. (TO DO: checks)
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{pred}{...}
#'   \item{Vg}{...},
#'   \item{Ve}{...}
#' }
#'
#' @export
nGE <- function(dat,
                Y,
                G,
                E,
                K) {
  ## Rename data columns for Y, G and E.
  dat <- renameYGE(dat = dat, Y = Y, G = G, E = E)
  dat <- dat[!is.na(dat$Y), ]
  dat <- droplevels(dat)
  if (is.list(K)) {
    for (i in 1:length(K)) {stopifnot(all(dat$G %in% colnames(K[[i]])))}
  } else {
    stopifnot(all(dat$G %in% colnames(K)))
  }
  genotypes <- levels(dat$G)
  environments <- levels(dat$E)
  nGeno <- nlevels(dat$G)
  nEnv <- nlevels(dat$E)
  if (is.list(K)) {
    ms <- setdiff(rownames(K[[1]]), genotypes)
  } else {
    ms <- setdiff(rownames(K), genotypes)
  }
  genotypes <- c(genotypes, ms)
  nGeno <- nGeno + length(ms)
  rownames(dat) <- paste0(dat$E, "_X_", dat$G)
  pred <- data.frame(G = rep(genotypes, nEnv),
                     E = rep(environments, each = nGeno),
                     value = NA,
                     PEV = NA)
  rownames(pred) <- paste0(pred$E, "_X_", pred$G)
  Vg <- Ve <- rep(NA, nEnv)
  dat$Ve <- rep(NA, nrow(dat))
  if (is.list(K)) {
    for (j in 1:nEnv) {
      envDat <- dat[dat$E == environments[j], ]
      envDat <- droplevels(envDat)
      obsGeno <- as.character(envDat$G[!is.na(envDat$Y)])
      if (length(obsGeno) > 0) {
        Kj <- K
        for (i in 1:length(Kj)) {
          Kj[[i]] <- Kj[[i]][obsGeno, obsGeno]
        }
        kinMod <- gaston::lmm.aireml(Y = envDat$Y, K = Kj, verbose = FALSE)
        pr <- rep(kinMod$BLUP_beta, nGeno)
        names(pr) <- paste0(environments[j], "_X_", genotypes)
        for (i in 1:length(K)) {
          pr <- pr + kinMod$tau[i] *
            as.numeric(K[[i]][genotypes, obsGeno] %*% kinMod$Py)
        }
        pred[names(pr), "value"] <- as.numeric(pr)
      }
    }
    return(list(pred = pred))
  } else {
    for (j in 1:nEnv) {
      envDat <- dat[dat$E == environments[j], ]
      envDat <- droplevels(envDat)
      kinMod  <- rrBLUP::kin.blup(data = envDat, geno = "G", pheno = "Y",
                                  K = K, PEV = TRUE)
      Vg[j] <- kinMod$Vg
      Ve[j] <- kinMod$Ve
      dat[dat$E == environments[j], "Ve"] <- kinMod$Ve
      names(kinMod$pred) <- paste0(environments[j], "_X_", names(kinMod$pred))
      pred[names(kinMod$pred), "value"] <- as.numeric(kinMod$pred)
      pred[names(kinMod$pred), "PEV"] <- as.numeric(kinMod$PEV)
    }
    return(list(pred = pred, Vg = Vg, Ve = Ve))
  }
}
