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
  nG <- nlevels(dat$G)
  nE <- nlevels(dat$E)
  if (is.list(K)) {
    ms <- setdiff(rownames(K[[1]]), levels(dat$G))
  } else {
    ms <- setdiff(rownames(K), levels(dat$G))
  }
  genotypes <- c(genotypes, ms)
  nG <- nG + length(ms)
  rownames(dat) <- paste0(dat$E, "_X_", dat$G)
  g <- data.frame(G = rep(genotypes, nE),
                  E = rep(environments, each = nG),
                  value = NA,
                  PEV = NA)
  rownames(g) <- paste0(g$E, "_X_", g$G)
  Vg <- Ve <- rep(NA, nE)
  dat$Ve <- rep(NA, nrow(dat))
  if (is.list(K)) {
    for (j in 1:nE) {
      dj <- dat[dat$E == environments[j], ]
      dj <- droplevels(dj)
      obs.geno <- as.character(dj$G[!is.na(dj$Y)])
      if (length(obs.geno) > 0) {
        Kj <- K
        for (i in 1:length(Kj)) {
          Kj[[i]] <- Kj[[i]][obs.geno, obs.geno]
        }
        r1 <- gaston::lmm.aireml(Y = dj$Y, K = Kj, verbose = FALSE)
        pr <- rep(r1$BLUP_beta, nG)
        names(pr) <- paste0(environments[j],'_X_',genotypes)
        for (i in 1:length(K)) {
          pr <- pr + r1$tau[i] *
            as.numeric(K[[i]][genotypes, obs.geno] %*% r1$Py)
        }
        g[names(pr), 3] <- as.numeric(pr)
      }
    }
    return(list(pred = g))
  } else {
    for (j in 1:nE) {
      dj <- dat[dat$E == environments[j],]
      dj <- droplevels(dj)
      b  <- rrBLUP::kin.blup(data = dj, geno = "G", pheno = "Y",
                             K = K, PEV = TRUE)
      Vg[j] <- b$Vg
      Ve[j] <- b$Ve
      dat[dat$E == environments[j], "Ve"] <- b$Ve
      names(b$pred) <- paste0(environments[j], "_X_", names(b$pred))
      g[names(b$pred), 3] <- as.numeric(b$pred)
      g[names(b$pred), 4] <- as.numeric(b$PEV)
    }
    return(list(pred = g, Vg = Vg, Ve = Ve))
  }
}
