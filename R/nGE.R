#' @title
#' this is the title
#'
#' @description
#' \loadmathjax
#' here comes the description
#'
#' @inheritParams GnE_glmnet
#'
#' @param K kinship matrix, or a list of such matrices. (TO DO: checks)
#'
#' @export
nGE <- function(dat,
                Y,
                G,
                E,
                K) {
  #                algorithm = c('',''),
  #                keepObserved = TRUE,
  #                predictAll = FALSE) {

  # from GnE_glmnet:
  ## Rename data columns for Y, G and E.
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

  if (!is.factor(dat$G)) {
    dat$G <- as.factor(dat$G)
  }
  if (!is.factor(dat$E)) {
    dat$E <- as.factor(dat$E)
  }

  dat <- dat[!is.na(dat$Y), ]
  dat <- droplevels(dat)

  if (is.list(K)) {
    for (i in 1:length(K)) {stopifnot(all(dat$G %in% colnames(K[[i]])))}
  } else {
    stopifnot(all(dat$G %in% colnames(K)))
  }

  genotypes    <- levels(dat$G)
  environments <- levels(dat$E)

  ng <- nlevels(dat$G)
  ne <- nlevels(dat$E)

  if (is.list(K)) {
    ms <- setdiff(rownames(K[[1]]), levels(dat$G))
  } else {
    ms <- setdiff(rownames(K), levels(dat$G))
  }

  genotypes <- c(genotypes, ms)
  ng        <- ng + length(ms)

  rownames(dat) <- paste0(dat$E,'_X_',dat$G)

  ###############

  g <- data.frame(G = rep(genotypes, ne),
                  E = rep(environments, each = ng),
                  value = rep(NA, ne * ng),
                  PEV   = rep(NA, ne * ng))

  rownames(g) <- paste0(g$E,'_X_',g$G)

  Vg <- rep(NA, ne); Ve <- rep(NA, ne)
  dat$Ve <- rep(NA, nrow(dat))

  if (is.list(K)) {

    for (j in 1:ne) {

      dj <- dat[dat$E == environments[j],]
      dj <- droplevels(dj)
      obs.geno <- as.character(dj$G[!is.na(dj$Y)])

      if (length(obs.geno) > 0) {

        K.j <- K
        for (i in 1:length(K.j)) {K.j[[i]] <- K.j[[i]][obs.geno, obs.geno]}
        r1 <- gaston::lmm.aireml(Y = dj$Y, K = K.j, verbose = F)
        pr <- rep(r1$BLUP_beta, ng)
        names(pr) <- paste0(environments[j],'_X_',genotypes)

        for (i in 1:length(K)) {
          pr <- pr + r1$tau[i] * as.numeric(K[[i]][genotypes, obs.geno] %*% r1$Py)
        }
        #dat[dat$E == environments[j],'Ve'] <- b$Ve
        #names(b$pred) <- paste0(environments[j],'_X_',names(b$pred))
        g[names(pr), 3] <- as.numeric(pr) # names(b$pred)
        #g[names(b$pred), 4] <- as.numeric(b$PEV)
      } else {

      }
    }

    return(list(pred = g))

  } else {

    for (j in 1:ne) {
      dj <- dat[dat$E == environments[j],]
      dj <- droplevels(dj)
      b  <- rrBLUP::kin.blup(data = dj, geno = 'G', pheno = "Y", K = K, PEV = TRUE)
      Vg[j] <- b$Vg; Ve[j] <- b$Ve
      dat[dat$E == environments[j],'Ve'] <- b$Ve
      names(b$pred) <- paste0(environments[j],'_X_',names(b$pred))
      g[names(b$pred), 3] <- as.numeric(b$pred) # names(b$pred)
      g[names(b$pred), 4] <- as.numeric(b$PEV)
    }

    # if (keepObserved==TRUE) {
    #   rn <- rownames(dat)[!is.na(dat[,"Y"])]
    #   g[rn, 3] <- dat[rn, "Y"]
    #   g[rn, 4] <- dat[rn, "Ve"]
    # }
    return(list(pred = g, Vg = Vg, Ve = Ve))
  }
}
