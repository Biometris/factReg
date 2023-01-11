evalAcc <- function(obs, pred, genotypes = NULL) {
  rownames(obs)  <- paste0(as.character(obs$E),'_X_',as.character(obs$G))
  rownames(pred) <- paste0(as.character(pred$E),'_X_',as.character(pred$G))
}
