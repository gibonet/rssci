
# Provo una versione che stima gli standard error della ecdf a più livelli
# in una volta sola
# Parto da una copia di sep4 (da sep4.R)

#' standard error of the empirical cumulative distribution 
#' function (of wages, for example) at level \code{probs} 
#' (the quantile level, comprised between 0 and 1)
#' 
#' @return a vector with the standard errors of the ecdf at level probs. probs
#' can be a vector of length greater than 1 (of values 0 < probs < 1)
#' 
#' @inheritParams statQuantile
#' @export
sep6 <- function(x, probs = 0.5, strata = NULL, psu = NULL, nh = NULL, 
                 th = NULL, Nh = NULL, mhi = NULL, thi = NULL, Mhi = NULL, 
                 weights = NULL,  crit = NULL){
  
  stra <- tapply(strata, psu, `[`, 1) 
  nh <- tapply(nh, strata, `[`, 1)
  th <- tapply(th, strata, `[`, 1)
  Nh <- nh / th
  mhi <- tapply(mhi, psu, `[`, 1)
  
  Mhi <- tapply(Mhi, psu, `[`, 1)
  thi <- tapply(thi, psu, `[`, 1)
  
  med <- computeQuantiles(x, weights, probs)
  
  med <- matrix(rep(med, each = length(x)), 
                nrow = length(x), ncol = length(med))
  
  probs <- matrix(rep(probs, each = length(x)),
                  nrow = length(x), ncol = length(probs))
  
  # for each SSU:
  zhij <- 1 * (x <= med)

  ej <- weights * (zhij - probs)

  ehi <- rowsum(ej, psu)  # Ev. reorder = FALSE

  NDhi <- matrix(rep(1, length(x)), nrow = length(x), byrow = FALSE)
  NDhi <- rowsum(NDhi, psu)
  NDhi <- matrix(rep(NDhi, ncol(med)), ncol = ncol(med), byrow = FALSE)
  
  ej_psu <- ehi / NDhi

  ehi2 <- rowsum(ej^2, psu)
  
  Bhi <- (ehi2 - NDhi * ej_psu^2) / (NDhi - 1)

  mhi <- matrix(rep(mhi, ncol(med)), ncol = ncol(med), byrow = FALSE)
  thi <- matrix(rep(thi, ncol(med)), ncol = ncol(med), byrow = FALSE)
  
  Bhi <- ((NDhi - 1) * Bhi + NDhi * (1 - NDhi / mhi) * (ehi / NDhi)^2) *
    (1 - thi) * mhi / (mhi - 1)
  
  Bhi[thi == 1] <- 0

  Bh <- rowsum(Bhi, stra, na.rm = TRUE)

  eh <- rowsum(ehi, stra)
  
  n_Ah <- matrix(rep(1, nrow(ehi)), nrow = nrow(ehi), byrow = FALSE)
  n_Ah <- rowsum(n_Ah, stra)
  n_Ah <- matrix(rep(n_Ah, ncol(med)), ncol = ncol(med), byrow = FALSE)
  
  ehi_stra <- rowsum(ehi, stra)
  Ah_stra <- ehi_stra / n_Ah
  
  ehi2_stra <- rowsum(ehi^2, stra)
  
  Ah <- (ehi2_stra - n_Ah * Ah_stra^2) / (n_Ah - 1)

  toth <- rowsum(NDhi, stra)
  dlh <- rowsum(NDhi - 1, stra)
  
  ne <- toth - dlh

  k <- (ne > 1)
  nh <- matrix(rep(nh, ncol(med)), ncol = ncol(med), byrow = FALSE)
  th <- matrix(rep(th, ncol(med)), ncol = ncol(med), byrow = FALSE)
  Ah[k] <- ((ne[k] - 1) * Ah[k] + ne[k] * (1 - ne[k] / nh[k])
            * (eh[k] / ne[k])^2) / (nh[k] - 1)

  Ah[ne == 1 & abs(th - 1) < 0.0000001] <- 0
  Ah <- Ah * nh * (1-th)
  V2sth <- ifelse (abs(th - 1) < 0.0000001, Bh, Ah + th *Bh)
  V2sth
  
  svhi <- rowsum(weights, psu)

  svh <- rowsum(svhi, stra)

  denom <- sum(svh)
  SV2st <- colSums(V2sth, na.rm = TRUE) / denom^2
  sep <- sqrt(SV2st)
  return(sep)
}

