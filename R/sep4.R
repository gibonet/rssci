
# Come sep3 ma sostituisco tutti i tapply(..., unique) con tapply(..., `[`, 1)

#' Estimates the standard error of the empirical cumulative distribution 
#' function (of wages, for example) at level \code{probs} 
#' (the quantile level, comprised between 0 and 1)
#' 
#' @inheritParams statQuantile
#' @export
sep4 <- function(x, probs = 0.5, strata = NULL, psu = NULL, nh = NULL, 
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
  
  # for each SSU:
  zhij <- 1 * (x <= med)
  ej <- weights * (zhij - probs)  
  
  ehi <- tapply(ej, psu, sum)
  Bhi <- tapply(ej, psu, stats::var)
  # number of SSU in the domain
  NDhi <- tapply(rep(1, length(x)), psu, sum)                  
  Bhi <- ((NDhi - 1) * Bhi + NDhi * (1 - NDhi / mhi) * (ehi / NDhi)^2) * 
    (1-thi) * mhi / (mhi - 1)
  # avoids NA if thi=1
  Bhi[thi == 1] <- 0                                             
  # sum of the weights in the PSU
  svhi <- tapply(weights, psu, sum)                            
  
  Bh <- tapply(Bhi, stra, function(v) {
    if (length(v) == 1) return(v)
    if (length(v) > 1) return(sum(v[!is.na(v)]))
  })
  
  eh <- tapply(ehi, stra, sum)
  Ah <- tapply(ehi, stra, stats::var)
  toth <- tapply(NDhi, stra, sum)                              
  dlh <- tapply(NDhi - 1, stra, sum)
  # number of SSU in the domain
  ne <- toth - dlh
  k <- (ne > 1)
  ne_cond <- ne[k]
  Ah[k] <- ((ne_cond - 1) * Ah[k] + ne_cond * (1 - ne_cond / nh[k])
                 * (eh[k] / ne_cond)^2) / (nh[k] - 1)
  # sum of the weights in the stratum
  svh <- tapply(svhi, stra, sum)
  
  # avoids NA if th = 1
  Ah[ne == 1 & abs(th - 1) < 0.0000001] <- 0
  Ah <- Ah *nh * (1-th)
  
  # B5 |-----------------------------------------------------------------------------------|
  # final computations
  # computes the variances of the strata
  V2sth <- ifelse (abs(th - 1) < 0.0000001, Bh, Ah + th *Bh)
  # total sum of weights
  denom <- sum(svh)
  # the global variance
  SV2st <- sum(V2sth[!is.na(V2sth)]) / denom^2
  # the standard error
  sep <- sqrt(SV2st)
  return(sep)
}
