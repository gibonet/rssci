# sep: estimated standard error of the empirical cumulative distribution 
# function (at level probs %in% ]0;1[)  

#' Estimates the standard error of the empirical cumulative distribution 
#' function (of wages, for example) at level \code{probs} 
#' (the quantile level, comprised between 0 and 1)
#' 
#' @inheritParams statQuantile
#' @export
sep <- function(x, probs = 0.5, strata = NULL, psu = NULL, nh = NULL, 
                th = NULL, Nh = NULL, mhi = NULL, thi = NULL, Mhi = NULL, 
                weights = NULL,  crit = NULL){
  
  # definition of a counting function
  compt <- function(variable) {return(length(unique(variable)))}
  
  # B1 |-----------------------------------------------------------------------------------|
  # if no stratum (or no psu), assumes there is only one stratum
  # (or psu)
  if (is.null(strata))
    strata <- rep(1, length(x))
  if (is.null(psu))
    psu <- rep(1, length(x))
  stra <- tapply(strata, psu, unique) 
  
  # B1.1 | {nh, Nh, th} sono 3 variabili in ogni caso sovrascritte !!!---------------------|
  # checks the data: PSU level
  ifelse (is.null(nh),
          # assumes domain corresponds to strata
          nh <- tapply(psu, strata, compt),
          nh <- tapply(nh, strata, unique)
  )
  ifelse (!is.null(Nh),
          {ifelse (!is.null(th),
                   {Nh <- tapply(Nh, strata, unique)
                   th <- tapply(th, strata, unique)},
                   {Nh <- tapply(Nh, strata, unique)
                   th <- nh/Nh}
          )},
          {ifelse (!is.null(th),
                   {th <- tapply(th, strata, unique)
                   Nh <- nh/th},
                   # assumes total sampling
                   {th <- rep(1, compt(stra))
                   Nh <- nh}
          )}
  )
  
  # B1.2 | {mhi, Mhi, thi} sono 3 variabili in ogni caso sovrascritte !!!------------------|
  # checks the data: SSU level
  ifelse (is.null(mhi),
          # assumes domain corresponds to strata
          # !!! questo comando non è corretto !!! vedi nota 4.
          # mhi <- tapply(x, psu, compt),
          # !!! sostituito con il seguente comando
          mhi <- tapply(x, psu, length),
          mhi <- tapply(mhi, psu, unique)
  )
  ifelse (!is.null(Mhi),
          {ifelse (!is.null(thi),
                   {Mhi <- tapply(Mhi, psu, unique)
                   thi <- tapply(thi, psu, unique)},
                   {Mhi <- tapply(Mhi, psu, unique)
                   thi <- mhi/Mhi}
          )},
          {ifelse (!is.null(thi),
                   {thi <- tapply(thi, psu, unique)
                   Mhi <- mhi/thi},
                   # assumes total sampling
                   {thi <- rep(1, compt(psu))
                   Mhi <- mhi}
          )}
  )
  
  # checks the data: weights
  if (is.null(weights))
    # computes sampling weights
    # !!! questo comando non è corretto !!! 
    # La lunghezza dei due vettori (th) e (thi) è diversa, codice non compilato
    weights <- 1/(th*thi)
  
  # B2 |-----------------------------------------------------------------------------------|
  #####################################################################
  # INSERIMENTO DELL'ARGOMENTO probs (SPERIMENTALE)
  # computation of the median
  med <- computeQuantiles(x, weights, probs)
  
  # for each SSU:
  zhij <- 1 * (x <= med)
  ej <- weights*(zhij - probs)  
  # sostituito .5 con probs (vedi pag. 22 del rapporto metodologico 
  # "Enquête suisse sur la structure des salaires 2000 / Plan d'échantillonage, 
  # pondération et méthode d'estimation pour le secteur privé", Monique Graf (2002))
  ##################################################################### 
  
  # B3 |-----------------------------------------------------------------------------------|
  # computation of the intra-PSU variance for each PSU:
  ehi <- tapply(ej, psu, sum)
  Bhi <- tapply(ej, psu, stats::var)
  # number of SSU in the domain
  NDhi <- tapply(rep(1, length(x)), psu, sum)                  
  Bhi <- ((NDhi-1)*Bhi+NDhi*(1-NDhi/mhi)*(ehi/NDhi)^2)*
    (1-thi)*mhi/(mhi-1)
  # avoids NA if thi=1
  Bhi[thi==1] <- 0                                             
  # sum of the weights in the PSU
  svhi <- tapply(weights, psu, sum)                            
  
  # sum of the intra-PSU variances
  # if there is only one psu in the stratum, and only one Bhi, returns 
  # this Bhi; otherwise, returns the sum of the Bhi in the stratum 
  # which are different from NA
  Bh <- tapply(Bhi, stra, function(v) {
    if (length(v)==1) return(v)
    if (length(v)>1) return(sum(v[!is.na(v)]))
  })
  
  # B4 |-----------------------------------------------------------------------------------|
  # computation of the inter-PSU variances for each stratum
  eh <- tapply(ehi, stra, sum)
  Ah <- tapply(ehi, stra, stats::var)
  toth <- tapply(NDhi, stra, sum)                              
  dlh <- tapply(NDhi-1, stra, sum)
  # number of SSU in the domain
  ne <- toth - dlh                                 
  Ah[ne>1] <- ((ne[ne>1]-1)*Ah[ne>1] + ne[ne>1]*(1-ne[ne>1]/nh[ne>1])
               *(eh[ne>1]/ne[ne>1])^2) / (nh[ne>1]-1)
  # sum of the weights in the stratum
  svh <- tapply(svhi, stra, sum)
  
  # checks if the imputation has to be done and, if it is necessary, 
  # does it
  if (!is.null(crit) & any(is.na(Ah)) & length(Ah[Ah!=NA])>0) {
    # computes a relative Ah variance
    Ahrel <- Ah/svh^2
    crit <- tapply(crit, strata, unique)
    # takes the mean of the relative variances for each value of crit
    mAhrel <- tapply(Ahrel[!is.na(Ahrel)], crit[!is.na(Ahrel)], mean)
    # proceeds to the imputation
    for (i in 1:length(Ah[ne==1])) {
      lab <- names(Ah[ne==1][i])
      lab <- crit[names(crit)==lab]
      # imputes the computed value if it exists
      if (lab %in% names(mAhrel))
        Ah[ne==1][i] <- svh[ne==1][i]^2*
        mAhrel[names(mAhrel)==lab]
    }
  }
  
  # avoids NA if th = 1
  Ah[ne==1 & abs(th-1)<0.0000001] <- 0
  Ah <- Ah*nh*(1-th)
  
  # B5 |-----------------------------------------------------------------------------------|
  # final computations
  # computes the variances of the strata
  V2sth <- ifelse (abs(th-1)<0.0000001, Bh, Ah+th*Bh)
  # total sum of weights
  denom <- sum(svh)
  # the global variance
  SV2st <- sum(V2sth[!is.na(V2sth)])/denom^2
  # the standard error
  sep <- sqrt(SV2st)
  return(sep)
}


#   #the 95% confidence interval
#   ###################################################################
#   # INSERIMENTO DELL'ARGOMENTO probs (SPERIMENTALE)
#   l.limit <- computeQuantiles(x, weights, probs - 1.96 * sep)
#   u.limit <- computeQuantiles(x, weights, probs + 1.96 * sep)
#   #the three coefficients of variation
#   cvs95 <- 100*max(med-l.limit, u.limit-med)/(1.96*med)
#   cl <- computeQuantiles(x, weights, probs-sep)
#   cu <- computeQuantiles(x, weights, probs+sep) 
#   cvs <- 100*max(med-cl, cu-med)/med
#   CVperc <- 100 * sep / probs  # QUA NON SO SE SOSTITURIRE 0.5 con probs 
#   (sì, vedi rapporto metodologico, in fondo a pag. 8)
#   #the number of strata, PSU and SSU
#   Nstrata <- compt(strata)
#   Npsu <- compt(psu)
#   Nssu <- length(x)
#   
#   return(cbind(l.limit, u.limit, quantile=med, cvs95, cvs, 
#                CVperc, Nstrata, Npsu, Nssu))
# }
