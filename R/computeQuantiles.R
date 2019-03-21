
#' calculates weighted quantiles
#' 
#' @param xx un vettore numerico (esempio: salari - MBLS nella RSS)
#' @param ww vettore di pesi positivi (esempio: nella RSS, GEWIBGRS)
#' @param qq vettore di numeri compresi tra 0 e 1, che indicano il o i livelli di quantile che si vogliono stimare (default 0.5: mediana)
#' @export
computeQuantiles <- function(xx, ww, qq = 0.5){
  
  # if no weights have been specified, returns non 
  # weighted quantiles
  if (missing(ww))
    return(stats::quantile(xx,probs=qq,na.rm=T))
  # |-----------------------------------------------------------------------------------|
  # otherwise, computes the partial sums of ww
  ord <- order(xx)
  # cum.w <- cumsum(ww[ord])[!is.na(xx)]/sum(ww[!is.na(xx)]) !vecchio codice errato!
  cum.w <- cumsum(ww[ord])[!is.na(xx[ord])]/sum(ww[!is.na(xx)])
  tmpS <- data.frame(matrix(rep(NA,2*length(qq)),nrow=2))
  tmpO <- data.frame(matrix(rep(NA,2*length(qq)),nrow=2))
  res <- c(rep(NA, length(qq)))
  # |-----------------------------------------------------------------------------------|
  # and computes each quantile
  for (i in 1:length(qq)) {
    # records the two sums directly greater than qq
    tmpS[i] <- cum.w[cum.w>=qq[i]][1:2]
    # and the corresponding orders (*)
    tmpO[i] <- ord[cum.w>=qq[i]][1:2]
    # if a sum is equal to qq, returns the mean of the two xx 
    # corresponding to (*), otherwise, the lowest
    res[i] <- (ifelse(abs(tmpS[1,i]-qq[i])<1e-010, 
                      mean(xx[tmpO[,i]]), xx[tmpO[1,i]]))
  }
  return(res)
}

