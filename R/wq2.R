
#' Empirical weighted quantile (without ordering)
#' 
#' This function has to be applied to already sorted data!!!
#' 
#' @param x A numeric vector
#' @param weights A vector of (positive) sample weights
#' @param probs a numeric vector with the desired quantile levels (default 0.5, the median)
#' @return The weighted quantile (a numeric vector)
#' @references Ferrez, J., Graf, M. (2007). Enquète suisse sur la structure des
#'  salaires. Programmes R pour l'intervalle de confiance de la médiane. 
#'  (Rapport de méthodes No. 338-0045). Neuchâtel: Office fédéral de statistique.
#'  
#' @examples
#' wq2(x = rnorm(100), weights = runif(100))
#' @export
wq2 <- function(x, weights, probs = c(0.5)){
  ord <- 1:length(x)
  cum.w <- cumsum(weights)/sum(weights)
  # tmpS <- tmpO <- data.frame(matrix(rep(NA, 2 * length(probs)), nrow = 2))
  res <- c(rep(NA, length(probs)))
  for (i in 1:length(probs)){
    k <- cum.w >= probs[i]
    # tmpS[i] <- cum.w[k][1:2]
    # tmpO[i] <- ord[k][1:2]
    # res[i] <- (ifelse(abs(tmpS[1, i] - probs[i]) < 1e-10, 
    #                   mean(x[tmpO[, i]]), x[tmpO[1, i]]))
    res[i] <- x[k][1]
  }
  return(res)
}

