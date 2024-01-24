
library(compiler)

##**Function for population projection until day 17 after fish death**
#ga= gamma which is tuning parameter
project <- function(pop_single, alive_single, ga) {

  # project parasite numbers beyond fish mortality
  n <- length(alive_single)
  k <- sum(alive_single == 1)
  if (k == n) return(pop_single)
  if (k == 0) return(matrix(0, nrow=4, ncol=n))
  if (k == 1) return(matrix(pop_single[,1], nrow=4, ncol=n))
  z <- log(colSums(pop_single[,1:k],na.rm=T))            
  al <- sum( (z[k] - z[1:(k-1)]) * ((k-1):1) * ga^((k-1):1),na.rm=T) / sum( ((k-1):1)^2 * ga^((k-1):1),na.rm=T)
  pop_single[,(k+1):n] <- pop_single[,k] %*% t( exp( (1:(n-k))*al ) )
  return(pop_single)
}

project_compiler=cmpfun(project) #converting function to byte-code compilation
