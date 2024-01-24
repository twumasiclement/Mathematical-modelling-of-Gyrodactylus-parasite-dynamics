
library(compiler)

#external scripts
#Script of functions for Galton-Watson & GMM estimation of the B-D-C parameters
source("MLE_catastrophe-script.r")
source("GMM-1st2nd-Steps-script.r")
source("BDC-GW-GMM-estimator-script.r")

#Script for population projection until day 17 after fish death to aid with ABC fitting
source("Project-Parasite-script.r")

na.one<- function (x) { #in case of 0/0 weight, set at 1
  #when population is zero we get NaN's in wa and/or wb, 
  #but wasserstein1d returns 0 anyway (including when wa=wb=1), which is what we want
  x[is.na(x)] <- 1
  return(x)
}



summary_func <- function(pop_single, alive_single, ga = 0.9,group_number,BDC_estimates) {
  #ga = 0.9 is tunning parameter
  pop_single<- project(pop_single=pop_single, alive_single=alive_single, ga=ga)
  S <- rep(NA, 17)
  S[1:9] <- log(1 + colSums(pop_single,na.rm=T)) #across time
  for (i in 1:4) { # when population is zero we get NA or NaN's in wa and/or wb, 
    #but wasserstein1d returns 0 anyway, which is what we want
    S[9+i] <- wasserstein1d(1:4, 1:4,      
                            wa=na.one(pop_single[,i]/sum(pop_single[,i]))+.001, 
                            wb=na.one(pop_single[,i+1]/sum(pop_single[,i+1]))+.001)
    
    
  }
  #finding the B-D-C parameter estimates depending on the parasite-fish group
  #using Galton-Watson approach and/or GMM estimator
  S[14:16]<- abs(BDC_estimates[group_number, ])# B-D-C parameter estimates b
  S[17] <- sum(alive_single == 1)  #time before death of fish
  return(S)
}

summary_compiler=cmpfun(summary_func)



na.inf.zero<- function(x){
  x[is.na(x)|is.finite(x)==FALSE]<- 0
  return(x)
}
