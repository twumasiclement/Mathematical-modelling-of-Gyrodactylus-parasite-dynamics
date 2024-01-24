

# Initial Prior and posterior distribution function
#prior distribution of model parameters (log scale)
prior<- function() {
  lb1<- runif(2, -4, 1) # birth of parasites (Gt3)
  logb11 <- max(lb1)      # birth rate for young parasites based on lb (Gt3)
  logb12<- min(lb1)      # birth rate for older parasites based on lb (Gt3)
  lb2<- runif(2, -4, 1) # birth of parasites (Gt)
  logb21 <- max(lb2)      # birth rate for young parasites based on lb (Gt)
  logb22<- min(lb2)      # birth rate for older parasites based on lb (Gt)
  lb3<- runif(2, -4, 1) # birth of parasites (Gb)
  logb31 <- max(lb3)      # birth rate for young parasites based on lb (Gb)
  logb32<- min(lb3)      # birth rate for older parasites based on lb (Gb) 
  
  ld1 <- runif(2, -5, 2) # death rates (Gt3)
  logd11 <- min(ld1)      # death rate without an immune response (Gt3)
  logd12 <- max(ld1)      # death rate with immune response (Gt3)
  ld2 <- runif(2, -5, 2) # death rates (Gt)
  logd21 <- min(ld2)      # death rate without an immune response (Gt)
  logd22 <- max(ld2)      # death rate with immune response (Gt) 
  ld3 <- runif(2, -5, 2) # death rates (Gb)
  logd31 <- min(ld3)      # death rate without an immune response (Gb)
  logd32 <- max(ld3)      # death rate with immune response (Gb) 
  
  logm<- runif(1, -4, 1) #movement rate
  logr <- runif(1, -10, 1)# immune response rate (base rate)
  logr1 <- runif(1, -10, 1)# immune response (adjustment for LA fish)
  logr2 <- runif(1, -10, 1)# immune response rate (adjustment for OS fish)
  logr3 <- runif(1, -10, 1)# immune response rate (adjustment for male fish)
  logs <- runif(1, -8, -2) #fish mortality rate (base rate)
  logs1 <- runif(1, -8, -2) #fish mortality (adjustment for male fish)
  loge1 <- runif(1, -8, log(2))  #rate of forward or backward movement/preference bias (Gt3)
  loge2 <- runif(1, -8, log(2))  #rate of forward or backward movement/preference bias (Gt)
  loge3 <- runif(1, -8, log(2))  #rate of forward or backward movement/preference bias (Gb)
  log_kappa <- runif(1, 4.5, 6.5)  #effective carrying capacity
  return(c(logb11, logb12,logb21, logb22,logb31, logb32,logd11, logd12,logd21, logd22,
           logd31, logd32, logm, logr,logr1,logr2,logr3,logs,logs1,loge1,loge2,loge3,log_kappa))
}



#Weighted-iterative ABC via
#Sequential Monte Carlo with Importance sampling
#Function for ABC calibration
ABC <- function(fork, pftn , n, w ) {
  # pftn is prior function
  # n is number of samples
  # w are weights across summary statistics for computing weighted distance
  dimS<-17 #dimension or number of ABC summary statistics
  number_of_parameters<- 23 #number of parameters to be estimated
  theta  <- matrix(nrow = n, ncol = number_of_parameters)# matrix of prior distributions
  theta_main<- list()
  #storing the summary stats across all fish and simulation realizations
  S_i <- NULL #saving saving stats of all fish for each simulation realizaton     
  #S is a matrix(nrow = n*total_fish, ncol = dimS) 
  d <- list()# weighted K-L distance
  SummaryStats_sim <- NULL
  output_without_errors<- list()

  
  for (i in 1:n) {
    print(paste("Simulation iteration=",i))
    theta[i, ] <- pftn()   
    output<- try(SimGroup_tauleap(theta1=theta[i, ],fish_sex=fishSex,fish_type=Fish_stock,
                              strain=Strain,fish_size=fishSize,error=0.002), silent = TRUE)
    
    #extract only simulations without errors
    if (!inherits(output, "try-error")) {
      k<- length(output_without_errors) + 1 #starting index (k=1) for main outputs
      theta_main[[k]]<- theta[k, ] 
      output_without_errors[[k]] <- output
    
    
    output_main<- output_without_errors[[k]]
    #B-D-C parameter estimates for the parasite-fish groups based on simulated data
    #for each simulation realisation
    BDC_estimates_sim<- GW_GMM_BDCestimator(X0=2,pop=output_main$pop_sim,
                                            output_main$alive_sim,group=parasite_fish)$BDC_estimates
    
    #Computing the summary stats for each group simulation realisation
    SummaryStats_sim[[k]] <- Summary_stats(pop=output_main$pop_sim,alive=output_main$alive_sim,
                                           BDC_estimates=BDC_estimates_sim) 
    
    #combining for all summary stats of parasite-fish groups for each simulation realisation
    SummaryStats_sim_combined<-do.call("rbind",SummaryStats_sim[[k]])
    
    
    #Combining the summary stats for the observed data for the parasite-fish groups
    SummaryStats_obs_combined<- do.call("rbind", summaries_obs)
    
    
    #Storing the summary stats and weighted distances (between summaries of observed and simulated data)
    S_i[[k]] <- SummaryStats_sim_combined
  
    d[[k]]<- w_distance(S_i[[k]],SummaryStats_obs_combined, weight=w)#weighted sum of squares distance
    #d[[k]] <- w_distance_avg(S_i[[k]],SummaryStats_obs_combined, weight=w)#weighted Euclidean distance
   
         }
    if (i %% 100 == 0) cat("fork", fork, "sample", i, "/", n, "\n")
  }
  
  
  S<-do.call("rbind",S_i) # summary stats matrix(nrow = n_complete*total_fish, ncol = dimS) 
 
  return(list(theta= do.call("rbind",theta_main), S=S, d=unlist(d)))
}

#Posterior function (draw from a new proposal distribution)
###Based on a multivariate normal kernel density


#Multivariate Normal kernel function given optimal bandwidth
#For peturbation
MultivNorm_rkernel<- function(Num,bandwidth_matrix){
  dim_k<- dim(bandwidth_matrix)[2]
  mean_vector<- rep(0,dim_k)
  return(tmvtnorm::rtmvnorm(n=1, mean=mean_vector, sigma=bandwidth_matrix, 
                lower=rep(-.1,dim_k),upper=rep(.1,dim_k), algorithm=c("gibbs")))
}





#Posterior function (draw from a new proposal distribution)
#Posterior function (draw from a new proposal distribution)
post <- function(samp=tha_post,importance_weight=weight,
                 optimal_bw_matrix=Sigma_optimal_t) {
  # new proposal from samp (previous prior samples)
  n <- dim(samp)[1]#dimension of proposal samples at t-1
  #randomly sampling from 1:n different sets of parameter values (importance sampling)
  #from previously accepted particles
  sample.particle<-sample(n, 1,prob=importance_weight)
  
  # Perturbation kernel: perturbing the particles for a new proposal  
  KDE_sampler<- samp[sample.particle, ]+MultivNorm_rkernel(Num=1,
                                          bandwidth_matrix=optimal_bw_matrix) 
  
  new_proposal<- KDE_sampler
  x<- new_proposal
  #birth rates (young and old parasites): birth rate of young>old
  x[1:2] <- sort(x[1:2], decreasing=TRUE);x[3:4] <- sort(x[3:4], decreasing=TRUE)
  x[5:6] <- sort(x[5:6], decreasing=TRUE) 
  #death rates (without and with immune response):death rate of without response< with response
  x[7:8] <- sort(x[7:8],decreasing=FALSE);  x[9:10] <- sort(x[9:10],decreasing=FALSE)
  x[11:12] <- sort(x[11:12],decreasing=FALSE)
  return(x)
}






