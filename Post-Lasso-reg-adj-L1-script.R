require("kedd")

#Gaussian kernel with bandwidth delta
guass_kernel<- function(dist, delta) {
  kern <- (1 / (sqrt(2 * pi) * delta)) * exp(- (dist^2) / (2 * delta^2))
  return(kern)
}


Post_Lasso_reg_adj<- function(post_distn,summary_obs) {
  #k #biasing parameter or lasso regression penalty parameter
  # post_dtn is the posterior sample
  # w are weights across summary statistics for computing weighted distance
  #storing the summary stats across all fish and simulation realizations
  S_i <- NULL #saving saving stats of all fish for each simulation realizaton 
  no_of_parameters<- 23
  posterior_mean_adj<- rep(NA,no_of_parameters) #storing adjusted posterior means
  #Combining the summary stats for the observed data for the parasite-fish groups
  m<- dim(post_distn)[1]# m=number of posterior samples
  d<- rep(0, m)# weighted distances given observed data
  p<- dim(summary_obs)[2] #dimension of summary statistics
  Unadj_dist<- post_distn
  SummaryStats_sim <- NULL
  X_Design_matrix<- matrix(NA, ncol=p,nrow=m) #design matrix
  # Weights based on Gaussian kernel for local-linear regression adjustment
  W<- matrix(0, ncol=m,nrow=m)
  X_bar=numeric(length=p) #saving weighted column means of design matrix 
  beta_lasso<-list()#store regression coefficients for each dependent variables
  
  for (i in 1:m) {
    theta<- as.vector(unlist(post_distn[i,]))
    
    output_sim<- SimGroup_tauleap(theta1=theta,fish_sex=fishSex,fish_type=Fish_stock,
                                  strain=Strain,fish_size=fishSize,error=0.002)
    #B-D-C parameter estimates for the parasite-fish groups based on simulated data
    #for each simulation realisation
    BDC_estimates_sim<- GW_GMM_BDCestimator(X0=2,pop=output_sim$pop_sim,
                                            output_sim$alive_sim,group=parasite_fish)$BDC_estimates
    
    #Computing the summary stats for each group simulation realisation
    SummaryStats_sim[[i]] <- Summary_stats(pop=output_sim$pop_sim,alive=output_sim$alive_sim,
                                           BDC_estimates=BDC_estimates_sim) 
    
    #combining for all summary stats of parasite-fish groups for each simulation realisation
    SummaryStats_sim_combined<-do.call("rbind",SummaryStats_sim[[i]])
    mean_diff<- apply(SummaryStats_sim_combined-summary_obs,2,mean,na.rm = TRUE)
    X_Design_matrix[i, ]<- mean_diff  #storing each row of design matrix X
    # Computing weights based on
    #Storing the average summary stats and distance (between summaries of observed and simulated data)
    S_i[[i]] <- SummaryStats_sim_combined
    w <-apply(S_i[[i]], 2, var, na.rm = TRUE)
    w<- w/sum(w) #normalised weight for computing distance
    d[i] <- w_distance(S1=S_i[[i]], S2=summary_obs, weight=w) 
  }
  
  distances<-na.inf.zero(d)
  #Adaptively choosing the bandwidth of the Gaussian kernel based on the distances
  #bandwidth<- KernSmooth::dpik(x=distances, "normal",scalest = "minim")
  bandwidth<- kedd::h.amise(x=distances, deriv.order =0,kernel = c("gaussian"))$h
  diag(W)<- guass_kernel(dist= distances,delta=bandwidth)
  theta_post<- as.matrix(post_distn)  
  weights<- diag(W)/sum(diag(W)) # (normalising) main diagonal of Weighting matrix
  
  
  #Transforming X and Y (posterior distribution and summary statistics)
  for(j in seq_along(posterior_mean_adj)){#For each jth model parameter, j=1,2,...23
    X<- X_Design_matrix
    Y<- theta_post[ ,j]
    
    #Step 1 (Centering X and Y)
    for (k in 1:p) X_bar[k]<- sum(weights*X[,k])
    
    for (k in 1:p) {
      X[, k]<- X[, k]-X_bar[k]
    }
    #finding the weighted mean of Y and centring
    Y_bar <- sum(weights*Y)
    Y<- Y- Y_bar
    
    
    #Step 2: scaling (centred X and Y) by weights
    
    for(k in 1:p) X[, k]<- sqrt(weights)*X[, k]
    Y<- sqrt(weights)*Y
    
    #Choose optimal value of k (the penalty paramters)
    # Using cross validation glmnet
    # Setting the range of lambda values
    options(warn = -1)
    #lambda_seq <- 10^seq(2, -2, by = -.1)
    #lasso_cv <- cv.glmnet(X, Y, alpha = 1, lambda =lambda_seq)
    lasso_cv <- cv.glmnet(X, Y, alpha = 1)
    # Best lambda value
    best_lambda <- lasso_cv$lambda.min
    k<-best_lambda
    #print(k)
    
    
    
    # calculate beta estimates corresponding to summary statistics X (standardised coefficients)
    
    # calculating beta estimates of predictors (without the intercept)
    best_model <- glmnet(X, Y, alpha = 1, lambda = k,intercept = FALSE)
    beta_lasso[[j]]<-   as.vector(coef(best_model))[-1]#no intercept
    alpha_estimate<- Y_bar - X_bar%*%beta_lasso[[j]]
    # calculate intercept estimates (adjusted posterior mean estimates)
    # Tranforming to original scale (due posterior samples were on log scale)
    posterior_mean_adj[j] <- exp(alpha_estimate)
    
    
    #Adjusting the posterior distribution
    Unadj_dist[,j]<- post_distn[, j]- X_Design_matrix%*%beta_lasso[[j]]
  }
  

  posterior_mean_unadj<- apply(exp(post_distn),2,mean)
  Posterior_mean_output<- data.frame(Adj_posterior_mean=posterior_mean_adj,
                      Unadj_posterior_mean=posterior_mean_unadj)
  
  
  return(list(X_Design_matrix=X,Posterior_mean_output=Posterior_mean_output,Adjusted_posterior_dist=Unadj_dist))
}
