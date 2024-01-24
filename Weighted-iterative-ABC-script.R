


Weighted_iterative_ABC<- function(N=1000,dimS=17,fish_total=Total_fish,
                                  numCores=numCores,ABC_time_steps=10,seed_num=1){

  # N= total number of samples
  n_cores <- numCores   #Run on  n cores
  n<- N/n_cores
  number_of_parameters<- 23 #number of parameters to be estimated
  #Storing importance weight for sequential sampling
  import_weights<- NULL
  #Storing weights corresponding to accepted samples
  w_accepted<- NULL
  dim_tha_post<-NULL #saving number of particles for each iteration
  S_i <- NULL #saving saving stats of all fish for each simulation realizaton  
  #ABC_time_steps= time for the algorithm to terminate
  eps<-NULL # storage for index of accepted particles
  #proportion of sample to retain during sequential sampling (20 steps)
  epsilon<- c(0.5,0.43,0.4,0.35,0.3,0.2,0.1,0.08,0.06,0.02)
  
  
  #if(N<1000){
  #  epsilon<- c(0.5,0.43,0.4,0.35,0.3,0.2,0.1,0.08,0.06,0.02)
  #}else if(N>=1000){
  #  epsilon<-c(0.5,0.3,0.2,0.1,0.08,0.07,0.06,0.03,0.02,0.01)
  #}
  
  
  d_i<-NULL;d<-NULL#storing weighted distances between summary stats at time t
  
  #For storing parameter values at time t
  theta_i<- NULL;theta<-NULL;theta_main<-NULL
  
  
  # for density plots (256 used here is 
  #the number of equally spaced points at which the density is to be estimated)
  x <- seq(from = -10, to = 7, length.out = 256)#range of prior distribution
  fx <- array(dim=c(ABC_time_steps+1, number_of_parameters, 256))
  
  #Setting seed number for reproducibility of results
  set.seed(seed_num)
  
  
  time0<- proc.time()
  for (t in 1:ABC_time_steps) {
    cat("ABC_time_steps", t, "\n")
    if (t == 1) {
      pftn <- prior
      #ABC_out <- mclapply(1:n_cores, ABC_without_w, pftn=pftn, n=n, mc.cores=n_cores)
      ABC_out <- mclapply(1:n_cores, ABC, pftn=pftn, n=n, w=w, mc.cores=n_cores)
      for (i in 1:n_cores) {
        theta_i[[i]] <- ABC_out[[i]]$theta
        S_i[[i]] <- ABC_out[[i]]$S
        d_i[[i]] <- ABC_out[[i]]$d 
      }
      
      theta[[1]]<- do.call("rbind",as.list(theta_i))# N by dimS matrix
      theta_main[[1]]<- na.omit(theta[[1]]) #remove NAs if any
      N_monte_carlo<- dim(theta_main[[1]])[1]
      print(paste("N_monte_carlo sample size at time t=",t,"","is",N_monte_carlo))
      
      
      
    }else{#if t>1
      #Calculating optimal value of the MVN kernel bandwidth matrix
      #N0= number of accepted particles
      #N1= total number of proposal samples
      #eps[[t]]= index of accepted samples
      #w_accepted[[t]]= weights corresponding to accepted particles
      #tha_post= accepted proposals at time t
      #eigenMatMult is a faster matrix multiplication function in C+
      #Optimal bandwith matrix for multivariate normal perturbation kernel
      #Sigma_optimal_t=Optimal bandwith matrix
      # using the Kullback--Leibler divergence minimisation approach
      #N1<-dim(theta[[t-1]])[1]
      #N0<- dim(tha_post)[1]
      #i <-  1:N1
      #k <- 1:N0
      #iandk<-expand.grid(i=i,k=k)
      #func_Sigmal_optimal<- function(rowId){
      #  i_index<- iandk$i[rowId];k_index<- iandk$k[rowId]
      #  A<-matrix(tha_post[k_index,]-theta[[t-1]][i_index, ])
      #  B<-t(matrix(tha_post[k_index,]-theta[[t-1]][i_index, ]))
        
       # return((import_weights[[t-1]][i_index]*w_accepted[[t-1]][k_index]*  
       #           (eigenMatMult(A, B, n_cores = n_cores))))
     # }  
     # Sigma_optimal_list<- lapply(1:nrow(iandk),func_Sigmal_optimal)
     # Sigma_optimal_t<- Reduce('+', Sigma_optimal_list)#matrix sum
      Sigma_optimal_t<- matrix(0,nrow=number_of_parameters,ncol=number_of_parameters)
      
      N1<-dim(theta_main[[t-1]])[1]
      N0<- dim(tha_post)[1]
      for(i in 1:N1) {
        for(k in 1:N0){
          Sigma_optimal_t<- Sigma_optimal_t+(import_weights[[t-1]][i]*w_accepted[[t-1]][k]*  
          (matrix(tha_post[k,]-theta_main[[t-1]][i, ])%*%t(matrix(tha_post[k,]-theta_main[[t-1]][i, ]))))
        }     
      }  
      
      
      
      
      
     ##Sampling from the Multivariate Normal Perturbation kernel
      weight<-w_accepted[[t-1]]
      pftn <- function() post(tha_post,weight,Sigma_optimal_t)
      #ABC_out <- mclapply(1:n_cores, ABC_without_w, pftn=pftn, n=n, mc.cores=n_cores)
      ABC_out <- mclapply(1:n_cores, ABC, pftn=pftn, n=n, w=w, mc.cores=n_cores)
      for (i in 1:n_cores) {
        theta_i[[i]] <- ABC_out[[i]]$theta
        S_i[[i]] <- ABC_out[[i]]$S
        d_i[[i]] <- ABC_out[[i]]$d 
      }
      #Combining theta at time t
      theta[[t]]<- as.matrix(do.call("rbind",as.list(theta_i)))# N by dimS=23 matrix
      theta_main[[t]]<- na.omit(theta[[t]]) #remove NAs if any
      print(paste("N_monte_carlo sample size at time t=",t,"","is",dim(theta_main[[t]])[1]))
      theta_main[[t]]<- theta_main[[t]][1:N_monte_carlo, ]
      
    
      
      
      #Re-weighting for importance sampling
      
      import_weights[[t]]<-rep(NA,length=N_monte_carlo)
      
    
      #Evaluating the perturbation kernel for each particle at time t
      dMVN_func<- function(i) mvtnorm::dmvnorm(x=theta_main[[t]][i, ], 
                                    mean = theta_main[[t-1]][i, ],sigma =Sigma_optimal_t)
      K_normal_kernel<- mclapply(1:dim(theta_main[[t]])[1],dMVN_func,mc.cores=n_cores)
      
      
      #### KDE of proposal probability distribution####
      #Estimating the optimal bandwidth
      density_proposals<- matrix(NA, nrow=length(import_weights[[t]]),ncol=number_of_parameters)
      N1<-length(unlist(K_normal_kernel))
      for(i in seq_along(import_weights[[t]])){
        
        density_proposals[i, ]<- ks::kde(x = theta_main[[t]][i, ],
                                         eval.points = theta_main[[t]][i, ])$estimate
        #KDE value for each proposal sample
        par.weight.numerator<-   mean(density_proposals[i, ])
        par.weight.denominator<- sum(import_weights[[t-1]][1:N1]*unlist(K_normal_kernel))
        import_weights[[t]][i]<- par.weight.numerator/par.weight.denominator                               
      }
      
      
      #normalizing weights
      import_weights[[t]]<- import_weights[[t]]/sum(import_weights[[t]]) 
      print(paste("Length of importance weight at time t=",t))
      print(length(import_weights[[t]]))
      
    }
    ##################################################################################      
    #Combining results from the ncores
    
    #Saving theta_main
    write.csv(theta_main[[t]], paste0("theta_t_",t,".csv")) 
    
    d[[t]]<-  do.call("c",as.list(d_i)) #length of N
    
    small_draws<- epsilon[t]*N #number of draw for posterior samples
    theta_dist<- cbind(theta_main[[t]],d[[t]])#adding the computed distance as extra column of theta matrix
    theta_dist<- na.omit(theta_dist) #omit realisation of NA distances
    eps[[t]]<- order(theta_dist[,24])[1:small_draws]#smallest distance index
    
    # choose posterior samples
    tha_post<-theta_dist[eps[[t]],][,-24]
    
    dim_tha_post[[t]]<- dim(tha_post)[1]
    #initialize importance weight for sequential sampling
    if(t==1)  import_weights[[1]]<- rep(1/N_monte_carlo,length=N_monte_carlo)
    #Weights corresponding to accepted proposal samples
    w_accepted[[t]]<- import_weights[[t]][eps[[t]]]
    w_accepted[[t]]<- w_accepted[[t]]/sum(w_accepted[[t]])#normalising accepted weights
    
    # update summary statistics weights (different from importance sampling weights) for weighted 
    #sum of squares distance between S_sim and S_obs (summary statistics)
    eps_dist_max <- sort(d[[t]])[small_draws]#max least distance
    S<-na.omit(do.call("rbind",S_i))#combining the summary stats [(N*fish_total) by 17 matrix]
    
    #current summary weights at time=t
    ## Summary weights based on inverse of variance
    w1inv <- apply(S[rep(d[[t]],fish_total)<=eps_dist_max , ], 2, var, na.rm = TRUE)#inverse weight
    #############Summary weights based on inverse of mean absolute deviation
     #w1inv <- apply(S[rep(d[[t]],fish_total)<=eps_dist_max , ], 2, madstat, na.rm = TRUE)#inverse weight
     w <- na.zero(2/(1/w + w1inv))# computing the geometric mean of current & previous summary weights
    
    
    # densities
    if (t == 1) {
      for (k in 1:number_of_parameters) { 
        fx[1,k,] <- density(theta_main[[1]][ ,k], from=-10, to=7, n=256)$y
        #saving the densities for each iteration
        write.csv(fx[1, ,], file = paste0("density_post_", 1,"_",N, ".csv"))
      }
    }
    for (k in 1:number_of_parameters) {
      fx[t+1,k,] <- density(tha_post[, k], from=-10, to=7, n=256)$y
      
      #saving the densities for each iteration
      write.csv(fx[t+1, ,], file = paste0("density_post_", t+1,"_",N,".csv"))
      
    }
    #saving importance weights 
    write.csv(import_weights[[t]], file = paste0("importance_weights_",t,"_",N,".csv")) 
    
    #accepted particles at each iteration
    write.csv(tha_post, file = paste0("theta_post_", t,"_",N, ".csv"))
    
    #saving weighted distance
    write.csv(d[[t]],file = paste0("weighted_distance_", t, "_",N,".csv"))
  }
  
  timef<- proc.time()-time0
  CPUtime<-sum(as.vector(timef)[-3])
  write.csv(CPUtime,file=paste0("CPUtime_", N, ".csv"))
  
  return(list(fx=fx,final_posterior=tha_post))
} #end of the weighted-iterative ABC algorithm

