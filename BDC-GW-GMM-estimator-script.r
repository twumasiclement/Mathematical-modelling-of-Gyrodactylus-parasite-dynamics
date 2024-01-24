
#Loading packages
library(transport) #For Wassertein distance computation
library(parallel) # For parallizing R codes
RNGkind("L'Ecuyer-CMRG") #Dealing with distinct seed value in R
#library(markovchain)
#library(diagram)
#library('latex2exp') #For adding LaTeX symbols to R plots
library(compiler)# byte code compilation
#library(data.table)
library("maxLik")#for maximum likelihood estimation/optimization

#Analytical probability of death due to catastrophe

#Function for constants and PGF
BDCconsts <- function(lambda, mu, rho,t) {
# Constants used in calculating distribution of BDC process at time t
  
  roots <- sort(Re(polyroot(c(mu, -(lambda+mu+rho), lambda))))
  v0 <- roots[1]
  v1 <- roots[2]
  sigma <- exp(-lambda*(v1 - v0)*t)
  k1 <- v0*v1*(1 - sigma)/(v1 - sigma*v0)
  k2 <- (v1*sigma - v0)/(v1 - sigma*v0)
  k3 <- (1 - sigma)/(v1 - sigma*v0)
  return(list(k1=k1, k2=k2, k3=k3, sigma=sigma, v0=v0, v1=v1))
}



#Function for the probability generating function G(z,t)
PGF_z<- function(lambda,mu,rho,t,z,m){
   #v0<-((lambda+mu+rho)-sqrt( ((lambda+mu+rho)^2)-4*mu*lambda))/(2*lambda)
   #v1<-((lambda+mu+rho)+sqrt( ((lambda+mu+rho)^2)-4*mu*lambda))/(2*lambda)
    constants=BDCconsts(lambda,mu,rho,t)
    v0<- constants$v0
    v1<- constants$v1
    sigma<- constants$sigma
    num<-(v0*v1*(1-sigma))+(z*(v1*sigma-v0))
    den<- v1-(sigma*v0)-(z*(1-sigma))
    return( (num/den)^m)
}
PGF_z_compiler<- cmpfun(PGF_z)
#PGF_z_compiler(lambda=0.513,mu=0.35,rho=0.003,t=1,z=1,m=2)

# Function of the Exact mean/1st moment of the BDC process 
First_moment<-function(b,d,c,t,m){
  #b,d,c are the birth,death and catastrophe rates; m=X0=2 and t=time
  roots <- sort(Re(polyroot(c(d, -(b+d+c), b))))
  v0 <- roots[1]
  v1 <- roots[2]
  sigma<-exp(-b*(v1-v0)*t)
  k1<-(v0*v1*(1-sigma))/(v1-(sigma*v0))
  k2<-((v1*sigma)-v0)/(v1-(sigma*v0))
  k3<-(1-sigma)/(v1-(sigma*v0))
  expectation=m*(((k1+k2)/(1-k3))^(m-1))*(k2+(k1*k3))*(1-k3)^-2
  return(expectation)#returns 1st moment
}



#Function of the 2nd moment of the BDC process
Second_moment<-function(b,d,c,t,m){
        roots <- sort(Re(polyroot(c(d, -(b+d+c), b))))
        v0 <- roots[1]
        v1 <- roots[2]
        sigma<-exp(-b*(v1-v0)*t)
        k1<- (v0*v1*(1-sigma))/(v1-(sigma*v0))
        k2<-((v1*sigma)-v0)/(v1-(sigma*v0))
        k3<-(1-sigma)/(v1-(sigma*v0))
       expectation<- m*(((k1+k2)/(1-k3))^(m-1))*(k2+(k1*k3))*(1-k3)^-2
    
       Second_derivative_pgf<-((2*m*k3*(k2+k1*k3))*((k1+k2)/(1-k3))^(m-1)*(1-k3)^-3 + 
                    m*(m-1)*(k2+k1*k3)^2*((k1+k2)/(1-k3))^(m-2)*(1-k3)^-4)
     
       Variance<- (Second_derivative_pgf+ expectation)-(expectation)^2
    
      Second_moment_results<- Variance + expectation^2
      return(Second_moment_results)#returns 2nd moment
}



#Function of the 3rd moment of the BDC process

Third_moment<-function(b,d,c,t,m){
       roots <- sort(Re(polyroot(c(d, -(b+d+c), b))))
       v0 <- roots[1]
       v1 <- roots[2]
       sigma<-exp(-b*(v1-v0)*t)
       k1<-(v0*v1*(1-sigma))/(v1-(sigma*v0))
       k2<-((v1*sigma)-v0)/(v1-(sigma*v0))
       k3<-(1-sigma)/(v1-(sigma*v0))
       expectation<- m*(((k1+k2)/(1-k3))^(m-1))*(k2+(k1*k3))*(1-k3)^-2
    
       Second_derivative_pgf<- ((2*m*k3*(k2+k1*k3))*((k1+k2)/(1-k3))^(m-1)*(1-k3)^-3 + 
                    m*(m-1)*(k2+k1*k3)^2*((k1+k2)/(1-k3))^(m-2)*(1-k3)^-4)
    
       Third_derivative_pgf<- 6*m*(k2+k1*k3)*(k3^2)*(((k1+k2)/(1-k3))^(m-1))*(1-k3)^(-4)+
                   6*m*(m-1)*((k2+k1*k3)^2)*k3*(((k1+k2)/(1-k3))^(m-2))*(1-k3)^(-5)+
                   m*(m-1)*(m-2)*((k2+k1*k3)^3)*(((k1+k2)/(1-k3))^(m-3))*(1-k3)^(-6)   
     
       Variance<- (Second_derivative_pgf+ expectation)-(expectation)^2
    
       Second_moment_results<- Variance + expectation^2
    
       Third_moment_results<- Third_derivative_pgf+(3*Second_moment_results)-(2*expectation)
       return(Third_moment_results)#returns 3rd moment
}


#Analytical probability of death due to catastrophe
#Estimating C(t)=P(catastrophe resulting in 0 population|host death)
Prob_catastrophe<- function(lambda,mu,rho,t,z=1,m=2){
             constant<- 1-PGF_z_compiler(lambda=lambda,mu=mu,rho=rho,t=t,z=z,m=m)
             #return the probability of catastrophic extinction
             return(constant)
    }


#Setting working directory
#setwd("/home/clement/Documents/Simulation_Data_Folder")

#External scripts

source("MLE_catastrophe-script.r")
source("GMM-1st2nd-Steps-script.r")

RestructureData_BDC<- function(pop,alive,group){
        #Inputs:pop=parasite population per region over time
                #alive= survival status over time
                # group=parasite-fish groups
    
       # to store parasite numbers over time as a dataframe for each parasite-fish
        ParasiteData_combined<- NULL 
       # to store  survival status as a dataframe for each parasite-fish
        SurvStatus_combined<- NULL 
        
       #Set NA in pop to state 0 denoting host death for the B-D-C estimation
       na.zero <- function (x) {
            x[is.na(x)] <- 0
            return(x)}

    
        for(pf in seq_along(group)){
            
            ParasiteData_combined[[pf]]<- matrix(NA,nrow=9, ncol=numF[[pf]])
            #Array for time steps fish was alive for each combination
            SurvStatus_combined[[pf]]<-  matrix(NA,nrow=9, ncol=numF[[pf]])
            for(i in 1:numF[[pf]]){
                #total parasites over time for each fish belonging to each parasite-fish group
                # state -1 in the BDC denote host death
                ParasiteData_combined[[pf]][,i]<- na.zero(apply(pop[[pf]][i,,],2,sum))
                SurvStatus_combined[[pf]][,i]<- alive[[pf]][i, ]
                }
        }
       return(list(PopTime_group=ParasiteData_combined,SurvTime_group=SurvStatus_combined))
}

GW_GMM_BDCestimator<-function(X0,pop,alive,group){
  #X0= initial parasites
 Parasite_data<- NULL; survival_data<- NULL
 #re-structuring the format of the data into the 9 parasite-fish groups
 data<- RestructureData_BDC(pop=pop,alive=alive,group=group)
                 
 time<-seq(1,17,by=2)
 # Parasite_data[[pf]][,fish_index]                             
 Prob_catastrophe_analytical=Prob_catastrophe_sample= matrix(0,nrow=length(time),ncol=length(group))
################################  Initialize GMM   ##############################
 #Computing catastrophic probability analytically and based on the sample data  
     
 #computing sample prob of catastrophe                   
                   
 time_index<- seq_along(time)
 for (pf in seq_along(group)){ 
    Parasite_data[[pf]]<- data$PopTime_group[[pf]]
    survival_data[[pf]]<- data$SurvTime_group[[pf]] 
    for(i in time_index){
       if(any(survival_data[[pf]][i,]==2)==TRUE){
         #print(paste("time=",time[i]))
          fish_dead_sim<-length(which(survival_data[[pf]][i, ]==2))
          #print( fish_dead_sim)
          Prob_catastrophe_sample[i,pf]<-fish_dead_sim/dim(survival_data[[pf]])[2]
            }
         }

     }
    

   #Let Zi_t be the population for fish i at time t
   # Let alive_status be the survival status of each fish
    Z=NULL; alive_status=NULL
    for(pf in seq_along(group)) {
        Z[[pf]]<-list()
        alive_status[[pf]]<-list()
                            }

    for(pf in seq_along(group)){
        for(k in 1:numF[[pf]]){     
           Z[[pf]][[k]]<- Parasite_data[[pf]][,k]       
           alive_status[[pf]][[k]]<- survival_data[[pf]][,k]  
                       }
               }
    
    
    #Computing the mean and variance for the Galton-Watson process based on fish survival
       # And for each k replicate
     mean_GW=NULL; var_GW=NULL;  mean_sum_num=NULL; mean_sum_den=NULL;var_sum=NULL 
                   
            #Computing the mean of GW process
     for(pf in seq_along(group)) mean_sum_num[[pf]]=mean_sum_den[[pf]]=0 #initial summation for the GW mean
         for(pf in seq_along(group)){
           for(k in 1:numF[[pf]]){
               if(all(survival_data[[pf]][,k]==1)==TRUE){
                   
                 mean_sum_num[[pf]]<-mean_sum_num[[pf]]+sum(Z[[pf]][[k]][1:9])# sum from t1 to t17
                 mean_sum_den[[pf]]<-mean_sum_den[[pf]]+sum(Z[[pf]][[k]][1:8])+X0 #sum from t0 to t15
                    
                                                            }
                                                  }
               mean_GW[[pf]]<-   one.ratio(mean_sum_num[[pf]]/mean_sum_den[[pf]])#if 0/0=1  
                                     }
    
       #computing the variance of GW process 
     for(pf in seq_along(group)) var_sum[[pf]]<-0 #initial summation for GW variance 
         for(pf in seq_along(group)){
            for(k in 1:numF[[pf]]){
               if(all(survival_data[[pf]][,k]==1)==TRUE){
                  
                  var_sum[[pf]]<-var_sum[[pf]]+ sum(Z[[pf]][[k]][1:9]*
                                (one.ratio(Z[[pf]][[k]][1:9]/c(X0,Z[[pf]][[k]][1:8])) -mean_GW[[pf]])^2)         
                                                    }
                                           }
                 var_GW[[pf]]<- var_sum[[pf]]/(numF[[pf]]*length(time))
                 #print(paste("Var_GW", var_GW[[k]]))
                 }
    
    #########    GMM estimation ####################
    birth_rate=NULL;death_rate=NULL; c_estimates<-NULL;delta_t=2; BDC_estimates=NULL
    GMM_resultsStep1=NULL; GMM_resultsStep2=NULL; weighting_matrix_cov=NULL; method=NULL 
            
       #Estimating the catastrophe rate using MLE when m>1 for GW estimation
    #time at death for each fish i and replicate/simulation run k
         t_death<-NULL; for(pf in seq_along(group))  t_death[[pf]]<-rep(NA,length=numF[[pf]])                
     for(pf in seq_along(group)){
           for(k in 1:numF[[pf]]){
               #time to death 
              t_death[[pf]][k]<-time[which(survival_data[[pf]][,k]==2)[1]]          
                           }
                    }
    
   
    
    for(pf in seq_along(group)){### begining of GW and GMM
       if(mean_GW[[pf]]>1){####  Consider GW if mean_GW>1
            method[[pf]]<-"GW estimation"
          birth_rate[[pf]]<-((log(mean_GW[[pf]])/(2*delta_t))*
                 (one.ratio(var_GW[[pf]]/(mean_GW[[pf]]*(mean_GW[[pf]]-1))) +1)) 
          death_rate[[pf]]<- ((log(mean_GW[[pf]])/(2*delta_t))*
                 (one.ratio(var_GW[[pf]]/(mean_GW[[pf]]*(mean_GW[[pf]]-1))) -1))
      
         #Computing MLE of catastrophe rate based on estimated birth and death rates              
          if(all(is.na(t_death[[pf]]))==FALSE){#if at least some fish are dead
         #estimates of the catastrophe rate
             c_estimates[[pf]]<-MLE_catastrophe_compiler(b_est=birth_rate[[pf]],
                           d_est<- death_rate[[pf]],dead_fish_time=na.omit(t_death[[pf]][k]))#as lists
                 }else if(all(is.na(t_death[[pf]]))==TRUE){ #if no fish is dead
             c_estimates[[pf]]<-0 
                     }
            
         BDC_estimates[[pf]]<-c(birth_rate[[pf]],death_rate[[pf]],c_estimates[[pf]])
    }else if(mean_GW[[pf]]<=1){ #Consider GMM
                 #print(paste("c_rate is empty=",all(is.na(t_death[[k]]))))
            method[[pf]]<-"GMM estimation"                  
           #First stage of GMM
           GMM_resultsStep1[[pf]]<- GMM_firstStep(prob_sample=Prob_catastrophe_sample[,pf],
                                                  x=as.data.frame(Parasite_data[[pf]]))
 
           weighting_matrix_cov[[pf]]<-Weight(x=as.data.frame(Parasite_data[[pf]]),
                                   prob_sample=Prob_catastrophe_sample[,pf],estimate1=GMM_resultsStep1[[pf]])
    
          #Second stage of GMM     
          GMM_resultsStep2[[pf]]<- GMM_2ndStep(prob_sample=Prob_catastrophe_sample[,pf],
                                               x=as.data.frame(Parasite_data[[pf]]),
                                   weighting_matrix=weighting_matrix_cov[[pf]])
           
          BDC_estimates[[pf]]<-  GMM_resultsStep2[[pf]]  
                                                       } ####GMM estimation ends
                            } #### end of GW and GMM
    
    
      
 BDC_estimates_df<-do.call("rbind", BDC_estimates) 
 return(list(BDC_estimates=BDC_estimates_df,method_used=unique(unlist(method))))       
                                  
        
}
