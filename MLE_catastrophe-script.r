
library(compiler)# byte code compilation
library("maxLik")#for maximum likelihood estimation/optimization

#Analytical probability of death due to catastrophe
na.zero<-function(x){
  x[is.na(x)]<-0
  return(x)
}
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

#Analytical probability of death due to catastrophe
#Estimating C(t)=P(catastrophe resulting in 0 population|host death)
Prob_catastrophe<- function(lambda,mu,rho,t,z=1,m=2){
             constant<- 1-PGF_z_compiler(lambda=lambda,mu=mu,rho=rho,t=t,z=z,m=m)
             #return the probability of catastrophic extinction
             return(constant)
    }

#Function for finding Maximum likelihood estimates for the catastrophe rate
#given the GW estimate of birth and death rates for the B-D-C model
MLE_catastrophe<-function(b_est,d_est,dead_fish_time){
          log_like<-0
      #LogLikelihood function to maximize
      Catastrophe_Loglik<-function(param){

         rho<-param[1]
 
      #log likelihood function for catastraphe rate
        for(i in dead_fish_time){ #sum across all dead fish for each group
            if(i>=3){#if the time to death>=3
               log_like<-log_like+na.zero(log(Prob_catastrophe(lambda=b_est,mu=d_est,rho=rho,t=i))-
               Prob_catastrophe(lambda=b_est,mu=d_est,rho=rho,t=(i-2)))
            }else{#if the time to death=1
               log_like<-log_like+na.zero(log(Prob_catastrophe(lambda=b_est,mu=d_est,rho=rho,t=i)-
               Prob_catastrophe(lambda=b_est,mu=d_est,rho=rho,t=0)))   
            }
            
                    }
        log_like
  }
  
    Catastrophe_Loglik_compiler=cmpfun(Catastrophe_Loglik)
  
  ## Inequality constraints:  rho>0
  
   estimates<-maxLik(logLik=Catastrophe_Loglik_compiler,start=c(rho= 1e-5)) 
  
  #return(print(summary(estimates))) #to print estimates & its standard errors and p-values 
  return(as.vector(estimates$estimate))#returning a vector of the estimates: birth, death & catastrophe
}

MLE_catastrophe_compiler<- cmpfun(MLE_catastrophe)
