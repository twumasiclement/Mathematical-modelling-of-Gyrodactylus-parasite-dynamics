
library(compiler)# byte code compilation
library("R.utils")

#Set the catastrophe state -1 to 0
zero.catastrophe <- function (x) {
    x[x<0] <- 0
    return(x)
}


#Set the ratio Z(i)/Z(i-1) to 1 if NA (due to case of 0/0 in Z(i)/Z(i-1))
one.ratio <- function (x) {
    x[is.na(x)|x==Inf|x==-Inf] <- 1
    return(x)
}



#function for sample moments
sample_mean_1st<- function(x) sum(x)/length(x)
sample_mean_2nd<- function(x) sum(x^2)/length(x)
sample_mean_3rd<- function(x) sum(x^3)/length(x)



time<-seq(1,17,by=2)

g_objectivefunc_firstStep <- function(x,prob_sample,fixed=c(FALSE,FALSE,FALSE)) {
        Prob_catastrophe_analytical<- rep(NA,length=length(time))
        params<-fixed
        function(p){
        params[!fixed]<-p
        #The three parameters to be optimized
        b1<-params[1]
        d1<-params[2]
        c1<-params[3]
 

        #Computing theoritical prob of catastrophe
        for(i in seq_along(time)){
            Prob_catastrophe_analytical[i]<- Prob_catastrophe(lambda=b1,mu=d1,rho=c1,t=time[i])                 
                                } 
     
        m1 <- First_moment(b=b1, d=d1,c=c1,t=seq(1,17,by=2),m=2)-apply(zero.catastrophe(x),1,sample_mean_1st)
        m2 <- Second_moment(b=b1, d=d1,c=c1,t=seq(1,17,by=2),m=2)-apply(zero.catastrophe(x),1,sample_mean_2nd)
        m3 <- Third_moment(b=b1, d=d1,c=c1,t=seq(1,17,by=2),m=2)-apply(zero.catastrophe(x),1,sample_mean_3rd)

       Catastrophe_Prob<- Prob_catastrophe_analytical- prob_sample
                                                                       
       gbar_theta<-c(mean(m1),mean(m2),mean(m3),mean(Catastrophe_Prob))
    
       Objective_func<- t(gbar_theta)%*%gbar_theta

                   } 
   
    }

#First step of GMM
GMM_firstStep<-function(prob_sample,x){
      objec_func<- g_objectivefunc_firstStep(x=x,prob_sample=prob_sample)
      initial<-c(2, 1, 0.001)# initial values to optimize over
      estimates<-constrOptim(initial, objec_func, NULL, 
                    ui=rbind(c(1,0,0),  # lambda>0 
                            c(0,1,0),    # mu >0
                            c(0,0,1) # rho > 0
                    ), 
      ci=c(0,0,0),method='Nelder-Mead')$par # the thresholds of the constraints (RHS)$par
   
      return(estimates)
               
  }
    

# Second step of GMM
#Second-step of the GMM optimization
#Function to calculating the weight matrix   
Weight<-function(x,prob_sample,estimate1){ 
 est_step1<- c(estimate1)
 Prob_catastrophe_analytical1<- rep(NA,length=length(time))
 #Computing theoretical prob of catastrophe
 for(i in seq_along(time)){
     Prob_catastrophe_analytical1[i]<- Prob_catastrophe(lambda=est_step1[1],
                    mu=est_step1[2],rho=est_step1[3],t=time[i])                 
                               } 
     
  m1 <- First_moment(b=est_step1[1],d=est_step1[2],c=est_step1[3],t=seq(1,17,by=2),m=2)-apply(zero.catastrophe(x),1,sample_mean_1st)
  m2 <- Second_moment(b=est_step1[1],d=est_step1[2],c=est_step1[3],t=seq(1,17,by=2),m=2)-apply(zero.catastrophe(x),1,sample_mean_2nd)
  m3 <- Third_moment(b=est_step1[1],d=est_step1[2],c=est_step1[3],t=seq(1,17,by=2),m=2)-apply(zero.catastrophe(x),1,sample_mean_3rd)
  Catastrophe_Prob<- Prob_catastrophe_analytical1- prob_sample
                                                                       
  g<-cbind(m1,m2,m3,Catastrophe_Prob)
     
  covariance_matrix<- cov(g)
  #Setting off-diagonals to 0 to obtain an invertible weighting (diagonal) matrix 
  #by assuming that the moment conditions are uncorrelated
  covariance_matrix[lower.tri(covariance_matrix)] <- 0
  covariance_matrix[upper.tri(covariance_matrix)] <- 0
    
  #Finding inverse for the covariance diagonal matrix 
  #NB: if the det(covariance_matrix) is close to 0, inverse(covariance_matrix)
    #using the solve() function in R doesn't converge.
     #weightmatrix<- solve(covariance_matrix)
     weightmatrix<- 1/covariance_matrix #finding reciprocal of entries
     weightmatrix[lower.tri(weightmatrix)] <- 0
     weightmatrix[upper.tri(weightmatrix)] <- 0
     weightmatrix
   
}
        

#Second optimization step                               
g_objectivefunc_2ndStep <- function(x,prob_sample,weighting_matrix,fixed=c(FALSE,FALSE,FALSE)) {
  Prob_catastrophe_analytical<- rep(NA,length=length(time))
  params<-fixed
  function(p){
        params[!fixed]<-p
        #The three parameters to be optimized
        b1<-params[1]
        d1<-params[2]
        c1<-params[3]
 

#Computing theoritical prob of catastrophe
 for(i in seq_along(time)){
     Prob_catastrophe_analytical[i]<- Prob_catastrophe(lambda=b1,mu=d1,rho=c1,t=time[i])                 
                            } 
     
  m1 <- First_moment(b=b1, d=d1,c=c1,t=seq(1,17,by=2),m=2)-apply(zero.catastrophe(x),1,sample_mean_1st)
  m2 <- Second_moment(b=b1, d=d1,c=c1,t=seq(1,17,by=2),m=2)-apply(zero.catastrophe(x),1,sample_mean_2nd)
  m3 <- Third_moment(b=b1, d=d1,c=c1,t=seq(1,17,by=2),m=2)-apply(zero.catastrophe(x),1,sample_mean_3rd)

  Catastrophe_Prob<- Prob_catastrophe_analytical- prob_sample
                                                                       
  gbar_theta<-c(mean(m1),mean(m2),mean(m3),mean(Catastrophe_Prob))
    
  Objective_func<- t(gbar_theta)%*%weighting_matrix%*%gbar_theta

                 } 
   
    } 

#second step of GMM
GMM_2ndStep<-function(prob_sample,x,weighting_matrix){
    objec_func<- g_objectivefunc_2ndStep(x=x,prob_sample=prob_sample,weighting_matrix=weighting_matrix)
    initial<-c(2, 1, 0.001)# initial values to optimize over
    estimates=constrOptim(initial, objec_func, NULL, 
                    ui=rbind(c(1,0,0),  # lambda>0 
                            c(0,1,0),    # mu >0
                            c(0,0,1) # rho > 0
            ),       
    ci=c(0,0,0),method='Nelder-Mead')$par # the thresholds of the constraints (RHS)$par
   
   return(estimates)
               
  }
    
