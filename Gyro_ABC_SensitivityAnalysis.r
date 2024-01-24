#Loading some relevant R packages
library(transport) #For Wasserstein distance computation
library(parallel) # For parallizing R codes
RNGkind("L'Ecuyer-CMRG") #Dealing with distinct seed value in R
#library(markovchain)
#library(diagram)
#library('latex2exp') #For adding LaTeX symbols to R plots
library(compiler)# byte code compilation
#library(data.table)
library("maxLik")#for maximum likelihood estimation/optimization
library("R.utils")
library("KernSmooth")
library("ks")
#install.packages("/Users/clementtwumasi/Desktop/kedd_1.0.3.tar.gz", repos = NULL, type="source")
library(kedd)
library("mvtnorm")
library(Rcpp)# for C+ matrix multiplication
library(tmvtnorm)
library("philentropy")
library("ie2misc")#for computing mean absolute deviation

#Other packages for post-processing analysis
library(MASS)
library("kader")
library("PerformanceAnalytics")
library(glmnet)
library(caret)
#library("kedd")
library(bayestestR)#for credible interval and HDI+ROPE decision
library("sjstats")#for HDI+ROPE decision
library(abind)
#library(sjmisc)
#library(rstanarm)
#library("HDInterval")
#library("LaplacesDemon")
library(dplyr)
library(ggplot2)
library(tidyverse)
library(lme4)
#library(merTools)
library(glmmTMB) # fitting generalised linear mixed models
library(bbmle) # general maximum likelihood estimation
#library(ggthemes)

#library(ggthemes)
library(showtext)
library("ggpubr")
theme_set(theme_bw()+
theme(panel.spacing=grid::unit(0,"lines")))
#library(sjPlot) #for plotting lmer and glmer mods
library(sjmisc) 
library(effects)
#library(sjstats) #use for r2 functions
library(jtools)
library(forcats)
library(rstatix)
library("sjPlot")
library("gplots")

options(repr.plot.width=8, repr.plot.height=8,repr.plot.res = 300) #Setting plot size

#Setting working directory
setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro")


#### External scripts for the simulation model and ABC #####

#Script of function for computing event rates
source("Computing-rates-script.r")

# Script of function for updating exact SSA
source("Update-exactSSA-script.r")

# Script of function for updating tau-leaping
source("Update-tauleaping-script.r")

# Script of function experimental descriptors 
# (fish type, strain, fish size, fish sex and areas of the 4 body regions)
source("Descriptors-Data-script.r")

# Script of function for simulating parasites only a single fish over time and across body regions
source("Simulation-single-fish-script.r")

# Script of function for simulating parasites for a group fish over time and across body regions
#Corresponding to the empirical data
source("Simulation-Group-fish-script.r")

#Script of functions for Galton-Watson & GMM estimation of the B-D-C parameters
source("MLE_catastrophe-script.r")
source("GMM-1st2nd-Steps-script.r")
source("BDC-GW-GMM-estimator-script.r")

#Script for population projection until day 17 after fish death to aid with ABC fitting
source("Project-Parasite-script.r")


#Script for computing Weighted Euclidean Distance function between observed and simulated summaries
source("Weighted-distance-script.R")

#Script for function used to compute the 17 summary statistics for ABC fitting:
#1.Log count across time (9)
#2.Wasserstein distance between the four body regions (4),
#3.Estimates of the B-D-C model parameters (3)
#4.Time before death (1)
source("summary-stats-ABC-script.r")

#Script of function for combining the 17 summary statistics
source("combine-summary-stats-script.R")

#Script wich contains functions of priors, ABC importance sampling
# and Function for ABC fitting
source("ABC-Importance-Sampling-Improved-script.R")

#Script for the Weighted-iterative ABC (with ABC-SMC and adaptive importance sampling)
source("Weighted-iterative-ABC-script.R")


#To detect the number of processors or cores of your computer 
numCores<-parallel::detectCores()
numCores


#Importing empirical data 
Combined_data <- read.csv(file="Parasite_Data.csv") 
#Importing data for area of the 8 body parts across 18 fish (measured in millimeters square)
Bodyparts_area<-read.csv(file="Area_Fish_bodyParts.csv")

#Experimental descriptors
Descriptors<- Experiment_descriptors(empirical_data=Combined_data)

fishSize <- Descriptors$fishSize #fish size
fishSex  <- Descriptors$fishSex #fish sex
Strain   <- Descriptors$Strain # parasite strain
Fish_stock<-Descriptors$Fish_stock #fish stock
numF    <-  Descriptors$numF # total fish for each parasite-fish group 
fishID  <-  Descriptors$fishID # fish IDs for each parasite-fish group 

Area_normalized<- Body_area(Area_data=Bodyparts_area)#body areas

#Parasite-fish groups
parasite_fish<-c("Gt3-OS","Gt3-LA","Gt3-UA","Gt-OS","Gt-LA","Gt-UA","Gb-OS","Gb-LA","Gb-UA") 


parameter_labels=c(expression(paste("b"[11])),expression(paste("b"[12])),
                   expression(paste("b"[21])),expression(paste("b"[22])), 
                   expression(paste("b"[31])),expression(paste("b"[32])),
                   expression(paste("d"[11])),expression(paste("d"[12])),
                   expression(paste("d"[21])),expression(paste("d"[22])), 
                   expression(paste("d"[31])),expression(paste("d"[32])),
                   expression(paste("m")),  expression(paste("r")),
                   expression(paste("r"[1])),expression(paste("r"[2])),
                   expression(paste("r"[3])),expression(paste("s")),
                   expression(paste("s"[1])),expression(paste(epsilon[1])),
                   expression(paste(epsilon[2])),expression(paste(epsilon[3])),
                   expression(paste(kappa))
)

fixed_parameter_values<- c(0.1796, 0.0572, 0.6463,  0.0351, 2.0180, 0.0237,  0.0054, 0.0967, 0.0327, 0.0710, 0.0538, 3.8354,
  0.0584, 0.4836, 0.0005, 0.0002, 0.0016, 0.0133, 0.0037, 0.0059, 1.4134, 1.0580, 118.0921)

fixed_parameter_values_log<- log(fixed_parameter_values)

cat("Fixed parameter values on Log scale:","\n")

data.frame(Parameters=paste(parameter_labels), 
           true_values_log=fixed_parameter_values_log,
          true_values_original=fixed_parameter_values)


##########. Simulating a pseudo-observed data at the posterior estimates ########
set.seed(12)#for reproducibility of results
pseudo_observed_data<-SimGroup_tauleap(theta1=fixed_parameter_values_log,fish_sex=fishSex,
                                       fish_type=Fish_stock,strain=Strain,fish_size=fishSize,error=0.002)

#view the list of outputs
names(pseudo_observed_data)

#view a pseudo-observed data
group_number<-8
pseudo_observed_data$pop_sim[[group_number]][9,,]#view simulated parasite data for 9th fish for group 8
pseudo_observed_data$alive_sim[[group_number]][9, ]
pseudo_observed_data$group[group_number]

#Storing simulated data in the required array structure
pop_pseudo_obs<-NULL;  alive_pseudo_obs<- NULL; 
#The Experimental descriptors of the observed data are the same as that of the pseudo-observed data
for (pf in 1:length(parasite_fish)){
  #Pseudo-observed data or matrix across 4 regions
  pop_pseudo_obs[[pf]] <- array(dim = c(numF[[pf]], 4, 9))
  #Array for time steps simulated fish was alive for each parasite-host combination
  alive_pseudo_obs[[pf]] <- array(dim = c(numF[[pf]], 9))     
  
  for (i in 1:numF[[pf]]) {
    pop_pseudo_obs[[pf]][i,,] <- pseudo_observed_data$pop_sim[[pf]][i,,]
    alive_pseudo_obs[[pf]][i, ] <- pseudo_observed_data$alive_sim[[pf]][i,  ]
  }
  
}

#B-D-C parameter estimates for the parasite-fish groups based on pseudo-observed data
BDC_estimates_obs<- GW_GMM_BDCestimator(X0=2,pop= pop_pseudo_obs,alive=alive_pseudo_obs,
                                        group=parasite_fish)$BDC_estimates
#BDC_estimates_obs



#Computing summary statistics for observed data across the parasite-fish groups
summaries_pseudo_obs<- Summary_stats(pop=pop_obs,alive=alive_obs,BDC_estimates=BDC_estimates_obs)
summaries_obs<-summaries_pseudo_obs

#View ABC summary statistics for the observed data
#summaries_obs[[1]]# for group 1 or Gt3-0S group

#Pseudo-observed summary statistics
S_obs<- do.call("rbind",summaries_obs)
dim(S_obs)

#Initial simulation inputs for A (parasite numbers) and  B (immune status)
A0 <- matrix(0, 4, 2)  
A0[1, 1] <- 2   #Intial parasites at the tail
B0 <- rep(1, 4)  #initial immune response at 4 body regions (no response)

#Transition matrix
J<- matrix(c(0,    1,    0,     0, 
             1/2,   0,    1/2,   0,
             0,    1/2,    0,   1/2,
             0,     0,     1,    0), 4, 4, byrow=TRUE)

#initial weight (inverse of summary statistics variances)

###### initial weights estimate
setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro")

#Importing saved results of the initial summary statistics weights
w=read.csv(file="w0.csv")#based on inverse variance
w<-as.vector(w[,-1])
print(w)


Total_fish<- dim(do.call("rbind",summaries_obs))[1]

#Setting working directory to save ABC results

sample_sizes<- 1500 
ABC_output<- NULL
options(warn=-1)

setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro/Posterior_Paramater identifiabity")

print(paste("At N =",sample_sizes))
ABC_output<- Weighted_iterative_ABC(N=sample_sizes,dimS=17,fish_total=Total_fish,
                                    numCores=numCores,ABC_time_steps=10,seed_num=1)

setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro/Posterior_Paramater identifiabity")

#Importing saved outputs from ABC-SMC fitting
ABC_time_steps<-10
N<- 1500
#Importing saved outputs   
distances_1500<- list();theta_1500<- list();theta_posterior_1500<- list()
  for(t in 1:ABC_time_steps){  
     distances_1500[[t]]<- read.csv(paste0("weighted_distance_",t,"_",N,".csv"))
     distances_1500[[t]]<-distances_1500[[t]][,-1]
     theta_1500[[t]]<- read.csv(paste0("theta_t_",t,".csv")) 
     theta_1500[[t]]<- theta_1500[[t]][,-1]  
     theta_posterior_1500[[t]]<- read.csv(paste0("theta_post_",t,"_",N,".csv"))
     theta_posterior_1500[[t]]<-theta_posterior_1500[[t]][,-1] 
    }


#Final posterior  
Final_posterior_1500<- read.csv(paste0("PosteriorFinal_",N, ".csv"))
Final_posterior_1500<-Final_posterior_1500[,-1]

# Importing prior, intermediate priors and posterior densities
#Importing posterior densities from Weighted-iterative ABC with Sequential Monte Carlo with importance sampling

density_post_1500<-NULL 
ABC_iterations<-ABC_time_steps+1

no.params<- length(parameter_labels)

for(j in 1:ABC_iterations) {
    density_post_1500[[j]]<-read.csv(file=paste0("density_post_",j,"_",N,".csv"))
    density_post_1500[[j]]<-density_post_1500[[j]][,-1]
          }

#Import the total CPU time to fit the complex model using the ABC-SMC sampler
N<- 1500 # Monte Carlo Sample Size
time_CPU_1500<- read.csv(file=paste0("CPUtime_",N,".csv"))
time_CPU_1500<- time_CPU_1500[,-1]
    

#Computational times
 print(paste("CPU time at N=", N,"",":",time_CPU_1500/10,"secs"))
 print(paste("CPU time at N=", N,"",":",time_CPU_1500/864000,"days"))


#Converting from seconds to days
library(lubridate)
print(paste("CPU time at N=", N,"",":",seconds_to_period(time_CPU_1500/10)))

dist_0.05quantile_1500<-rep(NA,ABC_time_steps)
dist_0.15quantile_1500<-rep(NA,ABC_time_steps)
dist_0.25quantile_1500<-rep(NA,ABC_time_steps)
dist_0.50quantile_1500<-rep(NA,ABC_time_steps)
dist_0.75quantile_1500<-rep(NA,ABC_time_steps)
dist_0.85quantile_1500<-rep(NA,ABC_time_steps)

for(t in 1:ABC_time_steps){
   dist_0.05quantile_1500[t]<- quantile(distances_1500[[t]],probs =0.05)
   dist_0.15quantile_1500[t]<- quantile(distances_1500[[t]],probs =0.15)
   dist_0.25quantile_1500[t]<- quantile(distances_1500[[t]],probs =0.25)
   dist_0.50quantile_1500[t]<- quantile(distances_1500[[t]],probs =0.5)
   dist_0.75quantile_1500[t]<- quantile(distances_1500[[t]],probs =0.75)
   dist_0.85quantile_1500[t]<- quantile(distances_1500[[t]],probs =0.85) 
}


plot(1:ABC_time_steps,dist_0.05quantile_1500,type="b",lwd=3,col="blue",
     ylim=c(0.55,.9),xlab="ABC-SMC time steps",pch=3
     ,main="", ylab="Quantiles of ABC-SMC distances")
axis(side=1, at = 1:ABC_time_steps, labels = paste(1:ABC_time_steps))
lines(1:ABC_time_steps,dist_0.15quantile_1500,type="b",lwd=3,col="green",pch=4)
lines(1:ABC_time_steps,dist_0.25quantile_1500,type="b",lwd=3,col="red",pch=13)
lines(1:ABC_time_steps,dist_0.50quantile_1500,type="b",lwd=3,col="black",pch=6)
lines(1:ABC_time_steps,dist_0.75quantile_1500,type="b",lwd=3,col="pink",pch=5)
lines(1:ABC_time_steps,dist_0.85quantile_1500,type="b",lwd=3,col="orange",pch=10)


legend(x=2,y=.9,legend=c("5th percentile","15th percentile",
                         "25th percentile","50th percentile","75th percentile","85th percentile"),
        col=c("blue","green","red","black","pink","orange"),pch=c(3,4,13,6,5,10),
       cex=1,box.lwd = 2,ncol =3,title="Percentiles:",bty="n")

#All model parameters
N<- 1500 # Monte Carlo Sample Size
x <- seq(from = -10, to = 7, length.out = 256)#range of prior distribution

par(mfrow=c(5,5), mar=c(4,2,1,0),font=2)
plot(NULL ,xaxt='n',yaxt='n',bty='n', xlab="", ylab="",xlim=0:1, ylim=0:1)
plot_colors <- c("blue","red","black")
legend(x=-0.1,y=1.2,c(paste0("First Prior 
(N=",N,")"),"Intermediate 
Priors" ,"Final Posterior"),
        col=c("blue","red","black"),bty="n",cex=1.05,box.lwd = .8,fill=c("blue","red","black"),horiz=FALSE)


for (k in 1:23) {
  plot(x, density_post_1500[[1]][k, ], type="l", ylim=c(0,5),
       xlab=parameter_labels[k],ylab="", col="blue",cex.lab=1.2,lwd=2.5,las=1) 
  
    
  for (j in 2:ABC_time_steps) {
    lines(x,density_post_1500[[j]][k, ], yaxt= "n",col="red",lwd=1,pch=4,ann=FALSE,yaxt="n")
  }
    
     lines(x, density_post_1500[[ABC_iterations]][k, ], col="black", lwd=1,ann=FALSE,yaxt="n")
   
}

#Estimating the model parameters from the unadjusted posterior samples and their corresponding credible intervals

  #computing posterior estimaates and their credible intervals (ETI)
  post_data_1500<-   Final_posterior_1500
  ci_eti_1500<- bayestestR::ci(as.data.frame(post_data_1500), ci = 0.95,method = "ETI") 
  ci_eti_1500<- exp(as.data.frame(ci_eti_1500[,3:4]))
   # estimating the  posterior  mean estimate
  Posterior_estimates_unadj_1500<- apply(exp(Final_posterior_1500),2,mean)
    

Posterior_est_unadj_data_1500<- data.frame(Posterior_estimate_unadj=Posterior_estimates_unadj_1500,
                                     Cred_Int_lower95_unadj= ci_eti_1500[,1],
                                     Cred_Int_upper95_unadj= ci_eti_1500[,2])


rownames(Posterior_est_unadj_data_1500)<- parameter_labels
#Print unadjusted posterior estimates with 95% credible intervals
Posterior_est_unadj_data_1500

#Save the unadjusted posterior estimate results

write.csv(Posterior_est_unadj_data_1500, "Posterior_estimates_unadj_data_1500.csv")

#Setting working directory
setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro")


#Function for posterior adjustment using Weighted Ridge and Lasso regression, respectively
source("Post-Ridge-reg-adj-L2-script.R")#ABC Post processing-Ridge/l2 correction
source("Post-Lasso-reg-adj-L1-script.R")#ABC Post processing-Lasso/L1 correction

set.seed(1) #for reproducibility
#Ridge or L2 posterior correction
Reg_adjust_L2<- function(i) Post_Ridge_reg_adj(post_distn=Final_posterior_1500,
                                    summary_obs= S_obs)

Output_adj_Ridge<- mclapply(1, Reg_adjust_L2,mc.cores=numCores)

#Based on Ridge regression
X_matrix<- as.data.frame(Output_adj_Ridge[[1]]$X_Design_matrix)

for (i in 1:dim(X_matrix)[2]) names(X_matrix)[i]<- paste0("S_",i)

chart.Correlation(X_matrix, histogram=TRUE, pch=19)

round(Output_adj_Ridge[[1]]$Posterior_mean_output,4)

set.seed(1234) #for reproducibility
#Lasso or L1 posterior correction
Reg_adjust_L1<- function(i) Post_Lasso_reg_adj(post_distn=Final_posterior_1500,
                                    summary_obs= S_obs)

Output_adj_Lasso<- mclapply(1,Reg_adjust_L1,mc.cores=numCores)

#Based on Lasso regression
X_matrix<- as.data.frame(Output_adj_Lasso[[1]]$X_Design_matrix)

for (i in 1:dim(X_matrix)[2]) names(X_matrix)[i]<- paste0("S_",i)

chart.Correlation(X_matrix, histogram=TRUE, pch=19)

round(Output_adj_Lasso[[1]]$Posterior_mean_output,4)

#Setting working directory
setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro/Posterior_Paramater identifiabity")




#Importing final posterior at different sample sizes of initial prior distn
Ridge_Adjusted_posterior<-NULL;Lasso_Adjusted_posterior<-NULL
draws<- c(1500)

for (n in draws){
    #On log scale
  Ridge_Adjusted_posterior[[n]]<- Output_adj_Ridge[[1]]$Adjusted_posterior_dist
  Lasso_Adjusted_posterior[[n]]<- Output_adj_Lasso[[1]]$Adjusted_posterior_dist
}

n<-draws
#Estimating the 95% Credible interval for the unadjusted posterior
ci_eti_Ridge_adj<- NULL;Posterior_estimates_Ridge_adj<- NULL
ci_eti_Lasso_adj<- NULL;Posterior_estimates_Lasso_adj<- NULL

 ####Computing credible intervals (ETI)---For Ridge adjustment##########
  Ridge_post_data<-  Ridge_Adjusted_posterior[[n]]
  #Estimated Credible Intervals On log scale
  ci_eti_Ridge_adj[[n]]<- as.data.frame(ci(as.data.frame(Ridge_post_data), ci = 0.95,method = "ETI"))[,3:4]
 #Estimated Credible Intervals On Original scale
  ci_eti_Ridge_adj[[n]]<- exp(ci_eti_Ridge_adj[[n]])
   # estimating the  posterior  mean estimate
  Posterior_estimates_Ridge_adj[[n]]<- Output_adj_Ridge[[1]]$Posterior_mean_output[,1]
    

Posterior_est_Ridge_adj_data<- data.frame(Ridge_Posterior_estimate_adj=Posterior_estimates_Ridge_adj[[1500]],
                                     Ridge_Cred_Int_lower95_adj= ci_eti_Ridge_adj[[1500]][,1],
                                     Ridge_Cred_Int_upper95_adj= ci_eti_Ridge_adj[[1500]][,2])


#Print Ridge adjusted and unadjusted posterior estimates with 95% credible intervals
Combined_posterior_Ridge_adj_unadj<- cbind(Posterior_est_Ridge_adj_data,Posterior_est_unadj_data_1500)


#Save output
write.csv(Combined_posterior_Ridge_adj_unadj,"Post_Ridge_Regression_output.csv")





 ####Computing credible intervals (ETI)---For Lasso adjustment##########
  Lasso_post_data<-  Lasso_Adjusted_posterior[[n]]
  #Estimated Credible Intervals On log scale
  ci_eti_Lasso_adj[[n]]<- as.data.frame(ci(as.data.frame(Lasso_post_data), ci = 0.95,method = "ETI"))[,3:4]
 #Estimated Credible Intervals On Original scale
  ci_eti_Lasso_adj[[n]]<- exp(ci_eti_Lasso_adj[[n]])
   # estimating the  posterior  mean estimate
  Posterior_estimates_Lasso_adj[[n]]<- Output_adj_Lasso[[1]]$Posterior_mean_output[,1]
    

Posterior_est_Lasso_adj_data<- data.frame(Lasso_Posterior_estimate_adj=Posterior_estimates_Lasso_adj[[1500]],
                                     Lasso_Cred_Int_lower95_adj= ci_eti_Lasso_adj[[1500]][,1],
                                     Lasso_Cred_Int_upper95_adj= ci_eti_Lasso_adj[[1500]][,2])


#Print Lasso adjusted and unadjusted posterior estimates with 95% credible intervals
Combined_posterior_Lasso_adj_unadj<- cbind(Posterior_est_Lasso_adj_data,Posterior_est_unadj_data_1500)


#Save output
write.csv(Combined_posterior_Lasso_adj_unadj,"Post_Lasso_Regression_output.csv")


#Saving Ridge and Lasso adjusted posterior distributions 
write.csv(Output_adj_Ridge[[1]]$Adjusted_posterior_dist,
                        paste0("Posterior_Ridge_adj_distn_log_",draws,".csv"))


write.csv(Output_adj_Lasso[[1]]$Adjusted_posterior_dist,
                        paste0("Posterior_Lasso_adj_distn_log_",draws,".csv"))

#Importing previously saved ABC Ridge post-processing results
#Setting working directory
setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro/Posterior_Paramater identifiabity")


Post_Ridge_Regression_output<- read.csv("Post_Ridge_Regression_output.csv")
Post_Ridge_Regression_output<- Post_Ridge_Regression_output[,-1]



Post_Lasso_Regression_output<- read.csv("Post_Lasso_Regression_output.csv")
Post_Lasso_Regression_output<- Post_Lasso_Regression_output[,-1]

#Setting working directory
setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro/Posterior_Paramater identifiabity")


#Importing adjusted and unadjusted posterior distributions
Posterior_Ridge_adj_distn_log<- read.csv("Posterior_Ridge_adj_distn_log_1500.csv")
Posterior_Ridge_adj_distn_log<- Posterior_Ridge_adj_distn_log[,-1]



Posterior_Lasso_adj_distn_log<- read.csv("Posterior_Lasso_adj_distn_log_1500.csv")
Posterior_Lasso_adj_distn_log<- Posterior_Lasso_adj_distn_log[,-1]



Posterior_Unadj__distn_log<- read.csv("PosteriorFinal_1500.csv")
Posterior_Unadj__distn_log<- Posterior_Unadj__distn_log[,-1]

#Estimating the kernel density of the adjusted posterior
number_of_parameters<- 23
density_post_Ridge_adjusted<- array(dim=c(number_of_parameters, 256))
density_post_Lasso_adjusted<- array(dim=c(number_of_parameters, 256))

for(k in 1:number_of_parameters){
  density_post_Ridge_adjusted[k, ]<- density(Posterior_Ridge_adj_distn_log[ ,k], from=-10, to=7, n=256)$y
  density_post_Lasso_adjusted[k, ]<- density(Posterior_Lasso_adj_distn_log[ ,k], from=-10, to=7, n=256)$y
}

par(mfrow=c(3,2),mar=c(4,4,1,1))

x <- seq(from = -10, to = 7, length.out = 256)#range of prior distribution

for (k in 1:6) {
  plot(x, density_post_1500[[1]][k, ], type="l", ylim=c(0,5),
       xlab=parameter_labels[k],ylab="Density", col="blue",cex.lab=1.2,lwd=2.5,las=1) 

 if(k==1) legend("topleft",c("First Prior","Intermediate Priors" ,"Unadjusted Posterior","Ridge Adjustment","Lasso Adjustment"),
         col=c("blue","red","black","green","orange"),bty="n",cex=1,box.lwd = 2,fill=c("blue","red","black","green","orange"))

    
  for (j in 2:10) {
    lines(x,density_post_1500[[j]][k, ], yaxt= "n",col="red",lwd=1,pch=4,ann=FALSE,yaxt="n")
  }
    
    
     lines(x, density_post_1500[[ABC_iterations]][k, ], col="black", lwd=1.5)
     lines(x,density_post_Ridge_adjusted[k, ], col="green",lty=3,lwd=2)
      lines(x,density_post_Lasso_adjusted[k, ], col="orange",lty=3,lwd=2)
   
    }

par(mfrow=c(3,2),mar=c(4,4,1,1))

x <- seq(from = -10, to = 7, length.out = 256)#range of prior distribution

for (k in 7:12) {
  plot(x, density_post_1500[[1]][k, ], type="l", ylim=c(0,5),
       xlab=parameter_labels[k],ylab="Density", col="blue",cex.lab=1.2,lwd=2.5,las=1) 

 if(k==7) legend("topright",c("First Prior","Intermediate Priors" ,"Unadjusted Posterior","Ridge Adjustment","Lasso Adjustment"),
         col=c("blue","red","black","green","orange"),bty="n",cex=1,box.lwd = 2,fill=c("blue","red","black","green","orange"))

    
  for (j in 2:10) {
    lines(x,density_post_1500[[j]][k, ], yaxt= "n",col="red",lwd=1,pch=4,ann=FALSE,yaxt="n")
  }
    
     lines(x, density_post_1500[[ABC_iterations]][k, ], col="black", lwd=1.5)
     lines(x,density_post_Ridge_adjusted[k, ], col="green",lty=3,lwd=2)
      lines(x,density_post_Lasso_adjusted[k, ], col="orange",lty=3,lwd=2)
   
    }

par(mfrow=c(3,2),mar=c(4,4,1,1))

x <- seq(from = -10, to = 7, length.out = 256)#range of prior distribution

for (k in 13:18) {
  plot(x, density_post_1500[[1]][k, ], type="l", ylim=c(0,5),
       xlab=parameter_labels[k],ylab="Density", col="blue",cex.lab=1.2,lwd=2.5,las=1) 

 if(k==14) legend("topleft",c("First Prior","Intermediate Priors" ,"Unadjusted Posterior","Ridge Adjustment","Lasso Adjustment"),
         col=c("blue","red","black","green","orange"),bty="n",cex=1,box.lwd = 2,fill=c("blue","red","black","green","orange"))

    
  for (j in 2:10) {
    lines(x,density_post_1500[[j]][k, ], yaxt= "n",col="red",lwd=1,pch=4,ann=FALSE,yaxt="n")
  }
    
     lines(x, density_post_1500[[ABC_iterations]][k, ], col="black", lwd=1.5)
     lines(x,density_post_Ridge_adjusted[k, ], col="green",lty=3,lwd=2)
      lines(x,density_post_Lasso_adjusted[k, ], col="orange",lty=3,lwd=2)
   
    }

par(mfrow=c(3,2),mar=c(4,4,1,1))

x <- seq(from = -10, to = 7, length.out = 256)#range of prior distribution

for (k in 19:23) {
  plot(x, density_post_1500[[1]][k, ], type="l", ylim=c(0,5),
       xlab=parameter_labels[k],ylab="Density", col="blue",cex.lab=1.2,lwd=2.5,las=1) 

 if(k==19) legend("topright",c("First Prior","Intermediate Priors" ,"Unadjusted Posterior","Ridge Adjustment","Lasso Adjustment"),
         col=c("blue","red","black","green","orange"),bty="n",cex=1,box.lwd = 2,fill=c("blue","red","black","green","orange"))

    
  for (j in 2:10) {
    lines(x,density_post_1500[[j]][k, ], yaxt= "n",col="red",lwd=1,pch=4,ann=FALSE,yaxt="n")
  }
    
     lines(x, density_post_1500[[ABC_iterations]][k, ], col="black", lwd=1.5)
     lines(x,density_post_Ridge_adjusted[k, ], col="green",lty=3,lwd=2)
      lines(x,density_post_Lasso_adjusted[k, ], col="orange",lty=3,lwd=2)
   
    }

Estimator_accuracy_func<- function(true_parameter_values_log, ridge_estimates_log, lasso_estimates_log,
                                  unadj_distn_log, ridge_distn_log, lasso_distn_log,
                                  parameter_labels){
    #number of posterior samples
    n<- dim(unadj_distn_log)[1]  
                                      
    #number of parameters
    n_parameters<- dim(unadj_distn_log)[2] 
                                      
    #bias point estimates based on the unadjusted posterior
    Unadjusted_estimates<- round(exp(apply(unadj_distn_log, 2, mean)), 4) #4 d.p
    bias_est_unadj<- log(Unadjusted_estimates)-true_parameter_values_log
    #Estimated Credible Intervals On log scale (based on the unadjusted posterior) 
    bias_distn_unadj<- matrix(NA, nrow=n, ncol=n_parameters)
    for(i in 1:30) bias_distn_unadj[i, ]<- as.numeric(unadj_distn_log[i, ])-true_parameter_values_log                              
    bias_95CredInt_unadj<-  as.data.frame(ci(as.data.frame(bias_distn_unadj), ci = 0.95,method = "ETI"))[,3:4]
    bias_unadj_df<-  cbind(round(true_parameter_values_log,4),round(log(Unadjusted_estimates),4)
                           ,round(bias_est_unadj,4),round(bias_95CredInt_unadj,4))                                
    names(bias_unadj_df)<- c("True_values_log", "Unadjusted_estimates_log","bias_est_unadj",
                             "bias_Lower95%CredInt_unadj","bias_Upper95%CredInt_unadj")
    rownames(bias_unadj_df)<- paste(parameter_labels) 
    
                                      
    #bias point estimates based on the ridge-adjusted posterior
    bias_est_ridge<- ridge_estimates_log-true_parameter_values_log
    #Estimated Credible Intervals On log scale (based on the ridge-adjusted posterior) 
    bias_distn_ridge<- matrix(NA, nrow=n, ncol=n_parameters)
    for(i in 1:30) bias_distn_ridge[i, ]<- as.numeric(ridge_distn_log[i, ])-true_parameter_values_log                              
    bias_95CredInt_ridge<-  as.data.frame(ci(as.data.frame(bias_distn_ridge), ci = 0.95,method = "ETI"))[,3:4]
    bias_ridge_df<-  cbind(round(true_parameter_values_log,4),round(ridge_estimates_log,4),
                           round(bias_est_ridge,4),round(bias_95CredInt_ridge,4))                                 
    names(bias_ridge_df)<- c("True_values_log", "Ridge_estimates_log","bias_est_ridge",
                             "bias_Lower95%CredInt_ridge","bias_Upper95%CredInt_ridge")
    rownames(bias_ridge_df)<- paste(parameter_labels)
       
    
    #bias point estimates based on the lasso-adjusted posterior
    bias_est_lasso<- lasso_estimates_log-true_parameter_values_log
     #Estimated Credible Intervals On log scale (based on the lasso-adjusted posterior) 
    bias_distn_lasso<- matrix(NA, nrow=n, ncol=n_parameters)
    for(i in 1:30) bias_distn_lasso[i, ]<- as.numeric(lasso_distn_log[i, ])-true_parameter_values_log                              
    bias_95CredInt_lasso<-  as.data.frame(ci(as.data.frame(bias_distn_lasso), ci = 0.95,method = "ETI"))[,3:4]
    bias_lasso_df<-  cbind(round(true_parameter_values_log,4),round(lasso_estimates_log,4),
                           round(bias_est_lasso,4),round(bias_95CredInt_lasso,4))                                 
    names(bias_lasso_df)<- c("True_values_log", "Lasso_estimates_log","bias_est_lasso",
                             "bias_Lower95%CredInt_lasso","bias_Upper95%CredInt_lasso")                                  
                                      
    rownames(bias_lasso_df)<- paste(parameter_labels)                                 
                                      
                                      
                                      
      return(list(bias_unadj_df=bias_unadj_df,bias_ridge_df=bias_ridge_df,bias_lasso_df=bias_lasso_df
                 ))                                
    
}




ABC_accuracy_measures<- Estimator_accuracy_func(true_parameter_values_log=fixed_parameter_values_log, 
                                ridge_estimates_log=log(Post_Ridge_Regression_output[, 1]), 
                                lasso_estimates_log=log(Post_Lasso_Regression_output[, 1]),
                                  unadj_distn_log=Final_posterior_1500,
                                ridge_distn_log=Posterior_Ridge_adj_distn_log, 
                                lasso_distn_log=Posterior_Lasso_adj_distn_log,
                                  parameter_labels=parameter_labels)

ABC_accuracy_measures_unadj<- ABC_accuracy_measures$bias_unadj_df
ABC_accuracy_measures_unadj$Var_est_unadj<- round(as.numeric(apply(Final_posterior_1500, 2, var)),4)
ABC_accuracy_measures_unadj$MSE_est_unadj<- round((ABC_accuracy_measures_unadj$Var_est_unadj+
                                            ABC_accuracy_measures$bias_unadj_df$bias_est_unadj^2),4)


cat("Pooled bias (based on unadjusted posterior estimate):", "\n")
apply(ABC_accuracy_measures$bias_unadj_df[, 3:5],2, mean)


cat("Estimated Variance (based on unadjusted posterior estimate):", "\n")
mean(ABC_accuracy_measures_unadj$Var_est_unadj)


cat("Estimated MSE (based on unadjusted posterior estimate):", "\n")
mean(ABC_accuracy_measures_unadj$MSE_est_unadj)

ABC_accuracy_measures_unadj

write.csv(ABC_accuracy_measures_unadj,"ABC_accuracy_measures_unadj.csv")

#Estimating the model parameters from the unadjusted posterior samples and their corresponding credible intervals

  #computing posterior estimaates and their credible intervals (ETI)
  post_data_1500<-   Final_posterior_1500
  ci_eti_1500<- bayestestR::ci(as.data.frame(post_data_1500), ci = 0.95,method = "ETI") 
  ci_eti_1500<- round(as.data.frame(ci_eti_1500[,3:4]),4)
   # estimating the  posterior  mean estimate
  Posterior_estimates_unadj_1500<- round(ABC_accuracy_measures_unadj[, 2],4)
    

Posterior_est_unadj_data_1500<- data.frame(Posterior_estimate_unadj=Posterior_estimates_unadj_1500,
                                     Cred_Int_lower95_unadj= ci_eti_1500[,1],
                                     Cred_Int_upper95_unadj= ci_eti_1500[,2])


rownames(Posterior_est_unadj_data_1500)<- parameter_labels
#Print unadjusted posterior estimates with 95% credible intervals
Posterior_est_unadj_data_1500

#Save the unadjusted posterior estimate results

write.csv(Posterior_est_unadj_data_1500, "Posterior_estimates_unadj_data_1500.csv")

ABC_accuracy_measures_ridge<- ABC_accuracy_measures$bias_ridge_df
ABC_accuracy_measures_ridge$Var_est_ridge<- round(as.numeric(apply(Posterior_Ridge_adj_distn_log, 2, var)),4)
ABC_accuracy_measures_ridge$MSE_est_ridge<- round((ABC_accuracy_measures_ridge$Var_est_ridge+
                                            ABC_accuracy_measures$bias_ridge_df$bias_est_ridge^2),4)


cat("Pooled bias (based on ridge-adjusted posterior estimate):", "\n")
apply(ABC_accuracy_measures$bias_ridge_df[, 3:5],2, mean)


cat("Estimated Variance (based on ridge-adjusted posterior estimate):", "\n")
mean(ABC_accuracy_measures_ridge$Var_est_ridge)


cat("Estimated MSE (based on ridge-adjusted posterior estimate):", "\n")
mean(ABC_accuracy_measures_ridge$MSE_est_ridge)


ABC_accuracy_measures_ridge
write.csv(ABC_accuracy_measures_ridge,"ABC_accuracy_measures_ridge.csv")

#Estimating the model parameters from the ridge-adjusted posterior samples and their corresponding credible intervals

  #computing posterior estimaates and their credible intervals (ETI)
  post_data_1500<-   Posterior_Ridge_adj_distn_log
  ci_eti_1500<- bayestestR::ci(as.data.frame(post_data_1500), ci = 0.95,method = "ETI") 
  ci_eti_1500<- round(as.data.frame(ci_eti_1500[,3:4]),4)
   # estimating the  posterior  mean estimate
  Posterior_estimates_ridge_1500<- round(log(Post_Ridge_Regression_output[, 1]),4)
    

Posterior_est_ridge_data_1500<- data.frame(Posterior_estimate_ridge=Posterior_estimates_ridge_1500,
                                     Cred_Int_lower95_ridge= ci_eti_1500[,1],
                                     Cred_Int_upper95_ridge= ci_eti_1500[,2])


rownames(Posterior_est_ridge_data_1500)<- parameter_labels
#Print unadjusted posterior estimates with 95% credible intervals
Posterior_est_ridge_data_1500

#Save the unadjusted posterior estimate results

write.csv(Posterior_est_ridge_data_1500, "Posterior_estimates_ridge_data_1500.csv")

ABC_accuracy_measures_lasso<- ABC_accuracy_measures$bias_lasso_df
ABC_accuracy_measures_lasso$Var_est_lasso<- round(as.numeric(apply(Posterior_Lasso_adj_distn_log, 2, var)),4)
ABC_accuracy_measures_lasso$MSE_est_lasso<- round((ABC_accuracy_measures_lasso$Var_est_lasso+
                                            ABC_accuracy_measures$bias_lasso_df$bias_est_lasso^2),4)


cat("Pooled bias (based on lasso-adjusted posterior estimate):", "\n")
apply(ABC_accuracy_measures$bias_lasso_df[, 3:5],2, mean)


cat("Estimated Variance (based on lasso-adjusted posterior estimate):", "\n")
mean(ABC_accuracy_measures_lasso$Var_est_lasso)


cat("Estimated MSE (based on lasso-adjusted posterior estimate):", "\n")
mean(ABC_accuracy_measures_lasso$MSE_est_lasso)

ABC_accuracy_measures_lasso
write.csv(ABC_accuracy_measures_lasso,"ABC_accuracy_measures_lasso.csv")

#Estimating the model parameters from the lasso-adjusted posterior samples and their corresponding credible intervals

  #computing posterior estimaates and their credible intervals (ETI)
  post_data_1500<-   Posterior_Lasso_adj_distn_log
  ci_eti_1500<- bayestestR::ci(as.data.frame(post_data_1500), ci = 0.95,method = "ETI") 
  ci_eti_1500<- round(as.data.frame(ci_eti_1500[,3:4]),4)
   # estimating the  posterior  mean estimate
  Posterior_estimates_lasso_1500<- round(log(Post_Lasso_Regression_output[, 1]),4)
    

Posterior_est_lasso_data_1500<- data.frame(Posterior_estimate_lasso=Posterior_estimates_lasso_1500,
                                     Cred_Int_lower95_lasso= ci_eti_1500[,1],
                                     Cred_Int_upper95_lasso= ci_eti_1500[,2])


rownames(Posterior_est_lasso_data_1500)<- parameter_labels
#Print unadjusted posterior estimates with 95% credible intervals
Posterior_est_lasso_data_1500

#Save the unadjusted posterior estimate results

write.csv(Posterior_est_lasso_data_1500, "Posterior_estimates_lasso_data_1500.csv")

#Setting working directory
setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro/Posterior_Paramater identifiabity")


ABC_accuracy_measures_unadj<- read.csv("ABC_accuracy_measures_unadj.csv")
ABC_accuracy_measures_ridge<- read.csv("ABC_accuracy_measures_ridge.csv")
ABC_accuracy_measures_lasso<- read.csv("ABC_accuracy_measures_lasso.csv")

parameter_labels=c(expression(paste("b"[11])),expression(paste("b"[12])),
                   expression(paste("b"[21])),expression(paste("b"[22])), 
                   expression(paste("b"[31])),expression(paste("b"[32])),
                   expression(paste("d"[11])),expression(paste("d"[12])),
                   expression(paste("d"[21])),expression(paste("d"[22])), 
                   expression(paste("d"[31])),expression(paste("d"[32])),
                   expression(paste("m")),  expression(paste("r")),
                   expression(paste("r"[1])),expression(paste("r"[2])),
                   expression(paste("r"[3])),expression(paste("s")),
                   expression(paste("s"[1])),expression(paste(epsilon[1])),
                   expression(paste(epsilon[2])),expression(paste(epsilon[3])),
                   expression(paste(kappa))
)


####Setting the plot size and resolultion (300 dpi)
options(repr.plot.width=8, repr.plot.height=8,repr.plot.res = 300) #Setting plot size

#o<-par(mar=c(0,4,2,2))#Run this before the code below

#nf<-layout(matrix(1:3, nrow=3,ncol=1))

# rest parameters
#par(o)#
#o<-par(mar=c(0,4,2,2))    


my_cols <- c("black","green","orange")

plot(ABC_accuracy_measures_unadj[, 4],type="p",col=my_cols[1],lwd=1,xaxt = "n",
     ylab="",xlab="",
     ylim=c(-3,6),pch=19,cex.lab=.9, cex.axis=.9)

arrows(x0=1:23, y0=ABC_accuracy_measures_unadj[, 5], x1=1:23, 
       y1=ABC_accuracy_measures_unadj[, 6], 
       code=3,  angle=90, length=0.05,col =my_cols[1],lwd=1)
 axis(1, at=1:23, labels=parameter_labels,font=3, hadj = .9,cex.axis=.8)

#abline(h=0, col="red",lwd=1.5,lty=2)

mtext(text = expression(paste("Bias(",hat(theta),  ") with 95% credible interval")),
      side = 2,
      line = 2.5,cex=1,font=2)

mtext(text = "Model parameters",
      side = 1,
      line = 2.5,cex=1,font=1)

legend(x=0.5,y=6.3,legend=c("Unadjusted posterior estimates",
                        "Ridge-adjusted posterior estimates", "Lasso-adjusted posterior estimates"),
        col=c("black","green","orange"),bty="n",
       cex=.85,box.lwd = .8,fill=c("black","green","orange"),
       ncol =1,title="Estimated bias based on:")


lines(ABC_accuracy_measures_ridge[, 4],type="p",col=my_cols[2],lwd=1,xaxt = "n",ylab="",xlab="",
    pch=19,cex.lab=1, cex.axis=1)

arrows(x0=1:23, y0=ABC_accuracy_measures_ridge[, 5], x1=1:23, 
       y1=ABC_accuracy_measures_ridge[, 6], 
       code=3,  angle=90, length=0.05,col =my_cols[2],lwd=1)



lines(ABC_accuracy_measures_lasso[, 4],type="p",col=my_cols[3],lwd=1,xaxt = "n",ylab="",xlab="",
    pch=19,cex.lab=1, cex.axis=1)

arrows(x0=1:23, y0=ABC_accuracy_measures_lasso[, 5], x1=1:23, 
       y1=ABC_accuracy_measures_lasso[, 6], 
       code=3,  angle=90, length=0.05,col =my_cols[3],lwd=1)


text(6.2,4.6, expression(paste("Pooled bias(",theta["unadj"],  ") = -0.025 (95% PI: -0.274, 0.219)")), 
                col=my_cols[1],font=2,cex=.85)
     
#abline(h=-0.025, col=my_cols[1],lwd=2,lty=2)

abline(h=0, col="grey",lwd=.8,lty=2)


text(6.2,4.3, expression(paste("Pooled bias(",theta["ridge"],  ") = -0.033 (95% PI: -0.267, 0.219)")), 
     col=my_cols[2],font=2,cex=.9)

#abline(h=-0.033, col=my_cols[2],lwd=1,lty=2)




text(6.2,4,expression(paste("Pooled bias(",theta["lasso"],  ") = 0.125 (95% PI: -0.096, 0.399)")),
     col=my_cols[3],font=2,cex=.9)
#abline(h=0.125,col=my_cols[3],lwd=1,lty=2)


#Variance
text(18, 6.1,expression(paste("Pooled Var(",hat(theta)["unadj"],  ") = 0.0209")),col=my_cols[1],font=2)
text(18, 5.7,expression(paste("Pooled Var(",hat(theta)["ridge"],  ") = 0.0203")),col=my_cols[2],font=2)
text(18, 5.3,expression(paste("Pooled Var(",hat(theta)["lasso"],  ") = 0.0253")),col=my_cols[3],font=2)


#MSE
text(18, 4.7,expression(paste("Pooled MSE(",hat(theta)["unadj"],  ") = 0.0783")),col=my_cols[1],font=2)
text(18, 4.3,expression(paste("Pooled MSE(",hat(theta)["ridge"],  ") = 0.0918")),col=my_cols[2],font=2)
text(18, 3.9,expression(paste("Pooled MSE(",hat(theta)["lasso"],  ") = 0.8087")),col=my_cols[3],font=2)





Unadj_posterior<- exp(Posterior_Unadj__distn_log)
Ridge_posterior<- exp(Posterior_Ridge_adj_distn_log)
Lasso_posterior<- exp(Posterior_Lasso_adj_distn_log)


Posterior_combined_df<- rbind(Unadj_posterior,Ridge_posterior,Lasso_posterior)
#Posterior_combined_df

Unadj_posterior_log<- log(Unadj_posterior)
Ridge_posterior_log<-log(Ridge_posterior)
Lasso_posterior_log<- log(Lasso_posterior)

library(vegan)


#The function computes dissimilarity indices that are useful for or popular with community ecologists
diss_indices <- vegdist(Posterior_combined_df)

#The ABC methods being compared
ABC_groups <- c( rep("Unadjusted ABC",dim(Unadj_posterior_log)[1]),rep("Ridge",dim(Ridge_posterior_log)[1]), 
            rep("Lasso",dim(Lasso_posterior_log)[1]))

ABC_groups<- as.factor(ABC_groups)
ABC_groups <- factor(ABC_groups, levels = c("Unadjusted ABC", "Ridge", "Lasso"))
#ABC_groups<- relevel(ABC_groups,ref="Unadjusted ABC")
levels(ABC_groups)



## Calculate multivariate dispersions
multivariate_dispersion_mod <- betadisper(d=diss_indices, group=ABC_groups,  
                                          bias.adjust = T)
multivariate_dispersion_mod 

## Perform test
anova(multivariate_dispersion_mod)[1, ]


permutest(multivariate_dispersion_mod, permutations = 99)

TukeyHSD(multivariate_dispersion_mod)

## with data ellipses instead of hulls
par(mfrow=c(2,1), mar=c(4,4,2,2))
my_cols <- c("black","green","orange")

plot(multivariate_dispersion_mod, ellipse = F, hull = T,label = TRUE, label.cex = .54,main="",col=my_cols, 
xlab="PCoA 1",ylab="PCoA 2",lwd=2, pch=c(6,8, 19),seg.lwd=.5,ylim=c(-0.1,.1)) # 1 sd data ellipse


text(-0.02, .1, "Multivariate homogeneity test of ABC posterior variances: p-value=0.999", col="black", font=2)
mtext(LETTERS[1], adj=0, line=1,font=2)



boxplot(multivariate_dispersion_mod,col=my_cols, xlab="ABC methods",ylab="Distance to posterior centroid")
mtext(paste(LETTERS[2]), adj=0, line=1,font=2)

print(paste("ABC unadjusted posterior estimates at N=",1500))
Parameter_estimates_unadj<- Post_Ridge_Regression_output[,4]
print(round(Parameter_estimates_unadj,4))

print(paste("Proposed Ridge adjusted posterior estimates at N=",1500))
Parameter_estimates_ridge<- Post_Ridge_Regression_output[,1]
print(round(Parameter_estimates_ridge,4))


print(paste("Proposed Lasso adjusted posterior estimates at N=",1500))
Parameter_estimates_lasso<- Post_Lasso_Regression_output[,1]
print(round(Parameter_estimates_lasso,4))

#A function to create dataframe of parasite numbers over time for either the simulated or observed data
Parasite_data_summary<- function(data_all_fish, observed_times= 1:9,
                                 Total_fish=sum(unlist(numF)),
                                 parasite_fish=parasite_fish,sim=TRUE,
                                 fishSex=fishSex,fishID=fishID,
                                 fishSize=fishSize,Strain=Strain,
                                 Fish_stock=Fish_stock,Group){
  #Group= the comparable groups observed data vs ABC methods
  
  sum_func<- function(x) sum(x, na.rm=TRUE)
  #mean_func<- function(x) mean(x, na.rm=TRUE)
  parasites_data<- data.frame(matrix(NA,nrow=Total_fish, ncol=length(observed_times)+7)) 
  fish_index<- 0
  
  for(pf in 1:9){
    for(i in 1:numF[[pf]]){
      fish_index <- fish_index +1
      #Total number of parasites at each timepoint per fish
      if(sim==TRUE) parasites_per_fish<- apply(data_all_fish$pop_sim[[pf]][ i, , ], 2, sum_func)
      else if(sim==FALSE) parasites_per_fish<- apply(data_all_fish[[pf]][ i, , ], 2, sum_func)
      
      parasites_data[fish_index , 1:9 ] <- as.vector( parasites_per_fish)
      parasites_data[fish_index , 10  ] <- parasite_fish[pf]
      parasites_data[fish_index , 11  ] <- fishSex[[pf]][i]
      parasites_data[fish_index , 12  ] <- fishID[[pf]][i]
      parasites_data[fish_index , 13  ] <- fishSize[[pf]][i]
      parasites_data[fish_index , 14  ] <- Strain[[pf]][i]
      parasites_data[fish_index , 15  ] <- Fish_stock[[pf]][i]
      parasites_data[fish_index,  16  ] <- paste0(Group)
      
    }
  }
  
  names(parasites_data)<- c(paste0("Parasites_t",observed_times), "Parafish_group","Sex",
                            "FishID","Fish_size","Parasite_strain","Fish_stock","Group")
  return(parasites_data)
}





#A function to combine observed and simulated datasets
Combine_sim_obs_data_func<- function(obs_data,sim_unadj,sim_ridge,sim_lasso, rep=rep, combine_all=TRUE){
  Combine_sim_obs_data<- NULL
  for(i in 1:rep){
    Combine_sim_obs_data[[i]]<- rbind(obs_data,sim_unadj[[i]],sim_ridge[[i]],sim_lasso[[i]] )
    Combine_sim_obs_data[[i]]$Replicate_index<-  paste0(i)
  }
  if(combine_all==FALSE) return(Combine_sim_obs_data)
  else if( combine_all==TRUE) return(do.call("rbind",Combine_sim_obs_data))
}

#Simulations based on the unadjusted posterior##########################
set.seed(1)
rep<-1 #Number of replicates/repetition
sim_func_unadj<- function(i) SimGroup_tauleap(theta1=log(Parameter_estimates_unadj),
            fish_sex=fishSex,fish_type=Fish_stock,strain=Strain,fish_size=fishSize,error=0.002)


output_sim_unadj1<- mclapply(1:rep, sim_func_unadj,mc.cores=numCores)

#Simulations based on the Ridge adjusted posterior#########################
set.seed(1)
rep<-1 #Number of replicates/repetition
sim_func_ridge<- function(i) SimGroup_tauleap(theta1=log(Parameter_estimates_ridge),
            fish_sex=fishSex,fish_type=Fish_stock,strain=Strain,fish_size=fishSize,error=0.002)

output_adj_ridge1<- mclapply(1:rep, sim_func_ridge,mc.cores=numCores)

#Simulations based on the Lasso adjusted posterior
set.seed(1)
rep<-1 #Number of replicates/repetition
sim_func_lasso<- function(i) SimGroup_tauleap(theta1=log(Parameter_estimates_lasso),
            fish_sex=fishSex,fish_type=Fish_stock,strain=Strain,fish_size=fishSize,error=0.002)

output_adj_lasso1<- mclapply(1:rep, sim_func_lasso,mc.cores=numCores)

Sim_data_unadj_reg<- NULL; Sim_data_ridge_reg<- NULL;Sim_data_lasso_reg<- NULL

for(i in 1:rep){
    #Simulations based on unadjusted posterior estimates
    Sim_data_unadj_reg[[i]]<- Parasite_data_summary(data_all_fish=output_sim_unadj1[[i]], 
                    observed_times= 1:9,Total_fish=sum(unlist(numF)),parasite_fish=parasite_fish,sim=TRUE,
                     fishSex=fishSex,fishID=fishID,fishSize=fishSize,Strain=Strain,Fish_stock=Fish_stock,
                    Group="Unadjusted ABC")
    
    Sim_data_ridge_reg[[i]]<- Parasite_data_summary(data_all_fish=output_adj_ridge1[[i]], 
                    observed_times= 1:9,Total_fish=sum(unlist(numF)),parasite_fish=parasite_fish,sim=TRUE,
                     fishSex=fishSex,fishID=fishID,fishSize=fishSize,Strain=Strain,Fish_stock=Fish_stock,
                    Group="Ridge adjustment")
    
    Sim_data_lasso_reg[[i]]<- Parasite_data_summary(data_all_fish=output_adj_lasso1[[i]], 
                    observed_times= 1:9,Total_fish=sum(unlist(numF)),parasite_fish=parasite_fish,sim=TRUE,
                     fishSex=fishSex,fishID=fishID,fishSize=fishSize,Strain=Strain,Fish_stock=Fish_stock,
                    Group="Lasso adjustment")
    
    }


#Observed data
Observed_data_reg<- Parasite_data_summary(data_all_fish=pop_pseudo_obs, observed_times= 1:9,
                    Total_fish=sum(unlist(numF)),parasite_fish=parasite_fish,sim=FALSE,
                     fishSex=fishSex,fishID=fishID,
                fishSize=fishSize,Strain=Strain,Fish_stock=Fish_stock,Group="Pseudo-observed data")

dim(Observed_data_reg)

head(Observed_data_reg)

#Combined data for mixed-effect regression model
Combine_sim_obs_data_reg<- Combine_sim_obs_data_func(obs_data=Observed_data_reg,sim_unadj= Sim_data_unadj_reg,
                    sim_ridge=Sim_data_ridge_reg,sim_lasso=Sim_data_lasso_reg, rep=rep, combine_all=TRUE)

dim(Combine_sim_obs_data_reg)

names(Combine_sim_obs_data_reg)

head(Combine_sim_obs_data_reg)

transform_structure<- function(data){
    data.timepoints<- NULL; times<- seq(1,17,by=2)
    
    for(t in 1:9){
        data.timepoints[[t]]<- data[, c(t,10:17)]
        data.timepoints[[t]][,1]<- log(data.timepoints[[t]][,1]+1) #on log scale
        names(data.timepoints[[t]])[1]<- "Parasites"
        data.timepoints[[t]]$Time<- time[t]
    }
    
    
    return(do.call("rbind",data.timepoints))
}

transdata_reg<- transform_structure(Combine_sim_obs_data_reg)
#Converting all character variables to factors
transdata_reg[sapply(transdata_reg, is.character)] <- lapply(transdata_reg[sapply(transdata_reg, is.character)], 
                                                           as.factor)



transdata_reg$Group<- relevel(transdata_reg$Group,ref="Pseudo-observed data")
transdata_reg$Sex<- relevel(transdata_reg$Sex,ref="Female")
transdata_reg$Parasite_strain<- relevel(transdata_reg$Parasite_strain,ref="Gt")
transdata_reg$Fish_stock<- relevel(transdata_reg$Fish_stock,ref="UA")


transdata_reg$Group<- factor(transdata_reg$Group,levels=c("Pseudo-observed data",
                        "Unadjusted ABC","Ridge adjustment","Lasso adjustment"))

levels(transdata_reg$Group)

transdata_reg_original<-transdata_reg
transdata_reg_original[,1]<- exp(transdata_reg_original[,1])#inverse log transformation



theme_set(
  theme_bw() +
    theme(legend.position = "top")
  )



p1<- transdata_reg%>%
  mutate(type=fct_reorder(as.factor(Time),Parasites),
        Group=fct_reorder(as.factor(Group),Parasites)) %>%
     ggplot(fill = factor(Group, levels=positions))+
geom_smooth(aes(x=Time, y=Parasites,color=Group,group=Group),
                          method = "loess",se = FALSE)+
labs(y="Parasite mean intensity") +
theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),

        panel.background = element_blank())+
scale_color_manual(values=c("black", "blue","green","orange"))+
xlim(1,17)+facet_wrap(~Parasite_strain+Fish_stock,  ncol=3)+
  guides(color = guide_legend(title = "Comparable groups:"))


p1

#transdata_reg

transdata_reg$Day<- as.factor(transdata_reg$Time)
levels(transdata_reg$Day)<- c("Day 1","Day 3","Day 5", "Day 7", "Day 9",
                              "Day 11","Day 13","Day 15","Day 17")


head(transdata_reg)

my_cols<- c("blue","black","green","orange")

# Visualize: Specify the comparisons you want
my_comparisons <- list( c("Pseudo-observed data", "Unadjusted ABC"),
                       c("Pseudo-observed data", "Ridge adjustment"),
             c("Pseudo-observed data", "Lasso adjustment"))

#subset(transdata_reg, Time==1)
bxp1<- ggboxplot(transdata_reg, x = "Group", y = "Parasites",
          fill = "Group")+scale_fill_manual(values=my_cols) #you can change "fill" to "color"

bxp1<- bxp1+stat_compare_means(comparisons = my_comparisons, hide.ns = F,
    label = c("p.signif"),palette = c('red'),label.x.npc ='center',#ref.group =".all."
  label.y.npc ='top')+ # Add pairwise comparisons p-value
  stat_compare_means(method = "kruskal.test", label.y = 15,label.x=1.5,size = 3)+     # Add global p-value
 theme(text = element_text(size = 13))+  
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
facet_wrap( ~ Day)

bxp1<- bxp1+labs(x="Groups",y="Parasite abundance (log-scale)",fill = "Groups:")+
geom_dotplot(binaxis = "y", stackdir = "center", binwidth = .05)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())




bxp1



library("FactoMineR")
library("factoextra")

data_obs<- Observed_data_reg[, 1:9]
data_unadj<- Sim_data_unadj_reg[[1]][, 1:9]
data_ridge<- Sim_data_ridge_reg[[1]][, 1:9]
data_lasso<- Sim_data_lasso_reg[[1]][, 1:9]

combined_data_parasite_numbers<- rbind(data_obs,data_unadj,data_ridge,data_lasso)
colnames(combined_data_parasite_numbers)<- c("Parasites_Day1","Parasites_Day3",
                "Parasites_Day5","Parasites_Day7","Parasites_Day9","Parasites_Day11",
                "Parasites_Day13","Parasites_Day15","Parasites_Day17")

group_ABC.methods<- c(Observed_data_reg$Group, Sim_data_unadj_reg[[1]]$Group,
                     Sim_data_ridge_reg[[1]]$Group, Sim_data_lasso_reg[[1]]$Group)


group_ABC.methods<- as.factor(group_ABC.methods)
group_ABC.methods<- factor(group_ABC.methods,levels=c("Pseudo-observed data",
                        "Unadjusted ABC","Ridge adjustment","Lasso adjustment"))


levels(group_ABC.methods)



group_Parafish_group<- c(Observed_data_reg$Parafish_group, Sim_data_unadj_reg[[1]]$Parafish_group,
                     Sim_data_ridge_reg[[1]]$Parafish_group, Sim_data_lasso_reg[[1]]$Parafish_group)



parasite_numbers.pca <- PCA(combined_data_parasite_numbers, graph = FALSE)
print(parasite_numbers.pca)

#The sum of all the eigenvalues give a total variance
get_eigenvalue(parasite_numbers.pca )#: Extract the eigenvalues/variances of principal components

var <- get_pca_var(parasite_numbers.pca)
var

# Coordinates
head(var$coord)
# Cos2: quality on the factore map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)

fviz_cos2(parasite_numbers.pca, choice = "var", axes = 1:2, sort.val="none")+
            labs(title="Cos2 of variables to Dim 1-2")



##Scree plot
#Proportion of variance explained by each principal component
#p1=fviz_eig(parasite_numbers.pca, addlabels = TRUE, ylim = c(0, 65),xlim=c(1,10),
   #         hjust = 0.3,main="A: Scree plot")+
#theme(axis.line = element_line(),
#        panel.grid.major = element_blank(),
 #       panel.grid.minor = element_blank(),

   #     panel.background = element_blank())


p1<- fviz_cos2(parasite_numbers.pca, choice = "var", axes = 1:2)+
            labs(title="A: Cos2 of variables to Dim 1-2")+
theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),

        panel.background = element_blank())

#Biplot
# Color by cos2 values: quality on the factor map
#Cos2 of variables to Dim 1-2
p2<- fviz_pca_biplot(parasite_numbers.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, # Avoid text overlapping
                   label = "var" , addEllipses=T , ellipse.level=0.95,title="B: PCA- Biplot")+
theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),

        panel.background = element_blank())



gridExtra::grid.arrange(p1,p2, ncol=1,nrow=2)

#Parasite distribution/numbers  - (Individual- PCA)
p1<- fviz_pca_ind(parasite_numbers.pca,
             # label = "none", # hide individual labels
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = group_ABC.methods, # color by groups
             palette = c("blue","black","green","orange"),
             addEllipses = TRUE, # Concentration ellipses
            #ellipse.type = "confidence",
             label = "var",
             legend.title = "Groups:",
             axes=c(1,2),
            repel = TRUE,
            title=""
             )+
theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),

        panel.background = element_blank())


p1+theme(legend.position="top")
