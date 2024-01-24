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

####Setting the plot size and resolultion (300 dpi)
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

#Script for computing weighted Euclidean distance (between simulated and observed summaries)
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
pop_obs <-  Descriptors$pop_obs # observed parasite numbers for each parasite-fish group 
alive_obs<- Descriptors$alive_obs # observed surviva status for each parasite-fish group 
Area_normalized<- Body_area(Area_data=Bodyparts_area)#body areas

#Parasite-fish groups
parasite_fish<-c("Gt3-OS","Gt3-LA","Gt3-UA","Gt-OS","Gt-LA","Gt-UA","Gb-OS","Gb-LA","Gb-UA") 

#B-D-C parameter estimates for the parasite-fish groups based on observed data
BDC_estimates_obs<- GW_GMM_BDCestimator(X0=2,pop=pop_obs,alive=alive_obs,
                                        group=parasite_fish)$BDC_estimates
BDC_estimates_obs

#Computing summary statistics for observed data across the parasite-fish groups
summaries_obs<-Summary_stats(pop=pop_obs,alive=alive_obs,BDC_estimates=BDC_estimates_obs)

#View ABC summary statistics for the observed data
#summaries_obs[[1]]# for group 1 or Gt3-0S group

S_obs<- do.call("rbind",summaries_obs)
dim(S_obs)

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

#Monte Carlo sample sizes of N=500

setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro/Updated ABC_results_N_500")


#Importing saved outputs from ABC-SMC fitting
ABC_time_steps<-10
N<- 500
#Importing saved outputs   
distances_500<- list();theta_500<- list();theta_posterior_500<- list()
  for(t in 1:ABC_time_steps){  
     distances_500[[t]]<- read.csv(paste0("weighted_distance_",t,"_",N,".csv"))
     distances_500[[t]]<-distances_500[[t]][,-1]
     theta_500[[t]]<- read.csv(paste0("theta_t_",t,".csv")) 
     theta_500[[t]]<- theta_500[[t]][,-1]  
     theta_posterior_500[[t]]<- read.csv(paste0("theta_post_",t,"_",N,".csv"))
     theta_posterior_500[[t]]<-theta_posterior_500[[t]][,-1] 
    }


#Final posterior  
Final_posterior_500<- read.csv(paste0("PosteriorFinal_",N, ".csv"))
Final_posterior_500<-Final_posterior_500[,-1]

# Importing prior, intermediate priors and posterior densities
#Importing posterior densities from Weighted-iterative ABC with Sequential Monte Carlo with importance sampling

density_post_500<-NULL 
ABC_iterations<-ABC_time_steps+1

no.params<- length(parameter_labels)

for(j in 1:ABC_iterations) {
    density_post_500[[j]]<-read.csv(file=paste0("density_post_",j,"_",N,".csv"))
    density_post_500[[j]]<-density_post_500[[j]][,-1]
          }

#Monte Carlo sample sizes of N=1000

setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro/Updated ABC_results_N_1000")


#Importing saved outputs from ABC-SMC fitting
ABC_time_steps<-10
N<- 1000
#Importing saved outputs   
distances_1000<- list();theta_1000<- list();theta_posterior_1000<- list()
  for(t in 1:ABC_time_steps){  
     distances_1000[[t]]<- read.csv(paste0("weighted_distance_",t,"_",N,".csv"))
     distances_1000[[t]]<-distances_1000[[t]][,-1]
     theta_1000[[t]]<- read.csv(paste0("theta_t_",t,".csv")) 
     theta_1000[[t]]<- theta_1000[[t]][,-1]  
     theta_posterior_1000[[t]]<- read.csv(paste0("theta_post_",t,"_",N,".csv"))
     theta_posterior_1000[[t]]<-theta_posterior_1000[[t]][,-1] 
    }


#Final posterior  
Final_posterior_1000<- read.csv(paste0("PosteriorFinal_",N, ".csv"))
Final_posterior_1000<-Final_posterior_1000[,-1]

# Importing prior, intermediate priors and posterior densities
#Importing posterior densities from Weighted-iterative ABC with Sequential Monte Carlo with importance sampling

density_post_1000<-NULL 
ABC_iterations<-ABC_time_steps+1

no.params<- length(parameter_labels)

for(j in 1:ABC_iterations) {
    density_post_1000[[j]]<-read.csv(file=paste0("density_post_",j,"_",N,".csv"))
    density_post_1000[[j]]<-density_post_1000[[j]][,-1]
          }

#Monte Carlo sample sizes of N=1500

setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro/Updated ABC_results_N_1500")



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

setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro/Updated ABC_results_N_500")


#Import the total CPU time to fit the complex model using the ABC-SMC sampler
N<- 500 # Monte Carlo Sample Size
time_CPU_500<- read.csv(file=paste0("CPUtime_",N,".csv"))
time_CPU_500<- time_CPU_500[,-1]
    

#Computational times
 print(paste("CPU time at N=", N,"",":",time_CPU_500/10,"secs"))
 print(paste("CPU time at N=", N,"",":",time_CPU_500/864000,"days"))


#Converting from seconds to days
library(lubridate)
print(paste("CPU time at N=", N,"",":",seconds_to_period(time_CPU_500/10)))

setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro/Updated ABC_results_N_1000")


#Import the total CPU time to fit the complex model using the ABC-SMC sampler
N<- 1000 # Monte Carlo Sample Size
time_CPU_1000<- read.csv(file=paste0("CPUtime_",N,".csv"))
time_CPU_1000<- time_CPU_1000[,-1]
    

#Computational times
 print(paste("CPU time at N=", N,"",":",time_CPU_1000/10,"secs"))
 print(paste("CPU time at N=", N,"",":",time_CPU_1000/864000,"days"))


#Converting from seconds to days
library(lubridate)
print(paste("CPU time at N=", N,"",":",seconds_to_period(time_CPU_1000/10)))

setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro/Updated ABC_results_N_1500")


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

#n=500
dist_0.05quantile_500<-rep(NA,ABC_time_steps)
dist_0.15quantile_500<-rep(NA,ABC_time_steps)
dist_0.25quantile_500<-rep(NA,ABC_time_steps)
dist_0.50quantile_500<-rep(NA,ABC_time_steps)
dist_0.75quantile_500<-rep(NA,ABC_time_steps)
dist_0.85quantile_500<-rep(NA,ABC_time_steps)

for(t in 1:ABC_time_steps){
   dist_0.05quantile_500[t]<- quantile(distances_500[[t]],probs =0.05)
   dist_0.15quantile_500[t]<- quantile(distances_500[[t]],probs =0.15)
   dist_0.25quantile_500[t]<- quantile(distances_500[[t]],probs =0.25)
   dist_0.50quantile_500[t]<- quantile(distances_500[[t]],probs =0.5)
   dist_0.75quantile_500[t]<- quantile(distances_500[[t]],probs =0.75)
   dist_0.85quantile_500[t]<- quantile(distances_500[[t]],probs =0.85) 
}


plot(1:ABC_time_steps,dist_0.05quantile_500,type="b",lwd=3,col="blue",
     ylim=c(0.6,.85),xlab="",pch=3,
    main="", ylab="")
#axis(side=1, at = 1:ABC_time_steps, labels = paste(1:ABC_time_steps))
lines(1:ABC_time_steps,dist_0.15quantile_500,type="b",lwd=3,col="green",pch=4)
lines(1:ABC_time_steps,dist_0.25quantile_500,type="b",lwd=3,col="red",pch=13)
lines(1:ABC_time_steps,dist_0.50quantile_500,type="b",lwd=3,col="black",pch=6)
lines(1:ABC_time_steps,dist_0.75quantile_500,type="b",lwd=3,col="pink",pch=5)
lines(1:ABC_time_steps,dist_0.85quantile_500,type="b",lwd=3,col="orange",pch=10)

#text(3,.85,"At Monte Carlo sample size of N=500",font=2,cex=1.5)


mtext(text = "ABC-SMC time steps",
      side = 1,
      line = 2.5,cex=1,font=2)


mtext(text = "Quantiles of ABC-SMC distances",
      side = 2,
      line = 2.5,cex=1,font=2)


legend(x=2,y=.85,legend=c("5th percentile","15th percentile",
                         "25th percentile","50th percentile","75th percentile","85th percentile"),
        col=c("blue","green","red","black","pink","orange"),pch=c(3,4,13,6,5,10),
       cex=1,box.lwd = 2,ncol =3,title="Percentiles:",bty="n")


dist_0.05quantile_1000<-rep(NA,ABC_time_steps)
dist_0.15quantile_1000<-rep(NA,ABC_time_steps)
dist_0.25quantile_1000<-rep(NA,ABC_time_steps)
dist_0.50quantile_1000<-rep(NA,ABC_time_steps)
dist_0.75quantile_1000<-rep(NA,ABC_time_steps)
dist_0.85quantile_1000<-rep(NA,ABC_time_steps)

for(t in 1:ABC_time_steps){
   dist_0.05quantile_1000[t]<- quantile(distances_1000[[t]],probs =0.05)
   dist_0.15quantile_1000[t]<- quantile(distances_1000[[t]],probs =0.15)
   dist_0.25quantile_1000[t]<- quantile(distances_1000[[t]],probs =0.25)
   dist_0.50quantile_1000[t]<- quantile(distances_1000[[t]],probs =0.5)
   dist_0.75quantile_1000[t]<- quantile(distances_1000[[t]],probs =0.75)
   dist_0.85quantile_1000[t]<- quantile(distances_1000[[t]],probs =0.85) 
}


plot(1:ABC_time_steps,dist_0.05quantile_1000,type="b",lwd=3,col="blue",
     ylim=c(0.55,.85),xlab="",pch=3,
     main="", ylab="")
#axis(side=1, at = 1:ABC_time_steps, labels = paste(1:ABC_time_steps))
lines(1:ABC_time_steps,dist_0.15quantile_1000,type="b",lwd=3,col="green",pch=4)
lines(1:ABC_time_steps,dist_0.25quantile_1000,type="b",lwd=3,col="red",pch=13)
lines(1:ABC_time_steps,dist_0.50quantile_1000,type="b",lwd=3,col="black",pch=6)
lines(1:ABC_time_steps,dist_0.75quantile_1000,type="b",lwd=3,col="pink",pch=5)
lines(1:ABC_time_steps,dist_0.85quantile_1000,type="b",lwd=3,col="orange",pch=10)

mtext(text = "ABC-SMC time steps",
      side = 1,
      line = 2.5,cex=1,font=2)


mtext(text = "Quantiles of ABC-SMC distances",
      side = 2,
      line = 2.5,cex=1,font=2)


legend(x=2,y=.85,legend=c("5th percentile","15th percentile",
                         "25th percentile","50th percentile","75th percentile","85th percentile"),
        col=c("blue","green","red","black","pink","orange"),pch=c(3,4,13,6,5,10),
       cex=1,box.lwd = 2,ncol =3,title="Percentiles:",bty="n")


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
     ylim=c(0.6,.85),xlab="",pch=3
     ,main="", ylab="", xaxt ="n")
axis(side=1, at = 1:ABC_time_steps, labels = paste(1:ABC_time_steps))
lines(1:ABC_time_steps,dist_0.15quantile_1500,type="b",lwd=3,col="green",pch=4)
lines(1:ABC_time_steps,dist_0.25quantile_1500,type="b",lwd=3,col="red",pch=13)
lines(1:ABC_time_steps,dist_0.50quantile_1500,type="b",lwd=3,col="black",pch=6)
lines(1:ABC_time_steps,dist_0.75quantile_1500,type="b",lwd=3,col="pink",pch=5)
lines(1:ABC_time_steps,dist_0.85quantile_1500,type="b",lwd=3,col="orange",pch=10)

mtext(text = "ABC-SMC time steps",
      side = 1,
      line = 2.5,cex=1,font=2)


mtext(text = "Quantiles of ABC-SMC distances",
      side = 2,
      line = 2.5,cex=1,font=2)


legend(x=2,y=.85,legend=c("5th percentile","15th percentile",
                         "25th percentile","50th percentile","75th percentile","85th percentile"),
        col=c("blue","green","red","black","pink","orange"),pch=c(3,4,13,6,5,10),
       cex=1,box.lwd = 2,ncol =3,title="Percentiles:",bty="n")


o<-par(mar=c(0,4,2,2))#Run this before the code below

nf<-layout(matrix(1:3, nrow=3,ncol=1))

# rest parameters
par(o)#
o<-par(mar=c(0,4,2,2))    


#n=500
dist_0.05quantile_500<-rep(NA,ABC_time_steps)
dist_0.15quantile_500<-rep(NA,ABC_time_steps)
dist_0.25quantile_500<-rep(NA,ABC_time_steps)
dist_0.50quantile_500<-rep(NA,ABC_time_steps)
dist_0.75quantile_500<-rep(NA,ABC_time_steps)
dist_0.85quantile_500<-rep(NA,ABC_time_steps)

for(t in 1:ABC_time_steps){
   dist_0.05quantile_500[t]<- quantile(distances_500[[t]],probs =0.05)
   dist_0.15quantile_500[t]<- quantile(distances_500[[t]],probs =0.15)
   dist_0.25quantile_500[t]<- quantile(distances_500[[t]],probs =0.25)
   dist_0.50quantile_500[t]<- quantile(distances_500[[t]],probs =0.5)
   dist_0.75quantile_500[t]<- quantile(distances_500[[t]],probs =0.75)
   dist_0.85quantile_500[t]<- quantile(distances_500[[t]],probs =0.85) 
}


plot(1:ABC_time_steps,dist_0.05quantile_500,type="b",lwd=3,col="blue",
     ylim=c(0.6,.85),xlab="",pch=3,
    main="", ylab="", xaxt ="n")
#axis(side=1, at = 1:ABC_time_steps, labels = paste(1:ABC_time_steps))
lines(1:ABC_time_steps,dist_0.15quantile_500,type="b",lwd=3,col="green",pch=4)
lines(1:ABC_time_steps,dist_0.25quantile_500,type="b",lwd=3,col="red",pch=13)
lines(1:ABC_time_steps,dist_0.50quantile_500,type="b",lwd=3,col="black",pch=6)
lines(1:ABC_time_steps,dist_0.75quantile_500,type="b",lwd=3,col="pink",pch=5)
lines(1:ABC_time_steps,dist_0.85quantile_500,type="b",lwd=3,col="orange",pch=10)

text(3,.85,"At Monte Carlo sample size of N=500",font=2,cex=1.5)


#N=1000
o<-par(mar=c(0,4,0,2))

dist_0.05quantile_1000<-rep(NA,ABC_time_steps)
dist_0.15quantile_1000<-rep(NA,ABC_time_steps)
dist_0.25quantile_1000<-rep(NA,ABC_time_steps)
dist_0.50quantile_1000<-rep(NA,ABC_time_steps)
dist_0.75quantile_1000<-rep(NA,ABC_time_steps)
dist_0.85quantile_1000<-rep(NA,ABC_time_steps)

for(t in 1:ABC_time_steps){
   dist_0.05quantile_1000[t]<- quantile(distances_1000[[t]],probs =0.05)
   dist_0.15quantile_1000[t]<- quantile(distances_1000[[t]],probs =0.15)
   dist_0.25quantile_1000[t]<- quantile(distances_1000[[t]],probs =0.25)
   dist_0.50quantile_1000[t]<- quantile(distances_1000[[t]],probs =0.5)
   dist_0.75quantile_1000[t]<- quantile(distances_1000[[t]],probs =0.75)
   dist_0.85quantile_1000[t]<- quantile(distances_1000[[t]],probs =0.85) 
}


plot(1:ABC_time_steps,dist_0.05quantile_1000,type="b",lwd=3,col="blue",
     ylim=c(0.55,.82),xlab="",pch=3,
     main="", ylab="", xaxt ="n")
#axis(side=1, at = 1:ABC_time_steps, labels = paste(1:ABC_time_steps))
lines(1:ABC_time_steps,dist_0.15quantile_1000,type="b",lwd=3,col="green",pch=4)
lines(1:ABC_time_steps,dist_0.25quantile_1000,type="b",lwd=3,col="red",pch=13)
lines(1:ABC_time_steps,dist_0.50quantile_1000,type="b",lwd=3,col="black",pch=6)
lines(1:ABC_time_steps,dist_0.75quantile_1000,type="b",lwd=3,col="pink",pch=5)
lines(1:ABC_time_steps,dist_0.85quantile_1000,type="b",lwd=3,col="orange",pch=10)

legend(x=4,y=.8,legend=c("5th percentile","15th percentile",
                         "25th percentile","50th percentile","75th percentile","85th percentile"),
        col=c("blue","green","red","black","pink","orange"),pch=c(3,4,13,6,5,10),
       cex=1.4,box.lwd = 2,ncol =3,title="Percentiles:",bty="n")

mtext(text = "Quantiles of ABC-SMC distances",
      side = 2,
      line = 2.5,cex=1,font=2)


text(3,.82,"At Monte Carlo sample size of N=1000",font=2,cex=1.5)

#n=1500
par(mar=c(4,4,0,2))

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
     ylim=c(0.55,.82),xlab="",pch=3
     ,main="", ylab="", xaxt ="n")
axis(side=1, at = 1:ABC_time_steps, labels = paste(1:ABC_time_steps))
lines(1:ABC_time_steps,dist_0.15quantile_1500,type="b",lwd=3,col="green",pch=4)
lines(1:ABC_time_steps,dist_0.25quantile_1500,type="b",lwd=3,col="red",pch=13)
lines(1:ABC_time_steps,dist_0.50quantile_1500,type="b",lwd=3,col="black",pch=6)
lines(1:ABC_time_steps,dist_0.75quantile_1500,type="b",lwd=3,col="pink",pch=5)
lines(1:ABC_time_steps,dist_0.85quantile_1500,type="b",lwd=3,col="orange",pch=10)

mtext(text = "ABC-SMC time steps",
      side = 1,
      line = 2.5,cex=1,font=2)

text(3,.82,"At Monte Carlo sample size of N=1500",font=2,cex=1.5)

#All model parameters
N<- 500 # Monte Carlo Sample Size
x <- seq(from = -10, to = 7, length.out = 256)#range of prior distribution

par(mfrow=c(5,5), mar=c(4,2,1,0),font=2)
plot(NULL ,xaxt='n',yaxt='n',bty='n', xlab="", ylab="",xlim=0:1, ylim=0:1)
plot_colors <- c("blue","red","black")
legend(x=-0.1,y=1.2,c(paste0("First Prior 
(N=",N,")"),"Intermediate 
Priors" ,"Final Posterior"),
        col=c("blue","red","black"),bty="n",cex=1.05,box.lwd = .8,fill=c("blue","red","black"),horiz=FALSE)


for (k in 1:23) {
  plot(x, density_post_500[[1]][k, ], type="l", ylim=c(0,5),
       xlab=parameter_labels[k],ylab="", col="blue",cex.lab=1.2,lwd=2.5,las=1) 
  
    
  for (j in 2:ABC_time_steps) {
    lines(x,density_post_500[[j]][k, ], yaxt= "n",col="red",lwd=1,pch=4,ann=FALSE,yaxt="n")
  }
    
     lines(x, density_post_500[[ABC_iterations]][k, ], col="black", lwd=1,ann=FALSE,yaxt="n")
   
}

#All model parameters
N<- 1000 # Monte Carlo Sample Size
x <- seq(from = -10, to = 7, length.out = 256)#range of prior distribution

par(mfrow=c(5,5), mar=c(4,2,1,0),font=2)
plot(NULL ,xaxt='n',yaxt='n',bty='n', xlab="", ylab="",xlim=0:1, ylim=0:1)
plot_colors <- c("blue","red","black")
legend(x=-0.1,y=1.2,c(paste0("First Prior 
(N=",N,")"),"Intermediate 
Priors" ,"Final Posterior"),
        col=c("blue","red","black"),bty="n",cex=1.05,box.lwd = .8,fill=c("blue","red","black"),horiz=FALSE)


for (k in 1:23) {
  plot(x, density_post_1000[[1]][k, ], type="l", ylim=c(0,5),
       xlab=parameter_labels[k],ylab="", col="blue",cex.lab=1.2,lwd=2.5,las=1) 
  
    
  for (j in 2:ABC_time_steps) {
    lines(x,density_post_1000[[j]][k, ], yaxt= "n",col="red",lwd=1,pch=4,ann=FALSE,yaxt="n")
  }
    
     lines(x, density_post_1000[[ABC_iterations]][k, ], col="black", lwd=1,ann=FALSE,yaxt="n")
   
}

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

setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro/Updated ABC_results_N_500")


#Estimating the model parameters from the unadjusted posterior samples and their corresponding credible intervals

  #computing posterior estimaates and their credible intervals (ETI)
  post_data_500<-   Final_posterior_500
  ci_eti_500<- bayestestR::ci(as.data.frame(post_data_500), ci = 0.95,method = "ETI") 
  ci_eti_500<- exp(as.data.frame(ci_eti_500[,3:4]))
   # estimating the  posterior  mean estimate
  Posterior_estimates_unadj_500<- apply(exp(Final_posterior_500),2,mean)
    

Posterior_est_unadj_data_500<- data.frame(Posterior_estimate_unadj=Posterior_estimates_unadj_500,
                                     Cred_Int_lower95_unadj= ci_eti_500[,1],
                                     Cred_Int_upper95_unadj= ci_eti_500[,2])


rownames(Posterior_est_unadj_data_500)<- parameter_labels
#Print unadjusted posterior estimates with 95% credible intervals
Posterior_est_unadj_data_500

#Save the unadjusted posterior estimate results

write.csv(Posterior_est_unadj_data_500, "Posterior_estimates_unadj_data_500.csv")

setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro/Updated ABC_results_N_1000")


#Estimating the model parameters from the unadjusted posterior samples and their corresponding credible intervals

  #computing posterior estimaates and their credible intervals (ETI)
  post_data_1000<-   Final_posterior_1000
  ci_eti_1000<- bayestestR::ci(as.data.frame(post_data_1000), ci = 0.95,method = "ETI") 
  ci_eti_1000<- exp(as.data.frame(ci_eti_1000[,3:4]))
   # estimating the  posterior  mean estimate
  Posterior_estimates_unadj_1000<- apply(exp(Final_posterior_1000),2,mean)
    

Posterior_est_unadj_data_1000<- data.frame(Posterior_estimate_unadj=Posterior_estimates_unadj_1000,
                                     Cred_Int_lower95_unadj= ci_eti_1000[,1],
                                     Cred_Int_upper95_unadj= ci_eti_1000[,2])


rownames(Posterior_est_unadj_data_1000)<- parameter_labels
#Print unadjusted posterior estimates with 95% credible intervals
Posterior_est_unadj_data_1000

#Save the unadjusted posterior estimate results

write.csv(Posterior_est_unadj_data_1000, "Posterior_estimates_unadj_data_1000.csv")

setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro/Updated ABC_results_N_1500")


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

#unadjusted posterior estimates
round(Posterior_est_unadj_data_1500, 4)

# Posterior estimates comparison at N=500, 1000 & 1500
Estimates_comparision<-  data.frame(Estimates_N500=round(Posterior_est_unadj_data_500[, 1], 4),
                                    Estimates_N1000=round(Posterior_est_unadj_data_1000[, 1], 4),
                                    Estimates_N1500=round(Posterior_est_unadj_data_1500[, 1], 4)
                                   )


rownames(Estimates_comparision)<- parameter_labels
Estimates_comparision

#options(warn=-1)

plot(log(Posterior_est_unadj_data_500[, 1]),type="p",col="blue",lwd=1,xaxt = "n",
     ylab="Posterior mean with 95% credible interval  (log scale)",xlab="Model parameters",
     ylim=c(-10,20),pch=19,cex.lab=1, cex.axis=1)

arrows(x0=1:23, y0=log(Posterior_est_unadj_data_500[, 2]), x1=1:23, y1=log(Posterior_est_unadj_data_500[, 3]), 
       code=3,  angle=90, length=0.05,col ="blue",lwd=1)
 axis(1, at=1:23, labels=parameter_labels,font=3, hadj = .9,cex.axis=.8)


legend(x=2,y=10,legend=c("N=500", "N=1000", "N=1500"),
        col=c("blue","green","red"),bty="n",
       cex=.85,box.lwd = .8,fill=c("blue","green","red"),
       ncol =1,title="Unadjusted posterior estimates at
Monte Carlo sample sizes:")

lines(log(Posterior_est_unadj_data_1000[, 1]),type="p",col="green",lwd=1,xaxt = "n",ylab="",xlab="",
     ylim=c(-9,10),pch=19,cex.lab=1, cex.axis=1)

arrows(x0=1:23, y0=log(Posterior_est_unadj_data_1000[, 2]), x1=1:23, y1=log(Posterior_est_unadj_data_1000[, 3]), 
       code=3,  angle=90, length=0.05,col ="green",lwd=1)


lines(log(Posterior_est_unadj_data_1500[, 1]),type="p",col="red",lwd=1,xaxt = "n",ylab="",xlab="",
     ylim=c(-9,10),pch=19,cex.lab=1, cex.axis=1)

arrows(x0=1:23, y0=log(Posterior_est_unadj_data_1500[, 2]), x1=1:23, y1=log(Posterior_est_unadj_data_1500[, 3]), 
       code=3,  angle=90, length=0.05,col ="red",lwd=1)

# calculate position of inset
plotdim <- par("plt")
xleft    = plotdim[2] - (plotdim[2] - plotdim[1]) * .35
xright   = plotdim[2]  #
ybottom  = plotdim[4] - (plotdim[4] - plotdim[3]) * 0.35  #
ytop     = plotdim[4]  #


# set position for inset
par(
  fig = c(xleft, xright, ybottom, ytop)
  , mar=c(0,0,0,0)
  , new=TRUE
  )



# add inset
plot(c(time_CPU_500,time_CPU_1000,time_CPU_1500),type = "b", col = c("black"), xlab = "N", ylab = "CPU time",
    xaxt = "n",lwd=3, yaxt="n" ,ylim=c(50000, 140000))

 axis(1, at=1:3, labels=c("500","1000","1500"),font=1, hadj = .7,cex.axis=.8)

# Specify the y-axis ticks
axis(2, at = c(50000, 90000, 130000), labels = c("5K", "90K", "130K"),cex=.5)

mtext(text = "N",
      side = 1,
      line = 2.5,cex=1,font=2)

mtext(text = "CPU time (in secs)",
      side = 2,
      line = 2.5,cex=1,font=2)

compute_bootstrap_se <- function(data, num_iterations = 1000) {
  n<- dim(data)[1]
  bootstrap_samples <- numeric(num_iterations)
  # Perform bootstrap resampling
  for (i in 1:num_iterations) {
    # Sample with replacement from the data
    bootstrap_sample <- sample(data, replace = TRUE)
      n<- length(bootstrap_sample)
    # Compute the mean of the bootstrap sample
    bootstrap_samples[i]<- sd(bootstrap_sample)/sqrt(n)
  }
  
  # Compute bootstrap standard error
  bootstrap_se <- mean(bootstrap_samples)
  Pred_int95<- quantile(bootstrap_samples, c(0.025,0.975))
  output<- data.frame(bootstrap_se=bootstrap_se,
                     Pred_int95_Lower=Pred_int95[1],
                     Pred_int95_Upper=Pred_int95[2]) 
    
   rownames(output)<-NULL 
  #Return the standard error estimate and its 95% prediction 
  return(output)
}



#Estimated standard errors at N=500
Std_errors_N500<- NULL


for(para in 1:dim(Final_posterior_500)[2]){
    Std_errors_N500[[para]]<- compute_bootstrap_se(Final_posterior_500[ ,para],num_iterations = 1000)
    
}


#Estimated standard errors at N=500
Std_errors_N500_df<- do.call("rbind",Std_errors_N500)
rownames(Std_errors_N500_df)<- parameter_labels
Std_errors_N500_df


cat("Pooled standard error at N=500:")
apply(Std_errors_N500_df,2,mean)

#Estimated standard errors at N=1000
Std_errors_N1000<- NULL


for(para in 1:dim(Final_posterior_1000)[2]){
    Std_errors_N1000[[para]]<- compute_bootstrap_se(Final_posterior_1000[ ,para],num_iterations = 1000)
    
}


#Estimated standard errors at N=1000
Std_errors_N1000_df<- do.call("rbind",Std_errors_N1000)
rownames(Std_errors_N1000_df)<- parameter_labels
Std_errors_N1000_df

cat("Pooled standard error at N=1000:")
apply(Std_errors_N1000_df,2,mean)

#Estimated standard errors at N=1000
Std_errors_N1500<- NULL


for(para in 1:dim(Final_posterior_1500)[2]){
    Std_errors_N1500[[para]]<- compute_bootstrap_se(Final_posterior_1500[ ,para],num_iterations = 1000)
    
}


#Estimated standard errors at N=1000
Std_errors_N1500_df<- do.call("rbind",Std_errors_N1500)
rownames(Std_errors_N1500_df)<- parameter_labels
Std_errors_N1500_df

cat("Pooled standard error at N=1500:")
apply(Std_errors_N1500_df,2,mean)

###### Plotting the bootstrap standard errors and their 95% corresponding prediction intervals########
#options(warn=-1)
#Posterior at N=500
plot(Std_errors_N500_df[, 1],type="p",col="blue",lwd=1,xaxt = "n",
     ylab="Bootstrap standard error of ABC posterior mean",xlab="Model parameters",
     ylim=c(0,.1),pch=19,cex.lab=1, cex.axis=1)

lines(Std_errors_N500_df[, 1],type="l",col="blue",lwd=2,xaxt = "n",
     ylab="Bootstrap standard error of the posterior mean",xlab="Model parameters",
     ylim=c(0,.1),pch=19,cex.lab=1, cex.axis=1)

arrows(x0=1:23, y0=Std_errors_N500_df[, 2], x1=1:23, y1=Std_errors_N500_df[, 3], 
       code=3,  angle=90, length=0.05,col ="blue",lwd=1)
 axis(1, at=1:23, labels=parameter_labels,font=3, hadj = .9,cex.axis=.8)


legend(x=16,y=.105,legend=c("N=500", "N=1000", "N=1500"),
        col=c("blue","green","red"),bty="n",
       cex=.85,box.lwd = .8,fill=c("blue","green","red"),
       ncol =1,title="Monte Carlo sample sizes:")


#Posterior at N=1000
lines(Std_errors_N1000_df[, 1],type="p",col="green",lwd=1,xaxt = "n",ylab="",xlab="",
     ylim=c(0,.08),pch=19,cex.lab=1, cex.axis=1)

lines(Std_errors_N1000_df[, 1],type="l",col="green",lwd=2,xaxt = "n",ylab="",xlab="",
     ylim=c(0,.08),pch=19,cex.lab=1, cex.axis=1)

arrows(x0=1:23, y0=Std_errors_N1000_df[, 2], x1=1:23, y1=Std_errors_N1000_df[, 3], 
       code=3,  angle=90, length=0.05,col ="green",lwd=1)


#Posterior at N=1500
lines(Std_errors_N1500_df[, 1],type="p",col="red",lwd=1,xaxt = "n",ylab="",xlab="",
     ylim=c(0,0.08),pch=19,cex.lab=1, cex.axis=1)

lines(Std_errors_N1500_df[, 1],type="l",col="red",lwd=2,xaxt = "n",ylab="",xlab="",
     ylim=c(0,0.08),pch=19,cex.lab=1, cex.axis=1)

arrows(x0=1:23, y0=Std_errors_N1500_df[, 2], x1=1:23, y1=Std_errors_N1500_df[, 3], 
       code=3,  angle=90, length=0.05,col ="red",lwd=1)




text(7.5,.1, "Pooled SE = 0.0355 (95% PI: 0.0204, 0.0341)", col="blue",font=1)
text(7.5,.095, "Pooled SE = 0.0274 (95% PI: 0.0198, 0.0481)", col="green",font=1)
text(7.5,.09, "Pooled SE = 0.0255 (95% PI: 0.0190, 0.0317)", col="red",font=1)

#Setting working directory
setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro")


#Function for posterior adjustment using Weighted Ridge and Lasso regression, respectively
source("Post-Ridge-reg-adj-L2-script.R")#ABC Post processing-Ridge/l2 correction
source("Post-Lasso-reg-adj-L1-script.R")#ABC Post processing-Lasso/L1 correction

#set.seed(12) #for reproducibility
set.seed(12) #remove if necessary
#Ridge or L2 posterior correction
Reg_adjust_L2<- function(i) Post_Ridge_reg_adj(post_distn=Final_posterior_1500,
                                    summary_obs= S_obs)

Output_adj_Ridge<- mclapply(1, Reg_adjust_L2,mc.cores=numCores)

Output_adj_Ridge[[1]]$Posterior_mean_output

#Lasso or L1 posterior correction
#set.seed(12) #for reproducibility
set.seed(12) #remove if necessary
Reg_adjust_L1<- function(i) Post_Lasso_reg_adj(post_distn=Final_posterior_1500,
                                    summary_obs= S_obs)

Output_adj_Lasso<- mclapply(1,Reg_adjust_L1,mc.cores=numCores)

Output_adj_Lasso[[1]]$Posterior_mean_output

#Output_adj_Ridge[[1]]$Adjusted_posterior_dist

#Based on Ridge regression
X_matrix<- as.data.frame(Output_adj_Ridge[[1]]$X_Design_matrix)

for (i in 1:dim(X_matrix)[2]) names(X_matrix)[i]<- paste0("S_",i)

chart.Correlation(X_matrix, histogram=TRUE, pch=19)

#Based on Lasso regression
X_matrix<- as.data.frame(Output_adj_Lasso[[1]]$X_Design_matrix)

for (i in 1:dim(X_matrix)[2]) names(X_matrix)[i]<- paste0("S_",i)

chart.Correlation(X_matrix, histogram=TRUE, pch=19)

#Setting working directory
setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro/Regression_adjust_SimulationModel_new")




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

setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro/Regression_adjust_SimulationModel_new")

#Saving Ridge and Lasso adjusted posterior distributions 
write.csv(Output_adj_Ridge[[1]]$Adjusted_posterior_dist,
                        paste0("Posterior_Ridge_adj_distn_log_",draws,".csv"))


write.csv(Output_adj_Lasso[[1]]$Adjusted_posterior_dist,
                        paste0("Posterior_Lasso_adj_distn_log_",draws,".csv"))

#Importing previously saved ABC Ridge post-processing results
#Setting working directory
setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro/Regression_adjust_SimulationModel_new")


Post_Ridge_Regression_output<- read.csv("Post_Ridge_Regression_output.csv")
Post_Ridge_Regression_output<- Post_Ridge_Regression_output[,-1]
rownames(Post_Ridge_Regression_output)<- parameter_labels


Post_Lasso_Regression_output<- read.csv("Post_Lasso_Regression_output.csv")
Post_Lasso_Regression_output<- Post_Lasso_Regression_output[,-1]
rownames(Post_Lasso_Regression_output)<- parameter_labels

#unadjusted posterior estimates
round(Posterior_est_unadj_data_1500, 6)


round(Post_Ridge_Regression_output[, 1:3],6)

round(Post_Lasso_Regression_output[, 1:3],6)

plot(log(Posterior_est_unadj_data_1500[,3]-Posterior_est_unadj_data_1500[,2]), col="black",type="b",ylim=c(-11,5))
lines(log(Post_Ridge_Regression_output[,3]-Post_Ridge_Regression_output[,2]), col="green",type="b")
lines(log(Post_Lasso_Regression_output[,3]-Post_Lasso_Regression_output[,2]), col="orange",type="b")



#Setting working directory
setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro/Regression_adjust_SimulationModel_new")


#Importing adjusted and unadjusted posterior distributions
Posterior_Ridge_adj_distn_log<- read.csv("Posterior_Ridge_adj_distn_log_1500.csv")
Posterior_Ridge_adj_distn_log<- Posterior_Ridge_adj_distn_log[,-1]



Posterior_Lasso_adj_distn_log<- read.csv("Posterior_Lasso_adj_distn_log_1500.csv")
Posterior_Lasso_adj_distn_log<- Posterior_Lasso_adj_distn_log[,-1]



Posterior_Unadj__distn_log<- read.csv("PosteriorFinal_1500.csv")
Posterior_Unadj__distn_log<- Posterior_Unadj__distn_log[,-1]



Unadj_posterior<- exp(Posterior_Unadj__distn_log)
Ridge_posterior<- exp(Posterior_Ridge_adj_distn_log)
Lasso_posterior<- exp(Posterior_Lasso_adj_distn_log)


Posterior_combined_df<- rbind(Unadj_posterior,Ridge_posterior,Lasso_posterior)

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


text(-0.02, .1, "Multivariate homogeneity test of ABC posterior variances: p-value=0.981", col="black", font=2)
mtext(LETTERS[1], adj=0, line=1,font=2)



boxplot(multivariate_dispersion_mod,col=my_cols, xlab="ABC methods",ylab="Distance to posterior centroid")
mtext(paste(LETTERS[2]), adj=0, line=1,font=2)

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

 if(k==7) legend("topleft",c("First Prior","Intermediate Priors" ,"Unadjusted Posterior","Ridge Adjustment","Lasso Adjustment"),
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

 if(k==19) legend("topleft",c("First Prior","Intermediate Priors" ,"Unadjusted Posterior","Ridge Adjustment","Lasso Adjustment"),
         col=c("blue","red","black","green","orange"),bty="n",cex=1,box.lwd = 2,fill=c("blue","red","black","green","orange"))

    
  for (j in 2:10) {
    lines(x,density_post_1500[[j]][k, ], yaxt= "n",col="red",lwd=1,pch=4,ann=FALSE,yaxt="n")
  }
    
     lines(x, density_post_1500[[ABC_iterations]][k, ], col="black", lwd=1.5)
     lines(x,density_post_Ridge_adjusted[k, ], col="green",lty=3,lwd=2)
      lines(x,density_post_Lasso_adjusted[k, ], col="orange",lty=3,lwd=2)
   
    }

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

set.seed(100820249)

rep<-1 #Number of replicates/repetition
sim_func_unadj<- function(i) SimGroup_tauleap(theta1=log(Parameter_estimates_unadj),
            fish_sex=fishSex,fish_type=Fish_stock,strain=Strain,fish_size=fishSize,error=0.002)


output_sim_unadj1<- mclapply(1:rep, sim_func_unadj,mc.cores=numCores)

#Simulations based on the Ridge adjusted posterior#########################

set.seed(100820249)


rep<-1 #Number of replicates/repetition
sim_func_ridge<- function(i) SimGroup_tauleap(theta1=log(Parameter_estimates_ridge),
            fish_sex=fishSex,fish_type=Fish_stock,strain=Strain,fish_size=fishSize,error=0.002)

output_sim_ridge1<- mclapply(1:rep, sim_func_ridge,mc.cores=numCores)

#Simulations based on the Lasso adjusted posterior
set.seed(100820249)



rep<-1 #Number of replicates/repetition
sim_func_lasso<- function(i) SimGroup_tauleap(theta1=log(Parameter_estimates_lasso),
            fish_sex=fishSex,fish_type=Fish_stock,strain=Strain,fish_size=fishSize,error=0.002)

output_sim_lasso1<- mclapply(1:rep, sim_func_lasso,mc.cores=numCores)

Sim_data_unadj_reg<- NULL; Sim_data_ridge_reg<- NULL;Sim_data_lasso_reg<- NULL

for(i in 1:rep){
    #Simulations based on unadjusted posterior estimates
    Sim_data_unadj_reg[[i]]<- Parasite_data_summary(data_all_fish=output_sim_unadj1[[i]], 
                    observed_times= 1:9,Total_fish=sum(unlist(numF)),parasite_fish=parasite_fish,sim=TRUE,
                     fishSex=fishSex,fishID=fishID,fishSize=fishSize,Strain=Strain,Fish_stock=Fish_stock,
                    Group="Unadjusted ABC")
    
    Sim_data_ridge_reg[[i]]<- Parasite_data_summary(data_all_fish=output_sim_ridge1[[i]], 
                    observed_times= 1:9,Total_fish=sum(unlist(numF)),parasite_fish=parasite_fish,sim=TRUE,
                     fishSex=fishSex,fishID=fishID,fishSize=fishSize,Strain=Strain,Fish_stock=Fish_stock,
                    Group="Ridge adjustment")
    
    Sim_data_lasso_reg[[i]]<- Parasite_data_summary(data_all_fish=output_sim_lasso1[[i]], 
                    observed_times= 1:9,Total_fish=sum(unlist(numF)),parasite_fish=parasite_fish,sim=TRUE,
                     fishSex=fishSex,fishID=fishID,fishSize=fishSize,Strain=Strain,Fish_stock=Fish_stock,
                    Group="Lasso adjustment")
    
    }


#Observed data
Observed_data_reg<- Parasite_data_summary(data_all_fish=pop_obs, observed_times= 1:9,
                    Total_fish=sum(unlist(numF)),parasite_fish=parasite_fish,sim=FALSE,
                     fishSex=fishSex,fishID=fishID,
                fishSize=fishSize,Strain=Strain,Fish_stock=Fish_stock,Group="Observed data")

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



transdata_reg$Group<- relevel(transdata_reg$Group,ref="Observed data")
transdata_reg$Sex<- relevel(transdata_reg$Sex,ref="Female")
transdata_reg$Parasite_strain<- relevel(transdata_reg$Parasite_strain,ref="Gt")
transdata_reg$Fish_stock<- relevel(transdata_reg$Fish_stock,ref="UA")

transdata_reg$Group<- factor(transdata_reg$Group,levels=c("Observed data",
                        "Unadjusted ABC","Ridge adjustment","Lasso adjustment"))

levels(transdata_reg$Group)

transdata_reg$Day<- as.factor(transdata_reg$Time)
levels(transdata_reg$Day)<- c("Day 1","Day 3","Day 5", "Day 7", "Day 9",
                              "Day 11","Day 13","Day 15","Day 17")


#head(transdata_reg)

p1<- transdata_reg %>%
  mutate(type=fct_reorder(as.factor(Time),Parasites),
        Group=fct_reorder(as.factor(Group),Parasites)) %>%
     ggplot()+geom_smooth(aes(x=Time, y=Parasites,color=Group,group=Group),method = "loess")+
labs(y="Parasite mean intensity") +
theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),

        panel.background = element_blank())


#p1+xlim(1,17)+facet_wrap(~Parasite_strain+Fish_stock,  ncol=3)+
 # guides(color = guide_legend(title = "Comparable groups"))

my_cols<- c("blue","black","green","orange")

# Visualize: Specify the comparisons you want
my_comparisons <- list( c("Observed data", "Unadjusted ABC"),
                       c("Observed data", "Ridge adjustment"),
             c("Observed data", "Lasso adjustment"))


#subset(transdata_reg, Time==1)
bxp1<- ggboxplot(transdata_reg, x = "Group", y = "Parasites",
          fill = "Group")+scale_fill_manual(values=my_cols) #you can change "fill" to "color"

bxp1<- bxp1+stat_compare_means(comparisons = my_comparisons, hide.ns = F,
    label = c("p.signif"),palette = c('red'),label.x.npc ='center',#ref.group =".all."
  label.y.npc ='top',# method = "wilcox",  # Use Wilcoxon test
                     p.adjust.method = "BH",  # Benjamini-Hochberg adjustment
                     alpha = 0.05)+ # Add pairwise comparisons p-value
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
group_ABC.methods<- factor(group_ABC.methods,levels=c("Observed data",
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

#Posterior distribution on original scale
Posterior_Lasso_adj_distn<- exp(Posterior_Lasso_adj_distn_log)

# Function to perform Region of Practical Equivalence (ROPE) and Highest Density Interval (HDI)
require(bayestestR)
ROPE_Cred_Int<- function(theta_distn_diff,parameter_labels, ci_percent=0.89){
    if(is.list(theta_distn_diff)==FALSE){
        sigma_d<- sd(theta_distn_diff)#standard deviation of differenced posterior samples
        output<- bayestestR::equivalence_test(theta_distn_diff, 
                #ROPE range is based on recommendation by Norman et al (2003)
                range =c(-.5*sigma_d,.5*sigma_d), ci = ci_percent,ci_method = "HDI")
    final_output<- cbind(parameter_labels,output);names(final_output)[1]<- "Parameter"
    return(final_output)
    }
    #theta_distn_diff= a list of posterior samples of differences of parameters of interest
    output<-list() #save ROPE+HDI results
    for(i in seq_along(parameter_labels)){
     sigma_d<- sd(theta_distn_diff[[i]])#standard deviation of differenced posterior samples
    output[[i]]<- bayestestR::equivalence_test(theta_distn_diff[[i]], 
                #ROPE range is based on recommendation by Norman et al (2003)
                range =c(-.5*sigma_d,.5*sigma_d), ci = ci_percent,ci_method = "HDI")
        }
    
     #Function returns ROPE interval, ROPE Percentage, ROPE equivalence decision
      # and the corresponding HDI 
    final_output<- do.call("rbind",output)
    final_output<- cbind(parameter_labels,final_output)
    names(final_output)[1]<- "Parameters"
    return(final_output)
}

#1a. Comparing significant difference between birth rates of young parasites across the three parasite strains
#Extracting corresponding adjusted posterior samples
b1_young_posterior<-  Posterior_Lasso_adj_distn[, c(1,3,5)]
b1_young_diff<- NULL

parameter_labels<- c("b_11-b_21","b_11-b_31","b_21-b_31")
dim<- length(parameter_labels)


g<-0#initialise number of groups
for(i in 1:(dim-1)){
    for(j in (dim-1):dim){
        if(i!=j){
            print(paste("i=",i,"","j=",j))
            g<- g+1
            #Compute posterior of the difference under H0
            b1_young_diff[[g]]<-b1_young_posterior[,i]-b1_young_posterior[,j]
        }
    }
      
 }       
 
#Performing ROPE+HDI test
ROPE_Cred_Int(theta_distn_diff=b1_young_diff,parameter_labels=parameter_labels, ci_percent=0.89)

#1b. Comparing significant difference between birth rates of old parasites across the three parasite strains
#Extracting corresponding adjusted posterior samples
b2_old_posterior<- Posterior_Lasso_adj_distn[, c(2,4,6)]
b2_old_diff<- NULL

parameter_labels<- c("b_12-b_22","b_12-b_32","b_22-b_32")
dim<- length(parameter_labels)


g<-0#initialise number of groups
for(i in 1:(dim-1)){
    for(j in (dim-1):dim){
        if(i!=j){
            print(paste("i=",i,"","j=",j))
            g<- g+1
            #Compute posterior of the difference under H0
            b2_old_diff[[g]]<-b2_old_posterior[,i]-b2_old_posterior[,j]
        }
    }
      
 }       
 
#Performing ROPE+HDI test
ROPE_Cred_Int(theta_distn_diff=b2_old_diff,parameter_labels=parameter_labels, ci_percent=0.89)

#1c. Comparing significant difference between birth rates of all parasite strains per age group 
#(between young and old parasites)
#Extracting corresponding adjusted posterior samples
b_posterior<- Posterior_Lasso_adj_distn[, 1:6]
b_diff<- NULL

parameter_labels<- c("b_11-b_12","b_21-b_22","b_31-b_32")
#Compute posterior of the difference under H0
b_diff[[1]]<- b_posterior[,1]-b_posterior[,2]
b_diff[[2]]<- b_posterior[,3]-b_posterior[,4]
b_diff[[3]]<- b_posterior[,5]-b_posterior[,6]
#Performing ROPE+HDI test
ROPE_Cred_Int(theta_distn_diff=b_diff,parameter_labels=parameter_labels, ci_percent=0.89)

#2a. Comparing significant difference between death rates in the absence of host immune response 
#across the three parasite strains

#Extracting corresponding adjusted posterior samples
d_nr_posterior<- Posterior_Lasso_adj_distn[, c(7,9,11)]
d_nr_diff<- NULL

parameter_labels<- c("d_11-d_21","d_11-d_31","d_21-d_31")
dim<- length(parameter_labels)


g<-0#initialise number of groups
for(i in 1:(dim-1)){
    for(j in (dim-1):dim){
        if(i!=j){
            print(paste("i=",i,"","j=",j))
            g<- g+1
            #Compute posterior of the difference under H0
            d_nr_diff[[g]]<-d_nr_posterior[,i]-d_nr_posterior[,j]
        }
    }
      
 }       
 
#Performing ROPE+HDI test
ROPE_Cred_Int(theta_distn_diff=d_nr_diff,parameter_labels=parameter_labels, ci_percent=0.89)

#2b. Comparing significant difference between death rates in the presence of host immune response across
#the three parasite strains

#Extracting corresponding adjusted posterior samples
d_r_posterior<- Posterior_Lasso_adj_distn[, c(8,10,12)]
d_r_diff<- NULL

parameter_labels<- c("d_12-d_22","d_12-d_32","d_22-d_32")
dim<- length(parameter_labels)


g<-0#initialise number of groups
for(i in 1:(dim-1)){
    for(j in (dim-1):dim){
        if(i!=j){
            print(paste("i=",i,"","j=",j))
            g<- g+1
            #Compute posterior of the difference under H0
            d_r_diff[[g]]<-d_r_posterior[,i]-d_r_posterior[,j]
        }
    }
      
 }       
 
#Performing ROPE+HDI test
ROPE_Cred_Int(theta_distn_diff=d_r_diff,parameter_labels=parameter_labels, ci_percent=0.89)

#2c. Comparing death rates of all parasite strains in the presence or absence of host immune response

#Extracting corresponding adjusted posterior samples
d_posterior<- Posterior_Lasso_adj_distn[, 7:12]
d_diff<- NULL

parameter_labels<- c("d_11-d_12","d_21-d_22","d_31-d_32")
#Compute posterior of the difference under H0
d_diff[[1]]<- d_posterior[,1]-d_posterior[,2]
d_diff[[2]]<- d_posterior[,3]-d_posterior[,4]
d_diff[[3]]<- d_posterior[,5]-d_posterior[,6]
#Performing ROPE+HDI test
ROPE_Cred_Int(theta_distn_diff=d_diff,parameter_labels=parameter_labels, ci_percent=0.89)

#3. Comparing the movement rate adjustment of all parasite strains in relation to microhabitat preference

#Extracting corresponding adjusted posterior samples
e_posterior<- Posterior_Lasso_adj_distn[, c(20,21,22)]
e_diff<- NULL

parameter_labels<-c("epsilon_1-epsilon_2","epsilon_1-epsilon_3","epsilon_2-epsilon_3")
dim<- length(parameter_labels)


g<-0#initialise number of groups
for(i in 1:(dim-1)){
    for(j in (dim-1):dim){
        if(i!=j){
            print(paste("i=",i,"","j=",j))
            g<- g+1
            #Compute posterior of the difference under H0
            e_diff[[g]]<-e_posterior[,i]-e_posterior[,j]
        }
    }
      
 }       
 
#Performing ROPE+HDI test
ROPE_Cred_Int(theta_distn_diff=e_diff,parameter_labels=parameter_labels, ci_percent=0.89)

#4. Comparing the immune response rate adjustments in relation to fish type and host sex as well as 
#effect of sex on host mortality

#Extracting corresponding adjusted posterior samples
r_posterior<- Posterior_Lasso_adj_distn[, c(15,16,17)]
r_diff<- NULL

parameter_labels<-c("r_1-r_2","r_1-r_3","r_2-r_3")
dim<- length(parameter_labels)


g<-0#initialise number of groups
for(i in 1:(dim-1)){
    for(j in (dim-1):dim){
        if(i!=j){
            print(paste("i=",i,"","j=",j))
            g<- g+1
            #Compute posterior of the difference under H0
            r_diff[[g]]<-r_posterior[,i]-r_posterior[,j]
        }
    }
      
 }       
 
#Performing ROPE+HDI test
ROPE_Cred_Int(theta_distn_diff=r_diff,parameter_labels=parameter_labels, ci_percent=0.89)

#unadjusted posterior estimates
round(Posterior_est_unadj_data_1500, 6)


#5. Comparing mortality rate of male fish and female fish

s1_posterior<- Posterior_Lasso_adj_distn[, 19]

parameter_labels<-c("s_1")
#Performing ROPE+HDI test
ROPE_Cred_Int(theta_distn_diff=s1_posterior,parameter_labels=parameter_labels, ci_percent=0.89)
