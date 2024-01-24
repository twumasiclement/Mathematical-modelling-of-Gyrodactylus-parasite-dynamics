######## The Main Script for Performing ABC fitting of the complex stochastic simulation model###


#Loading packages
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
library(kedd)
library("mvtnorm")
library(Rcpp)# for C+ matrix multiplication
library(tmvtnorm)
library("ie2misc")#for computing mean absolute deviation
#options(repr.plot.width=8, repr.plot.height=8,repr.plot.res = 300) #Setting plot size

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
summaries_obs[[1]]# for group 1 or Gt3-0S group




#Initial simulation inputs for A (parasite numbers) and  B (immune status)
A0 <- matrix(0, 4, 2)  
A0[1, 1] <- 2   #Intial parasites at the tail
B0 <- rep(1, 4)  #initial immune response at 4 body regions (no response)

#Transition matrix
J<- matrix(c(0,    1,    0,     0, 
             1/2,   0,    1/2,   0,
             0,    1/2,    0,   1/2,
             0,     0,     1,    0), 4, 4, byrow=TRUE)


#### NB: # Summary statistics weights are computed based on inverse of MAD as recommended by Prangle
# initial weights estimate
#dimS<-17 #length of summary statistics for ABC
#n0 <- 100 #change to n0 <- 1000 or 500
#saving summary statistic for each group simulation realisation
#for computing intial weights for ABC fitting
#SummaryStats_sim <- NULL;SummaryStats_sim_combined<-NULL 

#for (i in 1:n0) {
#theta<- prior()
# output<- SimGroup_tauleap(theta1=theta,fish_sex=fishSex,fish_type=Fish_stock,
 #                          strain=Strain,fish_size=fishSize,error=0.01)

#B-D-C parameter estimates for the parasite-fish groups based on simulated data
#for each simulation realisation
#BDC_estimates_sim<- GW_GMM_BDCestimator(X0=2,pop=output$pop_sim,
 #                          output$alive_sim,group=parasite_fish)$BDC_estimates

#Computing the summary stats for each group simulation realisation
#SummaryStats_sim[[i]] <- Summary_stats(pop=output$pop_sim,alive=output$alive_sim,
 #                               BDC_estimates=BDC_estimates_sim) 
#combining for all summary stats of parasite-fish groups for each simulation realisation
 #SummaryStats_sim_combined[[i]]<-do.call("rbind",SummaryStats_sim[[i]])

 #if (i %% 5 == 0) cat(i, "/", n0, "\n") 
 # }

#S0<- do.call("rbind",SummaryStats_sim_combined) #dimension is rows=(n0*total_fish) by cols=17 
#initial weight (inverse of summary statistics)

###### initial weights estimate
setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro")

#w=read.csv(file="initial_weight_w0.csv") #based inverse MAD
w=read.csv(file="w0.csv")#based on inverse variance
w<-as.vector(w[,-1])
print(w)


Total_fish<- dim(do.call("rbind",summaries_obs))[1]

#Setting working directory to save ABC results


sample_sizes<- 1500 #c(500, 1000,1500,2000,2500, 3000, 3500)
setwd("/Users/clementtwumasi/Desktop/Bulletin of Mathematical Biology/Modified ABC-SMC fitting_Gyro/Updated ABC_results_N_1500")
options(warn=-1)


print(paste("At N =",sample_sizes))
ABC_output<- Weighted_iterative_ABC(N=sample_sizes,dimS=17,fish_total=Total_fish,
                                    numCores=numCores,ABC_time_steps=10,seed_num=1)

################################################################################
ABC_time_steps<-10
fx<- ABC_output$fx

#Saving density of final posterior
write.csv(fx[ABC_time_steps+1,,],paste0("final_posterior_density_",sample_sizes,".csv"))


#saving posterior distribution of the 23 model parameters (final posterior)
write.csv(ABC_output$final_posterior,paste0("PosteriorFinal_",sample_sizes,".csv"))







############################################################################
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


N<-sample_sizes

x <- seq(from = -10, to = 7, length.out = 256)#range of prior distribution
# plot posteriors
pdf(paste0("PosteriorPlot_first8_",N,".pdf")) 
par(mfrow=c(4,2), mar=c(4,4,1,1))
for (k in 1:8) {
  plot(x, fx[1,k,], type="l", ylim=c(0,7),xlab=parameter_labels[k],ylab="Density", col="blue")
  for (j in 1:(ABC_time_steps-1)) {
    lines(x, fx[j+1,k,], col="red")
  }
  lines(x, fx[ABC_time_steps+1,k,], col="black", lwd=2)
  legend("topleft",c("First Prior","Intermediate Priors" ,"Final Posterior"),
         col=c("blue","red","black"),bty="n",cex=1,box.lwd = 2,fill=c("blue","red","black"))
}
dev.off()


# plot posteriors
pdf(paste0("PosteriorPlot_second8_",N,".pdf")) 
par(mfrow=c(4,2), mar=c(4,4,1,1))
for (k in 9:16) {
  plot(x, fx[1,k,], type="l", ylim=c(0,7),xlab=parameter_labels[k],ylab="Density", col="blue")
  for (j in 1:(ABC_time_steps-1)) {
    lines(x, fx[j+1,k,], col="red")
  }
  lines(x, fx[ABC_time_steps+1,k,], col="black", lwd=2)
  legend("topleft",c("First Prior","Intermediate Priors" ,"Final Posterior"),
         col=c("blue","red","black"),bty="n",cex=1,box.lwd = 2,fill=c("blue","red","black"))
}
dev.off()


# plot posteriors
pdf(paste0("PosteriorPlot_last7_",N,".pdf")) 
par(mfrow=c(4,2), mar=c(4,4,1,1))
for (k in 17:23) {
  plot(x, fx[1,k,], type="l", ylim=c(0,7),xlab=parameter_labels[k],ylab="Density", col="blue")
  for (j in 1:(ABC_time_steps-1)) {
    lines(x, fx[j+1,k,], col="red")
  }
  lines(x, fx[ABC_time_steps+1,k,], col="black", lwd=2)
  legend("topleft",c("First Prior","Intermediate Priors" ,"Final Posterior"),
         col=c("blue","red","black"),bty="n",cex=1,box.lwd = 2,fill=c("blue","red","black"))
}
dev.off()






