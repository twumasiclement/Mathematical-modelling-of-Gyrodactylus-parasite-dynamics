
#Loading packages
library(transport) #For Wassertein distance computation
library(parallel) # For parallizing R codes
RNGkind("L'Ecuyer-CMRG") #Dealing with distinct seed value in R
#library(markovchain)
#library(diagram)
library('latex2exp') #For adding LaTeX symbols to R plots
library(compiler)# byte code compilation

options(repr.plot.width=8, repr.plot.height=8,repr.plot.res = 300) #Setting plot size


#### External scripts for the simulation model #####

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


#Specific parameter values to test the stochastic simulation algorithm
fixed_par<- function() {
  b11 <- 0.668    # birth rate for young parasites (Gt3)
  b12<- 0.018    # birth rate for older parasites (Gt3)
  b21 <- 0.668    # birth rate for young parasites (Gt)
  b22<- 0.018    # birth rate for older parasites (Gt)  
  b31 <- 0.668    # birth rate for young parasites (Gb)
  b32<- 0.018    # birth rate for older parasites (Gb)  
    
  d11 <- 0.008   # death rate without an immune response (Gt3)
  d12 <-  0.071   # death rate with immune response (Gt3)
  d21 <- 0.008   # death rate without an immune response (Gt)
  d22 <-  0.071   # death rate with immune response (Gt)
  d31 <- 0.008   # death rate without an immune response (Gb)
  d32 <-  0.071   # death rate with immune response (Gb)
    
  m<- 0.083#movement rate
  r <- 0.001# immune response rate (base rate)
  r1 <- 0.351# immune response (adjustment for LA fish)
  r2 <- 0.196 # immune response rate (adjustment for OS fish)
  r3 <- 0.994 # immune response rate (adjustment for male fish)
  s <-0.009 #fish mortality rate (base rate)
  s1 <- 0.041#fish mortality (adjustment for male fish)
  e1 <- 0.545 #rate of forward or backward movement/preference bias (Gt3)
  e2 <- 0.333 #rate of forward or backward movement/preference bias (Gt)
  e3 <- 0.001 #rate of forward or backward movement/preference bias (Gb)
  kappa <- 182  #effective carrying capacity
  return(c(b11,b12,b21,b22,b31,b32,d11,d12,d21,d22,d31,d32,m,r,r1,r2,r3,s,s1,e1,e2,e3,kappa))
}


fixed_ParaValues<- fixed_par()
fixed_ParaValues

fixed_ParaValues_log<-log(fixed_ParaValues)
fixed_ParaValues_log#log scales
length(fixed_ParaValues_log)

#Initial simulation inputs for A (parasite numbers) and  B (immune status)
A0 <- matrix(0, 4, 2)  
A0[1, 1] <- 2   #Intial parasites at the tail
B0 <- rep(1, 4)  #initial immune response at 4 body regions (1=no response, 2=response)

#Transition matrix
J<- matrix(c(0,    1,    0,     0, 
             1/2,   0,    1/2,   0,
             0,    1/2,    0,   1/2,
             0,     0,     1,    0), 4, 4, byrow=TRUE)
        

#To simulate group of fish for each parasite-fish combination as observed in the empirical data
SimGroup_tauleap<- function(theta1,fish_sex,fish_type,strain,
                            fish_size,error){
    #Inputs: theta1= parameter values from prior distribution
    #fish_sex= sex of fish 
    #fish_type= type of fish
    #strain= parasite strain
    #fish_size= fish size
    #error= error bound of tau-leaping
    pop_sim<-NULL; alive_sim<- NULL;exploded_sim<-NULL;results=NULL;group=NULL
   
  for(pf in 1:9){ 
     pop_sim[[pf]]<- array(dim = c(numF[[pf]], 4, 9))
      #Array for time steps fish was alive for each combination
     alive_sim[[pf]]<-array(dim = c(numF[[pf]], 9))
      #Array for time steps parasites>pop_max for each combination
     exploded_sim[[pf]]<-array(dim = c(numF[[pf]], 9)) 
      
      for(i in 1:numF[[pf]]){
          results[[pf]]<-sim_tauleap_singlefish_compiler(A0=A0, B0=B0,J=J, 
                         b1=matrix(exp(theta1[1:2]), 2, 2),b2=matrix(exp(theta1[3:4]), 2, 2),
                         b3=matrix(exp(theta1[5:6]), 2, 2),d1=matrix(exp(theta1[7:8]), 2, 2, byrow=TRUE),       
                         d2=matrix(exp(theta1[9:10]), 2, 2,byrow=TRUE),d3=matrix(exp(theta1[11:12]),2, 2,byrow=TRUE),       
                         m=matrix(exp(theta1[13]), 2, 2),r=exp(theta1[14]),r1=exp(theta1[15]),r2=exp(theta1[16]), 
                         r3=exp(theta1[17]),s=exp(theta1[18]),s1=exp(theta1[19]),e1=exp(theta1[20]),
                         e2=exp(theta1[21]),e3=exp(theta1[22]),kappa=exp(theta1[23]),f=Area_normalized,
                         a=fish_size[[pf]][i],fish_sex=fish_sex[[pf]][i],fish_type=fish_type[[pf]][i],
                         strain=strain[[pf]][i],error=error)
          #print(results)
         pop_sim[[pf]][i, ,]<- results[[pf]]$pop
         alive_sim[[pf]][i, ]<-  results[[pf]]$alive
         exploded_sim[[pf]][i, ]<- results[[pf]]$exploded
         group[[pf]]<-results[[pf]]$parasite_fish
                           }   
                 }
    
    #Output: returns pop_sim (parasite pop at each location and time)
    #alive_sim: survival status of fish
    # exploded_sim: explosion status (whether parasite numbers>pop_max=10000)
    # group: the host-parasite groups being simulated
    
  return(list(pop_sim=pop_sim,alive_sim=alive_sim,exploded_sim=exploded_sim,group= unlist(group)))
    }

#Testing the group simulator
results<-SimGroup_tauleap(theta1=fixed_ParaValues_log,fish_sex=fishSex,fish_type=Fish_stock,strain=Strain,
                            fish_size=fishSize,error=0.01)
group_number<-1
results$pop_sim[[group_number]][8,,]#view simulated parasite data for 8th fish for group 1
results$alive_sim[[group_number]][8, ]
results$group[group_number]
