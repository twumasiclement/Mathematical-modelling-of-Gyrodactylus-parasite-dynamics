
#Setting working directory
#setwd("/home/clement/Documents/Simulation_Data_Folder")

library(compiler)# byte code compilation

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

na.zero<-function(x){
  x[is.na(x)]<-0
  return(x)
}


#We have assumed that the rate of host immune response depends on the parasites at each j=1,2,3,4 body location.
#Suppose we want to make host immune response dependent on fish sex (with 2 levels) and fish type (with three levels); we use the following (tabulated) rates depending on the sex and fish type of interested: 

#gender\fish type  |  UA        |  LA        |  OS  |
#-------------------------------------------------------------

#female            |r          |r+r1         |r+r2  |

#male              |r+r3        |r+r1+r3|     r+r2+r3|

#then, whichever rate you are using (depending on fish sex and fish type), you multiply that by the population.

#NB: #r = rate a single parasite increases immune state (base rate of immune response): 
 #r1= immune response rate adjustment for LA fish (LA fish considered as reference)
 # r2= immune response rate adjustment for OS fish  (OS fish considered as reference)
 # r3= immune response rate adjustment for Male fish (Male fish considered as reference)
 

#Also, host mortatility rate must depend on total parasite numbers, fish sex and fish size.

#gender          rate
#----------------------------------
#female            s (base rate)
#male              s+s1 (s1 is the adjustment for male fish)

#female fish is kept as reference. fish size is only considered in the population carrying capacity.
#Finally, microhabitat preference rate depends on parasite strain
#############################################################################



##### Description of variables#######
# State variables

  # A[j,k] gives the number of parasites at location j, age k, where

  #   j = 1 for Tail population

  #   j = 2 for Lower region population

  #   j = 3 for Upper region population

  #   j = 4 for head population

  #   k = 1 for young parasites (yet to give birth)

  #   k = 2 for old parasites (have given birth)

  # B[j] = immune response at location j (1 for no response; 2 for a response)

  # X = state of fish (1 for alive; 2 for dead)

 

  # Base simulation parameters

  # b1[k,el] = birth rate for parasites age k, when immune state is el (for Gt3)
  # b2[k,el] = birth rate for parasites age k, when immune state is el (for Gt)
  # b3[k,el] = birth rate for parasites age k, when immune state is el (for Gb)
    
  # d1[k,el] = death rate for parasites age k, when immune state is el (for Gt3)
  # d2[k,el] = death rate for parasites age k, when immune state is el (for Gt)
  # d3[k,el] = death rate for parasites age k, when immune state is el (for Gb)

  # m[k,el] = movement rate for parasites age k, when immune state is el

 # e = the adjustment to the movement rate for forward/backward movement

  # r = rate a single parasite increases immune state (base rate of immune response)

  # kappa = effective carrying capacity per unit area

  # area of each body part

  # s = rate a single parasite causes fish mortality
  # a=fish size
 

  # Additional simulation parameters

  # r1 = immune response rate adjustment for male fish (female fish considered as reference)

  # r2 = immune response rate adjustment for LA fish  (UA fish considered as reference)

  # r3 = immune response rate adjustment for OS fish  (UA fish considered as reference)

  # e1, e2, e3 = movement rate adjustment depending on parasite type

 #s = must depend on total parasite numbers, fish sex and fish size
 # s1 = host mortality rate with adjustment for male fish 
 

  # Experiment descriptors

  # fish_type (1 for UA, 2 for LA & 3 for OS) 

  # Parasite type (Gt3, Gt & Gb)

  # fish_sex (1 for female fish & 2 for male fish)

  # f= area of each body part (depends on size and gender)

sim_tauleap_singlefish<- function(A0, B0,J, b1,b2,b3,d1,d2,d3, m, r,r1,r2,r3,s,s1,e1,e2,
                                  e3,kappa,f,a,fish_sex,fish_type,strain,error){   

 #Inputs: inital conditions, parameter values, fish sex, 
 #fish type, parasite strain and error bound
 #f=body area (dependent on fish sex and size)
 #parasite_fish=c("Gt3-OS","Gt3-LA","Gt3-UA","Gt-OS","Gt-LA","Gt-UA","Gb-OS","Gb-LA","Gb-UA") 
     
 #strain-parasite type to be simulated 
 parasite_fish<-paste(strain,"-",fish_type)    

 pop=NULL;alive =NULL;exploded=NULL; Leap_sizes=NULL;A<-A0; B<- B0 

 save_ti <- c(1, 3, 5, 7, 9, 11, 13, 15, 17) #observed discrete times
 save_TF <- rep(FALSE, length(save_ti))
         
 ti<- 0 # time
 # parasite pop at each location (rows) and timepoint (cols)
  pop[[parasite_fish]] <- matrix(NA, 4, length(save_ti))
 # host fish status at each time point
  alive[[parasite_fish]] <- rep(2, length(save_ti)) 
    
  pop_ti <- rowSums(A)
  
  # host survival status (alive=1; dead=2)
  alive_ti <- 1 
  exploded[[parasite_fish]] <- FALSE   
    
  pop_max <- 10000 # stop the simulation if total population exceeds this limit
  X <- 1 # fish starts out alive
     
 while(sum(save_TF) < length(save_ti)){ ### beginning of main while loop
     
  #### Computing the rates ######
        computed_rates<- compute_rates_compiler(A=A, B=B, b1=b1,b2=b2,b3=b3, d1=d1,
                                     d2=d2,d3=d3, m=m,r=r,r1=r1,r2=r2,r3=r3,s=s,
                                     s1=s1,e1=e1,e2=e2,e3=e3,kappa=kappa,f=f,a=a,
                                     fish_sex=fish_sex,fish_type=fish_type,strain=strain)
     
        laB<-computed_rates$laB;laD<-computed_rates$laD;laM_forward<-computed_rates$laM_forward;
        laM_backward<-computed_rates$backward;laI<-computed_rates$laI;laX<-computed_rates$laX;
        la<-computed_rates$la;QB<-computed_rates$QB;QD<-computed_rates$QD;
        QM_forward<-computed_rates$QM_forward;QM_backward<-computed_rates$QM_backward;QI<-computed_rates$QI
     
###### Determining the  condition to switch between the exact SSA and the Tau-leaping algorithm #####     
     ####estimate of the leap size (tau)
      #Computing tau  based on the B-D-C leap size formula (Gillespie & Petzold, 2003)
        #selecting birth and deaths rates depending on parasite strain for leap size
       if(strain=="Gt3"){b_selected<-b1; d_selected<-d1}
       if(strain=="Gt"){b_selected<- b2; d_selected<-d2}
       if(strain=="Gb"){b_selected<- b3; d_selected<-d3}
     
        b_avg<- mean(b_selected[,1]) #finding average birth rate for young & old parasites 
        d_avg<-mean(d_selected[, 1])#finding average death rate with or without host immune response
        Leap_sizes[[1]]<- (error*(b_avg+d_avg))/(abs(b_avg-d_avg)*max(b_avg,d_avg))
        Leap_sizes[[2]]<- sum(A)*(error*(b_avg+d_avg))^2/((b_avg+d_avg)* max(b_avg^2,d_avg^2))

        tau<- na.zero(min(Leap_sizes[[1]],Leap_sizes[[2]])) # the leap size: time increment for tau-leaping 
        leap_condition<- na.zero((1/(10*la)))  #la = total rate
         
     
        if (sum(pop_ti) > pop_max) {
            exploded[[parasite_fish]] <- TRUE
            break
                                   }
     
        if (alive_ti == 2)  break
    
      #Running Tau-leaping if tau >leap_condition
        if(tau >leap_condition){ #Execute tau-leaping
              out<-taulealping_update_event_compiler(A=A,B=B,J=J,X=X,laB=laB,laD=laD,
                   laM_forward=laM_forward,laM_backward=laM_backward,laI=laI,laX=laX,la=la,
                   QB=QB,QD=QD,QM_forward=QM_forward,QM_backward=QM_backward,QI=QI,tau=tau) 
              X<- out$X
              A<- out$A
              B<- out$B
              time_increment=tau 
   # print(paste(rowSums(A),"",ti,"Results from Tauleaping",tau,"Condition",leap_condition))   
                               } #end of tau-leaping 

        else if(tau <=leap_condition){#Execute exact SSA if tau <=leap_condition
            
              out <- SSA_update_event_compiler(A=A,B=B,J=J,X=X,laB=laB,laD=laD,
                     laM_forward=laM_forward,laM_backward=laM_backward,laI=laI,
                     laX=laX,la=la,QB=QB,QD=QD,QM_forward=QM_forward,QM_backward=QM_backward,QI=QI)

            
              time_increment<- out$t_incr_SSA #time increment for SSA
              X<- out$X
              A<- out$A
              B<- out$B
         
    #print(paste(rowSums(A),"",ti,"Results from SSA since tau=",tau,"Condition",leap_condition))            
                                   } #end of exact SSA
              
  ti <- ti +time_increment #updating time ti
  save_new <- which((ti >= save_ti) & !save_TF)
            
         for (i in save_new) {
              pop[[parasite_fish]][,i] <- pop_ti
              alive[[parasite_fish]][i] <- alive_ti
                             }
              save_TF <- (ti >= save_ti)
     
              pop_ti <- rowSums(A)
              alive_ti <- X 
     
     #break if parasite number is negative at any body region
         if(any(pop_ti<0) ==TRUE) break  
     
            } #end main while loop  
 
  #Output: returns pop (parasite pop at each location and time)
  #alive: survival status of fish
  # exploded: explosion status (whether parasite numbers>pop_max=10000)
  # parasite_fish: the host-parasite group being simulated
  return(list(pop = pop[[parasite_fish]], alive = alive[[parasite_fish]],
            exploded = exploded[[parasite_fish]],parasite_fish=parasite_fish))
   
    
  }



sim_tauleap_singlefish_compiler<- cmpfun(sim_tauleap_singlefish) #Converting to  byte code compilation
