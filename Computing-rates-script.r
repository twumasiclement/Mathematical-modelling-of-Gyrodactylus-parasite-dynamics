
library(compiler)# byte code compilation

na.zero<-function(x){
  x[is.na(x)]<-0
  return(x)
}

# Function for computing rates based on fish sex, fish type and parasite strain
compute_rates_func<- function(A, B, b1,b2,b3, d1,d2,d3, m, r,r1,r2,r3,s,s1,e1,e2,e3,kappa,f,a,
                              fish_sex,fish_type,strain){
  
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
  
  
  
  # Additional simulation parameters
  
  # r1 = immune response rate adjustment for LA fish  (UA fish considered as reference)
  # r2 = immune response rate adjustment for OS fish  (UA fish considered as reference)
  # r3 = immune response rate adjustment for male fish (female fish considered as reference)
  
  # e1, e2, e3 = movement rate adjustment depending on parasite type (Gt3,Gt, Gb respectively)
  
  #s = must depend on total parasite numbers, fish sex and fish size
  # s1 = host mortality rate with adjustment for male fish 
  
  
  # Experiment descriptors
  
  # fish_type (1 for UA, 2 for LA & 3 for OS) 
  
  # Parasite type (Gt3, Gt & Gb)
  
  # fish_sex (1 for female fish & 2 for male fish)
  
  # f= area of each body part (depends on size and gender)
  # a= fish size
  
  r_selected<-0;s_selected<-0 #initialise
  r_matrix<- matrix(c(r,r+r1,r+r2,r+r3,r+r1+r3,r+r2+r3),nrow=2,ncol=3,byrow=T)
  #r_selected= be the selected rate based on adjustments for fish sex and fish type 
  #selecting the immune response rate combination to select depending on fish sex and fish type
  if(fish_sex=="F" & fish_type=="UA") {r_selected<-  r_matrix[1,1]}# base rate 
  if(fish_sex=="F" & fish_type=="LA") {r_selected<-  r_matrix[1,2]}# adjustment for LA fish
  if(fish_sex=="F" & fish_type=="OS") {r_selected<-  r_matrix[1,3]}# adjustment for OS fish
  if(fish_sex=="M" & fish_type=="UA") {r_selected<-  r_matrix[2,1]}# adjustment for male fish
  if(fish_sex=="M" & fish_type=="LA") {r_selected<-  r_matrix[2,2]}# adjustment for male fish & LA fish
  if(fish_sex=="M" & fish_type=="OS") {r_selected<-  r_matrix[2,3]}# adjustment for male fish & OS fish
  
  # selecting which host mortality rate combination & body areas to select fish sex
  if(fish_sex=="F"){
    s_selected<- s #base host mortality rate 
    #f=body_area
    f<-as.vector(f[,1])# body areas for female fish (tail, lower region, upper region & head)
  }
  
  if(fish_sex=="M"){
    s_selected<- s+s1 #host mortality rate with adjustment for male fish 
    #f=body_area
    f<-as.vector(f[,2])## body areas for male fish
  }    
  
  #selecting microhabitat preference rate depending on parasite strain
  if(strain=="Gt3"){e_selected<- e1}
  if(strain=="Gt"){e_selected<- e2}
  if(strain=="Gb"){e_selected<- e3}
  
  #selecting birth and deaths rates depending on parasite strain
  if(strain=="Gt3"){b_selected<-b1; d_selected<-d1}
  if(strain=="Gt"){b_selected<- b2; d_selected<- d2}
  if(strain=="Gb"){b_selected<- b3; d_selected<- d3}
  
  # birth rates; death rates; movement rates; immune response
  QB <- matrix(0, 4, 2) # QB[k,j] = birth rate for parasites location j age k 
  QD <- matrix(0, 4, 2) # QD[k,j] = death rate for parasites location j age k 
  QM_forward <- matrix(0, 4, 2) # QM[k,j] = movement rate for parasites location j age k 
  QM_backward <- matrix(0, 4, 2) 
  QI <- rep(0, 4) # QI[j] = rate at which location j increases immune response
  for (j in 1:4) {
    QI[j] <- sum(A[j, ]) * r_selected
    for (k in 1:2) {
      QB[j, k] <- A[j, k] * (1-(A[j, k]/(f[j]*a*kappa)))*b_selected[k, B[j]]
      QD[j, k] <- A[j, k] * (1-(A[j, k]/(f[j]*a*kappa)))*d_selected[k, B[j]]
      QM_forward[j, k] <- A[j, k] *m[k, B[j]]*e_selected
      QM_backward[j, k] <- A[j, k] *m[k, B[j]]*(1-e_selected)  
    }
  }
  
  # total rates (for the higher-level events)
  laB <- sum(QB) #total birth rate
  laD <- sum(QD) #total death rate
  laM_forward <- sum(QM_forward) # total rate for forward movement
  laM_backward <- sum(QM_backward) # total rate for backward movement
  laI <- sum(QI)# total rate of immuune response
  laX <- sum(A) * s_selected # host fish mortality rate
  la <- na.zero(abs(laB + laD + laM_forward+laM_backward+ laI + laX)) #overall total
  
  #Returns rates in relation to birth, death, movement, immune response, host mortality
  # and total rate (la)
  return(list(laB=laB,laD=laD,laM_forward=laM_forward,
              laM_backward=laM_backward,laI=laI,laX=laX,
              la=la,QB=QB,QD=QD,QM_forward=QM_forward,
              QM_backward=QM_backward,QI=QI))
  
}

compute_rates_compiler<- cmpfun(compute_rates_func)
