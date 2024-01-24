
library(compiler)# byte code compilation

# Function for updating exact SSA
#For for updating simulation events across the 4 body regions (Tail, Lower region, Upper region, Head)
SSA_update_event <-function(A,B,J,X,laB,laD,laM_forward,laM_backward,laI,laX,la,
                           QB,QD,QM_forward,QM_backward,QI) {
    
     #Inputs: 
  # A[j,k] gives the number of parasites at location j, age k, where
  # B[j] = immune response at location j (1 for no response; 2 for a response)
  #J is transition matrix 
  #X is survival status (1= alive, 2=dead)
  #And all rates in relation to birth, death, movement, immune response, host mortality
  # and total rate (la)

 if (la == 0) {
    return(list(A = A, B = B, t_incr_SSA = Inf, X = X)) # zero population
          }
    
  U <- runif(1, 0, la)  #uniform random number/generator
    
  if (U < laB) {# birth
    i <- sample(8, 1, prob = abs(QB))
    j <- ((i-1) %% 4) + 1 # location
    k <- ((i-1) %/% 4) + 1 # age
    if (k == 1) {
      A[j, 2] <- A[j, 2] + 1
    } else {
      A[j, 1] <- A[j, 1] + 1
    }
  } else if (U < sum(c(laB,laD))) {# death
    i <- sample(8, 1, prob = abs(QD))
    j <- ((i-1) %% 4) + 1 # location
    k <- ((i-1) %/% 4) + 1 # age
    A[j, k] <- A[j, k] - 1
  } else if (U < sum(c(laB,laD,laM_forward)) ) {# forward movement
      
    i <- sample(8, 1, prob = abs(QM_forward))
    j <- ((i-1) %% 4) + 1 # location
    k <- ((i-1) %/% 4) + 1 # age
    j_new <- sample(4, 1, prob =abs(J[j,]))# new location
    A[j, k] <- A[j, k] - 1
    A[j_new, k] <- A[j_new, k] + 1
  } else if (U < sum(c(laX,laI,laB,laD,laM_forward)) ){# backward movement
    i <- sample(8, 1, prob = abs(QM_backward))
    j <- ((i-1) %% 4) + 1 # location
    k <- ((i-1) %/% 4) + 1 # age
    j_new <- sample(4, 1, prob =abs(J[j,]))# new location
    A[j, k] <- A[j, k] - 1
    A[j_new, k] <- A[j_new, k] + 1
  }else if (U < sum(c(laB,laD,laM_forward,laM_backward,laI)) ){# immune response
      
    i <- sample(4, 1, prob = abs(QI))
    B[i] <- 2
  } else {# fish death
    X <- 2
  }
   t_incr_SSA <- rexp(1, la) # time increment for exact SSA

    #Output: returns A[j,k] gives the number of parasites at location j, age k, where
   #  # B[j] = immune response at location j (1 for no response; 2 for a response)
   # t_incr_SSA= time increment for exact SSA
   # X survival status
  return(list(A = A, B = B,t_incr_SSA=t_incr_SSA, X = X))
}

SSA_update_event_compiler<-cmpfun(SSA_update_event) #Converting function in C
