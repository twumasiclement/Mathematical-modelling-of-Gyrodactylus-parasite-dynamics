
library(compiler)# byte code compilation

#Function for updating tau-leaping
taulealping_update_event <-function(A,B,J,X,laB,laD,laM_forward,laM_backward,laI,laX,la,
                                   QB,QD,QM_forward,QM_backward,QI,tau) {
    
  #Inputs: 
  # A[j,k] gives the number of parasites at location j, age k, where
  # B[j] = immune response at location j (1 for no response; 2 for a response)
  #J is transition matrix 
  #X is survival status (1= alive, 2=dead)
  # tau is the leap size
  #And all rates in relation to birth, death, movement, immune response, host mortality
  # and total rate (la)
       
             U <- runif(1, 0, la)
             if(U<laX)  X <-2 #Fish mortality
             else if(U<sum(c(laX,laI))){  #Immune response
                 
                j <- sample(4, 1, prob = abs(QI))
                B[j] <- 2
           } else if (U<sum(c(laX,laI,laB,laD,laM_forward))){ #birth, death or forward movement
                i <- sample(8, 1, prob = abs(QB+ QD+QM_forward))
                j <- ((i-1) %% 4) + 1 # current location
                k <- ((i-1) %/% 4) + 1 # age
                j_new <- sample(4, 1, prob =abs(J[j,]))# new location
                A[j,k]<- A[j,k] + rpois(1,abs(laB*tau))-rpois(1,abs(laD*tau))-rpois(1,abs(laM_forward*tau))               
                A[j_new,k]<-  A[j_new,k]+rpois(1,abs(laB*tau))-rpois(1,abs(laD*tau))+rpois(1,abs(laM_forward*tau))
                 
          }  else if (U< sum(c(laX,laI,laB,laD,laM_forward,laM_backward))){ #brith, death or backward movement
                i <- sample(8, 1, prob = abs(QB+ QD+QM_backward))
                j <- ((i-1) %% 4) + 1 # current location
                k <- ((i-1) %/% 4) + 1 # age
                j_new <- sample(4, 1, prob =abs(J[j,]))# new location
                A[j,k]<- A[j,k]+rpois(1,abs(laB*tau))-rpois(1,abs(laD*tau))-rpois(1,abs(laM_backward*tau))
                A[j_new,k]<-A[j_new,k]+rpois(1,abs(laB*tau))-rpois(1,abs(laD*tau))+rpois(1,abs(laM_backward*tau)) 
             } 
      
    #Output: returns A[j,k] gives the number of parasites at location j, age k, where
    #   B[j] = immune response at location j (1 for no response; 2 for a response)
    # X=survival status
       #print(paste(A,ti,"Results from Tauleaping",tau,"Condition",leap_condition)) 
         return(list(A = A, B = B, X = X))
    
    }

taulealping_update_event_compiler<-cmpfun(taulealping_update_event)
