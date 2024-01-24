#Function to store summary stats across all fish for ABC
Summary_stats<-function(pop,alive,BDC_estimates){
  parasite_fish<-c("Gt3-OS","Gt3-LA","Gt3-UA","Gt-OS","Gt-LA",
                   "Gt-UA","Gb-OS","Gb-LA","Gb-UA") 
  dimS <- 17 #length of summary statistics
  summaries=NULL
  for (pf in 1:length(parasite_fish)){
    summaries[[pf]] <- array(dim=c(numF[[pf]],dimS))
    for (k in 1:numF[[pf]]) {
      summaries[[pf]][k,] <- summary_compiler(pop_single=pop[[pf]][k,,], 
                                              alive_single=alive[[pf]][k,],group_number=pf,
                                              BDC_estimates=BDC_estimates)
    }
  }
  return(summaries) #returns the 17 summary stats
}