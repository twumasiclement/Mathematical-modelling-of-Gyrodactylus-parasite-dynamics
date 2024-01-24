
# Function for computing weighted distance between simulated and observed summary statistics
#Weighted sums of squares across all fish between simulation and observed data
#where S1 and S2 are matrices of summary stats for the entire fish (between observed & simulation)

w_distance <- function(S1, S2, weight)  {
  #weight is a vector/one-dimensional array (summary stats weights)
  normalised_func<- function(x) x/sum(x,na.rm=TRUE)
  weight<- normalised_func(weight)
  n<- dim(S1)[1]
  Squared_diff_mat<- (S1-S2)^2 #squared difference between matrix S1 & S2
  #Multiplying vector to weights to each rows of  Squared_diff_mat
  Weighted_sq_diff<- lapply(1:dim(S1)[1],
                            function(k) weight*Squared_diff_mat[k, ])
  WSS<- do.call("sum",Weighted_sq_diff)#total weighted distances (WSS)
  return(sqrt(WSS/n))#return weighted sum of square distance 
}


distance_compiler=cmpfun(w_distance)







# Function for computing weighted Euclidean distance between simulated and observed summary statistics
#where S1 and S2 are matrices of summary stats for the entire fish (between observed & simulation)

w_distance_avg <- function(S1, S2, weight)  {
  #weight is a vector/one-dimensional array (summary stats weights)
  mean_func<- function(x) mean(x, na.rm=TRUE)
  normalised_func<- function(x) x/sum(x,na.rm=TRUE)
  
  S1_vector<- apply(S1, 2, mean_func)
  S2_vector<- apply(S2, 2, mean_func)
  dimS<- length(S1_vector)
  weight<- normalised_func(weight)
  Squared_diff<- (S1_vector- S2_vector)^2
  Weighted_sq_diff<-   weight*Squared_diff
  Weighted_distance<- sqrt(sum(Weighted_sq_diff)/dimS)

  
  #return weighted sum of square distance (scaled by 1/dimS)
  return(Weighted_distance)
}


