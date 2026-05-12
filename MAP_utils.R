################################################################################
# Helper functions for Markovian Arrival Processes (MAP) #######################
################################################################################

# Make sure to source EM_algorithms for matrix_exponential() function
# Rcpp::sourceCpp("EM_algorithms.cpp")

# Load package matrixdist for phase-type objects
library(matrixdist)

# Set initial parameters for MAP(alpha, C, D)
MAP <- function(p, structure = NULL){
  PH0 <- matrixdist::ph(structure = ifelse(is.null(structure), 
                                           "general", 
                                           structure), 
                        dimension = p)
  alpha <- PH0@pars$alpha
  D <- matrix(ncol = p, nrow = p)
  for(i in 1:p){
    for(j in 1:p){
      D[i,j] <- rexp(1)
    }
  }
  C <- PH0@pars$S
  diag(C) <- rep(0,p)
  diag(C) <- -rowSums(C+D)
  return(list(alpha = alpha, C = C, D = D))
}

# Computes the distribution of the underlying Markov jump process of a 
# MAP(alpha, C, D) after k arrivals from the ordered sample x
MAP_cond_dist <- function(x, alpha, C, D, k){
  p <- length(alpha)
  v <- alpha
  if(k>1){
    for(i in 1:(k-1)){
      v <- v %*% matrix_exponential(C * x[i]) %*% D
      v <- v / sum(v) 
    }
  }
  v <- v / sum(v)
  return(v)
}
