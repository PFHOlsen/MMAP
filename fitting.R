################################################################################
# Statistical fitting by EM algorithms #########################################
################################################################################

# MAP fitting ##################################################################

# Set MAP dimension
p <- 4

# Set number of iterations for EM algorithm
N_iter <- 500

set.seed(1)
# Create MAP object with general structure
MAP0 <- MAP(p = p, structure = "general")
# Extract parameters
alpha_MAP <- MAP0$alpha
C_MAP <- MAP0$C
D_MAP <- MAP0$D
# Store log-likelihood
LL_MAP_vec <- numeric(0)
# Run EM algorithm for MAPS - updates alpha_MAP, C_MAP, D_MAP
for(i in 1:N_iter){
  EM_MAP(x = x, 
         alpha = alpha_MAP, 
         C = C_MAP, 
         D = D_MAP, 
         benchmark = FALSE) # Set TRUE to see runtime for each part of E-step
  LL_MAP_vec[i] <- sum(log(D_weights(x = x, 
                                     alpha = alpha_MAP, 
                                     C = C_MAP, D = D_MAP)))
  cat("\r", paste0("Iteration: ", i, 
                   ", loglikelihood: ", round(LL_MAP_vec[i],3)))
}

# Final log-likelihood
tail(LL_MAP_vec,1)

# MMAP fitting #################################################################

# Set number of dimensions - here 3 corresponding to the areas EOA, NOA, SOA.
n <- 3

# Create list of n pxp matrices with non-zero elements 
M_ones <- matrix(rep(1,p*p), ncol = p)
R_list <- list(M_ones/n,M_ones/n,M_ones/n)

# Create list with arrival times for each area
x_list <- list(x_EOA, x_NOA, x_SOA)

# Run EM algorithm for MMAPs - updates elements in R_list
for(i in 1:10){
  R_list <- EM_MMAP(alpha = alpha_MAP, 
                    C = C_MAP, 
                    D = D_MAP, 
                    x_list = x_list, 
                    R_list = R_list)
  cat("\r", paste0("Iteration: ", i))
}

# Extract weight matrices
R1 <- R_list[[1]]
R2 <- R_list[[2]]
R3 <- R_list[[3]]

# Compute marginal representations

# Marginal representation for EOA
C_1 <- C_k(k = 0, C = C_MAP, D = D_MAP, R_list = R_list)
D_1 <- D_k(k = 0, D = D_MAP, R_list = R_list)
# Marginal representation for N_NOA
C_2 <- C_k(k = 1, C = C_MAP, D = D_MAP, R_list = R_list)
D_2 <- D_k(k = 1, D = D_MAP, R_list = R_list)
# Marginal representation for N_NOA
C_3 <- C_k(k = 2, C = C_MAP, D = D_MAP, R_list = R_list)
D_3 <- D_k(k = 2, D = D_MAP, R_list = R_list)

# Check that D_MAP was split correctly
all.equal(D_1 + D_2 + D_3, D_MAP)
