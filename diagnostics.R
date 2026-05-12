################################################################################
# Diagnostics ##################################################################
################################################################################

# Uniform pseudo-samples #######################################################

# Construct uniform pseudo-samples for aggregate MAP
U_full_MAP <- numeric(0)
for(i in 1:N){
  PHic <- matrixdist::ph(alpha = MAP_cond_dist(x = x, 
                                               alpha = alpha_MAP, 
                                               C = C_MAP, 
                                               D = D_MAP, 
                                               k = i),
                         S = C_MAP)
  U_full_MAP[i] <- matrixdist::cdf(PHic, x[i])
}
# Construct uniform pseudo-samples for EOA
U_EOA <- numeric(0)
for(i in 1:N_EOA){
  PHic <- matrixdist::ph(alpha = MAP_cond_dist(x = x_EOA, 
                                               alpha = alpha_MAP, 
                                               C = C_1, 
                                               D = D_1, 
                                               k = i),
                         S = C_1)
  U_EOA[i] <- matrixdist::cdf(PHic, x_EOA[i])
}
# Construct uniform pseudo-samples for NOA
U_NOA <- numeric(0)
for(i in 1:N_NOA){
  PHic <- matrixdist::ph(alpha = MAP_cond_dist(x = x_NOA, 
                                               alpha = alpha_MAP, 
                                               C = C_2, 
                                               D = D_2, 
                                               k = i),
                         S = C_2)
  U_NOA[i] <- matrixdist::cdf(PHic, x_NOA[i])
}
# Construct uniform pseudo-samples for SOA
U_SOA <- numeric(0)
for(i in 1:N_SOA){
  PHic <- matrixdist::ph(alpha = MAP_cond_dist(x = x_SOA, 
                                               alpha = alpha_MAP, 
                                               C = C_3, 
                                               D = D_3, 
                                               k = i),
                         S = C_3)
  U_SOA[i] <- matrixdist::cdf(PHic, x_SOA[i])
}

# Figures ######################################################################

# QQ-plot for aggregate MAP
plot(U_sq_full, U_sq_full, type = "l", col = "black", 
     xlab = "Empirical quantiles", ylab = "Theoretical quantiles",
     main = "QQ-plot - Full processes")
lines(U_sq_full, sort(U_full_MAP), type ="l", col = "#122947", lwd = 2)
legend("topleft", 
       legend = c("MAP"), 
       col = c("#122947"), 
       bty = "n",
       lwd = c(3,3),
       lty = c(1,1),
       pch = c(NA,NA),
       cex = 1.15, 
       text.col = "black", 
       horiz = FALSE , 
       inset = c(0.05))

# QQ-plot for marginal MAPs
plot(U_sq_full, U_sq_full, type = "l", col = "black", 
     xlab = "Empirical quantiles", ylab = "Theoretical quantiles",
     main = "QQ-plot - Marginal processes")

lines(U_sq_NOA, sort(U_NOA), type ="l", col = "#4b8325", lwd = 2, lty = "dotted")
lines(U_sq_EOA, sort(U_EOA), type ="l", col = "#D81B60", lwd = 2, lty = "dotted")
lines(U_sq_SOA, sort(U_SOA), type ="l", col = "#ffbd38", lwd = 2, lty = "dotted")

legend("topleft", 
       legend = c("NOA", "EOA", "SOA"), 
       col = c("#4b8325", "#D81B60", "#ffbd38"), 
       bty = "n",
       pch = c(15,15,15),
       cex = 1.15,
       pt.cex = 2,
       text.col = "black", 
       horiz = FALSE , 
       inset = c(0.05))
