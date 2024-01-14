# Importing required libraries
library("MASS")
library("scipy")

C_to_F <- function(C) {
  F = (C * 9/5) + 32
  return F
}

# Function to calculate density of oil using Standing method
calculate_rho_o_Standing_method <- function(gamma_o, Rs, gamma_g, T) {
  rho_o <- (gamma_o * (62.4 / 32.2)) / (1 + (5.615 * Rs * gamma_g / (T + 460)))
  return(rho_o)
}

gamma_o <- 0.8 # specific gravity of oil
T <- 40 # temperature in C
Rs <- 0 # Standard gas-oil ratio in scf/STB
gamma_g <- 0 # Specific gravity of gas
T <- C_to_F(T)
rho_o <- calculate_rho_o_Standing_method(gamma_o, Rs, gamma_g, T)

cat(rho_o, 'lbm/ft^3')