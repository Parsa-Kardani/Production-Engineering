# Importing required libraries
library("MASS")
library("scipy")

# Function to calculate specific gravity of oil using API gravity
calculate_gamma_o_using_api <- function(api) {
  gamma_o <- (141.5 / (api + 131.5))
  return(gamma_o)
}


# Function to calculate density of oil using Standing method
calculate_rho_o_Standing_method <- function(gamma_o, Rs, gamma_g, T) {
  rho_o <- (gamma_o * (62.4 / 32.2)) / (1 + (5.615 * Rs * gamma_g / (T + 460)))
  return(rho_o)
}

api <- 25 # API gravity
T <- 100 # Temperature in F
Rs <- 0 # Standard gas-oil ratio in scf/STB
gamma_g <- 0 # Specific gravity of gas

gamma_o <- calculate_gamma_o_using_api(api)
rho_o <- calculate_rho_o_Standing_method(gamma_o, Rs, gamma_g, T)

cat(rho_o, 'lbm/ft^3')