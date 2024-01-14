# Importing required libraries
library("MASS")
library("scipy")

# Function to calculate specific gravity of oil using API gravity
calculate_gamma_o_using_api <- function(api) {
  gamma_o <- (141.5 / api) + 131.5
  return(gamma_o)
}

# Function to calculate density of oil using Standing method
calculate_rho_o_Standing_method <- function(gamma_o, Rs, gamma_g, T) {
  rho_o <- (gamma_o * (62.4 / 32.2)) / (1 + (5.615 * Rs * gamma_g / (T + 460)))
  return(rho_o)
}

# Function to calculate mu_o using Standing method
calculate_mu_o_Standing_method <- function(api, Rs, p, pb, t) {
  A <- 10^(0.43 + (8.33 / api))
  mu_od <- (0.32 + 1.8 * 10^7 / (api^4.53)) * (360/(t + 200))^A
  c <- 8.62 * 10^(-5) * Rs
  d <- 1.1 * 10^(-3) * Rs
  e <- 3.74 * 10^(-3) * Rs
  b <- 0.68 / 10^c + 0.25 / 10^d + 0.062 / 10^e
  a <- Rs * (2.2 * 10^(-7) * Rs - 7.4 * 10^(-4))
  mu_ob <- 10^a * mu_od^b
  mu_o <- mu_ob + 0.001 * (p - pb) * (0.024 * mu_ob^1.6 + 0.38 * mu_ob^0.56)
  return(c(mu_o, mu_ob))
}

# API gravity
api <- 35
# Temperature in F
T <- 120
# Standard gas-oil ratio in scf/STB
Rs <- 800
# Pressure in psia
P <- 3000
# Bubble point pressure in psia
Pb <- 2500

# Specific gravity of gas
gamma_g <- 0.77
# Calculate gamma_o using API
gamma_o <- calculate_gamma_o_using_api(api)
# Calculate rho_o using Standing method
rho_o <- calculate_rho_o_Standing_method(gamma_o, Rs, gamma_g, T)
# Calculate mu_o and mu_ob using Standing method
mu_o <- calculate_mu_o_Standing_method(api, Rs, P, Pb, T)
mu_ob <- calculate_mu_ob_Standing_method(api, Rs, Pb, T)

cat("Above bubble point pressure hence Rs = GOR and since Rs for both is the same rho_o Would be the same as well.\n")
cat("The density of oil is: ", rho_o, " lbm/ft^3\n")
cat("The viscosity of oil at bubble point pressure is: ", mu_ob, " cp\n")
cat("The viscosity of oil at ", P, " psia is: ", mu_o, " cp\n")
