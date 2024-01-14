# Importing required libraries
library("MASS")
library("scipy")

C_to_F <- function(C) {
  F = (C * 9/5) + 32
  return F
}

mPa_to_psia <- function(mPa) {
  psia <- mPa * 145.038
  return(psia)
}

calculate_api = function(gamma_o) {
  api = 141.5 / gamma_o - 131.5
  return(api)
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

# Convert temperature from Celsius to Fahrenheit
C_to_F <- function(T_C) {
  T_F <- (T_C * 9/5) + 32
  return(T_F)
}

# Convert pressure from mPa to psia
mPa_to_psia <- function(P_mPa) {
  P_psia <- P_mPa * 0.145038
  return(P_psia)
}

# Calculate API gravity
calculate_api <- function(gamma_o) {
  API <- (141.5 / gamma_o) - 131.5
  return(API)
}

# Calculate oil density using Standing method
calculate_rho_o_Standing_method <- function(gamma_o, Rs, gamma_g, T) {
  rho_o <- (gamma_o * 62.4) / (5.615 * (1 + (Rs * gamma_g * T) / 460))
  return(rho_o)
}

# Calculate oil viscosity using Standing method
calculate_mu_o_Standing_method <- function(api, Rs, P, Pb, T) {
  mu_ob <- 10 ^ (3.0324 - 0.02023 * api)
  mu_o <- mu_ob * exp(3.5 * log((P / Pb) * (1 + 0.0395 * (Rs * T) ^ 0.5) - 1))
  return(list(mu_o, mu_ob))
}

T <- 50 # Temperature in C
Rs <- 4000 # Standard gas-oil ratio in sm^3/m^3
P <- 20 # Pressure in mPa
Pb <- 15 # Bubble point pressure in mPa
gamma_g <- 0.77 # Specific gravity of gas
gamma_o <- 0.8 # Specific gravity of oil

T <- C_to_F(T)
P <- mPa_to_psia(P)
Pb <- mPa_to_psia(Pb)
Rs <- Rs / 5.615 # convert sm^3/m^3 to scf/STB

api <- calculate_api(gamma_o)
rho_o <- calculate_rho_o_Standing_method(gamma_o, Rs, gamma_g, T)
mu_o <- calculate_mu_o_Standing_method(api, Rs, P, Pb, T)
mu_ob <- calculate_mu_o_Standing_method(api, Rs, P, Pb, T)

cat("Above bubble point pressure hence Rs = GOR and since Rs for both is the same rho_o Would be the same as well.\n")
cat("The density of oil is: ", rho_o, " lbm/ft^3\n")
cat("The viscosity of oil at bubble point pressure is: ", mu_ob, " cp\n")
cat("The viscosity of oil at ", P, " psia is: ", mu_o, " cp\n")

