# Importing required libraries
library("MASS")
library("scipy")
library(gaspr)

mPa_to_psia <- function(mPa) {
  psia <- mPa * 145.038
  return(psia)
}

# Convert temperature from Celsius to Fahrenheit
C_to_F <- function(T_C) {
  T_F <- (T_C * 9/5) + 32
  return(T_F)
}

#' Calculate Ppc using Ahmed method
#' 
#' This function calculates Ppc using the Ahmed method.
#' 
#' @param gamma_g The gamma_g value
#' @param yN2 The mole fraction of N2
#' @param yCO2 The mole fraction of CO2
#' @param yH2S The mole fraction of H2S
#' @return The calculated Ppc
calculate_Ppc_ahmed_method <- function(gamma_g, yN2, yCO2, yH2S) {
  Ppc <- 678 - 50 * (gamma_g - 0.5) - 206.7 * yN2 + 440 * yCO2 + 606.7 * yH2S
  return(Ppc)
}

#' Calculate Tpc using Ahmed method
#' 
#' This function calculates Tpc using the Ahmed method.
#' 
#' @param gamma_g The gamma_g value
#' @param yN2 The mole fraction of N2
#' @param yCO2 The mole fraction of CO2
#' @param yH2S The mole fraction of H2S
#' @return The calculated Tpc
calculate_Tpc_ahmed_method <- function(gamma_g, yN2, yCO2, yH2S) {
  Tpc <- 326 + 315.7 * (gamma_g - 0.5) - 240 * yN2 - 83.3 * yCO2 + 133.3 * yH2S
  return(Tpc)
}

calculate_zfactor_brill_and_beggs_method <- function(Tpr, Ppr) {
  F <- 0.3106 - 0.49 * Tpr + 0.1824 * Tpr ^ 2
  E <- 9 * (Tpr - 1)
  D <- 10 ^ F
  C <- 0.132 - 0.32 * log10(Tpr)
  B <- (0.62 - 0.23 * Tpr) * Ppr + (0.066 / (Tpr - 0.86) - 0.037) * Ppr ^ 2 + 0.32 * Ppr ^ 6 / 10 ^ E
  A <- 1.39 * (Tpr - 0.92) ^ 0.5 - 0.36 * Tpr - 0.10
  z <- A + (1 - A) / exp(B) + C * Ppr ^ D
  return(z)
}

calculate_zfactor_hall_yarborough_method <- function(Tpr, Ppr) {
  tr <- 1 / Tpr
  A <- 0.06125 * tr * exp(-1.2 * (1 - tr) ^ 2)
  B <- tr * (14.76 - 9.76 * tr + 4.58 * tr ^ 2)
  C <- tr * (90.7 - 242.2 * tr + 42.4 * tr ^ 2)
  D <- 2.18 + 2.82 * tr
  f <- function(Y) {
    return((Y + Y ^ 2 + Y ^ 3 - Y ^ 4) / (1 - Y) ^ 3 - A * Ppr - B * Y ^ 2 + C * Y ^ D)
  }
  sol <- uniroot(f, interval = c(0, 0.99999))
  Y <- sol$root
  z <- A * Ppr / Y
  return(z)
}

# Define constants
gamma_g <- 0.65  # specific gravity of gas
T <- 80  # temperature in deg C
P1 <- 5  # pressure in mPa
P2 <- 10  # pressure in mPa
P3 <- 50  # pressure in mPa

# Convert units
T <- C_to_F(T)
P1 <- mPa_to_psia(P1)
P2 <- mPa_to_psia(P2)
P3 <- mPa_to_psia(P3)

# Calculate pseudocritical properties
Ppc <- calculate_Ppc_ahmed_method(gamma_g, 0, 0, 0)
Tpc <- calculate_Tpc_ahmed_method(gamma_g, 0, 0, 0)
Tpr <- (T + 459.67) / Tpc
Ppr1 <- P1 / Ppc
Ppr2 <- P2 / Ppc
Ppr3 <- P3 / Ppc

# Calculate z factors using different methods
z1_1 <- calculate_zfactor_brill_and_beggs_method(Tpr, Ppr1)
z1_2 <- calculate_zfactor_hall_yarborough_method(Tpr, Ppr1)

z2_1 <- calculate_zfactor_brill_and_beggs_method(Tpr, Ppr2)
z2_2 <- calculate_zfactor_hall_yarborough_method(Tpr, Ppr2)

z3_1 <- calculate_zfactor_brill_and_beggs_method(Tpr, Ppr3)
z3_2 <- calculate_zfactor_hall_yarborough_method(Tpr, Ppr3)

# Print results
cat("at P = 5 mPa, z for Brill and Beggs method = ", z1_1, "\n")
cat("at P = 5 mPa, z for Hall Yarborough method = ", z1_2, "\n")
cat("at P = 10 mPa, z for Brill and Beggs method = ", z2_1, "\n")
cat("at P = 10 mPa, z for Hall Yarborough method = ", z2_2, "\n")
cat("at P = 50 mPa, z for Brill and Beggs method = ", z3_1, "\n")
cat("at P = 50 mPa, z for Hall Yarborough method = ", z3_2, "\n")