# Importing required libraries
library("MASS")
library("scipy")

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

calculate_mu_carr_et_al_method <- function(gamma_g, yN2, yCO2, yH2S, t, Ppr, Tpr) {
  mu1_hc <- 8.188 * 10 ** (-3) - 6.15 * 10 ** (-3) * log10(gamma_g) + (1.709 * 10 ** (-5) - 2.062 * 10 ** (-6) * gamma_g) * t
  mu1_n2 <- (9.59 * 10 **(-3) + 8.48 * 10 ** (-3) * log(gamma_g)) * yN2
  mu1_c2 <- (6.24 * 10 ** (-3) + 9.08 * 10 ** (-3) * log(gamma_g)) * yCO2
  mu1_h2s <- (3.73 * 10 ** (-3) + 8.49 * 10 ** (-3) * log(gamma_g)) * yH2S
  mu1 <- mu1_hc + mu1_n2 + mu1_c2 + mu1_h2s
  a0 <- -2.46211820
  a1 <- 2.97054714
  a2 <- -0.28626405
  a3 <- 0.00805420
  a4 <- 2.80860949
  a5 <- -3.49803305
  a6 <- 0.36037302
  a7 <- -0.01044324
  a8 <- -0.79338568
  a9 <- 1.39643306
  a10 <- -0.14914493
  a11 <- 0.00441016
  a12 <- 0.08393872
  a13 <- -0.18640885
  a14 <- 0.02033679
  a15 <- -0.00060958
  mu_r <- a0 + a1 * Ppr + a2 * Ppr ** 2 + a3 * Ppr ** 3 + Tpr * (a4 + a5 * Ppr + a6 * Ppr ** 2 + a7 * Ppr ** 3) + Tpr ** 2 * (a8 + a9 * Ppr + a10 * Ppr ** 2 + a11 * Ppr ** 3) + Tpr ** 3 * (a12 + a13 * Ppr + a14 * Ppr ** 2 + a15 * Ppr ** 3)
  mu_g <- mu1 * exp(mu_r) / Tpr
  return(mu_g)
}


gamma_g <- 0.7 # Specific gravity of gas
T <- 90 # Temperature in C
P1 <- 1 # Pressure in mPa
P2 <- 5 # Pressure in mPa
P3 <- 10 # Pressure in mPa
P4 <- 50 # Pressure in mPa

T <- C_to_F(T)
P1 <- mPa_to_psia(P1)
P2 <- mPa_to_psia(P2)
P3 <- mPa_to_psia(P3)
P4 <- mPa_to_psia(P4)

Ppc <- calculate_Ppc_ahmed_method(gamma_g, 0, 0, 0)
Tpc <- calculate_Tpc_ahmed_method(gamma_g, 0, 0, 0)
Tpr <- (T + 459.67) / Tpc
Ppr1 <- P1 / Ppc
Ppr2 <- P2 / Ppc
Ppr3 <- P3 / Ppc
Ppr4 <- P4 / Ppc

mu_g1 <- calculate_mu_carr_et_al_method(gamma_g, 0, 0, 0, T, Ppr1, Tpr)
mu_g2 <- calculate_mu_carr_et_al_method(gamma_g, 0, 0, 0, T, Ppr2, Tpr)
mu_g3 <- calculate_mu_carr_et_al_method(gamma_g, 0, 0, 0, T, Ppr3, Tpr)
mu_g4 <- calculate_mu_carr_et_al_method(gamma_g, 0, 0, 0, T, Ppr4, Tpr)

cat("gas viscosity at 1 mPa is", mu_g1, "cp\n")
cat("gas viscosity at 5 mPa is", mu_g2, "cp\n")
cat("gas viscosity at 10 mPa is", mu_g3, "cp\n")
cat("gas viscosity at 50 mPa is", mu_g4, "cp\n")
