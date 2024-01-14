# Importing required libraries
library("MASS")
library("scipy")
library(ggplot2)

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

# Convert temperature from Celsius to Fahrenheit
C_to_F <- function(T_C) {
  T_F <- (T_C * 9/5) + 32
  return(T_F)
}

gamma_g <- 0.65  # specific gravity of gas
T <- 110  # temperature in C
T <- C_to_F(T)
Ppc <- calculate_Ppc_ahmed_method(gamma_g, 0, 0, 0)
Ppr <- seq(14, 8000) / Ppc

plot(seq(14, 8000), Ppr, type = "l", xlab = "Pressure (psi)", ylab = "Ppr", main = "Pseudo-Reduced Pressure vs Pressure")