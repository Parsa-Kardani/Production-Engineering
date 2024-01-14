# Importing required libraries
library("MASS")
library("scipy")
library(ggplot2)

vogels_equation_for_q <- function(J_star, p, pb, pwf) {
  if (pwf >= pb) {
    q <- J_star * (p - pwf)
  } else {
    q <- J_star * (p - pb) + J_star * pb / 1.8 * (1 - 0.2 * (pwf / pb) - 0.8 * (pwf / pb) ^ 2)
  }
  return(q)
}

# Function to calculate Inflow Performance Relationship (IPR) for single liquid phase radial pseudo steady state
IPR_for_single_liquid_phase_radial_pseudo_steady_state <- function(k, h, Bo, mu_o, rw, re, s) {
  J_star <- k * h / (141.2 * Bo * mu_o * (log(re / rw) - 0.75 + s))
  return(J_star)
}

# Define input parameters
phi <- 0.2
k <- 80  # md
h <- 55  # ft
P <- 4500  # psia
pb <- 4500  # psia
Bo <- 1.1  # rb/stb
mu_o <- 1.8  # cp
ct <- 0.000013  # psi^-1
A <- 640  # acres
re <- 2980  # ft
rw <- 0.328  # ft
s <- 2

# Calculate J_star using IPR_for_single_liquid_phase_radial_pseudo_steady_state function
J_star <- IPR_for_single_liquid_phase_radial_pseudo_steady_state(k, h, bo, mu_o, rw, re, P, S)

# Calculate pwf and q values
pwf <- seq(P, 0, by = -1)
q <- sapply(pwf, function(i) vogels_equation_for_q(IPR_for_single_liquid_phase_radial_pseudo_steady_state(k, h, Bo, mu_o, rw, re, s), P, pb, i))

# Create plot
df <- data.frame(q = q, pwf = pwf)
ggplot(df, aes(x = q, y = pwf)) +
  geom_line() +
  xlab('q (stb/d)') +
  ylab('Pwf (psi)')