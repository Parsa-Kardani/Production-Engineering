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
phi <- 0.25
k <- 100 # md
h <- 55 # ft
P <- 5000 # psia
pb <- 3000 # psia
Bo <- 1.2 # rb/stb
mu_o <- 1.8 # cp
ct <- 0.000013 # psi^-1
A <- 640 # acres
re <- 2980 # ft
rw <- 0.328 # ft
s <- 5.5

# Calculate J_star using relevant function
J_star <- IPR_for_single_liquid_phase_radial_pseudo_steady_state(k, h, Bo, mu_o, rw, re, s)

# Generate pwf values
pwf <- seq(0, P, by = 1)

# Calculate q using Vogel's equation
q <- sapply(pwf, function(p) vogels_equation_for_q(J_star, P, pb, p))

# Create plot
df <- data.frame(q = q, pwf = pwf)
ggplot(df, aes(x = q, y = pwf)) +
  geom_line() +
  xlab('q (stb/d)') +
  ylab('Pwf (psi)')