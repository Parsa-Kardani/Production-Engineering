# Importing required libraries
library("MASS")
library("scipy")
library(ggplot2)


# Function to calculate Inflow Performance Relationship (IPR) for single liquid phase radial transient
IPR_for_single_liquid_phase_radial_transient <- function(k, h, phi, ct, Bo, mu_o, rw, t, s) {
  J_star <- k * h / (162.6 * Bo * mu_o * (log10(t) + log10(k / (phi * mu_o * ct * rw^2)) - 3.23 + 0.87 * s))
  return(J_star)
}

# Function to calculate Inflow Performance Relationship (IPR) for single liquid phase radial steady state
IPR_for_single_liquid_phase_radial_steady_state <- function(k, h, Bo, mu_o, rw, re, s) {
  J_star <- k * h / (141.2 * Bo * mu_o * (log(re / rw) + s))
  return(J_star)
}

# Function to calculate Inflow Performance Relationship (IPR) for single liquid phase radial pseudo steady state
IPR_for_single_liquid_phase_radial_pseudo_steady_state <- function(k, h, Bo, mu_o, rw, re, s) {
  J_star <- k * h / (141.2 * Bo * mu_o * (log(re / rw) - 0.75 + s))
  return(J_star)
}

vogels_equation_for_q <- function(J_star, p, pb, pwf) {
  if (pwf >= pb) {
    q <- J_star * (p - pwf)
  } else {
    q <- J_star * (p - pb) + J_star * pb / 1.8 * (1 - 0.2 * (pwf / pb) - 0.8 * (pwf / pb) ^ 2)
  }
  return(q)
}

# Define input parameters
phi <- 0.25
k <- 10 # md
h <- 50 # ft
P <- 5000 # psia
pb <- 100 # psia
bo <- 1.2 # rb/stb
mu_o <- 1.5 # cp
ct <- 0.0000125 # psi^-1
rw <- 0.328 # ft
S <- 5
t <- 30 * 24

# Calculate J_star1 using IPR_for_single_liquid_phase_radial_transient function
J_star1 <- IPR_for_single_liquid_phase_radial_transient(k, h, phi, ct, bo, mu_o, rw, t, S)

# Calculate q1 using vogels_equation_for_q function
pwf1 <- seq(P + 1)
q1 <- sapply(pwf1, function(i) vogels_equation_for_q(J_star1, P, pb, i))

# Calculate J_star2 using IPR_for_single_liquid_phase_radial_steady_state function
J_star2 <- IPR_for_single_liquid_phase_radial_steady_state(k, h, bo, mu_o, rw, re , P, S)

# Calculate q2 using vogels_equation_for_q function
pwf2 <- seq(P + 1)
q2 <- sapply(pwf2, function(i) vogels_equation_for_q(J_star2, P, pb, i))

# Calculate J_star3 using IPR_for_single_liquid_phase_radial_pseudo_steady_state function
J_star3 <- IPR_for_single_liquid_phase_radial_pseudo_steady_state(k, h, bo, mu_o, rw, re , P, S)

# Calculate q3 using vogels_equation_for_q function
pwf3 <- seq(P + 1)
q3 <- sapply(pwf3, function(i) vogels_equation_for_q(J_star3, P, pb, i))

# Create plots
par(mfrow=c(1,3))
plot(q1, pwf1, xlab='q (stb/d)', ylab='Pwf (psi)', main='transient flow')
plot(q2, pwf2, xlab='q (stb/d)', ylab='Pwf (psi)', main='steady state flow')
plot(q3, pwf3, xlab='q (stb/d)', ylab='Pwf (psi)', main='pseudo steady state flow')


