# Importing required libraries
library("MASS")
library("scipy")
library(ggplot2)

# Function to calculate future IPR using Vogel's method
future_IPR_vogels_method <- function(J_star_p, Bo_p, mu_o_p, kro_p, Bo_f, mu_o_f, kro_f) {
  J_star_f <- J_star_p * (kro_f / (Bo_f * mu_o_f)) / (kro_p / (Bo_p * mu_o_p))
  return(J_star_f)
}

# Function to calculate future production rate using Vogel's equation
vogels_equation_for_q_future <- function(J_star_f, Pf, pwf) {
  q <- J_star_f * Pf / 1.8 * (1 - 0.2 * (pwf / Pf) - 0.8 * (pwf / Pf) ^ 2)
  return(q)
}

vogels_equation_for_q <- function(J_star, p, pb, pwf) {
  if (pwf >= pb) {
    q <- J_star * (p - pwf)
  } else {
    q <- J_star * (p - pb) + J_star * pb / 1.8 * (1 - 0.2 * (pwf / pb) - 0.8 * (pwf / pb) ^ 2)
  }
  return(q)
}

# Define the parameters
Pp <- 2200 # psia
Pf <- 1500 # psia
J_star_p <- 1.25 # stb/d/psi
mu_o_p <- 3.55 # cp
mu_o_f <- 3.85 # cp
Bo_p <- 1.2 # rb/stb
Bo_f <- 1.15 # rb/stb
kr_p <- 0.82 # md
kr_f <- 0.65 # md

# Calculate future IPR using Vogel's method
J_star_f <- future_IPR_vogels_method(J_star_p, Bo_p, mu_o_p, kr_p, Bo_f, mu_o_f, kr_f)

# Create lists for pwf_p and pwf_f
pwf_p <- 0:Pp
pwf_f <- 0:Pf

# Calculate q_p and q_f
q_p <- sapply(pwf_p, function(i) vogels_equation_for_q(J_star_p, Pp, Pp, i))
q_f <- sapply(pwf_f, function(i) vogels_equation_for_q_future(J_star_f, Pf, i))

# Plot the results
plot(q_p, pwf_p, type='l', xlab='q (stb/d)', ylab='Pwf (psi)', col='blue', main='Present vs Future IPR')
lines(q_f, pwf_f, col='red')
legend('topright', legend=c('present', 'future'), col=c('blue', 'red'), lty=1)