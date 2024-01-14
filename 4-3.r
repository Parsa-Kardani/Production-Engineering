# Importing required libraries
library("MASS")
library("scipy")

mass_for_1_stb <- function(gamma_o, gamma_w, gamma_g, wor, gor) {
  M <- 350.17 * (gamma_o + wor * gamma_w) + 0.0765 * gor * gamma_g
  return(M)
}

d_rho_v <- function(M, D, qo) {
  d_rho_v <- 1.4737 * 10 ^ -5 * M * qo * 12 / D 
  return(d_rho_v)
}

f2F <- function(d_rho_v) {
  f2F <- 10^(1.444 - 2.5 * log10(d_rho_v))
  return(f2F * 4)
}

guo_ghalambor <- function(theta, L, D, qg, gamma_g, qo, gamma_o, qw, gamma_w, qs, gamma_s, T_avg, A, fm, P_head) {
  theta <- (theta * pi) / 180
  cos_theta <- cos(theta)
  a <- ((0.0765 * qg * gamma_g) + (350 * qo * gamma_o) + (350 * qw * gamma_w) + (62.4 * qs * gamma_s)) / (4.07 * T_avg * qg) * cos_theta
  b <- (5.615 * qo + 5.615 * qw + qs) / (4.07 * T_avg * qg)
  c <- 0.00678 * T_avg * qg / A
  d <- 0.0016666 * (5.615 * qo + 5.615 * qw + qs) / A
  e <- fm / (2 * 32.17 * D) / cos_theta
  m <- c * d * e / (cos_theta + d^2 * e)
  n <- c^2 * e * cos_theta / (cos_theta + d^2 * e)^2
  
  P_head <- P_head * 144
  
  equation <- function(p) {
    return(b * p + (1 - 2 * b * m) / 2 * log(abs(((p + m)^2 + n) / ((P_head + m)^2 + n))) - (m + b * n - b * m^2) / sqrt(n) * (atan((p + m) / sqrt(n)) - atan((P_head + m) / sqrt(n))) - a * (1 + d^2 * e) * L)
  }
  
  p_solution <- uniroot(equation, c(0, 100))$root / 144  # Starting guess for p is 100
  return(as.numeric(p_solution))
}

# Define the given parameters
L <- 8000 # ft
theta <- 50 # degree
D <- 1.995 / 12 # ft
qg <- 500000 # scf/d
gamma_g <- 0.75
qo <- 2000 # stb/d
gamma_o <- 0.85
qw <- 500 # stb/d
gamma_w <- 1.05
qs <- 4 # ft3/d
gamma_s <- 2.65
T_head <- 100 # F
T_bottom <- 170 # F
P_head <- 500 # psia

# Calculate the average temperature in Rankine
T_avg <- (T_head + T_bottom) / 2 + 460

# Calculate the cross-sectional area of the wellbore
A <- 3.14159 * (D * 12) ^ 2 / 4

# Calculate the mass for 1 stb
M <- mass_for_1_stb(gamma_o, gamma_w, gamma_g, qw / qo, qg / qo)

# Calculate the change in density with respect to the change in velocity
d_rho_v <- d_rho_v(M, D * 12, qo)

# Calculate the friction factor
fm <- f2F(d_rho_v) * 4

# Calculate the pressure at the bottom of the well
p <- guo_ghalambor(theta, L, D, qg, gamma_g, qo, gamma_o, qw, gamma_w, qs, gamma_s, T_avg, A, fm, P_head)

# Print the result
cat("The pressure at the bottom of the well is", format(p, digits = 2), "psia\n")
