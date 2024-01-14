# Importing required libraries
library("MASS")
library("scipy")

# Function to calculate specific gravity of oil using API gravity
calculate_gamma_o_using_api <- function(api) {
  gamma_o <- (141.5 / (api + 131.5))
  return(gamma_o)
}

calculate_rho_using_gamma_o <- function(gamma_o) {
  rho <- 62.4 * gamma_o
  return(rho)
}

dz_inclinced_wellbore <- function(L, alpha) {
  alpha <- (pi/180) * alpha
  dz <- L * cos(alpha)
  return(dz)
}

u_in_tubing <- function(q, D) {
  u <- 4 * q * 5.615 / (86400 * pi * D^2) 
  return(u)
}

reynolds_number <- function(rho, q, D, mu) {
  Re <- 1.48 * q * rho / (mu * D)
  return(Re)
}

friction_factor <- function(Re, delta, d) {
  if (Re < 2000) {
    f_F <- 16 / Re
  }
  if (Re > 2100) {
    eps <- delta / d
    f_F <- 1 / (-4 * log10(eps / 3.7065 - 5.0452 / Re * log10(eps^1.1098 / 2.8257 + (7.149 / Re)^0.8981)))^2
  }
  return(f_F)
}

pressure_drop_in_single_phase_flow <- function(rho, dz, du, u, L, D, f_F) {
  gc <- 32.174
  g <- gc
  dp <- g / gc * rho * dz + rho * du^2 / (2 * gc) + rho * L / D * f_F * u^2 /gc
  dp <- dp / 144
  return(dp)
}

# Input parameters
q <- 1000 # bbl/day
api <- 16 # degree API
mu <- 5 # cp
alpha <- 3 # degree
L <- 1000 # ft
d <- 2.259 # in
D <- d / 12 # ft
delta <- 0.001

# Calculate gamma_o
gamma_o <- calculate_gamma_o_using_api(api)

# Calculate rho
rho <- calculate_rho_using_gamma_o(gamma_o)

# Calculate dz
dz <- dz_inclinced_wellbore(L, alpha)

# Calculate u
u <- u_in_tubing(q, D)

# Calculate Reynolds number
Re <- reynolds_number(rho, q, d, mu)

# Calculate friction factor
f_F <- friction_factor(Re, delta, d)

# Calculate pressure drop
delta_p <- pressure_drop_in_single_phase_flow(rho, dz, 0, u, L, D, f_F)

# Print the result
cat("Pressure drop over 1000 ft of tubing is", round(delta_p, 2), "psi\n")