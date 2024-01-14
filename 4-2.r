# Importing required libraries
library("MASS")
library("scipy")

mass_for_1_stb <- function(gamma_o, gamma_w, gamma_g, wor, gor) {
  M <- 350.17 * (gamma_o + wor * gamma_w) + 0.0765 * gor * gamma_g
  return(M)
}

calculate_Rs <- function(gamma_g, api, p, t) {
  Rs <- gamma_g * (p / 18 * 10^(0.0125 * api) / 10^(0.00091 * t)) ^ 1.2048
  return(Rs)
}

calculate_Bo_Standing_method <- function(Rs, gamma_o, gamma_g, t) {
  Bo <- 0.9759 + 0.00012 * (Rs * ((gamma_g/gamma_o) ^ 0.5) + 1.25 * t) ^ 1.2
  return(Bo)
}

volume_for_1_stb <- function(Bo, wor, bw, gor, Rs, p, T, z) {
  V <- 5.615 * (Bo + wor * bw) + (gor - Rs) * (14.7 / p) * ((T + 460) / 520) * z
  return(V)
}

d_rho_v <- function(M, D, qo) {
  d_rho_v <- 1.4737 * 10^-5 * M * qo * 12 / D 
  return(d_rho_v)
}

f2F <- function(d_rho_v) {
  f2F <- 10^(1.444 - 2.5 * log10(d_rho_v))
  return(f2F * 4)
}

calculate_gamma_o_using_api <- function(api) {
  gamma_o <- (141.5/(api + 131.5))
  return(gamma_o)
}

calculate_Ppc_ahmed_method <- function(gamma_g, yN2, yCO2, yH2S) {
  Ppc <- 678 - 50 * (gamma_g - 0.5) - 206.7 * yN2 + 440 * yCO2 + 606.7 * yH2S
  return(Ppc)
}

calculate_Tpc_ahmed_method <- function(gamma_g, yN2, yCO2, yH2S) {
  Tpc <- 326 + 315.7 * (gamma_g - 0.5) - 240 * yN2 - 83.3 * yCO2 + 133.3 * yH2S
  return(Tpc)
}

friction_term <- function(f2F, qo, M, D) {
  k <- f2F * qo^2 * M^2 / (7.4137 * 10^10 * D^5)
  return(k)
}

calculate_zfactor_brill_and_beggs_method <- function(Tpr, Ppr) {
  F <- 0.3106 - 0.49 * Tpr + 0.1824 * Tpr^2
  E <- 9 * (Tpr -1)
  D <- 10^F
  C <- 0.132 - 0.32 * log10(Tpr)
  B <- (0.62 - 0.23 * Tpr) * Ppr + (0.066 / (Tpr - 0.86) - 0.037) * Ppr^2 + 0.32 * Ppr^6 / 10^E
  A <- 1.39 * (Tpr - 0.92)^0.5 - 0.36 * Tpr - 0.10
  z <- A + (1 - A) / exp(B) + C * Ppr^D
  return(z)
}

calculating_rho_at_any_point <- function(api, p, t, gamma_g, wor, bw, gor, M) {
  gamma_o <- calculate_gamma_o_using_api(api)
  Rs <- calculate_Rs(gamma_g, api, p, t)
  Bo <- calculate_Bo_Standing_method(Rs, gamma_o, gamma_g, t)
  Ppc <- calculate_Ppc_ahmed_method(gamma_g, 0, 0, 0)
  Tpc <- calculate_Tpc_ahmed_method(gamma_g, 0, 0, 0)
  Tpr <- (t + 460) / Tpc
  Ppr <- p / Ppc
  z <- calculate_zfactor_brill_and_beggs_method(Tpr, Ppr)
  V <- volume_for_1_stb(Bo, wor, bw, gor, Rs, p, t, z)
  rho <- M / V
  return(rho)
}

finding_bhp <- function(P_wellhead, rho_wellhead, k, L, api, bottomhole_temperature, gamma_g, wor, bw, gor, M) {
  pbh <- P_wellhead
  error_h <- 1
  for (i in 1:10) {
    rho_bottomhole <- calculating_rho_at_any_point(api, pbh, bottomhole_temperature, gamma_g, wor, bw, gor, M)    
    rho_avg <- (rho_wellhead + rho_bottomhole) / 2
    pbh <- P_wellhead + (rho_avg + k / rho_avg) * L / 144
    error_h <- 144 * (pbh - P_wellhead) / (rho_avg + k / rho_avg) - L
  }
  return(pbh)
}


# Given values
D <- 1.66
P_wellhead <- 300 # psia
production_rate <- 2000 # STB/day
gas_liquid_ratio <- 800 # scf/STB
wc <- 0.30
api <- 40 # api degree
gamma_w <- 1.05
gamma_g <- 0.7
bw <- 1.2 # bbl/STB
wellhead_temperature <- 100 # F
L <- 8000 # ft
bottomhole_temperature <- 170 # F

# Calculations
gor <- gas_liquid_ratio / (1 - wc)
wor <- wc / (1 - wc)
qo <- production_rate * (1 - wc)

gamma_o <- calculate_gamma_o_using_api(api)
M <- mass_for_1_stb(gamma_o, gamma_w, gamma_g, wor, gor)
Rs_wellhead <- calculate_Rs(gamma_g, api, P_wellhead, wellhead_temperature)
Bo_wellhead <- calculate_Bo_Standing_method(Rs_wellhead, gamma_o, gamma_g, wellhead_temperature)
Ppc <- calculate_Ppc_ahmed_method(gamma_g, 0, 0, 0)
Tpc <- calculate_Tpc_ahmed_method(gamma_g, 0, 0, 0)
Tpr <- (wellhead_temperature + 460) / Tpc
Ppr <- P_wellhead / Ppc
z_wellhead <- calculate_zfactor_brill_and_beggs_method(Tpr, Ppr)
V_wellhead <- volume_for_1_stb(Bo_wellhead, wor, bw, gor, Rs_wellhead, P_wellhead, wellhead_temperature, z_wellhead)
rho_wellhead <- M / V_wellhead
d_rho_v <- d_rho_v(M, D, qo)
f2F <- f2F(d_rho_v)
k <- friction_term(f2F, qo, M, D / 12)
bhp <- finding_bhp(P_wellhead, rho_wellhead, k, L, api, bottomhole_temperature, gamma_g, wor, bw, gor, M)

cat("The bottomhole using Poettmann-Carpenter method is:", round(bhp, 2), "psia\n")