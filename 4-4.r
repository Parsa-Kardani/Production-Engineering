# Importing required libraries
library("MASS")
library("scipy")
library(ggplot2)

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
  mu_r <- a0 + a1 * Ppr + a2 * Ppr ** 2 + a3 * Ppr ** 3 + Tpr *(a4 + a5 * Ppr + a6 * Ppr ** 2 + a7 * Ppr ** 3) + Tpr ** 2 * (a8 + a9 * Ppr + a10 * Ppr ** 2 + a11 * Ppr ** 3) + Tpr ** 3 * (a12 + a13 * Ppr + a14 * Ppr ** 2 + a15 * Ppr ** 3)
  mu_g <- mu1 * exp(mu_r) / Tpr
  return(mu_g)
}


Hagedorn_Brown_Correlation <- function(L, p_head, T_head, T_bot, dz, gamma_g, gamma_l, ift, D, m_t, qg, A, mu_l, u_sl) {
  p <- c(p_head)
  depth <- c(0)
  Tpc <- calculate_Tpc_ahmed_method(gamma_g, 0, 0, 0)
  Ppc <- calculate_Ppc_ahmed_method(gamma_g, 0, 0, 0)
  for (i in 1:(L / dz)) {
    eps <- 0.0006
    T <- depth[i] * (T_bot - T_head) / L + T_head 
    Tpr <- (T + 460) / Tpc
    Ppr <- p[i] / Ppc
    z <- calculate_zfactor_brill_and_beggs_method(Tpr, Ppr)
    mu_g <- calculate_mu_carr_et_al_method(gamma_g, 0, 0, 0, T, Ppr, Tpr)
    u_sg <- 1 / A * qg * z * (460 + T) / (460+60) * (14.7 / p[i]) /86400
    N_vl <- 1.938 * u_sl * (62.4 * gamma_l / ift) ** 0.25
    N_vg <- 1.938 * u_sg * (62.4 * gamma_l / ift) ** 0.25
    N_d <- 120.872 * D * sqrt(62.4 * gamma_l / ift)
    N_l <- 0.15726 * mu_l * (1 / (62.4 * gamma_l * ift ** 3)) ** 0.25
    x1 <- log10(N_l) + 3
    Y <- - 2.69851 + 0.15840954 * x1 - 0.55099756 * x1 ** 2 + 0.54784917 * x1 ** 3 - 0.12194578 * x1 ** 4
    Cn_l <- 10 ** Y
    x2 <- N_vl / N_vg ** 0.575 * (p[i] / 14.7) ** 0.1 * Cn_l / N_d
    y_L_psi <- - 0.10306578 + 0.617774 * (log10(x2) + 6) - 0.632946 * (log10(x2) + 6) ** 2 + 0.29598 * (log10(x2) + 6) ** 3 - 0.0401 * (log10(x2) + 6) ** 4 
    x3 <- 0.012
    # x3 <- N_vg * N_l ** 0.38 / N_d ** 2.14
    psi <- 0.91162574 - 4.82175636 * x3 + 1232.25036621 * x3 ** 2 - 22253.57617 * x3 ** 3 + 116174.28125 * x3 ** 4
    y_L <- psi * y_L_psi
    Re <- 0.022 * m_t / ((D * 12) * mu_l ** y_L * mu_g ** (1 - y_L))
    f <- 1 / (- 4 * log10(eps / 3.7065 - 5.0452 / Re * log10(eps ** 1.1098 / 2.8257 + (7.149 / Re) ** 0.8981))) ** 2
    rho_g <- 28.97 * gamma_g * p[i] / z / 10.73 / (460 + T)
    rho_avg <- y_L * gamma_l * 62.4 + (1 - y_L) * rho_g
    p_gradient <- 1 / 144 * (rho_avg + f * m_t ** 2 / 7.413 / 10000000000 / D ** 5 / rho_avg)
    p <- c(p, p[i] + p_gradient * dz)
    depth <- c(depth, depth[i] + dz)
  }
  return(list(p = p, depth = depth))
}

# Define input parameters
L <- 6000 # ft
D <- 1.995 / 12 # ft
api <- 30 # API
mu <- 2 # cp
glr <- 500 # scf/bbl
gamma_g <- 0.65
p_head <- 100 # psia
T_head <- 80 # F
T_bot <- 140 # F
liquid_rate <- 1500 # bbl/day
water_cut <- 0.2 # fraction
ift <- 30 # dynes/cm
gamma_w <- 1.05 

# Calculate intermediate variables
qw <- liquid_rate * water_cut
gamma_o <- calculate_gamma_o_using_api(api)  # Assuming a function calculate_gamma_o_using_api exists
gamma_l <- ((liquid_rate - qw) * gamma_o + qw * gamma_w) / liquid_rate
qg <- glr * liquid_rate
A <- 3.141592653 * (D) ^ 2 / 4
dz <- 100
u_sl <- liquid_rate * 5.615 / 86400 / A
mu_l <- (mu * (liquid_rate - qw) + 0.5 * qw) / liquid_rate
m_t <- gamma_l * 62.4 * liquid_rate * 5.615 + 0.0765 * gamma_g * qg


# Call Hagedorn_Brown_Correlation function to get pressure and depth
result <- Hagedorn_Brown_Correlation(L, p_head, T_head, T_bot, dz, gamma_g, gamma_l, ift, D, m_t, qg, A, mu_l, u_sl)
p <- result$p
depth <- result$depth


# Plot the results
ggplot(data = data.frame(p = p, depth = depth), aes(x = p, y = depth)) +
  geom_line() +
  scale_y_reverse() +
  xlab("Pressure (psia)") +
  ylab("Depth (ft)") +
  theme_minimal()


