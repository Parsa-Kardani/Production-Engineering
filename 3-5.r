# Importing required libraries
library("MASS")
library("scipy")
library(ggplot2)

qmax_for_test_points_two_phase <- function(q1, pwf1, P) {
  qmax <- q1 / (1 - 0.2 * (pwf1 / P) - 0.8 * (pwf1 / P) ^ 2)
  return(qmax)
}

J_star_using_qmax <- function(qmax, p) {
  J_star <- qmax / p * 1.8
  return(J_star)
}

fetkovichs_equation <- function(q1, pwf1, q2, pwf2, P, pwf) {
  n <- log10(q1 / q2) / log10((P ^ 2 - pwf1 ^ 2) / (P ^ 2 - pwf2 ^ 2))
  C <- q1 / (P ^ 2 - pwf1 ^ 2) ^ n
  q <- C * (P ^ 2 - pwf ^ 2) ^ n
  return(q)
}

# Define the given parameters
P <- 3500 # psia
pwf1 <- 2500 # psia
q1 <- 600 # stb/d
pwf2 <- 1500 # psia
q2 <- 900 # stb/d

J_star <- J_star_using_qmax(qmax, p)
qmax <-qmax_for_test_points_two_phase(q1, pwf1, P)

# Create a sequence of pwf values
pwf <- seq(0, P, by = 1)

# Initialize empty lists to store calculated q values
q_fetkovich <- numeric(length(pwf))
q_vogel <- numeric(length(pwf))

# Calculate q values using Fetkovich's and Vogel's equations
for (i in 1:length(pwf)) {
  q_fetkovich[i] <- fetkovichs_equation(q1, pwf1, q2, pwf2, P, pwf[i])
  q_vogel[i] <- vogels_equation_for_q(J_star, P, P, pwf[i])
}

# Create a data frame for plotting
df <- data.frame(pwf = pwf, q_fetkovich = q_fetkovich, q_vogel = q_vogel)

# Plot the results
ggplot(df, aes(x = q_fetkovich, y = pwf)) +
  geom_line(color = "blue", size = 1, linetype = "solid") +
  geom_line(aes(x = q_vogel), color = "red", size = 1, linetype = "solid") +
  xlab("q (stb/d)") +
  ylab("Pwf (psi)") +
  scale_color_manual(values = c("blue", "red"), labels = c("fetkovich", "vogel")) +
  theme_minimal()