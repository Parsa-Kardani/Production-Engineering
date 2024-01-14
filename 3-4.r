# Importing required libraries
library("MASS")
library("scipy")
library(ggplot2)

IPR_using_test_points <- function(P, pb, pwf1, q1) {
  if (pwf1 >= pb) {
    J_star <- q1 / (P - pwf1)
  } else {
    J_star <- q1 / ((P - pb) + pb / 1.8 * (1 - 0.2 * (pwf1 / pb) - 0.8 * (pwf1 / pb) ^ 2))
  }
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

P <- 5500 # psia
pb <- 3500 # psia
pwf1_A <- 4000 # psia
q1_A <- 400 # stb/d
pwf1_B <- 2000 # psia
q1_B <- 1000 # stb/d

J_star_A <- IPR_using_test_points(P, pb, pwf1_A, q1_A)
J_star_B <- IPR_using_test_points(P, pb, pwf1_B, q1_B)

pwf <- seq(0, P, by = 1)
q_A <- sapply(pwf, function(i) vogels_equation_for_q(J_star_A, P, pb, i))
q_B <- sapply(pwf, function(i) vogels_equation_for_q(J_star_B, P, pb, i))

plot(q_A, pwf, type = "l", xlab = "q (stb/d)", ylab = "Pwf (psi)", col = "blue", ylim = c(0, P))
lines(q_B, pwf, col = "red")
legend("topright", legend = c("well A", "well B"), col = c("blue", "red"), lty = 1)