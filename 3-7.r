# Importing required libraries
library("MASS")
library("scipy")
library(ggplot2)

# Function to calculate J_prime
J_prime <- function(J_prime_i, pi, pe) {
  J_prime <- J_prime_i * (pe / pi)
  return(J_prime)
}

# Function to apply Fetkovich's method for future well performance
fetkovichs_method_for_future_well_performance <- function(pe, pwf, J_prime) {
  q <- J_prime * (pe^2 - pwf^2)
  return(q)
}

# Define constants
pi <- 3000 # psia
J_prime_0 <- 4 * 10 ** - 4 # stb/d/psi^2
pe_1 <- 2500 # psia
pe_2 <- 2000 # psia
pe_3 <- 1500 # psia
pe_4 <- 1000 # psia

# Create pwf vectors
pwf_1 <- 0:pe_1
pwf_2 <- 0:pe_2
pwf_3 <- 0:pe_3
pwf_4 <- 0:pe_4

# Calculate J_prime values
J_prime_1 <- J_prime(J_prime_0, pi, pe_1)
J_prime_2 <- J_prime(J_prime_0, pi, pe_2)
J_prime_3 <- J_prime(J_prime_0, pi, pe_3)
J_prime_4 <- J_prime(J_prime_0, pi, pe_4)

# Calculate q values
q_1 <- sapply(pwf_1, function(i) fetkovichs_method_for_future_well_performance(pe_1, i, J_prime_1))
q_2 <- sapply(pwf_2, function(i) fetkovichs_method_for_future_well_performance(pe_2, i, J_prime_2))
q_3 <- sapply(pwf_3, function(i) fetkovichs_method_for_future_well_performance(pe_3, i, J_prime_3))
q_4 <- sapply(pwf_4, function(i) fetkovichs_method_for_future_well_performance(pe_4, i, J_prime_4))

# Create plot
df <- data.frame(q = c(q_1, q_2, q_3, q_4), pwf = c(pwf_1, pwf_2, pwf_3, pwf_4), pe = rep(c(pe_1, pe_2, pe_3, pe_4), each = length(pwf_1)))
ggplot(df, aes(x = q, y = pwf, color = factor(pe))) + geom_line() + xlab("q (stb/d)") + ylab("Pwf (psi)") + labs(color = "pe") + theme_minimal()

