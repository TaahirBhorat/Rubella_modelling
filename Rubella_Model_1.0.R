library(deSolve)
library(ggplot2)
library(reshape2)

# Function to define the rubella model
rubella_model <- function(time, state, parameters) {
  n_age <- parameters$n_age
  n_compartments <- 5  # M, S, I, R, V
  
  #  compartment states by age group
  M <- state[1:n_age]
  S <- state[(n_age + 1):(2 * n_age)]
  I <- state[(2 * n_age + 1):(3 * n_age)]
  R <- state[(3 * n_age + 1):(4 * n_age)]
  V <- state[(4 * n_age + 1):(5 * n_age)]
  
  with(as.list(parameters), {
    # Initialize FOI per each age group
    phi <- numeric(n_age)
    for (a in 1:n_age) {
      # FOI calculation
      phi[a] <- 1 - exp(-sum(beta[a, ] * I) / sum(M + S + I + R + V))
    }
    
    dM <- numeric(n_age)
    dS <- numeric(n_age)
    dI <- numeric(n_age)
    dR <- numeric(n_age)
    dV <- numeric(n_age)
    
    for (a in 1:n_age) {
      # Aging transitions
      aging_in_M <- ifelse(a > 1, aging_rate[a - 1] * M[a - 1], 0)
      aging_in_S <- ifelse(a > 1, aging_rate[a - 1] * S[a - 1], 0)
      aging_in_I <- ifelse(a > 1, aging_rate[a - 1] * I[a - 1], 0)
      aging_in_R <- ifelse(a > 1, aging_rate[a - 1] * R[a - 1], 0)
      aging_in_V <- ifelse(a > 1, aging_rate[a - 1] * V[a - 1], 0)
      
      # Compartment Transitions
      dM[a] <- -da[a] * M[a] - mu[a] * M[a] + aging_in_M - aging_rate[a] * M[a]
      dS[a] <- da[a] * M[a] - phi[a] * (1 - va[a]) * S[a] - mu[a] * S[a] + aging_in_S - aging_rate[a] * S[a] - va[a] * S[a]
      dI[a] <- phi[a] * (1 - va[a]) * S[a] - (gamm[a] + mu[a]) * I[a] + aging_in_I - aging_rate[a] * I[a]
      dR[a] <- gamm[a] * I[a] - mu[a] * R[a] + aging_in_R - aging_rate[a] * R[a]
      dV[a] <- va[a] * S[a] - mu[a] * V[a] + aging_in_V - aging_rate[a] * V[a]
    }
    
    list(c(dM, dS, dI, dR, dV))
  })
}

# Initialize parameters with aging and survival rates
initialize_parameters <- function(n_age) {
  list(
    beta = matrix(runif(n_age * n_age, min = 0.01, max = 0.05), nrow = n_age, ncol = n_age),
    da = rep(0.1, n_age),
    va = rep(0.05, n_age),
    gamm = rep(0.1, n_age),
    mu = rep(0.01, n_age),
    aging_rate = rep(1/12, n_age),
    n_age = n_age
  )
}

# Initialize states with compartments and age stratification
initialize_state <- function(n_age) {
  initial_state <- c(
    M = rep(50, n_age),
    S = rep(1000, n_age),
    I = rep(5, n_age),
    R = rep(100, n_age),
    V = rep(200, n_age)
  )
  return(initial_state)
}

# 4 age group test
n_age <- 4
parameters <- initialize_parameters(n_age)
initial_state <- initialize_state(n_age)

# Daily for a year
times <- seq(0, 365, by = 1)
results <- ode(y = initial_state, times = times, func = rubella_model, parms = parameters)

results_df <- as.data.frame(results)

# Melting and plotting
results_long <- melt(results_df, id.vars = "time")

ggplot(data = results_long, aes(x = time, y = value, color = variable)) +
  geom_line() +
  labs(title = "Rubella Model with Accurate Transitions", x = "Time (days)", y = "Population", color = "Compartment") +
  theme_minimal()

