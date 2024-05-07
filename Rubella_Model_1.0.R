library(deSolve)

# Define the rubella model function
rubella_model <- function(time, state, parameters) {
  n_age <- length(state) / 5  # Assuming 5 compartments
  M <- state[1:n_age]
  S <- state[(n_age + 1):(2 * n_age)]
  I <- state[(2 * n_age + 1):(3 * n_age)]
  R <- state[(3 * n_age + 1):(4 * n_age)]
  V <- state[(4 * n_age + 1):(5 * n_age)]
  
  with(as.list(parameters), {
    # Calculate the force of infection per age group
    lambda <- numeric(n_age)
    for (a in 1:n_age) {
      # Compute the total force of infection for age group 'a'
      lambda[a] <- sum(beta[a, ] * I) / sum(S + M + I + R + V)
    }
    
    # Initialize the differential change vectors
    dM <- numeric(n_age)
    dS <- numeric(n_age)
    dI <- numeric(n_age)
    dR <- numeric(n_age)
    dV <- numeric(n_age)
    
    # Equations for each age group
    for (a in 1:n_age) {
      dM[a] <- -mu[a] * M[a] - lambda[a] * M[a] + aging_rate[a] * (ifelse(a > 1, M[a - 1], 0))
      dS[a] <- -mu[a] * S[a] - lambda[a] * S[a] + aging_rate[a] * (ifelse(a > 1, S[a - 1], 0)) - vaccination_rate[a] * S[a]
      dI[a] <- lambda[a] * (S[a] + M[a]) - (gamma[a] + mu[a]) * I[a] + aging_rate[a] * (ifelse(a > 1, I[a - 1], 0))
      dR[a] <- gamma[a] * I[a] - mu[a] * R[a] + aging_rate[a] * (ifelse(a > 1, R[a - 1], 0))
      dV[a] <- -mu[a] * V[a] + aging_rate[a] * (ifelse(a > 1, V[a - 1], 0))
    }
    
    # Return list of all derivatives
    list(c(dM, dS, dI, dR, dV))
  })
}

# Parameters 
parameters <- list(
  beta = matrix(runif(9, min = 0.01, max = 0.05), nrow = 3, ncol = 3),  # Example contact matrix
  gamma = rep(0.1, 3),  # Recovery rates
  mu = rep(0.01, 3),    # Natural mortality rates
  aging_rate = rep(1/12, 3),  # Aging rates per month (example values)
  vaccination_rate = rep(0.001, 3)  # Vaccination rates per age group (example)
)

# Initial state conditions per age group
initial_state <- c(rep(50, 3), rep(1000, 3), rep(5, 3), rep(100, 3), rep(200, 3))

# Time points to solve over
times <- seq(0, 365, by = 1)

# Solve the model
results <- ode(y = initial_state, times = times, func = rubella_model, parms = parameters)

# Convert results to a data frame for visualization
results_df <- as.data.frame(results)


# Plot all compartments 

library(ggplot2)
results_long <- reshape2::melt(results_df, id.vars = "time")

ggplot(data = results_long, aes(x = time, y = value, color = variable)) +
  geom_line() +
  labs(title = "Age-Stratified Rubella Model", x = "Time (days)", y = "Population", color = "Compartment") +
  theme_minimal()

