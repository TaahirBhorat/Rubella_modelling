library(deSolve)
library(ggplot2)
library(reshape2)


rubella_model <- function(time, state, parameters) {
  n_age <- parameters$n_age
  n_compartments <- 5  # M, S, I, R, V
  

  M <- state[1:n_age]
  S <- state[(n_age + 1):(2 * n_age)]
  I <- state[(2 * n_age + 1):(3 * n_age)]
  R <- state[(3 * n_age + 1):(4 * n_age)]
  V <- state[(4 * n_age + 1):(5 * n_age)]
  
  with(as.list(parameters), {
    # Seasonal transmission function
    seasonality_factor <- 1 + seasonal_amp * cos(2 * pi * (time %% 365) / 365)
    
    # Force of infection calculation per age group
    phi <- numeric(n_age)
    for (a in 1:n_age) {
      phi[a] <- 1 - exp(-sum(beta[a, ] * I) * seasonality_factor / sum(M + S + I + R + V))
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
      
      dM[a] <- -da[a] * M[a] - mu[a] * M[a] + aging_in_M - aging_rate[a] * M[a]
      dS[a] <- da[a] * M[a] - phi[a] * S[a] - mu[a] * S[a] + aging_in_S - aging_rate[a] * S[a] - va[a] * S[a]
      dI[a] <- phi[a] * S[a] - (gamma[a] + mu[a]) * I[a] + aging_in_I - aging_rate[a] * I[a]
      dR[a] <- gamma[a] * I[a] - mu[a] * R[a] + aging_in_R - aging_rate[a] * R[a]
      dV[a] <- va[a] * S[a] - mu[a] * V[a] + aging_in_V - aging_rate[a] * V[a]
    }
    
    list(c(dM, dS, dI, dR, dV))
  })
}


initialize_parameters <- function(n_age) {

  beta <- matrix(runif(n_age * n_age, min = 0.05, max = 0.15), nrow = n_age, ncol = n_age)  # Transmission rates
  aging_rate <- rep(1/12, n_age)  # Aging rate (example)
  da <- rep(0.05, n_age)  # Loss of maternal immunity
  va <- rep(0.03, n_age)  # Vaccination rate
  gamma <- rep(0.1, n_age)  # Recovery rate
  mu <- rep(0.01, n_age)  # Natural death rate
  seasonal_amp <- 0.2  # Seasonal amplitude
  
  return(list(
    beta = beta,
    da = da,
    va = va,
    gamma = gamma,
    mu = mu,
    aging_rate = aging_rate,
    seasonal_amp = seasonal_amp,
    n_age = n_age
  ))
}

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

n_age <- 37
parameters <- initialize_parameters(n_age)
initial_state <- initialize_state(n_age)

# Simulate for 30 years with 17-day intervals
times <- seq(0, 30 * 365, by = 17)
results <- ode(y = initial_state, times = times, func = rubella_model, parms = parameters)

results_df <- as.data.frame(results)

results_df$total_infected <- rowSums(results_df[, grepl("I", colnames(results_df))])
ggplot(results_df, aes(x = time, y = total_infected)) +
  geom_line() +
  labs(title = "Total Rubella Cases Over Time", x = "Time (days)", y = "Total Infected") +
  theme_minimal()

average_age <- sapply(1:nrow(results_df), function(row) {
  infected <- results_df[row, grepl("I", names(results_df))]
  total_cases <- sum(infected)
  if (total_cases > 0) {
    sum(infected * 1:n_age) / total_cases
  } else {
    NA
  }
})

average_age_df <- data.frame(time = results_df$time, average_age = average_age)
ggplot(average_age_df, aes(x = time, y = average_age)) +
  geom_line(na.rm = TRUE) +
  labs(title = "Average Age of Rubella Cases Over Time", x = "Time (days)", y = "Average Age") +
  theme_minimal()
