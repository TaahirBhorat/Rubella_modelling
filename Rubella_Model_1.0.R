library(deSolve)

rubella_model <- function(t, state, parameters) {
  # Extract parameters and states
  N_age <- parameters$N_age
  beta <- parameters$beta
  gamma <- parameters$gamma
  mu <- parameters$mu
  delta <- parameters$delta
  aging_rate <- parameters$aging_rate
  
  # Initialize rate of change for each compartment
  dM <- numeric(N_age)
  dS <- numeric(N_age)
  dI <- numeric(N_age)
  dR <- numeric(N_age)
  dV <- numeric(N_age)
  
  # Dynamics for each age class
  for (a in 1:N_age) {
    # Aging out of the current class to the next one
    aging_out <- ifelse(a < N_age, aging_rate[a], 0)
    aging_in <- ifelse(a > 1, aging_rate[a - 1], 0)
    
    # Force of infection
    force_of_infection <- sum(beta[a, ] * state['I'] / sum(state['N']))
    
    # Model equations
    dM[a] <- -mu[a] * state['M'][a] - force_of_infection * state['M'][a] + aging_in * state['M'][a-1] - aging_out * state['M'][a]
    dS[a] <- -mu[a] * state['S'][a] - force_of_infection * state['S'][a] + delta[a] * state['R'][a] + aging_in * state['S'][a-1] - aging_out * state['S'][a]
    dI[a] <- force_of_infection * (state['S'][a] + state['M'][a]) - (gamma[a] + mu[a]) * state['I'][a] + aging_in * state['I'][a-1] - aging_out * state['I'][a]
    dR[a] <- gamma[a] * state['I'][a] - (mu[a] + delta[a]) * state['R'][a] + aging_in * state['R'][a-1] - aging_out * state['R'][a]
    dV[a] <- -mu[a] * state['V'][a] + aging_in * state['V'][a-1] - aging_out * state['V'][a]
  }
  
  # Return the rate of change
  list(c(dM, dS, dI, dR, dV))
}
