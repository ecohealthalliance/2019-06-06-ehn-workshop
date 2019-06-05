# Working script for this workshop
# Find it at https://github.com/ecohealthalliance/2019-06-06-ehn-workshop
#
# Launch the work environment at http://mybinder.org/v2/gh/ecohealthalliance/2019-06-06-ehn-workshop/master?urlpath=rstudio



state_0 = c(S=2000, I = 1)
state_0
params = c(beta = 0.001)

binomial_chain_step = function(state, parameters) {
  S = state["S"]
  I = state["I"]
  beta = parameters["beta"]
  I_next = rbinom(n = 1, size = S, prob = 1 - exp(-beta*I))
  S_next = S - I_next
  state_next = c(S = S_next, I = I_next)
  return(state_next)
}

binomial_chain_step(state_0, params)

binomial_chain_simulate = function(state, parameters, n_steps) {
  output = matrix(nrow = n_steps + 1, ncol = 3)
  output[1,1] = 0
  output[1,2] = state["S"]
  output[1,3] = state["I"]
  colnames(output) <- c("time", "S", "I")
  for (step in 1:n_steps) {
    output[step + 1, 1] = step
    output[step + 1, 2:3] = binomial_chain_step(output[step, 2:3], parameters)
  }
  return(output)
}
params = c(beta=0.001)
binomial_chain_simulate(state_0, params, 10)
results <- binomial_chain_simulate(state_0, params, 100)
plot(I~time, data=results, type="l")
params
state_0
params = c(beta=0.01)
results2 <- binomial_chain_simulate(state_0, params, 100)
plot(I~time, data=results2, type="l")

n_sims = 100
n_steps = 30

state_0 = c(S=2000, I=200)
params = c(beta=0.05)
results_all = list()
for (k in 1:n_sims) {
  results_all[[k]] = binomial_chain_simulate(state_0, params, n_steps)
}

plot(c(0,30), c(0,2000), type="n", xlab="time", ylab="I")
for (k in 1:n_sims) {
  lines(I~time, data=results_all[[k]], type="l", col="blue")
}

## Two-species model

state2_0 = c(S1 = 2000, I1 = 1, S2 = 2000, I2 = 0)
parameters2 = c(beta1 = 0.001, beta2 = 0.001, beta12 = 0.0001)

binomial_chain_step2 = function(state, parameters) {
  S1 = state["S1"]
  S2 = state["S2"]
  I1 = state["I1"]
  I2 = state["I2"]
  beta1 = parameters["beta1"]
  beta2 = parameters["beta2"]
  beta12 = parameters["beta12"]
  I1_next = rbinom(n = 1, size = S1, prob= 1 -exp(-beta1 * I1))
  S1_next = S1 - I1_next
  I2_next = rbinom(n = 1, size = S2, prob = 1 - exp(-(beta2 * I2 + beta12 * I1)))
  S2_next = S2 - I2_next
  return(c(S1_next, I1_next, S2_next, I2_next))
}
binomial_chain_step2(state2_0, parameters2)

binomial_chain_simulate2 = function(state, parameters, n_steps) {
  output = matrix(nrow = n_steps + 1, ncol = 5)
  colnames(output) <- c("time", "S1", "I1", "S2", "I2")
  output[1,1] <- 0
  output[1, 2:5] <- state
  for (step in 1:n_steps) {
    output[step + 1, 1] <- step
    output[step +1, 2:5] <- binomial_chain_step2(output[step, 2:5], parameters)
  }
  return(output)
}

binomial_chain_simulate2(state2_0, parameters2, 100)

n_sims = 1000
n_steps = 30

state2_0 = c(S1 = 2000, I1 = 1, S2 = 2000, I2 = 0)
parameters2 = c(beta1 = 0.0005, beta2 = 0.005, beta12 = 0.00001)
results_all2 = list()
for (k in 1:n_sims) {
  results_all2[[k]] = binomial_chain_simulate2(state2_0, parameters2, n_steps)
}

plot(c(0,30), c(0,1500), type="n", xlab="time", ylab="I")
for (k in 1:n_sims) {
  lines(I1~time, data=results_all2[[k]], type="l", col="blue")
  lines(I2~time, data=results_all2[[k]], type="l", col="red")
}


