# Setting the model parameters
## ------------------------------------------------------------------------

# Do three things here: Create an object (state_0),
# Make an object with two variables in it, and assign that objecct to state_0
state_0 <- c(S=2000, I=1)
state_0

my_parms <- c(beta=0.001)
my_parms

# Defining the model function
## ------------------------------------------------------------------------
binomial_chain_step = function(state, parameters) {
  S = state["S"]
  I = state["I"]
  beta = parameters["beta"]
  I_next = rbinom(n = 1, size = S, prob = 1 - exp(-beta*I))
  S_next = S - I_next
  state_next = c(S = S_next, I = I_next)
  return(state_next)
}

# Now let's run this function, using our state at parameters as inputs.
## ------------------------------------------------------------------------
binomial_chain_step(state_0, my_parms)

# Now we write another function that runs this many times a row in a _loop_,
# and store the results each time. It needs another input - how many steps to run?
## ------------------------------------------------------------------------
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

# Rub a simulation! Do it again, but store the results
## ------------------------------------------------------------------------
binomial_chain_simulate(state_0, my_parms, 10)
results <- binomial_chain_simulate(state_0, my_parms, 100)

#' Plot these resultss
## ------------------------------------------------------------------------
plot(I~time, data=results, type="l")


# Now do the simulation many times
## ------------------------------------------------------------------------
n_sims <- 100
n_step <- 30
results_all <- list()
for (k in 1:n_sims) {
  results_all[[k]] <- binomial_chain_simulate(state_0,my_parms, n_step)
}


# Plot the many simulations
## ------------------------------------------------------------------------
plot(c(0,20),c(0,400),type="n",xlab="time",ylab="I")
for (k in 1:n_sims) {
  lines(I~time, data = results_all[[k]], type="l")
}

##### Exercise! #####
#
# Explore the dynamics of the system for different values of
# $\beta$, as well as different initial values of $S$ and $I$.
#'


# Now extend this model to a *two species* version to look at spillover scenarios.
# Now we have to classes each of susceptible and infected individuals: $S1, I1, S2, I2$
## ------------------------------------------------------------------------
state2_0 = c(S1 = 2000, I1 = 1, S2 = 2000, I2 = 0)
parameters2 = c(beta1 = 0.001, beta2 = 0.001, beta12 = 0.0001)

# Define the two species model step
## ------------------------------------------------------------------------
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

# Run this once to test
## ------------------------------------------------------------------------
binomial_chain_step2(state2_0, parameters2)

# Now make this a whole simulation sequence for two-species
## ------------------------------------------------------------------------
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

# Run a two-species simulation
## ------------------------------------------------------------------------
binomial_chain_simulate2(state2_0, parameters2, 100)
results2 <- binomial_chain_simulate2(state2_0, parameters2, 100)

# Plot a two-species simulation
## ------------------------------------------------------------------------
plot(I1~time, data = results2, type="l", col="blue")
lines(I2~time, data= results2, type="l", col="red")

# Now make a loop for many two-species simulations
## ------------------------------------------------------------------------
n_sims = 1000
n_steps = 30

state2_0 = c(S1 = 2000, I1 = 1, S2 = 2000, I2 = 0)
parameters2 = c(beta1 = 0.0005, beta2 = 0.005, beta12 = 0.00001)
results_all2 = list()
for (k in 1:n_sims) {
  results_all2[[k]] = binomial_chain_simulate2(state2_0, parameters2, n_steps)
}

# Plot the many simulations
## ------------------------------------------------------------------------
plot(c(0,30), c(0,1500), type="n", xlab="time", ylab="I")
for (k in 1:n_sims) {
  lines(I1~time, data=results_all2[[k]], type="l", col="blue")
  lines(I2~time, data=results_all2[[k]], type="l", col="red")
}


##### Exercise! #####
#
# 1. Simulate a disease, such as rabies, where spillover is common from the reservoir host to
#    the spillover host, but there is little spread in the spillover hosts.
# 2.  Simulate a disease, such as influenza, where spillover from wild hosts
#     is extremely rare, but the disease is highly contagious in spillover
#     hosts
# 3. Make the disease in (2) less contagious in the wild host.  What type of
#    patterns do we see?
#

##### Susceptible-Infected-Recovered (SIR) Models
#
# SIR Models
#

## ------------------------------------------------------------------------
# (When running thhis part script on your  owncomputer, you may need to to run
# install.packages("deSolve")
# to install this package before proceeding
library(deSolve)


# Write a function to define the differential equation
## ------------------------------------------------------------------------
sir_diff_eqs <- function (time, state, parameters) {
  ## first extract the state variables
  S <- state["S"]
  I <- state["I"]
  R <- state["R"]
  N <- S+I+R
  ## now extract the parameters
  beta <-  parameters["beta"]
  gamma <- parameters["gamma"]
  mu <-    parameters["mu"]
  B <-     parameters["B"]
  ## now code the model equations
  dSdt <- B-beta*S*I/N-mu*S
  dIdt <- beta*S*I/N-(mu+gamma)*I
  dRdt <- gamma*I-mu*R
  ## combine results into a single vector
  derivs <- list(c(dSdt,dIdt,dRdt))
  ## return result as a list!
  return(derivs)
}

# Write a function to calculate $R_0$.
## ------------------------------------------------------------------------
R0 <- function(parameters) {
  R0 = parameters["beta"] /
         (parameters["mu"]+parameters["gamma"])
}

# Define the times at which we want solutions, assign some
# values to the parameters, and specify the **initial conditions**,
# **i.e.**, the values of the state variables $S$, $I$, and $R$ at the
# beginning of the simulation:
## ------------------------------------------------------------------------
times <- seq(0,30,by=1/120)
parameters  <- c(B=1/70,mu=1/70,N=1,beta=400,gamma=365/14)
state_0 <- c(S=1-0.001-0.9,I=0.001,R=0.9)


# Now we can simulate a model trajectory with the `ode` command:
## ------------------------------------------------------------------------
out <- ode(state_0,times,sir_diff_eqs,parameters)

# Plot the results
## ------------------------------------------------------------------------
plot(R~time,data=out,type='l', ylim=c(0,1), ylab="Population")
lines(I~time,data=out,type='l',col='red'); par(new=TRUE)
lines(S~time,data=out,type='l',  col='blue')

# Plot the variables against each other, instead of time, e.g., in "phase space"
## ------------------------------------------------------------------------
plot(I~S,data=out) #,type='b',log='xy',yaxt='n',xlab='S',cex=0.5)

#### Exercise ####

# Explore the dynamics of the system for different values of the
# `beta` and `B` parameters by simulating and plotting trajectories
# as time series and in phase space (e.g., I vs. S).  How  do the $\beta$, $B$,
# and $R0$ relate to the type of trajectories you get?

#  Big Challenges (Pick one, if there is time, do as group)
#
#  1.  What if you have a disease that has a latent period in the host before
#  it starts infecting other hosts?  Change the SIR model to an SEIR model
#  with four compartments (Susceptible, EXPOSED, Infectious, Recovered) and
#  plot the dynamics of all four.
#  2.  What if there are seasonal dynamics to the disease? Change the model
#  so that the value of either $\B$ cycles up and down, representing seasonality
#  of reproduction, or $\beta$ cycles up and down, representing seasonality in
#  contact.  (Hint: there are `sin()` and `cos()` functions, and `time` is
#  a variable in your differential equation function already).  Plot several
#  parameterizations of this model and compare to the version already created
