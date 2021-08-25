#### Biodiversity modeling course 2021

## Day 3 - Distributions 

curve(dnorm(x, mean = 0, sd = 0.05), 
      xlim = c(-3, 3))

dnorm(0.01, mean = 0, sd = 0.05)

curve(dbeta(x, 3, 7))

curve(dbinom(x, size = 20, p = 0.2),
      n = 21, 
      type = "p",
      xlim = c(0, 20))

curve(dbinom(x, size = 20, p = 0.8),
      n = 21, 
      type = "p",
      xlim = c(0, 20),
      col = "red", add = TRUE)

curve(dbinom(x, size = 20, p = 0.1),
      n = 21, 
      type = "p",
      xlim = c(0, 20),
      col = "blue", add =TRUE)

library(extraDistr)

curve(dbbinom(x, size = 1000, alpha = 5, beta = 13),
      n = 1001, 
      type = "p",
      xlim = c(0, 1000))

curve(dbbinom(x, size = 1000, alpha = 2, beta = 5),
      n = 1001, 
      type = "p",
      xlim = c(0, 1000),
      col = "red", add = TRUE)

curve(dbbinom(x, size = 1000, alpha = 14, beta = 10),
      n = 1001, 
      type = "p",
      xlim = c(0, 1000),
      col = "blue", add = TRUE)



## Maximum likelihood

X <- c(4, 2, 10, 5, 8, 4)

# Maximize log-likelihood
ll <- function(X, u, sd) {
  sum(dnorm(X, u, sd, log = TRUE))
}

ll(X, 5, 1)
ll(X, 5.5, 1)
ll(X, 5.5, 2.95)

# Exercise 2.1
X <- c(0, 1, 1, 0, 1, 1, 0, 0, 0, 0)

ll <- function(X, p) {
  sum(dbern(X, p, log = TRUE))
}

library(purrr)
P <- seq(from = 0.01, to = 0.99, by = 0.01)
ll_P <- P %>% map_dbl(ll, X = X) 

plot(P, ll_P)

P_maxll <- P[ll_P == max(ll_P)]


### Exercise 2.2

## Pseudo-code of a likelihood function 

## Define function f(obs, covariates, parameters) 
## Load obs, covariates, parameters
## Compute mu = f(covariates, parameters)
## Compute p_obs = pdf(obs, mu, parameters)
## Compute ll = sum(log(p_obs))

library(dplyr)
library(ggplot2)
library(extraDistr)

# Load data
sutton <- readr::read_csv2("data/sutton.csv")

acsa <- sutton$acsa
y <- sutton$y

# Plot data
sutton %>%
  ggplot(aes(x = y, y = acsa)) +
  geom_point(col = "red")

## Linear regression 

# Define a set of parameter values
a <- seq(from = 14, to = 16, by = 0.1)
b <- seq(from = -0.1, to = 0, by = 0.001)
s <- seq(from = 6, to = 7, by = 0.1)

params = expand.grid(a = a, b = b, s = s)

# Each combination of parameter values has its own log-likelihood
params$ll = 0

# Function to compute the log-likelihood of given parameter values
ll <- function(x, mu, s) {
  sum(dnorm(x, mu, s, log = TRUE))
}

# Compute log-likelihood for all combinations of parameter values
for (i in 1:nrow(params)) {
  # Compute the expected value for each elevation value given the parameters
  mu_i <- params[i, "a"] + params[i, "b"] * y
  # Compute log-likelihood of the parameters
  params[i, "ll"] <- ll(acsa, mu_i, params[i, "s"])
}

# Set of parameter values that maximize the log-likelihood
params_max <- params[which.max(params$ll),]

# Plot data with fitted parameter values (max log-likelihood)
sutton %>%
  ggplot(aes(x = y, y = acsa)) +
  geom_point(col = "red") + 
  geom_abline(intercept = params_max$a, slope = params_max$b)


## Logistic regression 

# Convert count data to binary data
acsa_bin <- rep(0, length(acsa))
for (i in 1:length(acsa)) {
  if (acsa[i] > 0) {
    acsa_bin[i] <- 1
  }
}

sutton$acsa_bin <- acsa_bin

# Plot binary data
sutton %>%
  ggplot(aes(x = y, y = acsa_bin)) +
  geom_point(col = "red")

sutton %>%
  ggplot(aes(x = y)) +
  geom_histogram() +
  facet_wrap(~acsa_bin)

# Define a set of parameter values
a <- seq(from = 0, to = 10, by = 0.1)
b <- seq(from = -1, to = 0, by = 0.001)

params = expand.grid(a = a, b = b)

# Each combination of parameter values has its own log-likelihood
params$ll = 0

# Function to compute the log-likelihood of given parameter values
ll <- function(x, p) {
  sum(dbern(x, p, log = TRUE))
}

# Compute log-likelihood for all combinations of parameter values
for (i in 1:nrow(params)) {
  # Compute the expected value for each elevation value given the parameters
  z_i <- params[i, "a"] + params[i, "b"] * y
  p_i <- 1 / (1 + exp(-z_i))
  # Compute log-likelihood of the parameters
  params[i, "ll"] <- ll(acsa_bin, p_i)
}

# Set of parameter values that maximize the log-likelihood
params_max <- params[which.max(params$ll),]

# Plot with fitted parameter values (maximum likelihood)
sutton$p_logist <- 1 / (1 + exp(-(params_max$a + params_max$b * y)))
  
sutton %>%
  ggplot(aes(x = y, y = acsa_bin)) +
  geom_point(col = "red") +
  geom_line(aes(x = y, y = p_logist))


## Poisson 

# Define a set of parameter values
a <- seq(from = 2, to = 4, by = 0.1)
b <- seq(from = -1, to = 0, by = 0.001)

params = expand.grid(a = a, b = b)

# Each combination of parameter values has its own log-likelihood
params$ll = 0

# Function to compute the log-likelihood of given parameter values
ll <- function(x, l) {
  sum(dpois(x, l, log = TRUE))
}

# Compute log-likelihood for all combinations of parameter values
for (i in 1:nrow(params)) {
  # Compute the expected value for each elevation value given the parameters
  z_i <- params[i, "a"] + params[i, "b"] * y
  l_i <- exp(z_i)
  # Compute log-likelihood of the parameters
  params[i, "ll"] <- ll(acsa, l_i)
}

# Set of parameter values that maximize the log-likelihood
params_max <- params[which.max(params$ll),]




#### Problem to solve ####

library(dplyr)

# Read and explore data
transit <- readr::read_delim("Data/transitions.txt", delim = " ", col_names = c("x","ID", "temp", "state1", "state2", "interval"), skip = 1)

transit_table <- table(transit$state1, transit$state2)

# Subset data to get transitions from state B only
transit_B <- transit[transit$state1 == "B",]

# Compute the frequency of each transition 
transit_BB <- sum(transit_B$state2 == "B")
transit_BM <- sum(transit_B$state2 == "M")
transit_BR <- sum(transit_B$state2 == "R")
transit_BT <- sum(transit_B$state2 == "T")

# Write likelihood function with 4 parameters (transition probabilities)
ll <- function(xB, xM, xR, xT, pB, pM, pR, pT) {
  dmultinom(x = c(xB, xM, xR, xT), 
            prob = c(pB, pM, pR, pT), 
            log = TRUE)
}

# Test function with candidate parameter values
ll(xB = transit_BB, xM = transit_BM, xR = transit_BM, xT =  transit_BT,
   pB = 0.25, pM = 0.25, pR = 0.25, pT = 0.25)

# Propose alternative candidate values and compare
ll(xB = transit_BB, xM = transit_BM, xR = transit_BM, xT =  transit_BT,
   pB = 0.5, pM = 0.3, pR = 0.2, pT = 0)

# Find maximum likelihood estimates using a grid search 
# Candidate parameters (must sum to 1 and all individual probabilities between 0 and 1)
pB <- seq(from = 0, to = 1, by = 0.01)
pM <- seq(from = 0, to = 1, by = 0.01)
pR <- seq(from = 0, to = 1, by = 0.01)

params = expand.grid(pB = pB, pM = pM, pR = pR)
params$pT = round(1 - params$pB - params$pM - params$pR, 2)
params <- params[params$pT >= 0,]

# Find maximum likelihood 
params$ll <- 0

for (i in 1:nrow(params)) {
  params[i, "ll"] <- ll(xB = transit_BB, 
                    xM = transit_BM, 
                    xR = transit_BM, 
                    xT =  transit_BT,
                    pB = params[i, "pB"], 
                    pM = params[i, "pM"], 
                    pR = params[i, "pR"], 
                    pT = params[i, "pT"])
}

params_max <- params[which.max(params$ll),]



## Redo the whole thing but with the entire dataset

# Observed transitions (col = year t, row = year t + 1)
transit_table <- table(transit$state2, transit$state1)

# Write likelihood function with 4 parameters (transition probabilities)
ll <- function(xB, xM, xR, xT, pB, pM, pR, pT) {
  dmultinom(x = c(xB, xM, xR, xT), 
            prob = c(pB, pM, pR, pT), 
            log = TRUE)
}

# Find maximum likelihood estimates using a grid search 
# Candidate parameters (must sum to 1 and all individual probabilities between 0 and 1)
pB <- seq(from = 0, to = 1, by = 0.01)
pM <- seq(from = 0, to = 1, by = 0.01)
pR <- seq(from = 0, to = 1, by = 0.01)

params = expand.grid(pB = pB, pM = pM, pR = pR)
params$pT = round(1 - params$pB - params$pM - params$pR, 2)
params <- params[params$pT >= 0,]

# Find maximum likelihood for each column independently 
params$ll <- 0

params_B <- params
params_M <- params
params_R <- params
params_T <- params

find_maxll <- function(params, X) {
    for (i in 1:nrow(params)) {
          params[i, "ll"] <- ll(xB = transit_table["B", X], 
                                xM = transit_table["M", X], 
                                xR = transit_table["R", X],
                                xT = transit_table["T", X],
                                pB = params[i, "pB"], 
                                pM = params[i, "pM"], 
                                pR = params[i, "pR"], 
                                pT = params[i, "pT"])
          }
    params_max <- params[which.max(params$ll),]
    return (params_max)
}

params_B_max <- find_maxll(params_B, "B")
params_M_max <- find_maxll(params_M, "M")
params_R_max <- find_maxll(params_R, "R")
params_T_max <- find_maxll(params_T, "T")

# Get transition matrix and compare with observed proportions
T <- matrix(c(params_B_max$pB, params_M_max$pB, params_R_max$pB, params_T_max$pB,
            params_B_max$pM, params_M_max$pM, params_R_max$pM, params_T_max$pM,
            params_B_max$pR, params_M_max$pR, params_R_max$pR, params_T_max$pR,
            params_B_max$pT, params_M_max$pT, params_R_max$pT, params_T_max$pT), 
            nrow = 4, 
            ncol = 4,
            byrow = TRUE)

transit_table_prop <- round(apply(transit_table, 2, function(x) {x/sum(x)}), 2)
