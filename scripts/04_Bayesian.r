#### Biodiversity modeling course 2021

## Day 7 - Bayesian approaches

## Grid search

# Generate data
set.seed(1859)
x_obs <- rbinom(10, 6, 0.6)
x_obs

# Make a vector of probabilities for p
p <- seq(from = 0.1, to = 0.99, length.out = 99)

# Think about the priors for alpha and phi
mu <- 0.5 # average probability
phi <- 4 # concentration around this value
curve(dbeta(x, mu * phi, (1-mu)*phi))

# Find the probability of each of these values in p (prior)
prior_dens <- dbeta(p, 2, 2)
plot(p, prior_dens)

# Find the likelihood for each of these values of p (likelihood)
likelihood_fct <- function(p) prod(dbinom(x = x_obs, size = 6, prob = p, log=FALSE))
likelihood_dens <- sapply(p, likelihood_fct)

# Multiply these two columns (prior and likelihood) together and normalize
norm <- sum(prior_dens * likelihood_dens)
post_p <- prior_dens * likelihood_dens / norm

plot(p, post_p)
points(p, prior_dens / 50, col = "red")


## Conjugate priors

post_conj <- function(x) dbeta(x, 2 + sum(x_obs), 2 + sum(6 - x_obs))
curve(post_conj(x))

post_conj_p <- post_conj(p)
factor <- sum(post_conj_p) 

curve(post_conj(x))
points(p, post_p*factor)




### Coates & Burton (1999) problem
### DOI:10.1139/cjfr-29-9-1374

library(purrr)
library(truncnorm)

# Number of observations for each species
N <- c(93, 77, 72, 91, 80)

# Parameter estimates for each species (a_mean, a_sd, s_mean, s_sd)
a_estimates <- c(198.5, 260.6, 125.7, 279.0, 506.8)
a_SE <- c(12.1, 23.3, 16.0, 28.6, 71.1)

s_estimates <- c(11.92, 7.56, 2.77, 4.39, 4.55)
s_SE <- c(1.62, 1.04, 0.44, 0.42, 0.37)

# Compute the standard deviations from the standard error for each species
a_sd <- a_SE * sqrt(N)
s_sd <- s_SE * sqrt(N)

# Range of values of L
L <- seq(from = 0, to = 100, length.out = 101)

# Generate fake data
generate_data <- function(L, a_mean, a_sd, s_mean, s_sd) {
    
    # Sample parameter values
    a <- rtruncnorm(n = 1, a = 0, b = Inf, mean = a_mean, sd = a_sd)
    s <- rtruncnorm(n = 1, a = 0, b = Inf, mean = s_mean, sd = s_sd)
    sigma <- rexp(n = 1, rate = 2)

    # Compute the average of the distribution (scientific model)
    mu <- (a*L)/((a/s) + L)
    
    # Generate a value of y
    y <- rtruncnorm(n = 1, a = 0, b = Inf, mean = mu, sd = sigma)
    return (y)
}

# Generate fake data (without sampling parameters)
generate_data_clean <- function(L, a, s, sigma) {
 
    # Compute the average of the distribution (scientific model)
    mu <- (a*L)/((a/s) + L)
    
    # Generate a value of y
    y <- rtruncnorm(n = 1, a = 0, b = Inf, mean = mu, sd = sigma)
    return (y)
}

# Simulate data for the first species
y_sim <- sapply(L, function(x) generate_data(L = x, a_mean = a_estimates[1], a_sd = a_sd[1], s_mean = s_estimates[1], s_sd = s_sd[1]))
plot(L, y_sim)

# Simulate data for the first species
y_sim_clean <- sapply(L, function(x) generate_data_clean(L = x, a = a_estimates[1], s = s_estimates[1], sigma = 10))
plot(L, y_sim_clean)

# Likelihood function
log_likelihood <- function(a, s, sigma) {
    mu <- (a*L)/((a/s) + L)
    sum(log(dtruncnorm(x = y_sim_clean, a = 0, b = Inf, mean = mu, sd = sigma)))
}

# Define a set of parameter values
a <- seq(from = 185, to = 215, by = 1)
s <- seq(from = 5, to = 15, by = 1)
sigma <- seq(from = 8, to = 15, by = 1)

params = expand.grid(a = a, s = s, sigma = sigma)

# Compute log-likelihood for all combinations of parameter values
params$ll <- 0

for (i in 1:nrow(params)) {
    params[i, "ll"] <- log_likelihood(params[i, "a"], params[i, "s"], params[i, "sigma"])
}

params[which.max(params$ll),]

# Effects of priors on the maximum likelihood estimates
params$prior_a <- dtruncnorm(params$a, a = 0, b = Inf, mean = a_estimates[1], sd = a_sd[1])
params$prior_s <- dtruncnorm(params$s, a = 0, b = Inf, mean = s_estimates[1], sd = s_sd[1])
params$prior_sigma <- rexp(params$sigma, rate = 2)
params$post <- with(params, exp(ll) * prior_a * prior_s * prior_sigma)

params[which.max(params$post),]
