#### Biodiversity modeling course 2021

## Day 9-10 - Sampling the posterior 

## Exercise 1 (Rejection sampling)

rej_sample <- function(n, mu, sigma) {
    # Define candidate function 
    c <- function(x) dunif(x, mu - 3 * sigma, mu + 3 * sigma)
    # Define target function
    t <- function(x) dnorm(x, mu, sigma)
    # Define constant M
    M <- 3
    # Instantiate a result vector 
    result <- rep(NA, n)
    # Rejection sampling 
    for (i in 1:n) {
        while(is.na(result[i])) {
        # Draw sample 
        x <- runif(1, mu - 3 * sigma, mu + 3 * sigma)
        # Compute acceptance probability
        p <- t(x) / (M * c(x))
        # Draw random value 
        rand <- runif(1, 0, 1)
        # Accept of reject x
        if (rand < p) result[i] <- x 
        }
    }
    return (result)
}

sample <- rej_sample(1000, 0, 1)
hist(sample)



## Metropolis-Hastings algorithm

metropolis <- function(N = 5000, X0 = 0, A) {
    # N: number of MCMC iterations
    # X0: starting value of the chain
    # A: step parameter
    # Object to store results
    chain <- numeric(N + 1)
    chain[1] <- X0

    # DEFINE candidate distribution
    cand_fn <- function(y, A) runif(1, y-A, y+A)  

    # DEFINE target distribution
    target_fn <- dnorm(x, 0, 1)  
 
    # DEFINE acceptance probability function
    accept_fn <- function(current.value, previous.value) target_fn(current.value)/target_fn(previous.value)

    for(step in 2:(N+1)) {
        # PROPOSE candidate
        current.value <- cand_fn(chain[step-1], A)
        # COMPUTE acceptance probability
        p <- accept_fn(current.value, chain[step-1])
        # DO rejection sampling
        rand <- runif(1)
        if(rand < p) chain[t] <- current.value
        else chain[t] <- chain[t 1]
    }
}



## Gibbs sampling 

mh.sampler <- function(previous, A, target, prior, other.params,
                        data, cand_fn = function(x, A) runif(1, x-A, x+A)) {
  # propose and select candidate values using metropolis–hastings algorith
  # cand_fn: the function used to draw samples
  # previous: the previous state in the chain
  # A: the tuning parameter
  # target : function(x) returning the target density at X
  candidate.value <- cand_fn(previous, A)
  p <- mh.acceptance(candidate.value, previous, target, prior,
                      other.params, data)
  U <- runif(1)
  if (U < p) {
    result <- candidate.value
  } else {
    result <- previous
  }
  return(result)
}

# Conditional mu
mu.conditional <- function(mu, mu.prior, sigma, data) {
  # log-likelihood
  lik <- exp(sum(dnorm(data, mu, sigma, log=TRUE)))
  # log prior
  prior <- exp(-(mu-mu.prior[1])^2 / (2*mu.prior[2]^2))
  # Return the conditional mu
  return(lik*prior)
}

# Conditional sigma
sigma.conditional <- function(sigma, tau.prior, mu, data) {
  # Compute the posterior
  # Make sure sigma is positive
  if (sigma <= 0) {
    return(-Inf)
  } else {
    # Likelihood
    lik <- exp(sum(dnorm(data, mu, sigma, log = TRUE)))
    # Prior
    prior <- exp((tau.prior[1]^2-1)*log(tau(sigma))-
                 tau(sigma)*tau.prior[2])
  }
  # Exponentiate back
  return(lik*prior)
}

# Acceptance
mh.acceptance <- function(current, previous, target, prior, other.params, data) {
  # target : function returning target density at parameter x
  num <- target(current, prior, other.params, data)
  denom <- target(previous, prior, other.params, data)
  if (denom == 0) return(-1)
  else return(num/denom)
}

# Run the whole thing
N <- 5000 # number of iterations of the mcmc
starting <- c(10,4) # starting values for mu and sigma
A <- c(1.4, 0.8) # tuning parameters for the candidate function
X <- c(21.4, 17.64, 18.31, 15.12, 14.40, 15, 19.59, 15.06, 15.71, 14.65) # data set
mu.prior <- c(14.5, 0.6) # prior parameters for the mean (normal distribution)
tau.prior <- c(1.28367, 3.068) # prior parameters for sigma (gamma distribution)
## Small function to get a parameter of the gamma
tau <- function(sigma) 1/(sigma^2)
# Object to store the results
chain <- matrix(nrow=N+1, ncol=2)
chain [1,] <- starting
# Run the chain
for(t in 2:(N+1)) {
    # Sample for mu, fixing sigma
    chain[t ,1] <- mh.sampler(previous = chain[t-1,1], A = A[1],
                               target = mu.conditional,
                               prior = mu.prior,
                               other.params = chain[t-1,2], data = X)
    # Sample for sigma, fixing mu
    chain[t ,2] <- mh.sampler(previous = chain[t-1,2], A = A[2],
                               target = sigma.conditional,
                               prior = tau.prior,
                               other.params = chain[t-1,1], data = X)
}





#### Posterior distribution of the hemlock example (species 1) ####
#### WORK IN PROGRESS

library(truncnorm)

## Simulate data 

# Generate fake data (without sampling parameters)
generate_data <- function(L, a, s, sigma) {
 
    # Compute the average of the distribution (scientific model)
    mu <- (a*L)/((a/s) + L)
    
    # Generate a value of y
    y <- rtruncnorm(n = 1, a = 0, b = Inf, mean = mu, sd = sigma)
    return (y)
}

# Values of L 
L <- runif(n = 1000, min = 0, max = 100)

# Simulate data 
y <- sapply(L, function(x) generate_data(L = x, a = 198.5, s = 11.92, sigma = 10))
plot(L, y)


## Gibbs sampling 

# Candidate function 
cand_fn = function(x, A) {
    # A: A tuning parameter 
    runif(1, x-A, x+A)
}

mh.sampler <- function(previous, A, target, prior, param1, param2,
                        L, data, cand_fn = cand_fn) {
  # propose and select candidate values using metropolis–hastings algorith
  # cand_fn: the function used to draw samples
  # previous: the previous state in the chain
  # A: the tuning parameter
  # target : function(x) returning the target density at X
  candidate.value <- cand_fn(previous, A)
  p <- mh.acceptance(candidate.value, previous, target, prior,
                      other.params, data)
  U <- runif(1)
  if (U < p) {
    result <- candidate.value
  } else {
    result <- previous
  }
  return(result)
}

# Conditional a
a.conditional <- function(a, a.prior, s, sigma, L, data) {
  # Compute mu
  mu <- (a*L)/((a/s) + L)
  # log-likelihood
  lik <- exp(sum(dnorm(data, mu, sigma, log=TRUE)))
  # log prior
  prior <- dnorm(a, a.prior, sigma, log=FALSE)
  # Return the conditional a
  return(lik*prior)
}

# Conditional s
s.conditional <- function(s, s.prior, a, sigma, L, data) {
  # Compute mu
  mu <- (a*L)/((a/s) + L)
  # log-likelihood
  lik <- exp(sum(dnorm(data, mu, sigma, log=TRUE)))
  # log prior
  prior <- dnorm(s, s.prior, sigma, log=FALSE)
  # Return the conditional s
  return(lik*prior)
}

# Conditional sigma
sigma.conditional <- function(sigma, sigma.prior, a, s, L, data) { 
    # Compute mu
    mu <- (a*L)/((a/s) + L)
    # Likelihood
    lik <- exp(sum(dnorm(data, mu, sigma, log=TRUE)))
    # Prior
    prior <- dexp(sigma, sigma.prior, log=FALSE)
    # Return the conditional sigma 
    return(lik*prior)
}

# Acceptance
mh.acceptance <- function(current, previous, target, prior, param1, param2, L, data) {
  # target : function returning target density at parameter x
  num <- target(current, prior, param1, param2, L, data)
  denom <- target(previous, prior, param1, param2, L, data)
  if (denom == 0) return(-1)
  else return(num/denom)
}

# Run the whole thing
N <- 5000 # number of iterations of the mcmc
starting <- c(150, 15, 12) # starting values for a, s and sigma
A <- c(1.5, 0.5, 0.5) # tuning parameters for the candidate function
X <- y # data set
a.prior <- c(200, 15) # prior parameters for a (normal distribution)
s.prior <- c(15, 2) # prior parameters for s (normal distribution)
sigma.prior <- 0.5 # prior parameters for sigma (exponential distribution)
# Object to store the results
chain <- matrix(nrow=N+1, ncol=3)
chain[1,] <- starting
# Run the chain
for(t in 2:(N+1)) {
    # Sample for a, fixing s and sigma
    chain[t,1] <- mh.sampler(previous = chain[t-1,1], 
                               A = A[1],
                               target = a.conditional,
                               prior = a.prior,
                               param1 = chain[t-1,2],
                               param2 = chain[t-1,3],
                               L = L,
                               data = X)
    # Sample for s, fixing a and sigma
    chain[t,2] <- mh.sampler(previous = chain[t-1,2], 
                               A = A[2],
                               target = s.conditional,
                               prior = s.prior,
                               param1 = chain[t-1,1],
                               param2 = chain[t-1,3],
                               L = L, 
                               data = X)
    # Sample for sigma, fixing a and s
    chain[t,3] <- mh.sampler(previous = chain[t-1,3], 
                               A = A[3],
                               target = sigma.conditional,
                               prior = sigma.prior,
                               param1 = chain[t-1,1],
                               param2 = chain[t-1,2],
                               L = L, 
                               data = X)
}