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

