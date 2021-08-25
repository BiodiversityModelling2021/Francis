#### Biodiversity modeling course 2021

## Day 5 - Optimization

## Exercise 3.1

library(hrbrthemes)

# Function to optimize
h <- function(x,y) {
    (x*sin(20*y)+y*sin(20*x))^2*cosh(sin(10*x)*x)+(x*cos(10*y)-y*sin(10*x))^2*cosh(cos(20*y)*y)
    }

# Propose a set of parameter values
X <- seq(from = -10, to = 10, by = 1)
Y <- seq(from = -10, to = 10, by = 1)

grid <- expand.grid(x = X, y = Y)
grid$h <- 0

# Find the value of h for all combinations of (x,y)
for (i in 1:nrow(grid)) {
    grid[i, "h"] <- h(x = grid[i, "x"], y = grid[i, "y"])
}
h(100, -100)
# Find local optimum 
optim <- grid[which.max(grid$h),]

# Plot h and parameters (x, y)
grid %>% ggplot(aes(x, y, fill= h)) + 
    geom_tile() +
    theme_ipsum() +
    theme(legend.position="none")


## Exercise 3.2 (simulated annealing)

# Simulate data
mu <- 5
sd <- 1
y <- rnorm(n = 1000, mean = mu, sd = sd)

# Log-likelihood function (function to optimize)
ll <- function(y, mu, sd) {
    sum(dnorm(x = y, mean = mu, sd = sd, log=TRUE))
}

# Initial values
mu_i <- 8
sd_i <- 3

# Temperature function and number of iterations
nsim = 10000
T <- exp(-100/nsim)

# Simulating algorithm 
for (i in 1:nsim) {

    # Draw a value for the average and compute the difference in log-likelihood from previous trial
    mu_test <- mu_i + runif(n = 1, min = -0.1, max = 0.1)
    diff_mu <- ll(y, mu_test, sd_i) - ll(y, mu_i, sd_i)
    
    # Accept if the difference is positive or with a probability p if it's negative
    if(diff_mu > 0) {
        mu_i <- mu_test
    }
    else {
        p <- exp(diff_mu/T)
        P <- runif(n = 1, min = 0, max = 1)
        if(P < p) {
            mu_i <- mu_test
        }
    }   

    # Draw a value for the standard deviation and compute the difference in log-likelihood from previous trial
    sd_test <- runif(n = 1, min = 0.01, max = 2)
    diff_sd <- ll(y, mu_i, sd_test) - ll(y, mu_i, sd_i)

    # Accept if the difference is positive or with a probability p if it's negative
    if(diff_sd > 0) {
        sd_i <- sd_test
    }
    else {
        p <- exp(diff_sd/T)
        P <- runif(n = 1, min = 0, max = 1)
        if(P < p) {
            sd_i <- sd_test
        }
    }

    # Update the temperature
    T <- exp(-100*(i+1)/nsim)  
}

mu_i
sd_i


## Use the previous algorithm to find optimum of function h

# Initial values
x_i <- 4
y_i <- 2

# Temperature function and number of iterations
nsim = 10000
T <- exp(-100/nsim)

# Simulating algorithm 
for (i in 1:nsim) {

    # Draw a value for x and compute the difference in log-likelihood from previous trial
    x_test <- runif(n = 1, min = -100, max = 100)
    diff_x <- h(x_test, y_i) - h(x_i, y_i)
    
    # Accept if the difference is positive or with a probability p if it's negative
    if(diff_x > 0) {
        x_i <- x_test
    }
    else {
        p <- exp(diff_x/T)
        P <- runif(n = 1, min = 0, max = 1)
        if(P < p) {
            x_i <- x_test
        }
    }   

    # Draw a value for y and compute the difference in log-likelihood from previous trial
    y_test <- runif(n = 1, min = -100, max = 0)
    diff_y <- h(x_i, y_test) - h(x_i, y_i)

    # Accept if the difference is positive or with a probability p if it's negative
    if(diff_y > 0) {
        y_i <- y_test
    }
    else {
        p <- exp(diff_y/T)
        P <- runif(n = 1, min = 0, max = 1)
        if(P < p) {
            y_i <- y_test
        }
    }

    # Update the temperature
    T <- exp(-100*(i+1)/nsim)  
}

x_i
y_i




#### Problem to solve (network stability) ####

library(tictoc)
library(tidyverse)
library(patchwork)

#1.Write a function to generate a random matrix for pre-specified S and C
random_matrix <- function(S, C, sigma) {
    L <- matrix(0, nr = S, nc = S)

    int <- matrix(rnorm(S^2, mean = 0,
                    sd = sigma), nr = S, nc = S)

    rand <- matrix(runif(S^2, 0, 1), nr =S, nc = S)

    L[rand<C] <- int[rand<C]

    return (L)
}

#2.Write a function to compute the real part of the largest eigen value of a matrix L

L_eigen <- function(L){
  real_values <- Re(eigen(L)$value)
  max(real_values[which(real_values != 0)])
}

#3 Write a function to perform simulated annealing with a matrix L as input and returning an optimized matrix L_opt as output

#"Candifate function" or here the funciton that switch aroud the eigen values in our matrix
c_x <- function(L) {
  # Get 2 indexes to permute
  ind_1 <- sample(which(L != 0), 1) # must not be zero
  ind_2 <- sample((1:length(L))[-ind_1], 1)
  # Permute indexes in new matrix
  L_opt <- L
  L_opt[ind_1] <- L[ind_2]
  L_opt[ind_2] <- L[ind_1]
  return(L_opt)
}

# Set conditions for the simulated annealing sequence
T_fn <- function(T0, alpha, step) T0*exp(alpha*step)
step = 1000
diff = 0.01
# exp(-diff/T_fn(T0, alpha, step))


# Initiate the algorithm
T0 <- 10
alpha = -0.001
nsteps <- 12000
# L0 <- random_matrix(S, C, sigma)

# Main loop
simulated_annealing <- function(L0, nsteps, T0, alpha) {
  eigens <- rep(0, nsteps)

  tic("Main loop"); for(step in 1:nsteps) {

    # Switch 2 values in the adjacency matrix
    L1 <- c_x(L0)

    # Evaluate the maximum eigenvalue
    L1_eigen <- L_eigen(L1)
    L0_eigen <- L_eigen(L0)

    # Compute the difference
    diff <- L1_eigen - L0_eigen

    # Accept if improvement (smaller largest eigenvalue)
    if(diff < 0) L0 <- L1

    # Accept wrong candidates
    else {
        p <- exp(-diff/T_fn(T0, alpha, step))
        if(runif(1)<p) L0 <- L1
        }

    # Record values to add later
    eigens[step] <- L0_eigen
  }; toc()

  # Plot the result
  eigens_plot <- tibble(eigens) %>%
    ggplot(aes(x = 1:length(eigens), y = eigens)) +
    geom_line() +
    labs(x = "Steps", y = "Minimum eigenvalue")

  # Collect results
  results <- list(
    eigens = tibble(eigens),
    L_opt = L0,
    plot = eigens_plot
  )

  return(results)
}

#4 Run the optimization procedure for a gradient of species richness and record L and L_opt
optimize_S <- function(S, nsteps, T0, alpha){
  # Define values
  C = 0.1
  sigma = 0.5

  # Generate random matrix
  L0 <- random_matrix(S, C, sigma)

  # Optimize
  results <- simulated_annealing(L0 = L0, nsteps = nsteps, T0 = T0, alpha = alpha)
}
res10 <- optimize_S(S = 10, nsteps = nsteps, T0 = T0, alpha = alpha)
res25 <- optimize_S(S = 25, nsteps = nsteps, T0 = T0, alpha = alpha)
res50 <- optimize_S(S = 50, nsteps = nsteps, T0 = T0, alpha = alpha)

res10$plot
res25$plot
res50$plot

Re(eigen(res25$L_opt)$value)

res <- map(
  c(10, 25, 50),
  function(x) optimize_S(S = x, nsteps = nsteps, T0 = T0, alpha = alpha)
)

res[[1]]$plot/
  res[[2]]$plot/
  res[[3]]$plot

#5 Evaluate if the "potential for stability" ($L-L_{opt}$) relates to complexity
res[[1]]$L_opt
