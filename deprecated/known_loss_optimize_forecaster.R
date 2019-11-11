#### DEFINING LOSS FUNCTION
####
####

# define hinge-like loss function
loss_function <- function(x,y){
  4 * pmax(x-y,0) + 1 * pmax(y-x,0)
}

# (x,y) grid
x_vector <- seq(from = 0, to = 1, length.out = 20)
y_vector <- seq(from = 0, to = 1, length.out = 20)
xy_grid <- expand.grid(x = x_vector, y = y_vector)

# compute loss for each grid point
xy_grid$loss <- NA
for(i in 1:nrow(xy_grid)){
  xy_grid$loss[i] <- loss_function(x = xy_grid$x[i], y = xy_grid$y[i])
}

# satellite view of loss function
library(ggplot2)
ggplot(data = xy_grid) +
  geom_tile(mapping = aes(x = x, y = y, fill = loss))

# mountain view of loss function
library(lattice)
wireframe(loss ~ x * y, data = xy_grid, 
          drape = TRUE, ylab = expression(tau), xlab = "Demand", zlab = "True Loss",
          col.regions = colorRampPalette(c("red", "blue"))(100))

### BUILD NON-LINEAR DEMAND TIME SERIES
###
###

# ----------------------------------------------------------- #
# Regime-Switching Markov                                     #
# ----------------------------------------------------------- #
#                                                             
# Let (D_1, \ldots, D_n) denote the demand time series.       
# We assume that at any time t,                               
#
#           D_t = \mu_{R_t} + \varepsilon_{t}
#
# where R_t \in \{1, \ldots, K\} is the regime at time t.
#
# The regimes \{ R_t \} follows a Markov-chain with unknown
# transition probabilities. 
#
# ----------------------------------------------------------- #

# number of regimes
K <- 3

# transition kernel
transition_kernel <- matrix(
  data = c(
    0.5, 0.4, 0.1,
    0.4, 0.5, 0.1,
    0.2, 0.3, 0.5
  ),
  nrow = K, ncol = K
)

# \mu vector
mu_vector <- c(10, 7, 2)

# noise level
sigma <- 1

# number of observations
n <- 10000
n_train <- floor(n/2)
                 
# prepare simulated data
sim <- data.frame('t' = 1:n, 'regime' = NA, 'demand' = NA)

# simulate regime
regime_prior <- rep(1/K, K)
sim$regime[1] <- sample(x = 1:K, size = 1, prob = regime_prior)
for(i in 2:n){
    sim$regime[i] <- sample(x = 1:K, size = 1, prob = transition_kernel[sim$regime[i-1],])
}

# simulate demand
sim$demand <- mu_vector[sim$regime] + rnorm(n = n, mean = 0, sd = sigma)

# plot training data
plot.ts(sim[1:50,])

### ESTIMATE REGIME-SWITCHING MODEL
###
###

### EM estimation of gaussian mixture with sigma = 1, K = 3
library(mixtools)
demand_em <- normalmixEM(x = sim$demand[1:n_train], sigma = 1, k = K, sd.constr = sigma)

# estimate of mu_vector
demand_em$mu

# demand appartenance to regime index
demand_em$posterior

#
posterior <- demand_em$posterior[,order(demand_em$mu,decreasing = TRUE)]

# kernel estimate of transition matrix
transition_est <- array(data = NA, dim = c(K,K,2))
for(i in 1:K){
  for(j in 1:K){
    transition_est[i,j,1] <- sum(posterior[1:(n_train-1),i] * posterior[2:n_train,j]) / sum(posterior[1:(n_train-1),i])
    transition_est[i,j,2] <- sum((sim$regime[1:(n_train-1)]==i)*(sim$regime[2:n_train]==j)) / sum((sim$regime[1:(n_train-1)]==i))
  }
}
transition_est[,,1] # classic count estimate
transition_est[,,2] # oracle
transition_kernel

### model-based cost optimization !
###
###

# now find demand value d such that $\tau^* (d) = tau_star_sq$.
tau_star <- function(d){
  tau_vector <- seq(from = tau_min, to = tau_max, length.out = 500)
  res <- sapply(X = tau_vector, FUN = function(tau){
    loss_function(x = tau / tau_norm, y = d / yy_norm)
  })
  tau_opt <- tau_vector[which(res == min(res))]
  return(tau_opt)
}
range(sim$demand[1:n_train])
demand_min <- min(sim$demand[1:n_train])
demand_max <- max(sim$demand[1:n_train])

# problem constants
yy_norm <- max(mu_vector) - min(mu_vector)
tau_min <- 0.2
tau_max <- 0.9
tau_norm <- tau_max - tau_min
tau_vector <- seq(from = tau_min, to = tau_max, length.out = 500)

#
demand_vector <- seq(from = demand_min, to = demand_max, length.out = 100)
tau_star_eval <- sapply(X = demand_vector, FUN = function(d){
  tau_star(d = d)
})

# simulate from posterior distribution (yy sample), in order to compute expectations
val <- 7
val_likelihood <- dnorm(x = val, mean = demand_em$mu, sd = demand_em$sigma)
val_likelihood <- val_likelihood / sum(val_likelihood)
next_val_likelihood <- val_likelihood %*% transition_est[,,2]
yy <- NULL
for(i in 1:K){
  temp <- rnorm(n = ceiling(1000 * next_val_likelihood[i]), 
                mean = demand_em$mu[i], sd = demand_em$sigma[i])
  yy <- c(yy, temp)  
}
range(yy)

# find tau^{**} optimal value
res <- sapply(X = tau_vector, FUN = function(tau){
  mean(loss_function(x = tau / tau_norm, y = yy / yy_norm))
})
plot(x = tau_vector, y = res, type = 'l')
tau_star_sq <- tau_vector[which(res == min(res))]

# find OK demand levels
temp <- abs(tau_star_eval - tau_star_sq)
ok_demand_levels <- which(temp == min(temp))

# plot finding
plot(x = demand_vector, y = tau_star_eval,type = 'l')
abline(h = tau_star_sq, col = 'red')


### SEMI-MARKOV MODELING
### LAWS ON PERIODS (STATIONARY LAW).
###

### DECISION-TREE OPTIMIZED BY TRUE LOSS

#### quantile gradient boosting : choice of quantile

### OPERA-CONTEXTUAL WEIGHTS
### genere scenario

# GENERATE SCENARIO
# ------------- #
#
# (features, target)
# construct K predictors
# optimal weights for target prediction
# conditionnal bagging
#
# pondération fixe 1/K : robuste
#
# utiliser AR confidence interval for sanity check 

bsts
