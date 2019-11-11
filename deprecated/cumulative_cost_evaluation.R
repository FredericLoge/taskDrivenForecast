set.seed(235667)
BBB <- 100
seeds <- floor(runif(n = BBB, min = 1, max = 1e06))
theta_par <- 0.00

val <- sapply(X = 1:BBB, FUN = function(s){
  
  set.seed(seeds[s])
  
### PROBLEM CONSTANTS FROM AFOREMENTIONED EQUATIONS ------------------------------------------------------------------------

# set of known constants in equations (1), (2), (3) and constraints on \tau
problem_constants <- list(
  'a' = 300, # slope of (1)
  'b' = 100, # intercept of (1)
  'c' = 0.9, # production efficiency of our units (2)
  'd' = 10,  # unit cost of input energy (3)
  'e' = 100,  # unit cost / penalty due to lack of supply (3)
  'f' = 20,   # unit reward per production output
  'tau_min' = 0.2, # minimal functioning level
  'tau_max' = 0.9  # maximal functioning level
)
pars <- problem_constants

# add to set problem_constants the bounds of production capacity
problem_constants$P_min <- with(problem_constants, c*(a*tau_min + b))
problem_constants$P_max <- with(problem_constants, c*(a*tau_max + b)) 


### NUMBER OF TIME UNITS ---------------------------------------------------------------------------------------------------

# number of time units
n <- 500


### SIMULATE DEMAND SIGNAL UNDER AR(p) -------------------------------------------------------------------------------------

# paramaters of the non-zero mean AR(p) model for the demand signal
demand_parameters <- list('theta' = c(theta_par), 'sigma2e' = 1)
demand_parameters$p <- length(demand_parameters$theta)

# simulate demand signal Y_{t} = \mu + \sum_{i \in [p]} \theta_i Y_{t-i} + e_t
# an AR(p) well normalized to fit our production capabilities ranges : we assume the unit is "well"-sized
y <- sim_ar_p(n = n, sigma2e = demand_parameters$sigma2e, theta = demand_parameters$theta)
y_low <-  - 1 * sqrt( 2 * log(2 * n / 0.01) / (1 - 0.8^2))
y_high <- (- y_low)
y <- with(problem_constants, P_min + (y - y_low) / (y_high - y_low) * (P_max - P_min))

### SOLVE TWO OPTIM PROBLEMS : "BEST" FORECASTERS --------------------------------------------------------------------------

# generate "forecast" matrix
mat <- gen_forecast_matrix(y = y, p = demand_parameters$p)

# select train-test split ## FUNCTION!
index_train <- 1:floor(n/2)
line_contains_na <- apply(mat[index_train,], 1, function(x) any(is.na(x)))
index_train <- index_train[!line_contains_na]
index_test <- (floor(n/2)+1):n
line_contains_na <- apply(mat[index_test,], 1, function(x) any(is.na(x)))
index_test <- index_test[!line_contains_na]

# find \theta^*_0 := argmin_{\theta} L_0(Y, X\theta), where
# L_0(y, y^) is the L2-norm of vector y-y^
theta_star_0 <- optim_L0(y = mat[index_train, 'y'], x = mat[index_train, -1])

# define grid search controls for optimization of L1 metric
# also search controls
theta_search_grid_controls <- list(
  'p' = demand_parameters$p+1,
  'lower_bound' = c(100, rep(-2, demand_parameters$p)),
  'upper_bound' = c(500, rep(+2, demand_parameters$p)),
  'nb_points' = rep(100, demand_parameters$p+1)
)
theta_search_grid_controls <- list(
  'p' = demand_parameters$p+1,
  'lower_bound' = c(-300, rep(-2, demand_parameters$p)),
  'upper_bound' = c(300, rep(2, demand_parameters$p)),
  'nb_points' = rep(100, demand_parameters$p+1)
)

# A DEPLACER DANS LA PARTIE VIZ -> AUTRE SCRIPT ?
# find \theta^*_1 := argmin_{\theta} L_1(Y, X\theta), where
# L_1(y, y^) := L(\tau(y^), y) is the cost defined in (3)
theta_star_1 <- optim_L1(y = mat[index_train, 'y'], x = mat[index_train, -1], search_grid_controls = theta_search_grid_controls)

# direct gradient optimization
opt_L1_res <- optim_grad_L1(y = mat[index_train, 'y'], x = mat[index_train, -1], 
                            lb = theta_search_grid_controls$lower_bound, up = theta_search_grid_controls$upper_bound, 
                            setup_pars = problem_constants, init_pars = theta_star_0$theta_hat)
opt_L1_res$par

### GRAPHIC REPRESENTATION OF THE LOSS FUNCTIONS ---------------------------------------------------------------------------

# label vectors, called in the creation of g just below
v1 <- c('L0 loss', 'L1 loss')
v2 <- c('Opt for L0 loss', 'Opt for L1 loss')

# represent uni-plot, normalized loss functions

### EVALUATION ON THE TEST SET ---------------------------------------------------------------------------------------------

# ADD FUNCTIONS COMPUTING LOSS FUNCTIONS

# nb model parameters
q <- 2

# reframe thetas
theta_stars <- rbind(
  theta_star_0$theta_hat,
  theta_star_1$theta_hat,
  rep(NA, q)
)
theta_stars$loss_test <- 0
theta_stars$L0loss_test <- 0

# get (x,y) test
ytest = mat[index_test, 'y']
xtest = mat[index_test, -1]

# compute test errors and add to theta_stars
y_hat_test <- list()
for(i in 1:(nrow(theta_stars))){
  y_hat <- xtest %*% as.numeric(theta_stars[i,1:q])
  y_hat_test[[i]] <- y_hat
  theta_stars$loss_test[i] <- eval_theta(D = ytest, Dhat = y_hat)
  theta_stars$L0loss_test[i] <- mean((ytest-y_hat)^2)
}
i = nrow(theta_stars)
y_hat <- ytest
y_hat_test[[i]] <- y_hat
theta_stars$loss_test[i] <- eval_theta(D = ytest, Dhat = y_hat)
theta_stars$L0loss_test[i] <- mean((ytest-y_hat)^2)

# # check model validity on train
# ytrain <- mat[index_train, 'y']
# xtrain <- mat[index_train, -1]
# y_hat_train <- xtrain %*% as.numeric(theta_stars[1, 1:q])
# hist(y_hat_train - ytrain)

### add evaluation regarding performance of the system
# Note that P <- c*(a*\tau + b), combining (1) and (2). As such, the range of demand which the
# production may handle is [c*(a*\tau_{min} + b) ; c*(a*\tau_{max} + b)].

# compute test errors and add to theta_stars
res_df <- list() 
for(i in 1:(0 + nrow(theta_stars))){
  temp <- data.frame('method' = c('L0_Loss', 'L1_Loss', 'ORACLE')[i], 'ytest' = ytest, 't' = 1:length(ytest))
  if(i == nrow(theta_stars)){
    y_hat <- ytest
  }else if(i == 1 + nrow(theta_stars)){
    y_hat <- as.numeric(predict(rf0, newdata = data.frame(xtest)))
  }else{
    y_hat <- xtest %*% as.numeric(theta_stars[i,1:q])
  }
  temp$y_hat <- y_hat
  ## as.numeric(xtest %*% as.numeric(theta_stars[i,1:q]))
  temp$tauSel <- optimal_tau(pars = pars, D = temp$y_hat)
  temp$costs <- cost_foo(pars = pars, D = ytest, tau = temp$tauSel)
  temp$cumsum_costs <- cumsum(temp$costs)
  res_df[[i]] <- temp
}
res_df <- do.call(rbind.data.frame, res_df)

return(res_df$cumsum_costs[res_df$t == max(res_df$t)])

})

val2 <- data.frame(
  'method' = rep(c('Agnostic', 'Global', 'Oracle'), each = 100),
  'value' = c(val[1,], val[2,], val[3,])
)
library(ggplot2)
ggplot(data=val2) +
  geom_density(mapping = aes(x = value/250, fill = method), alpha = 0.3) +
  ggtitle(label = 'Cumulative cost density (/250), per method, evaluated on 100 simulated trajectories.') + 
  labs(col = expression(Customer~demand~italic(d))) +
  xlab(label = '') +
  ylab(label = '') +
  ggtitle(label = expression(Density~of~frac(1,"T'-T")~sum(tilde(L)~group("(",list(tau^"*"~(hat(f)~(D[t^"-"])),D[t]),")"), t==T+1, "T'"))) +
  labs(fill = 'Method') +
  theme(
    axis.text = element_text(size = 20),
    axis.title.x = element_text(size = 20, margin = margin(t = 20)),
    axis.title.y = element_text(size = 20, margin = margin(r = 20), angle = 0, vjust = 0.5),
    plot.title = element_text(size = 30, margin = margin(b = 20))
  )

val3 <- data.frame(
  'method' = rep(c('Agnostic', 'Global'), each = 100),
  'value' = c(val[1,]-val[3,], val[2,]-val[3,])
)
ggplot(data=val3) +
  geom_density(mapping = aes(x = value/250, fill = method), alpha = 0.3) +
  ggtitle(label = 'Cumulative cost density (/250), per method, evaluated on 100 simulated trajectories.') + 
  labs(col = expression(Customer~demand~italic(d))) +
  xlab(label = '') +
  ylab(label = '') +
  ggtitle(label = 'Cumulative Regret') + 
#  ggtitle(label = expression(Density~of~frac(1,"T'-T")~sum(tilde(L)~group("(",list(tau^"*"~(hat(f)~(D[t^"-"])),D[t]),")") - tilde(L)~group("(",list(tau^"*"~(D[t]),D[t]),")"), t==T+1, "T'"))) +
  labs(fill = 'Method') +
  theme(
    axis.text = element_text(size = 20),
    axis.title.x = element_text(size = 20, margin = margin(t = 20)),
    axis.title.y = element_text(size = 20, margin = margin(r = 20), angle = 0, vjust = 0.5),
    plot.title = element_text(size = 30, margin = margin(b = 20))
  )


str(val)
boxplot(t(val))
boxplot(val[2,] / val[1,])
