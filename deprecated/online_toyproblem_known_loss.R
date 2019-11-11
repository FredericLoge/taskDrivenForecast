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

with(problem_constants, c(b * (1-0.8), 0.8))

### NUMBER OF TIME UNITS ---------------------------------------------------------------------------------------------------

# number of time units
n <- 500


### SIMULATE DEMAND SIGNAL UNDER AR(p) -------------------------------------------------------------------------------------

# paramaters of the non-zero mean AR(p) model for the demand signal
demand_parameters <- list('theta' = c(0.8), 'sigma2e' = 1)
demand_parameters$p <- length(demand_parameters$theta)

# simulate demand signal Y_{t} = \mu + \sum_{i \in [p]} \theta_i Y_{t-i} + e_t
# an AR(p) well normalized to fit our production capabilities ranges : we assume the unit is "well"-sized
set.seed(54567)
y <- sim_ar_p(n = n, sigma2e = demand_parameters$sigma2e, theta = demand_parameters$theta)
y_low <- - 1 * sqrt( 2 * log(2 * n / 0.01) / (1 - 0.8^2))
y_high <- (- y_low)
y <- with(problem_constants, P_min + (y - y_low) / (y_high - y_low) * (P_max - P_min))

# theta0 <- 0
# theta1 <- 0.8
# alpha0 <- with(problem_constants, P_min + ( - y_low) / (y_high - y_low) * (P_max - P_min))
# alpha1 <- with(problem_constants,  1 / (y_high - y_low) * (P_max - P_min))
# theta1
# alpha0 + alpha1 * theta0 - alpha0 * theta1
# alpha1
# P_min + (y - y_low) 

# viz demand signal, autocorrelation function on y and y differentiated -> validating simulation
plot.ts(y)
acf(y)
acf(diff(y))

### SOLVE TWO OPTIM PROBLEMS : "BEST" FORECASTERS --------------------------------------------------------------------------

# generate "forecast" matrix
mat <- gen_forecast_matrix(y = y, p = demand_parameters$p)

# select train-test split ## FUNCTION!
index_train <- 1:5 # 1:floor(n/2)
line_contains_na <- apply(mat[index_train,], 1, function(x) any(is.na(x)))
index_train <- index_train[!line_contains_na]
# index_test <- (floor(n/2)+1):n
# line_contains_na <- apply(mat[index_test,], 1, function(x) any(is.na(x)))
# index_test <- index_test[!line_contains_na]

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

#
reset_times <- (1.8)^(c(2:15))
mat
res <- data.frame(
  loss_test_E0 = rep(NA, 500),
  loss_test_E1 = rep(NA, 500)
)
for(time_index in 6:500){
  if(min(abs(time_index - reset_times)) < 1){
    index_train <- 1:(time_index-1)
    line_contains_na <- apply(mat[index_train,], 1, function(x) any(is.na(x)))
    index_train <- index_train[!line_contains_na]
    theta_star_0 <- optim_L0(y = mat[index_train, 'y'], x = mat[index_train, -1])
    theta_star_1 <- optim_L1(y = mat[index_train, 'y'], x = mat[index_train, -1], search_grid_controls = theta_search_grid_controls)
  }
  res$loss_test_E0[time_index] <- eval_theta(D = mat[time_index,1], Dhat = sum(mat[time_index,-1] * (theta_star_0$theta_hat)))
  res$loss_test_E1[time_index] <- eval_theta(D = mat[time_index,1], Dhat = sum(mat[time_index,-1] * (theta_star_1$theta_hat)))
}
colSums(res, na.rm = TRUE)  /(500 - 6)

### GRAPHIC REPRESENTATION OF THE LOSS FUNCTIONS ---------------------------------------------------------------------------

# re-create bi-plots in order to visualize well - only makes sense in the situation of 2 variable estimates
require(lattice)
levelplot(x = loss ~ dim1 + dim2, data = theta_star_1$search_grid)
levelplot(x = L0loss ~ dim1 + dim2, data = theta_star_1$search_grid)

# re-create single plot, taking the minimal value observed under given point-dimension
uni_plot <- aggregate(x = theta_star_1$search_grid[, c('loss', 'L0loss')], 
                      by = list('dim' = theta_star_1$search_grid$dim2), 
                      FUN = min)

# normalize uni-plot
uni_plot$loss <- uni_plot$loss / max(uni_plot$loss)
uni_plot$L0loss <- uni_plot$L0loss / max(uni_plot$L0loss)

# identify best points for each loss function
index_best_loss <- which(uni_plot$loss == min(uni_plot$loss))
index_best_L0loss <- which(uni_plot$L0loss == min(uni_plot$L0loss))

# label vectors, called in the creation of g just below
v1 <- c('L0 loss', 'L1 loss')
v2 <- c('Opt for L0 loss', 'Opt for L1 loss')

# represent uni-plot, normalized loss functions
library(ggplot2)
g <- ggplot(data = uni_plot, aes(x = dim)) +
  # loss functions over chosen dimension
  geom_line(mapping = aes(x = dim, y = loss, group = 'loss', colour = v1[2])) +
  geom_line(mapping = aes(x = dim, y = L0loss, group = 'loss', colour = v1[1])) + 
  # optimal points for each metric
  geom_point(data = uni_plot[index_best_loss,], mapping = aes(x = dim, y = loss, colour = v1[2], shape = v2[2]), cex = 3) +
  geom_point(data = uni_plot[index_best_loss,], mapping = aes(x = dim, y = L0loss, colour = v1[2], shape = v2[2]), cex = 3) +
  geom_point(data = uni_plot[index_best_L0loss,], mapping = aes(x = dim, y = loss, colour = v1[1], shape = v2[1]), cex = 3) +
  geom_point(data = uni_plot[index_best_L0loss,], mapping = aes(x = dim, y = L0loss, colour = v1[1], shape = v2[1]), cex = 3) +
  # aesthetics
  ggtitle(label = 'Loss functions depending on model parameters') +
  theme(legend.justification = "top", text = element_text(size = 12)) +
  xlab(label = '') + ylab(label = '') + labs(colour = '')

# focus on graphic region
g
g + ylim(0, 0.5) + xlim(99, 105)
g + ylim(0, 0.1) + xlim(0, 1)

# fit random forest
library(randomForest)
df_train_rf <-  data.frame(mat[index_train,])
rf0 <- randomForest(y~x1, data = df_train_rf)
rf0

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
# add ORACLE
i = nrow(theta_stars)
y_hat <- ytest
y_hat_test[[i]] <- y_hat
theta_stars$loss_test[i] <- eval_theta(D = ytest, Dhat = y_hat)
theta_stars$L0loss_test[i] <- mean((ytest-y_hat)^2)


# time series plot
plot.ts(ytest, ylim = c(150, 330))
lines(y_hat_test[[1]], col = 'red')
lines(y_hat_test[[2]], col = 'blue')
legend('bottomright', lty = 1, legend = c('True signal', 'OLS', 'Proper training'), col = c('black', 'red', 'blue'))

# error histogram plot - clearly shows how we have a tendency to overestimate true demand values, in order to cover the risk properly !
par(mfrow = c(2,1))
hist(ytest - y_hat_test[[1]], breaks = seq(-100, +100, 10))
hist(ytest - y_hat_test[[2]], breaks = seq(-100, +100, 10))
par(mfrow = c(1,1)) # RESET !

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
str(res_df)

res_df$cumsum_costs[res_df$t == max(res_df$t)]

# plot of Y predictions vs realized
library(ggplot2)
g1a <- ggplot(data = res_df) +
  geom_line(mapping = aes(x = t, y = ytest, col = 'True Demand'), colour = 'black') +
  geom_line(mapping = aes(x = t, y = y_hat, col = method)) +
  ggtitle('True Demand and Predictions')
g1b <- ggplot(data = res_df, mapping = aes(x = ytest, y = y_hat, col = method)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  geom_abline(slope = 1, intercept = 0) +
  ggtitle('True Demand versus Predictions')
# plot of chosen control \tau
g2a <- ggplot(data = res_df) +
  geom_line(mapping = aes(x = t, y = tauSel, col = method)) +
  ggtitle(label = 'Choices of tau over test period, temporal view, per method.')
g2b <- ggplot(data = res_df) +
  geom_boxplot(mapping = aes(x =  method, y = tauSel, fill = method)) +
  ggtitle(label = 'Choices of tau over test period, per method.')
# plot of costs
g3a <- ggplot(data = res_df) +
  geom_line(mapping = aes(x = t, y = costs, col = method)) +
  ggtitle(label = 'Costs realized over test period, temporal view, per method.')
g3b <- ggplot(data = res_df) +
  geom_density(mapping = aes(x = costs, fill = method), alpha = 0.5, col = NA) +
  ggtitle(label = 'Costs density over test period, per method.')
# cumulated costs
g4 <- ggplot(data = res_df) +
  geom_line(mapping = aes(x = t, y = cumsum_costs / max(abs(cumsum_costs)), col = method)) +
  ggtitle(label = paste0('Cumulated costs - normalized to 1 - nominal : ', floor(max(res_df$cumsum_costs)), ', per method.'))

library(gridExtra)  
grid.arrange(g1a, g1b, g2a, g2b, g3a, g3b, g4, layout_matrix = rbind(c(1,3,5,5),c(2,4,6,7)))
grid.arrange(g1a, g1b, g2a, g2b)
grid.arrange(g3a, g3b, g4, layout_matrix = rbind(c(1,1),c(2,3)))
