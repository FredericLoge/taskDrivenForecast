# Consider the following control problem.
#
# \tau, representing the level of use of the machines in the prod unit is our only control variable.
# \tau is bounded in [\tau_{min} ; \tau_{max}].
#
# Q (the elec quantity needed) and P (the production in output) are given by : 
#   Q <- a*\tau + b     (1)
#   P <- c*Q            (2)
# with (a,b,c) known quantities. (a,b) represent the conversion between the use level and the electric
# needs, (c) represents the production efficiency of our units.
#
# We define the cost function :
#   L(\tau, D) = d*(P-D)_{+} + e*(P-D)_{-}, (3)
# where (x)_{+} = max(0,x), (x)_{-} = max(0,-x) and (d,e) known quantities. d is the unit cost of vented
# molecules (P > D), e is the unit cost / penalty due to lack of supply. 
#
# To complete this control problem, we add the following assumptions :
# - at each decision epoch, we decide of the value of \tau for the next time step
# - the optimal decision would be based on the actual demand on this time step, but since we don't have
#   access to it, we will forecast this demand based on prior demands.
#
# Note that P <- c*(a*\tau + b), combining (1) and (2). As such, the range of demand which the
# production may handle is [c*(a*\tau_{min} + b) ; c*(a*\tau_{max} + b)].

# Here is the content of this R file :
# - problem constants from equations above
# - number of time units
# - simulate demand signal under AR(p)
# - solve two optim problems -> "best" forecasters
# - graphic representation of the loss functions
# - evaluation on test set

# needs installation prior :p
library(taskDrivenRandomForest)
source('data_sim.R')
source('custom_metrics.R')
source('task_driven_lm.R')

### PROBLEM CONSTANTS FROM AFOREMENTIONED EQUATIONS -------------------------------------------------

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

### SIMULATE DATA ---------------------------------------------------------------------------------

# number of time units
n <- 5000
SIGNAL_TYPE = 2
NOISE_TYPE = 2

# simulate signal 
signal <- sim_signal(n = n, signal_type = SIGNAL_TYPE)

# simulate noise 
noise <- sim_noise(n = n, noise_type = 2)

# add signal and noise, weigthed to satisfy some SNR ratio
snr <- 1
k <- sqrt(var(signal) / (var(noise) * snr))
target <- signal + k * noise

# viz target values
plot.ts(target[1:100])
acf(target)
acf(diff(target))

# retrieve pertinent features, depending on chosen signal and noise types
features <- get_features(y = target, signal_type = SIGNAL_TYPE, noise_type = NOISE_TYPE)

# viz features
str(features)
head(features)

# build proper dataframe to work with
mat <- cbind.data.frame(target, features)
mat$contains_na <- rowSums(is.na(mat) > 0)
mat$is_train <- TRUE
mat$is_train[(floor(n/2)+1):n] <- FALSE
mat_col_index_target <- 1
mat_col_index_features <- 1 + 1:ncol(features)

# erase NAs
mat <- na.omit(mat)

### LINEAR MODELS ---------------------------------------------------------------------------------

# find \theta^*_0 := argmin_{\theta} L_0(Y, X\theta), where
# L_0(y, y^) is the L2-norm of vector y-y^
opt_L0_res <- optim_L0(y = mat[mat$is_train, mat_col_index_target], x = mat[mat$is_train, mat_col_index_features])
agnostic_lm_param <- opt_L0_res$theta_hat
  
# define grid search controls for optimization of L1 metric
# also search controls
theta_search_grid_controls_scenario1 <- list(
  'p' = length(mat_col_index_features),
  'lower_bound' = c(0, rep(-200, length(mat_col_index_features)-1)),
  'upper_bound' = c(500, rep(+200, length(mat_col_index_features)-1)),
  'nb_points' = rep(100, length(mat_col_index_features))
)
theta_search_grid_controls_scenario2 <- list(
  'p' = length(mat_col_index_features),
  'lower_bound' = c(0, rep(-1, length(mat_col_index_features)-1)),
  'upper_bound' = c(500, rep(+1, length(mat_col_index_features)-1)),
  'nb_points' = rep(100, length(mat_col_index_features))
)
theta_search_grid_controls <- theta_search_grid_controls_scenario2

# direct gradient optimization
opt_L1_res <- optim_grad_L1(y = mat[mat$is_train, mat_col_index_target], x = mat[mat$is_train, mat_col_index_features], 
                            lb = theta_search_grid_controls$lower_bound, up = theta_search_grid_controls$upper_bound, 
                            setup_pars = problem_constants, init_pars = opt_L0_res$theta_hat)
goal_driven_lm_param <- opt_L1_res$par

### RANDOM FORESTS ---------------------------------------------------------------------------------

source('custom_metrics.R')
custom_loss <- function(y, y_pred){
  eval_theta(D = y, Dhat = y_pred)
}

# random forest hyper-parameters
N_TREES <- 10
MAX_DEPTH <- 8

# fit single CART tree
library(rpart)
agnostic_rpart <- rpart(target ~ ., data = mat[mat$is_train, c(mat_col_index_target, mat_col_index_features)], control = rpart.control(maxdepth = 8))

# fit agnostic random forest
library(randomForest)
agnostic_rf <- randomForest(target ~ ., data = mat[mat$is_train, c(mat_col_index_target, mat_col_index_features)], mtry = 2, ntree = N_TREES)
agnostic_rf

# redo predictions based on initial partition proposed by classic variance splitting
my_forest <- recalibrate_forest_predictions(rfobj = agnostic_rf, 
                                            x_train = mat[mat$is_train, mat_col_index_features], 
                                            y_train = mat[mat$is_train, mat_col_index_target], 
                                            customized_loss_foo = custom_loss)

# fit customized random forest
goal_driven_rf <- build_rf(y = mat[mat$is_train, mat_col_index_target], 
                           x = mat[mat$is_train, mat_col_index_features],
                           customized_loss_foo = custom_loss, nb_points = 200, nb_points_y = 50, min_data_size = 10, 
                           n_trees = N_TREES, max_depth = MAX_DEPTH, bootstrap_prop = 2/3)
# check that we have the right number of trees built
length(goal_driven_rf)
# look at second tree built
str(goal_driven_rf[[2]], 1)

### EVALUATION ON THE TEST SET ---------------------------------------------------------------------------------------------

#
agnostic_rpart_pred <- predict(agnostic_rpart, mat[mat$is_train == FALSE, mat_col_index_features])

# randomForest baseline, interesting to amke comparisons
agnostic_rf_pred <- predict(agnostic_rf, mat[mat$is_train == FALSE, mat_col_index_features])

# recalibrated randomForest, compute predictions for each tree and aggregate with mean
agnostic_to_goal_rf_pred <- predict_from_new_forest(myforest = my_forest, 
                                                    x_test = mat[mat$is_train == F, mat_col_index_features]) 
agnostic_to_goal_rf_pred_mean <- rowMeans(agnostic_to_goal_rf_pred)

# goal driven splits, but classic prediction metric
goal_driven_rf_classic_pred <- t(apply(X = mat[mat$is_train == FALSE, mat_col_index_features], MARGIN = 1, FUN = function(x){
  predict_from_rf(rf = goal_driven_rf, x_vector = x, nb_points_y = 50, customized_loss_foo = sum_squared_errors)
}))
goal_driven_rf_classic_pred_mean <- rowMeans(goal_driven_rf_classic_pred)

# goal driven splits and goal driven prediction
goal_driven_rf_pred <- t(apply(X = mat[mat$is_train == FALSE, mat_col_index_features], MARGIN = 1, FUN = function(x){
  predict_from_rf(rf = goal_driven_rf, x_vector = x, nb_points_y = 50, customized_loss_foo = custom_loss)
}))
goal_driven_rf_pred_mean <- rowMeans(goal_driven_rf_pred)
## apply(goal_driven_rf_pred, 1, sd)

# agnostic linear model prediction
agnostic_lm_pred <- as.numeric( as.matrix( mat[mat$is_train == FALSE, mat_col_index_features] ) %*% agnostic_lm_param )

# task driven linear model prediction
goal_driven_lm_pred <- as.numeric( as.matrix( mat[mat$is_train == FALSE, mat_col_index_features] ) %*% goal_driven_lm_param )

# combine all the predictions in a big dataframe
ypred_test <- cbind.data.frame(
  "oracle" = mat[mat$is_train == F, mat_col_index_target],
  "agnostic_lm" = agnostic_lm_pred,
  "goal_driven_lm" = goal_driven_lm_pred,
  "agnostic_rpart" = agnostic_rpart_pred,
  "agnostic_rf" = agnostic_rf_pred,
  "agnostic_to_goal_rf" = agnostic_to_goal_rf_pred_mean,
  "goal_driven_rf_classic" = goal_driven_rf_classic_pred_mean,
  "goal_driven_rf" = goal_driven_rf_pred_mean
)

# correlation matrix between predictions -> it may be perfectly correlated, and yet there is something to gain - see below !
cor(ypred_test)

# plot test set (pred vs reality)
color_vector <- rainbow(n = ncol(ypred_test))
plot(x = ypred_test$oracle, y = ypred_test$oracle)
for(pred_index in 2:ncol(ypred_test)){
  points(x = ypred_test$oracle, y = ypred_test[,pred_index], col = color_vector[pred_index])
}
legend("bottom", legend = colnames(ypred_test), col = color_vector, lty = 1, lwd = 2)

# nicer plot
library(ggplot2)
temp <- data.table::melt(ypred_test, id.vars = "oracle")
ggplot(data = temp, mapping = aes(x = oracle, y = value, col = variable)) + 
  geom_point(alpha = 0.3)  +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = 'loess')

# compute losses
rmse_loss <- numeric(ncol(ypred_test))
names(rmse_loss) <- colnames(ypred_test)
goal_loss <- rmse_loss
for(pred_index in 1:ncol(ypred_test)){
  goal_loss[pred_index] <- eval_theta(D = ypred_test$oracle, Dhat = ypred_test[,pred_index])
  rmse_loss[pred_index] <- mean((ypred_test$oracle - ypred_test[,pred_index])^2)
}

# quick look
sort(goal_loss)
sort(rmse_loss)

# raw numbers plots
barplot(sort(goal_loss), col = 'blue')
barplot(sort(rmse_loss), col = 'red')

# normalized plots
normalize <- function(x){ (x - min(x)) / (max(x) - min(x)) }
scores <- cbind(goal_loss, rmse_loss)
scores <- scores[order(scores[,1]),]
par(mfrow = c(2,1))
barplot(normalize(scores[,1]), col = 'blue', main = "Task loss")
barplot(normalize(scores[,2]), col = 'red', main = "RMSE")
par(mfrow = c(1,1))
