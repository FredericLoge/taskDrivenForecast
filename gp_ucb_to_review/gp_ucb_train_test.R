# define observation function
get_noisy_loss <- function(tau, d, s2){
  fx <- cost_foo(pars = problem_constants, D = d, tau = tau)
  noise <- rnorm(n = 1, mean = 0, sd = sqrt(s2))
  return(fx + noise)
}

# define kernel
# only choice for now: gaussian kernel 
# with common bandwidth parameter
my_kernel <- function(x0, x1, type = 'exp', l, scales){
  exp( -(2*l^2)^{-1} * sum(((x0 - x1)/scales)^2) )
}
KERNEL_VALUE_EQUALS <- my_kernel(x0 = rep(1,2), x1 = rep(1,2), l = 1, scales = 1)

# observe training sample, uniformly sampled
Ntrain <- 5000
XmatTrain <- array(data = NA, dim = c(Ntrain, 2))
LvecTrain <- array(data = NA, dim = c(Ntrain, 1))
for(i in 1:Ntrain){
  ei <- rnorm(n = 1, mean = 0, sd = 12)
  if(i==1){
    XmatTrain[i,2] <- true_theta_par[1] + true_theta_par[2] * ei
  }else{
    XmatTrain[i,2] <- true_theta_par[1] + true_theta_par[2] * XmatTrain[i-1,2] + ei
  }
  XmatTrain[i,1] <- runif(n = 1, min = 0.2, max = 0.9)
  LvecTrain[i,1] <- get_noisy_loss(tau = XmatTrain[i,1], d = XmatTrain[i,2], s2 = 1e-04)
}

# prepare nice kernel interpolation
kernel_bandwith <- 0.1
Xmin <- with(problem_constants, c(tau_min, P_min))
Xmax <- with(problem_constants, c(tau_max, P_max))
Kmat <- array(data = NA, dim = c(Ntrain,Ntrain))
for(i in 1:Ntrain){
  Kmat[i,i] <- 1
  temp <- apply(XmatTrain[-i,], 1, function(x){
    my_kernel(x0 = x, x1 = XmatTrain[i,], l = kernel_bandwith, scales = sqrt(Xmax-Xmin))
  })
  Kmat[i,-i] <- Kmat[-i,i] <- temp
}
summary(c(Kmat))
sigma2 <- 1
normi <- solve(Kmat + diag(Ntrain) * sigma2)
normi_times_LvecTrain <- normi %*% LvecTrain

#
compute_distance_vector <- function(new_x){
  apply(X = XmatTrain, MARGIN = 1, FUN = function(x){
    my_kernel(x0 = new_x, x1 = x, l = kernel_bandwith, scales = sqrt(Xmax-Xmin))
  })
}

#
NB_GRID_1 <- 20
estimated_grid <- expand.grid(
  'tau' =  with(problem_constants, seq(from = tau_min, to = tau_max, length.out = NB_GRID_1)),
  'demand' = with(problem_constants, seq(from = P_min, to = P_max, length.out = NB_GRID_1))
)
estimated_grid$mu_estimate <- NA
estimated_grid$var_estimate <- NA
for(i in 1:nrow(estimated_grid)){
  cat('\t', i, '/', nrow(estimated_grid))
  kTi <- compute_distance_vector(new_x = estimated_grid[i,1:2])
  estimated_grid$mu_estimate[i] <- as.numeric( t(kTi) %*% normi_times_LvecTrain )
  estimated_grid$var_estimate[i] <- KERNEL_VALUE_EQUALS - t(kTi) %*% normi %*% kTi
  if(estimated_grid$var_estimate[i] < 0) stop('Negative Variance ...')
}
levels_tau <- estimated_grid$tau[1:NB_GRID_1]
g0 <- ggplot(data = estimated_grid) +
  geom_tile(mapping = aes(x = demand, y = tau, fill = mu_estimate)) + 
  scale_fill_gradient(low = 'blue', high = 'red')
g1 <- ggplot(data = estimated_grid) +
  geom_tile(mapping = aes(x = demand, y = tau, fill = var_estimate)) + 
  scale_fill_gradient(low = 'blue', high = 'red')
library(gridExtra)
grid.arrange(g0,g1)

estimated_grid$true_loss <- NA
for(i in 1:nrow(estimated_grid)){
  cat('\t', i, '/', nrow(estimated_grid))
  estimated_grid$true_loss[i] <- cost_foo(pars = pars, 
                                          D = estimated_grid$demand[i],
                                          tau = estimated_grid$tau[i])
}


z <- reshape2::dcast(data = estimated_grid, formula = tau~demand, value.var = "mu_estimate")
str(z)
rownames(z) <- z$tau
z <- z[,-1]
colnames(z) <- paste0('x', 1:ncol(z))
z <- as.matrix(z)
str(z)
library(lattice)
wireframe(z, drape = TRUE, ylab = expression(tau), xlab = "Demand", zlab = "Loss estimate",
          col.regions = colorRampPalette(c("red", "blue"))(100))

estimated_grid$mu_estimate_opt <- FALSE
for(d in unique(estimated_grid$demand)){
  cond <- which(estimated_grid$demand == d)
  index <- which(estimated_grid$mu_estimate[cond] == min(estimated_grid$mu_estimate[cond]))
  cond <- cond[index]
  estimated_grid$mu_estimate_opt[cond] <- TRUE
}
table(estimated_grid$mu_estimate_opt)
ggplot(data = estimated_grid) +
  geom_tile(mapping = aes(y = tau, x = demand, fill = mu_estimate)) +
  geom_point(data = estimated_grid[estimated_grid$mu_estimate_opt,], mapping = aes(y = tau, x = demand), col  ='red', cex = 10) +
  scale_fill_gradient(low = 'blue', high = 'red') +
  xlab('Demand') +
  ylab(expression(tau)) +
  labs(fill = 'Loss estimate') +
  theme(legend.position = 'top', legend.key.size = unit(x = 2, units = 'cm'), text = element_text(size = 20),
        axis.title.x = element_text(size = 25, angle = 0, vjust = 1/2, margin = margin(r = 20, t = 20)),
        axis.title.y = element_text(size = 25, angle = 0, vjust = 1/2, margin = margin(r = 20, t = 20))
        )

estimated_grid$ucb <- estimated_grid$mu_estimate + sqrt(estimated_grid$var_estimate)
estimated_grid$ucb_opt <- FALSE
for(d in unique(estimated_grid$demand)){
  cond <- which(estimated_grid$demand == d)
  index <- which(estimated_grid$ucb [cond] == min(estimated_grid$ucb[cond]))
  cond <- cond[index]
  estimated_grid$ucb_opt[cond] <- TRUE
}
table(estimated_grid$ucb_opt)
ggplot(data = estimated_grid) +
  geom_tile(mapping = aes(y = tau, x = demand, fill = ucb)) +
  geom_point(data = estimated_grid[estimated_grid$ucb_opt,], mapping = aes(y = tau, x = demand), col  ='red', cex = 10) +
  scale_fill_gradient(low = 'blue', high = 'red') +
  xlab('Demand') +
  ylab(expression(tau)) +
  labs(fill = 'UCB estimate') +
  theme(legend.position = 'top', legend.key.size = unit(x = 2, units = 'cm'), text = element_text(size = 20),
        axis.title.x = element_text(size = 25, angle = 0, vjust = 1/2, margin = margin(r = 20, t = 20)),
        axis.title.y = element_text(size = 25, angle = 0, vjust = 1/2, margin = margin(r = 20, t = 20))
  )

estimated_grid$loss_opt <- FALSE
for(d in unique(estimated_grid$demand)){
  cond <- which(estimated_grid$demand == d)
  index <- which(estimated_grid$true_loss[cond] == min(estimated_grid$true_loss[cond]))
  cond <- cond[index]
  estimated_grid$loss_opt[cond] <- TRUE
}
table(estimated_grid$loss_opt)
ggplot(data = estimated_grid) +
  geom_tile(mapping = aes(y = tau, x = demand, fill = true_loss)) +
  geom_point(data = estimated_grid[estimated_grid$loss_opt,], mapping = aes(y = tau, x = demand), col  ='red', cex = 10) +
  scale_fill_gradient(low = 'blue', high = 'red') +
  xlab('Demand') +
  ylab(expression(tau)) +
  labs(fill = 'True Loss') +
  theme(legend.position = 'top', legend.key.size = unit(x = 2, units = 'cm'), text = element_text(size = 20),
        axis.title.x = element_text(size = 25, angle = 0, vjust = 1/2, margin = margin(r = 20, t = 20)),
        axis.title.y = element_text(size = 25, angle = 0, vjust = 1/2, margin = margin(r = 20, t = 20))
  )

plot(XmatTrain)

range(estimated_grid$demand)
XmatTrainDF <- data.frame(XmatTrain)
str(XmatTrainDF)
XmatTrainDF <- XmatTrainDF[XmatTrainDF$X1]
ggplot(data = XmatTrainDF,mapping = aes(x = X2, y = X1)) +
  geom_bin2d() +
  xlab('Demand') +
  ylab(expression(tau)) +
  labs(fill = 'Obs Count') +
  theme(legend.position = 'top', legend.key.size = unit(x = 2, units = 'cm'), text = element_text(size = 20),
        axis.title.x = element_text(size = 25, angle = 0, vjust = 1/2, margin = margin(r = 20, t = 20)),
        axis.title.y = element_text(size = 25, angle = 0, vjust = 1/2, margin = margin(r = 20, t = 20))
  ) +
  xlim(range(estimated_grid$demand)) +
  ylim(range(estimated_grid$tau))

  
#
new_foo <- function(theta, alpha = 35){
  losses <- sapply(seq(2, Ntrain, 20), function(jj){
    PredDemand <- theta[1] + theta[2] * XmatTrain[jj-1,2]
    temp <- (abs(estimated_grid$demand - PredDemand))
    crit <- estimated_grid[temp == min(temp),]
    tau_star <- crit$tau[which.min(crit$mu_estimate - 4 * alpha * sqrt(crit$var_estimate))]
    ### tau_star <- levels_tau[as.numeric(which.min(crit))]
    kTi <- compute_distance_vector(new_x = c(tau_star, XmatTrain[jj,2]))
    lossi <- as.numeric( t(kTi) %*% normi_times_LvecTrain )
    lossi <- lossi - 4 * alpha * sqrt( t(kTi) %*% normi %*% kTi )
    return(lossi) 
  })
  return(mean(losses))
}
# theta_est <- optim(par = c(0,0), fn = new_foo)
# theta_est <- optim(par = c(0,0), fn = new_foo, method = 'L-BFGS-B', lower = c(0,0), upper = c(100, 1.5))
# theta_est$par
theta_prop <- expand.grid(
  't0' = seq(from = 0, to = 100, length.out = 10),
  't1' = seq(from = 0, to = 1.5, length.out = 10)
)
theta_prop <- expand.grid(
  't0' = seq(from = 0, to = 100, length.out = 5),
  't1' = seq(from = 0.5, to = 1, length.out = 5)
)
theta_prop$val <- NA
### ORACLE new function
# new_foo <- function(theta){
#   losses <- sapply(2:Ntrain, function(jj){
#     PredDemand <- theta[1] + theta[2] * XmatTrain[jj-1,2]
#     tau_star <- optim(par = 0.5, fn = function(tau){
#       get_loss(tau = tau, d = PredDemand)
#     }, lower = 0, upper = 1, method = 'L-BFGS-B')
#     cost_foo(pars = problem_constants, D = XmatTrain[jj,2], tau = tau_star$par)
#   })
#   return(mean(losses))
# }
for(i in 1:nrow(theta_prop)){
  cat('\t', i, '/', nrow(theta_prop))
  theta_prop$val[i] <- new_foo(theta = as.numeric(theta_prop[i,1:2]), alpha = 1000) 
}
ggplot(data = theta_prop) +
  geom_tile(mapping = aes(x = t0, y = t1, fill = val)) + 
  geom_text(mapping = aes(x = t0, y = t1, label = round(val))) + 
  scale_fill_gradient(low = 'blue', high = 'red')
theta_prop[order(theta_prop$val)[1:10],]

#
xx <- cbind(1, XmatTrain[-nrow(XmatTrain),2])
theta_est_1 <- as.numeric(solve(t(xx) %*% xx) %*% t(xx) %*% XmatTrain[-1,2])
index <- 0
index <- index + 1
theta_est_2 <- as.numeric( theta_prop[(order(theta_prop$val))[index], 1:2] )

# we observe that with standard method proposed, 
# we aren't able to propose a better solution
theta_label <- function(x){
  apply(x, 1, function(xx){
    paste0('(', round(xx[1]),"; ", round(xx[2],2), ')')
  })
}

####
#### this does not work.
####

#
theta_est <- theta_est_2 # c(70, 0.87)
alpha=1000
set.seed(38765)
Ntest <- 500
XmatTest <- array(data = NA, dim = c(Ntest, 2))
LvecTest <- array(data = NA, dim = c(Ntest, 2))
for(i in 1:Ntest){
  ei <- rnorm(n = 1, mean = 0, sd = 12)
  if(i==1){
    XmatTest[i,2] <- true_theta_par[1] + true_theta_par[2] * XmatTrain[nrow(XmatTrain),2] + ei
    PredDemand <- theta_est[1] + theta_est[2] * XmatTrain[nrow(XmatTrain),2]
  }else{
    XmatTest[i,2] <- true_theta_par[1] + true_theta_par[2] * XmatTest[i-1,2] + ei
    PredDemand <- theta_est[1] + theta_est[2] * XmatTest[i-1,2]
  }
  temp <- (abs(estimated_grid$demand - PredDemand))
  crit <- estimated_grid[temp == min(temp),]
  tau_star <- crit$tau[which.min(crit$mu_estimate - alpha * sqrt(crit$var_estimate))]
  ### tau_star <- crit$tau[which.min(crit$mu_estimate)]
  # tau_star <- optim(par = 0.5, fn = function(tau){
  #   get_loss(tau = tau, d = PredDemand)
  # }, lower = 0, upper = 1, method = 'L-BFGS-B')
  # temp <- (abs(estimated_grid$demand - PredDemand))
  # crit <- estimated_grid$mu_estimate[temp == min(temp)]
  # XmatTest[i,1] <- levels_tau[as.numeric(which.min(crit))]
  XmatTest[i,1] <- tau_star
  LvecTest[i,1] <- get_loss(tau = XmatTest[i,1], d = XmatTest[i,2])
  opt <- optim(par = 0.5, fn = function(tau){
    get_loss(tau = tau, d = XmatTest[i,2])    
  }, lower = 0, upper = 1, method = 'L-BFGS-B')
  LvecTest[i,2] <- get_loss(tau = opt$par, d = XmatTest[i,2])
}
# plot.ts(cumsum(LvecTest[,1] - LvecTest[,2]), lwd = 3, col = 'black')
lines(cumsum(LvecTest[,1] - LvecTest[,2]), lwd = 3, col = 'red')
legend('topleft', legend = paste0(c('Agnostic ', 'Global    '), theta_label(rbind(theta_est_1, theta_est_2))), 
       lty = 1, lwd = 3, cex = 2, col = c('black', 'red'))
title(main = 'Cumulative regret on test set, after N = 5000 training samples, with unknown loss function.')


plot(LvecTest)
abline(0,1)
# kTi <- compute_distance_vector(new_x)
# mui <- as.numeric( t(kTi) %*% normi_times_LvecTrain )
# sigma2i <- KERNEL_VALUE_EQUALS - diag( t(kTi) %*% normi %*% kTi )

trainData <- data.frame(XmatTrain)
g0 <- ggplot(data = estimated_grid) +
  geom_tile(mapping = aes(x = demand, y = tau, fill = mu_estimate)) +
  geom_point(data = trainData, mapping = aes(x = X2, y = X1)) +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') +
  xlim(min(estimated_grid$demand), max(estimated_grid$demand))
g1 <- ggplot(data = pars_grid) +
  geom_tile(mapping = aes(x = d, y = tau, fill = val)) +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red')
grid.arrange(g0,g1)
