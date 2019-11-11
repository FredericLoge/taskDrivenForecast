set.seed(2854)

# define observation function
get_loss <- function(tau, d){
  fx <- cost_foo(pars = problem_constants, D = d, tau = tau)
  return(fx)
}
add_noise <- function(s2){
  return(rnorm(n = 1, mean = 0, sd = sqrt(s2)))
}

# define kernel
# only choice for now: gaussian kernel 
# with common bandwidth parameter
my_kernel <- function(x0, x1, type = 'exp', l, scales){
  exp( -(2*l^2)^{-1} * sum(((x0 - x1)/scales)^2) )
}
KERNEL_VALUE_EQUALS <- my_kernel(x0 = rep(1,2), x1 = rep(1,2), l = 1, scales = 1)

###
###
###

# define bounds of (\tau, demand)
Xmin <- with(problem_constants, c(tau_min, P_min))
Xmax <- with(problem_constants, c(tau_max, P_max))

# number of time steps
N <- 1000

# dimension of x (\tau, demand)
p <- 2

# prepare matrices
Kmat <- array(data = NA, dim = c(N,N))
lossVec <- array(data = NA, dim = c(N,2))
Xmat <- array(data = NA, dim = c(N, p))

# choose observation error level for the loss function
sigma2 <- 1^2

# theta parameter for the demand 
true_theta_par <- c(50, 0.8)

# observation error the demand
d_sigma2 <- 12

# kernel bandwidth
kernel_bandwith <- 0.1

# number of candidates sampled
B <- 20

# parameters of UCB
alpha <- sqrt(1000)
delta <- 0.0005

# random initial point
x0 <- runif(n = 2, min = Xmin, max = Xmax)
Xmat[1,] <- x0
lossVec[1,1] <- get_loss(tau = x0[1], d = x0[2])
lossVec[1,2] <- lossVec[1,1] + add_noise(s2 = sigma2)
Kmat[1,1] <- 1

## keep_computing_lambda <- TRUE

#
method_f_hat = 'NOT_AGNOSTIC'

#
losses <- array(data = NA, dim = c(N,2))
theta_est <- array(data = NA, dim = c(N,p))

#
theta_search_grid_controls <- list(
  'p' = demand_parameters$p+1,
  'lower_bound' = c(0, rep(0, demand_parameters$p)),
  'upper_bound' = c(100, rep(1.5, demand_parameters$p)),
  'nb_points' = rep(100, demand_parameters$p+1)
)

# minimal time
T0 <- 6

# status update
pb <- txtProgressBar(min = 1, max = N, initial = 1)

#
init_par <- c(0,0)
# setup_theta <- expand.grid(
#   'theta0' = seq(from = 0, to = 100, length.out = 5),
#   'theta1' = seq(from = 0, to = 1.5, length.out = 5)
# )
setup_theta <- data.frame(
  'theta0' = runif(n = 25, min = theta_search_grid_controls$lower_bound[1], max = theta_search_grid_controls$upper_bound[1]),
  'theta1' = runif(n = 25, min = theta_search_grid_controls$lower_bound[2], max = theta_search_grid_controls$upper_bound[2])
)

#
timepoint <- array(data = 0, dim = c(N,1))
starttime <- Sys.time()
timepoint[1] <- starttime

for(i in 2:N){
  
  # status update
  setTxtProgressBar(pb, value = i)
  
  # (X_T + \sigma^2 I)^{-1}
  normi <- solve(Kmat[1:(i-1), 1:(i-1)] + diag(i-1) * sigma2)
  normi_times_lossVec <- normi %*% lossVec[1:(i-1),2]
  
  # recalibrate \hat{f}_t
  if(method_f_hat == 'AGNOSTIC'){
    if(i <= T0){
      pred_dt <- max(0,Xmat[i-1,2])
    }else{
      ## add intercept !!! 
      xx <- cbind(1, Xmat[1:(i-2),2])
      init_par <- as.numeric( solve(t(xx) %*% xx) %*% t(xx) %*% Xmat[2:(i-1),2] )
      pred_dt <- max(0,init_par[1] + init_par[2] * Xmat[i-1,2])
    }
  }else{
    if(i <= T0){
      pred_dt <- max(0,Xmat[i-1,2])
    }else{
      perform_update <- (i %in% (1 + T0 * ceiling(1.5^{0:13})))
      if(perform_update){
        setup_theta <- data.frame(
          'theta0' = runif(n = 25, min = theta_search_grid_controls$lower_bound[1], max = theta_search_grid_controls$upper_bound[1]),
          'theta1' = runif(n = 25, min = theta_search_grid_controls$lower_bound[2], max = theta_search_grid_controls$upper_bound[2])
        )
        setup_theta_i <- setup_theta
        setup_theta_i$loss <- NA
        for(counter in 1:nrow(setup_theta_i)){
          setup_theta_i$loss[counter] <- foo_to_optim_complicated(theta_vec = as.numeric(setup_theta_i[counter,1:2]))
        }
        init_par <- as.numeric(setup_theta_i[which.min(setup_theta_i$loss),1:2])
        # setup_theta <- data.frame(
        #   'theta0' = runif(n = 100, min = theta_search_grid_controls$lower_bound[1], max = theta_search_grid_controls$upper_bound[1]),
        #   'theta1' = runif(n = 100, min = theta_search_grid_controls$lower_bound[2], max = theta_search_grid_controls$upper_bound[2])
        # )
        # setup_theta_dist <- apply(setup_theta, 1, function(x) sum(((x - init_par)/c(100,1))^2))
        # setup_theta[order(setup_theta_dist)[1:10],]
      }
      # predict Dt
      pred_dt <- max(0, init_par[1] + init_par[2] * Xmat[i-1,2])
    }
  }
  
  # sample B new x = (tau, dthat) candidates
  newx <- cbind(
    runif(n = B, min = Xmin[1], max = Xmax[1]),
    pred_dt
  )
  
  # compute kernel quantities
  kTi <- array(data = NA, dim = c(i-1, B))
  if(i == 2){
    kTi[1,] <- apply(X = newx, MARGIN = 1, FUN = function(x){ my_kernel(x0 = x, x1 = Xmat[1,], l = kernel_bandwith, scales = Xmax-Xmin) })
  }else{
    for(j in 1:B){
      kTi[,j] <- apply(X = Xmat[1:(i-1),], MARGIN = 1, FUN = function(x){ 
        my_kernel(x0 = newx[j,], x1 = x, l = kernel_bandwith, scales = Xmax-Xmin) 
      })
    }
  }
  
  # compute estimates of mean and variance
  mui <- as.numeric( t(kTi) %*% normi_times_lossVec )
  sigma2i <- KERNEL_VALUE_EQUALS - sapply(1:B, function(j){ t(kTi[,j]) %*% normi %*% kTi[,j] })
  stopifnot(sigma2i > 0)
  
  # compute UCB criterion : remember that we minimize here =)
  ucbi <- mui - alpha * sqrt(log((i-1) / delta) * sigma2i)
  
  # select based on UCB criterion
  sel <- which.min( ucbi )
  
  # if(keep_computing_lambda){
  #   crit <- sapply(1:B, function(j){
  #     Kj <- Kmat[1:i,1:i]
  #     Kj[i,i] <- 1
  #     Kj[1:(i-1),i] <- Kj[i,1:(i-1)] <- kTi[,j]
  #     eigen(x = Kj, only.values = TRUE)$values[i]
  #   })
  #   weight <- 10
  #   if(max(crit) < 0) keep_computing_lambda <- FALSE
  # }else{
  #   weight <- 0
  # }
  # sel <- which.max(mui / min(abs(mui)) + weight * crit / min(crit))
  
  # update x records, observe noisy observation and update kernel matrix
  Xmat[i,1] <- newx[sel,1]
  Xmat[i,2] <- max(0 ,true_theta_par[1] + true_theta_par[2] * Xmat[i-1,2] + rnorm(n = 1, mean = 0, sd = d_sigma2))
  lossVec[i,1] <- get_loss(tau = Xmat[i,1], d = Xmat[i,2])
  lossVec[i,2] <- lossVec[i,1] + add_noise(s2 = sigma2)
  Kmat[i,i] <- 1
  if(i == 2){
    temp <- my_kernel(x0 = Xmat[2,], x1 = Xmat[1,], l = kernel_bandwith, scales = Xmax-Xmin)
  }else{
    temp <- apply(X = Xmat[1:(i-1),], MARGIN = 1, FUN = function(x){
      my_kernel(x0 = Xmat[i,], x1 = x, l = kernel_bandwith, scales = Xmax-Xmin)
    })
  }
  # temp <- array(data = NA, dim = c(i-1))
  # for(k in 1:(i-1)){
  #   temp[k] <- my_kernel(x0 = Xmat[i,], x1 = Xmat[k,], l = kernel_bandwith, scales = Xmax-Xmin)
  # }
  Kmat[1:(i-1),i] <- Kmat[i,1:(i-1)] <- temp
  
  # record regret
  losses[i,1] <- lossVec[i,1]
  opt <- optim(par = 0.5, fn = function(tau){
    get_loss(tau = tau, d = Xmat[i,2])    
  }, lower = 0, upper = 1, method = 'L-BFGS-B')
  losses[i,2] <- get_loss(tau = opt$par, d = Xmat[i,2]) 

  # record theta estimated
  theta_est[i,] <- init_par
  
  #
  timepoint[i] <- as.numeric(difftime(time2 = starttime, time1 = Sys.time(), units = 'mins'))
  cat('\n', i, '\t', round(timepoint[i]-timepoint[i-1],2), '\t', round(init_par,2))
  
}

###
###
###

# represent regret information
plot.ts(cumsum(c(0,losses[-1,1])))
lines(cumsum(c(0,losses[-1,2])), lwd = 3)
plot.ts((cumsum(c(0,losses[-1,1]) - c(0,losses[-1,2]))), lwd = 3)
lines(cumsum(c(0,losses[-1,1]) - c(0,losses[-1,2])), lwd = 3, col = 'red')

plot.ts((cumsum(c(0,losses[2:(i-1),1]) - c(0,losses[2:(i-1),2]))), lwd = 3)

losses_not_agnostic <- readRDS(file = 'history_loss_theta_NOT_AGNOSTIC_00.RDS')
lines(cumsum(c(0,losses_not_agnostic$losses[-1,1]) - c(0,losses_not_agnostic$losses[-1,2])), lwd = 3, col = 'red')

losses_agnostic <- readRDS(file = 'history_loss_theta_AGNOSTIC.RDS')
lines(cumsum(c(0,losses_agnostic$losses[-1,1]) - c(0,losses_agnostic$losses[-1,2])), lwd = 3, col = 'red')

saveRDS(object = list('losses' = losses, 'theta_est' = theta_est),
        file = 'history_loss_theta_NOT_AGNOSTIC_00.RDS')

saveRDS(object = list('losses' = losses, 'theta_est' = theta_est),
        file = 'history_loss_theta_NOT_AGNOSTIC.RDS')
## history_loss_agnostic <- losses
## history_theta_agnostic <- theta_est

plot.ts(theta_est)
###
###
###
setup_theta <- data.frame(
  'theta0' = runif(n = 250, min = theta_search_grid_controls$lower_bound[1], max = theta_search_grid_controls$upper_bound[1]),
  'theta1' = runif(n = 250, min = theta_search_grid_controls$lower_bound[2], max = theta_search_grid_controls$upper_bound[2])
)
setup_theta_i <- setup_theta
setup_theta_i$loss <- NA
for(counter in 1:nrow(setup_theta_i)){
  setup_theta_i$loss[counter] <- foo_to_optim_complicated(theta_vec = as.numeric(setup_theta_i[counter,1:2]))
}
hist(setup_theta_i$loss)
setup_theta_i[order(setup_theta_i$loss),]

# sample B new x = (tau, dthat) candidates
Bmy <- 200
B2 <- floor(sqrt(Bmy))
newx <- expand.grid(
  seq(from = Xmin[1], to = Xmax[1], length.out = B2),
  seq(from = Xmin[2], to = Xmax[2], length.out = B2)
)
B3 <- nrow(newx)

# compute kernel quantities
N_stop <- i-1
kTi <- array(data = NA, dim = c(N_stop, B3))
for(j in 1:B3){
    kTi[,j] <- apply(X = Xmat[1:(i-1),], MARGIN = 1, FUN = function(x){
      my_kernel(x0 = newx[j,], x1 = x, l = kernel_bandwith, scales = c(Xmax-Xmin))
      })
}
normi <- solve(Kmat[1:N_stop, 1:N_stop] + diag(N_stop) * sigma2)

# compute estimates of mean and variance
mui <- as.numeric( t(kTi) %*% normi %*% lossVec[1:N_stop,2])
sigma2i <- KERNEL_VALUE_EQUALS - sapply(1:B3, function(j){ t(kTi[,j]) %*% normi %*% kTi[,j] })
stopifnot(sigma2i > 0)
ucbi <- mui - (alpha) * sqrt(log((i-1) / delta) * sigma2i)
dfff <- data.frame('tau' = newx[,1], 'demand' = newx[,2], 'val' = ucbi, 'sumw' = colSums(kTi))
dfff$opt <- FALSE
levels_d <- unique(dfff$demand)
for(i in 1:length(levels_d)){
  index <- which.min(dfff$val[dfff$demand == levels_d[i]])
  dfff$opt[which(dfff$demand == levels_d[i])[index]] <- TRUE
}
# ggplot(data = dfff) +
#   geom_point(mapping = aes(x = tau, y = demand, col = val))
XmatDF <- data.frame(Xmat)
colnames(XmatDF) <- c('tau', 'demand')
g0 <- ggplot(data = dfff) +
  geom_tile(mapping = aes(x = demand, y = tau, fill = val)) + 
  scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue') +
  geom_point(data = XmatDF, mapping = aes(x = demand, y = tau)) +
  geom_point(data = dfff[dfff$opt,], mapping = aes(x = demand, y = tau), col = 'red', pch = 20, cex = 10) +
  geom_point(data = pars_grid[pars_grid$opt,], mapping = aes(x = d, y = tau), col = 'blue', pch = 20, cex = 10) 
 #+ geom_label(mapping = aes(x = demand, y = tau, label = round(sumw,0)))
g1 <- ggplot(data = pars_grid) +
  geom_tile(mapping = aes(x = d, y = tau, fill = val)) +
  scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue') +
  geom_point(data = dfff[dfff$opt,], mapping = aes(x = demand, y = tau), col = 'red', pch = 20, cex = 10) +
  geom_point(data = pars_grid[pars_grid$opt,], mapping = aes(x = d, y = tau), col = 'blue', pch = 20, cex = 10)
library(gridExtra)
grid.arrange(g0,g1)

### A bunch of plots ...
###
###

#
XmatDF <- data.frame(Xmat)
XmatDF$ID <- 1:nrow(XmatDF)
XmatDF$VAL <- yvec
str(XmatDF)

#
library(ggplot2)
ggplot(data = XmatDF[1:i,]) +
  geom_point(mapping = aes(x = X1, y = X2, col = ID)) +
  scale_color_continuous(low = 'white', high = 'blue')

#
plot(x = XmatDF$X1, y = XmatDF$X2, pch = 20, cex = 0.5)
lines(x = XmatDF$X1, y = XmatDF$X2, pch = 20, cex = 0.5)
XmatDF[i-1,]
pars_grid[which.max(pars_grid$val),]
ggplot(data = pars_grid) +
  geom_tile(mapping = aes(x = tau, y = d, fill = val))
ggplot(data = XmatDF) +
  geom_point(mapping = aes(x = X1, y = X2, col = VAL))
