foo_to_optim_complicated <- function(theta_vec){
  losses <- numeric(i-2)
  for(jjj in 1:(i-2)){
    
    # predict new demand
    pred_dt <- max(0, theta_vec[1] + theta_vec[2] * Xmat[jjj,2])
    
    # sample B new x = (tau, dthat) candidates
    BBBB <- 20
    newx <- cbind(runif(n = BBBB, min = Xmin[1], max = Xmax[1]), pred_dt)
    
    # compute kernel quantities
    kTi <- array(data = NA, dim = c(i-1, BBBB))
    for(j in 1:BBBB){
      kTi[,j] <- apply(X = Xmat[1:(i-1),], MARGIN = 1, FUN = function(x){ 
        my_kernel(x0 = newx[j,], x1 = x, l = kernel_bandwith, scales = Xmax-Xmin) 
      })
    }
    # for(k in 1:(i-1)){
    #   for(j in 1:BBBB){
    #     kTi[k,j] <- my_kernel(x0 = newx[j,], x1 = Xmat[k,], l = kernel_bandwith, scales = Xmax-Xmin)
    #   }
    # }
    
    # compute estimates of mean and variance
    mui <- as.numeric( t(kTi) %*% normi_times_lossVec )
    sigma2i <- KERNEL_VALUE_EQUALS - sapply(1:BBBB, function(j){ t(kTi[,j]) %*% normi %*% kTi[,j] })
    
    # choose \tau for the predicted value
    opt_taui <- which.min(mui - alpha * sqrt(log((i-1) / delta) * sigma2i))
    
    # compute new kernel quantities
    kTi <- apply(X = Xmat[1:(i-1),], MARGIN = 1, FUN = function(x){ 
      my_kernel(x0 = c(newx[opt_taui,1], Xmat[jjj+1,2]), x1 = x, l = kernel_bandwith, scales = Xmax-Xmin) 
    })
    # kTi <- array(data = NA, dim = c(i-1, 1))
    # for(k in 1:(i-1)){
    #   kTi[k,1] <- my_kernel(x0 = c(newx[opt_taui,1], Xmat[jjj+1,2]), 
    #                         x1 = Xmat[k,], l = kernel_bandwith, scales = Xmax-Xmin)
    # }
    
    # predict loss of i
    lossi <- as.numeric( t(kTi) %*% normi_times_lossVec )
    
    losses[jjj] <- lossi
    
  }
  mean(losses) + 1e05 * (any(theta_search_grid_controls$lower_bound > theta_vec) +
                           any(theta_search_grid_controls$upper_bound < theta_vec))
}

#######
#######
#######
#######
#######

setup_theta <- data.frame(
  'theta0' = runif(n = 25, min = theta_search_grid_controls$lower_bound[1], max = theta_search_grid_controls$upper_bound[1]),
  'theta1' = runif(n = 25, min = theta_search_grid_controls$lower_bound[2], max = theta_search_grid_controls$upper_bound[2])
)

foo_to_optim <- function(theta_vec){
  losses <- numeric(i-2)
  for(jjj in 1:(i-2)){
    
    # predict new demand
    pred_dt <- max(0, theta_vec[1] + theta_vec[2] * Xmat[jjj,2])
    
    # sample B new x = (tau, dthat) candidates
    BBBB <- 20
    newx <- cbind(seq(from = Xmin[1], to = Xmax[1], length.out = BBBB), pred_dt)
    
    # compute kernel quantities
    kTi <- array(data = NA, dim = c(i-1, BBBB))
    for(k in 1:(i-1)){
      for(j in 1:BBBB){
        kTi[k,j] <- my_kernel(x0 = newx[j,], x1 = Xmat[k,], l = kernel_bandwith, scales = Xmax-Xmin)
      }
    }
    
    # compute estimates of mean and variance
    mui <- as.numeric( t(kTi) %*% normi %*% lossVec[1:(i-1),2] )
    sigma2i <- KERNEL_VALUE_EQUALS - diag( t(kTi) %*% normi %*% kTi )
    
    # choose \tau for the predicted value
    opt_taui <- which.min(mui - alpha * sqrt(log((i-1) / delta) * sigma2i))
    
    # compute new kernel quantities
    kTi <- array(data = NA, dim = c(i-1, 1))
    for(k in 1:(i-1)){
      kTi[k,1] <- my_kernel(x0 = c(newx[opt_taui,1], Xmat[jjj+1,2]), 
                            x1 = Xmat[k,], l = kernel_bandwith, scales = Xmax-Xmin)
    }
    
    # predict loss of i
    lossi <- as.numeric( t(kTi) %*% normi %*% lossVec[1:(i-1),2] )
    
    losses[jjj] <- lossi
    
  }
  mean(losses) + 1e05 * (any(theta_search_grid_controls$lower_bound > theta_vec) +
                           any(theta_search_grid_controls$upper_bound < theta_vec))
}

setup_theta_i <- setup_theta
setup_theta_i$loss <- NA
for(counter in 1:nrow(setup_theta_i)){
  setup_theta_i$loss[counter] <- foo_to_optim(theta_vec = as.numeric(setup_theta_i[counter,1:2]))
}
# 
# hist(setup_theta_i$loss)
# ggplot(data = setup_theta_i) +
#   geom_tile(mapping = aes(x = theta0, y = theta1, fill = loss)) +
#   geom_label(mapping = aes(x = theta0, y = theta1, label = loss))

setup_theta_i[which.min(setup_theta_i$loss),]

# optim(par = c(83, 0), fn = foo_to_optim, method = 'L-BFGS-B', 
#       lower = theta_search_grid_controls$lower_bound, upper = theta_search_grid_controls$upper_bound)
# optim(par = c(83, 0), fn = foo_to_optim, method = 'L-BFGS-B', control = list(parscale = c(200,1), ndeps = c(1e-1, 1e-1)), 
#       lower = c(0, 0), upper = c(200,1))

