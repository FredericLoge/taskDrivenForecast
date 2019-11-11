optim_L0 <- function(y, x){
  
  x <- as.matrix(x)
  
  # theta estimate
  theta_hat <- solve(t(x) %*% x) %*% t(x) %*% y
  
  # error term variance estimate
  sigma2e_hat <- mean( (y - x %*% theta_hat)^2 )

  # variance estimate of theta hat
  theta_var_hat <- sigma2e_hat * solve(t(x) %*% x)
  
  # convert theta_hat to a nicer vector w/ names
  theta_hat_ <- as.numeric(theta_hat)
  names(theta_hat_) <- rownames(theta_hat)
  
  # return estimates
  list(
    'theta_hat' = theta_hat_,
    'sigma2e_hat' = sigma2e_hat,
    'theta_var_hat' = theta_var_hat
  ) 
  
}


#' Optimize over L1 function by grad evaluation
optim_grad_L1 <- function(y, x, lb, up, init_pars, setup_pars){
  
  x <- as.matrix(x)
  
  # function to optimize
  of <- function(pars, y, x, setup_pars){
    tauSel <- optimal_tau(pars = setup_pars, D = x %*% pars)
    costs <- cost_foo(pars = setup_pars, D = y, tau = tauSel)
    return(mean(costs))
  }
  
  # optimization
  opt <- optim(par = init_pars, fn = of, method = 'L-BFGS-B', lower = lb, upper = up, y = y, x = x, setup_pars = setup_pars)
  return(opt)
  
}


#' Optimize over L1 function by deterministic grid search
optim_L1 <- function(y, x, search_grid_controls){
  
  # extract p, reused several times
  p <- search_grid_controls$p
  
  # initialize search grid with first element search controls
  search_grid <- seq(from = search_grid_controls$lower_bound[1], 
                     to = search_grid_controls$upper_bound[1], 
                     length.out = search_grid_controls$nb_points[1])
  
  # extend search grid for each other element
  if(p > 1){
    for(i in 2:p){
      new_grid <- seq(from = search_grid_controls$lower_bound[i], 
                      to = search_grid_controls$upper_bound[i], 
                      length.out = search_grid_controls$nb_points[i])
      search_grid <- expand.grid(search_grid, new_grid)      
    }
  }else{ # or adapt to a correct format
    search_grid <- data.frame(search_grid)
  }
  
  # add nice column names
  colnames(search_grid) <- paste0('dim', 1:p)

  # for all configurations proposed, compute loss on (y,x) sample, solving LP problem
  search_grid$loss <- NA
  search_grid$L0loss <- NA
  for(i in 1:nrow(search_grid)){
    y_hat <- x %*% as.numeric(search_grid[i,1:p])
    search_grid$loss[i] <- eval_theta(D = y, Dhat = y_hat)
    search_grid$L0loss[i] <- mean((y-y_hat)^2)
  }
  
  # identify set of optimum thetas
  minimizers_row_index <- which(x = search_grid$loss == min(search_grid$loss))
  minimizers <- search_grid[minimizers_row_index, 1:p]
  
  # return search grid, optimal thetas
  l <- list('search_grid' = search_grid, 'theta_hat' = minimizers)
  return(l)
  
}

