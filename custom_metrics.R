#' @title Cost function of optim problem
#' @param pars unit model parameters
#' @param D real demand 
#' @param tau control chosen
cost_foo <- function(pars, D, tau){
  I <- pars$a*tau + pars$b
  P <- pars$c*I
  cost <- pars$d*I + pars$e*pmax(0, D-P) - pars$f*pmin(D,P)
  return(cost)
}

#' @title Find optimal tau based on a reliable demand information
#' @param pars unit model parameters
#' @param D real demand 
optimal_tau <- function(pars, D){
  # design simple search grid
  tau_grid <- seq(from = pars$tau_min, to = pars$tau_max, length.out = 30)
  # compute cost function
  ot <- cost_foo(pars = pars, D = D, tau = tau_grid)
  # return tau minimizing cost function
  return(tau_grid[which.min(ot)])
}

#' @title Evaluate the loss function over a selected theta
#' @param D
#' @param Dhat
eval_theta <- function(D, Dhat){
  tauSel <- sapply(X = Dhat, FUN = function(Dhat){
    otau <- optimal_tau(pars = pars, D = Dhat)
  })
  costs <- cost_foo(pars = pars, D = D, tau = tauSel)
  return(mean(costs))
}

# return optimal tau for demand vector D
optimal_tau <- function(pars, D){
  
  # we try to adjust production to demand as best we can :
  ##tau_0 <- (D - pars$c*pars$b) / (pars$c*pars$a)
  tau_0 <- with(pars, (D - c*b) / (c*a))
  
  # if the optimal tau is out of bounds, put it on the bounds
  tau <- pmin(tau_0, pars$tau_max)
  tau <- pmax(tau, pars$tau_min)
  
  # return optimal tau control
  return(tau)
  
}

# eval theta
eval_theta <- function(D, Dhat){
  tauSel <- optimal_tau(pars = pars, D = Dhat)
  costs <- cost_foo(pars = pars, D = D, tau = tauSel)
  return(mean(costs))
}