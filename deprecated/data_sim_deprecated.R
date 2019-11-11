

#' @title Simulate demand series, assuming AR(p)
#' @param n number of observations
#' @param sigma2e error variance
#' @param theta ar coefficients
sim_ar_p <- function(n, sigma2e, theta){
  
  # storing simulation in vec s
  s <- numeric(n)
  
  # get p from vec beta
  p <- length(theta)
  
  # simulate under AR(p)
  for(i in 2:n){
    if(i <= p){
      s[i] <- sum(s[(i-1):1] * theta[1:(i-1)])
    }else{
      s[i] <- sum(s[(i-1):(i-p)] * theta)
    }
    s[i] <- s[i] + rnorm(n = 1, mean = 0, sd = sqrt(sigma2e))
  }
  
  # return sample
  return(s)
  
}

# generate "forecast" matrix for demand signal
gen_forecast_matrix <- function(y, p){
  
  # generate matrix mat with rows : (Y_{t-1}, ..., Y_{t-p}) 
  mat <- sapply(X = seq(from = 1, to = p, by = 1), FUN = function(i){
    c(rep(NA, i), y[1:(n-i)])
  })
  
  # add on left side of mat : (Y_{t}, 1), actual signal and intercept
  mat <- cbind(y, 1, mat)
  
  # add nice column names
  colnames(mat) <- c('y', paste0('x', 0:p))
  
  # return matrix
  return(mat)
}

