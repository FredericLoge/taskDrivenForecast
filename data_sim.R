sim_ar_p_2 <- function(n, sigma2e, theta){
  
  # storing simulation in vec s
  s <- numeric(n)
  
  # get p from vec beta
  p <- length(theta)
  
  # simulate under AR(p)
  for(i in 2:n){
    s[i] <- sign(s[i-1]) * sqrt(abs(s[i-1])) * theta
    s[i] <- s[i] + rnorm(n = 1, mean = 0, sd = sqrt(sigma2e))
  }
  
  # return sample
  return(s)
  
}

sim_ar_1 <- function(n, sigma2e, theta){
  
  # storing simulation in vec s
  s <- numeric(n)
  
  # get p from vec beta
  p <- length(theta)
  
  # simulate under AR(p)
  for(i in 2:n){
    s[i] <- s[i-1] * theta[1] + rnorm(n = 1, mean = 0, sd = sqrt(sigma2e))
  }

  # normalize  
  s_range <- sqrt( 2 * log(2 * n / 0.01) / (1 - theta^2))
  s <- with(problem_constants, P_min + (s + s_range) / (2 * s_range) * (P_max - P_min))
  
  # return sample
  return(s)
  
}

sim_ar_1_features <- function(y){
  mat <- cbind("intercept" = 1,  "y_lag_1" = c(rep(NA, 1), y[1:(length(y)-1)]))
  return(mat)
}

daily_constant_scenario <- function(n){
  s <- rep(x = c(0, 10, -10, 5, 0, 5, 4), ceiling(n/7))
  s <- s[1:n]
  s <- with(problem_constants, (s - min(s)) / (max(s) - min(s)) * (P_max - P_min) + P_min)
  return(s)
}

daily_constant_scenario_features <- function(n){
  mat <- cbind(1, FactoMineR::tab.disjonctif(tab = rep(1:7, ceiling(n/7))[1:n])[,-1])
  colnames(mat) <- c('intercept', paste0('day', 2:7))
  return(mat)
}

rasymnorm <- function(n = 1, mean = 0, sd_left = 1, sd_right = sd_left){
  is_left <- runif(n = n, min = -1/2, max = +1/2) > 0
  o <- numeric(n)
  o[is_left == T] <- - abs(rnorm(n = sum(is_left), mean = 0, sd = sd_left))
  o[is_left == F] <- + abs(rnorm(n = n - sum(is_left), mean = 0, sd = sd_right))
  o <- o + mean
  return(o)
}
# sam <- rasymnorm(n = 100000, mean = 10, sd_left = 5, sd_right = 20)
# hist(sam)

rasymexpnorm <- function(n = 1, mean = 0, sd_left = 1, sd_right = sd_left){
  is_left <- runif(n = n, min = -1/2, max = +1/2) > 0
  o <- numeric(n)
  o[is_left == T] <- - 2^(abs(rnorm(n = sum(is_left), mean = 0, sd = sd_left)))
  o[is_left == F] <- + 2^(abs(rnorm(n = n - sum(is_left), mean = 0, sd = sd_right)))
  o <- o + mean
  return(o)
}
# sam <- rasymexpnorm(n = 100000, mean = 10, sd_left = 1, sd_right = 1)
# hist(sam)

rasymweibull <- function(n = 1, mean = 0, scale_left = 1, scale_right = scale_left){
  is_left <- runif(n = n, min = -1/2, max = +1/2) > 0
  v <- NULL
  LAMBDA <- 1
  K <- 1
  while(length(v) < n){
    newv <- rweibull(n = n, shape = K, scale = LAMBDA)
    mode_weibull <- LAMBDA * (1 - 1/K)^(1/K)
    ok_cond <- (newv > mode_weibull)
    newv <- newv[ok_cond]
    v <- c(v, newv - mode_weibull)
  }
  v <- v[1:n]
  o <- numeric(n)
  o[is_left == T] <- - scale_left * v[is_left == T]
  o[is_left == F] <- + scale_right * v[is_left == F]
  o <- o + mean
  return(o)
}
# sam <- rasymweibull(n = 100000, mean = 0, scale_left = 5, scale_right = 20)
# hist(sam, breaks = 30)

###
sim_noise <- function(n, noise_type){
  switch(noise_type,
         "1" = rnorm(n = n, mean = 0, sd = 5),
         "2" = rasymnorm(n = n, mean = 0, sd_left = 5, sd_right = 20),
         "3" = rasymexpnorm(n = n, mean = 0, sd_left = 1, sd_right = 1)
  )
}

# sam <- sim_scenario_3(n = 1000, noise_type = 3, snr = 1)
# plot.ts(sam)

###
sim_signal <- function(n, signal_type){
  switch(signal_type,
         "1" = daily_constant_scenario(n),
         "2" = sim_ar_1(n = n, sigma2e = 5, theta = 0.8)
  )
}

###
get_features <- function(y, signal_type, noise_type){
  # noise_type will be ignored for now
  switch(signal_type,
         "1" = daily_constant_scenario_features(n),
         "2" = sim_ar_1_features(y = y)
  )
  
}
