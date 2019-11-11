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

###
###
###

#
find_maxima <- FALSE

# define bounds of (\tau, demand)
Xmin <- with(problem_constants, c(tau_min, P_min))
Xmax <- with(problem_constants, c(tau_max, P_max))

# number of time steps
N <- 100

# dimension of x (\tau, demand)
p <- 2

# prepare matrices
Kmat <- array(data = NA, dim = c(N,N))
yvec <- array(data = NA, dim = c(N,1))
Xmat <- array(data = NA, dim = c(N, p))

# choose observation error level
sigma2 <- 0.02^2

# kernel bandwidth
kernel_bandwith <- 1

# number of candidates sampled
B <- 10

# parameters of UCB
alpha <- 35
delta <- 0.0005

# random initial point
x0 <- runif(n = 2, min = Xmin, max = Xmax)
Xmat[1,] <- x0
yvec[1,1] <- get_noisy_loss(tau = x0[1], d = x0[2], s2 = sigma2) 
Kmat[1,1] <- 1

## keep_computing_lambda <- TRUE

# status update
pb <- txtProgressBar(min = 1, max = N, initial = 1)

#
for(i in 2:N){
  
  # status update
  setTxtProgressBar(pb, value = i)
  
  # (X_T + \sigma^2 I)^{-1}
  normi <- solve(Kmat[1:(i-1), 1:(i-1)] + diag(i-1) * sigma2)
  
  # sample B new x candidates
  newx <- cbind(
    runif(n = B, min = Xmin[1], max = Xmax[1]),
    runif(n = B, min = Xmin[2], max = Xmax[2])
  )
  
  # compute kernel quantities
  kTi <- array(data = NA, dim = c(i-1, B))
  for(k in 1:(i-1)){
    for(j in 1:B){
      kTi[k,j] <- my_kernel(x0 = newx[j,], x1 = Xmat[k,], l = kernel_bandwith, scales = Xmax-Xmin)
    }
  }
  
  # compute estimates of mean and variance
  mui <- as.numeric( t(kTi) %*% normi %*% yvec[1:(i-1),1] )
  sigma2i <- KERNEL_VALUE_EQUALS - diag( t(kTi) %*% normi %*% kTi )
  
  if(find_maxima){
    # compute UCB criterion
    ucbi <- mui + alpha * sqrt(log((i-1) / delta) * sigma2i)
    # select based on UCB criterion
    sel <- which.max( ucbi )
  }else{
    # compute LCB criterion
    ucbi <- (- mui) + alpha * sqrt(log((i-1) / delta) * sigma2i)
    # select based on UCB criterion
    sel <- which.min( ucbi )
  }
  
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
  Xmat[i,] <- newx[sel,]
  yvec[i,1] <- get_noisy_loss(tau = newx[sel,1], d = newx[sel,2], s2 = sigma2) 
  Kmat[i,i] <- 1
  Kmat[1:(i-1),i] <- Kmat[i,1:(i-1)] <- kTi[,sel]
  
}

# represent regret information
plot.ts(cumsum(rep(1,N)), col = 'black', lwd = 4)
if(find_maxima == TRUE){
  lines(cumsum(1 - yvec / max(pars_grid$val)), col = 'blue')
}else{
  lines(cumsum(1 - yvec / min(pars_grid$val)), col = 'blue')
}
legend('topleft', col = 'blue', legend = 'Regret incurred', lty = 1)

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

