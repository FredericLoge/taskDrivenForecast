muu <- runif(n = 5, min = 0, max = 10)
muu <- sort(muu, decreasing = TRUE)
a = 1
astar = 1:5
v <- sapply(1:5, function(a){
  sum(muu[astar < a] - muu[a]) + 1/2 * sum(muu[astar > a] - muu[a])
})
sum(v)
sum(muu * ((5+1)/2 - astar))

x <- runif(n = 100000, min = -100, max = 100)
s2 <- 1
mu_vector <- seq(from = 0, to = 10, length.out = 100)
res <- sapply(X = mu_vector, FUN = function(mu){
  v1 <- mean(exp(- x^2 / (2*s2)) * exp( - (x-mu)^2 / 2*s2))
  v2 <- sqrt(s2 * pi) * sqrt(exp( - (mu)^2 / 2*s2))
  return(c(v1, v2))
})
res <- t(res)
mean(res[,1] < (res[,2]))
 plot.ts(sqrt(res[,2]))
lines(res[,1], col = 'red')

index <- 70:100
plot.ts(sqrt(res[index,2]))
lines(res[index,1], col = 'red')

# # model parameters
# fx <- function(x, i, noisy = FALSE, s2){
#   # set up Sound-to-Noise (variance) ratio
#   SNR <- 5
#   #
#   fx <- sin(x*(2*pi))
#   #
#   ux <- sin(x*(2*pi) - i*pi/12) / sqrt(SNR)
#   #
#   if(noisy){
#     if(missing(s2)) stop('Variance par is missing ...')
#     noise <- rnorm(n = length(x), mean = 0, sd = s2)
#     return( fx + ux + noise )
#   }else{
#     return( fx + ux )
#   }
# }
# 
# # model parameters
# fx <- function(x, i, noisy = FALSE, s2){
#   #
#   fx <- sin(x*(2*pi))
#   #
#   ux <- rep(0, length(x))
#   ux[i == 1] <- 1*exp(- 1/2 * (0.15)^(-2) * (x[i == 1] - 0.15)^2)
#   ux[i == 2] <- -1*(x[i == 2] < 0.5 & x[i == 2] > 0.35)
#   ux[i == 3] <- 1*(x[i == 3] < 0.8 & x[i == 3] > 0.2)
#   #
#   if(noisy){
#     if(missing(s2)) stop('Variance par is missing ...')
#     noise <- rnorm(n = length(x), mean = 0, sd = s2)
#     return( fx + ux + noise )
#   }else{
#     return( fx + ux )
#   }
# }
# 
# norm_fx <- sapply(1:3, function(i){
#   sum(abs(fx(x = x, i = i, noisy = FALSE))) / length(x)
# })
# 
# x <- seq(from = 0, to = 1, length.out = 100)
# YLIM <- c(-1, 1)
# par(mfrow = c(m,1))
# for(i in 1:m){
#   plot(x = x, y = get_obs(x = x, i = i, noisy = FALSE), type = 'l', col = 'black', lwd = 2, ylim = YLIM)
#   lines(x = x, y = get_obs(x = x, i = 0, noisy = FALSE), type = 'l', col = 'red', lwd = 2, lty = 2)
# }
# par(mfrow = c(1,1))

# bounds of variable X
x_bounds <- c(0,1)

# number of individuals
m <- 9

# generate smooth function by generating some random points and fitting a spline
gen_smooth_function <- function(seed){
  set.seed(seed)
  K <- 20
  x <- seq(0, 1, length.out = K)
  ## sort(runif(n = K, min = 0, max = 1))
  y <- runif(n = K, min = - sqrt(12)/2, max = sqrt(12)/2)
  ll <- loess(formula = y~x)
  return(ll)
}

# x-base coordinates
x <- seq(0, 1, length.out = 100)

# generate base function f and individual functions u_i
fx <- gen_smooth_function(seed = 18543)
ux <- list()
ux_mean <- rep(0, m)
set.seed(3743) ; seeds <- floor(runif(n = m, min = 1753, max = 298753))
for(i in 1:m){
  ux[[i]] <- gen_smooth_function(seed = seeds[i])
  ux_mean[i] <- mean(predict(object = ux[[i]], newdata = x))
}

# sound to "noise" ratio, handling the balance Variance(f) / Variance(u_i)
SNR <- 10

# observation function
get_obs <- function(x, i, noisy = FALSE, s2, only.fx = FALSE, only.ux = FALSE){
  #
  fx_obs <- predict(object = fx, newdata = x)
  if(i == 0 | only.fx) return(fx_obs)
  #
  ux_obs <- predict(object = ux[[i]], newdata = x) - ux_mean[i]
  ux_obs <- ux_obs / sqrt(SNR)
  if( only.ux ) return(ux_obs)
  #
  if(noisy){
    if(missing(s2)) stop('Variance par is missing ...')
    noise <- rnorm(n = length(x), mean = 0, sd = sqrt(s2))
    return( fx_obs + ux_obs + noise )
  }else{
    return( fx_obs + ux_obs )
  }
} 

#
library(ggplot2)

# plot generated functions
x <- seq(from = 0, to = 1, length.out = 100)
temp <- lapply(X = 0:m, FUN = function(i){
  y <- get_obs(x = x, i = i, noisy = FALSE)
  cbind(x, y, i)
})
temp <- do.call(rbind.data.frame, temp)
colnames(temp) <- c('x', 'y', 'index')
temp$index_label <- paste0('Ind #', temp$index)
temp$index_label[temp$index == 0] <- 'Base function'
temp$index_label <- factor(temp$index_label)
ggplot(data = NULL) +
  geom_line(data = temp[temp$index != 0,], mapping = aes(x = x , y = y, group = index_label, col = index_label), lty = 2) + 
  geom_line(data = temp[temp$index == 0,], mapping = aes(x = x , y = y), col = 'black', lty = 1)  
  
#
Ni <- sort(floor(exp(runif(n = m, min = 2, max = 5))))
N <- sum(Ni)
df <- data.frame(
  'index' = 1,
  'x' = runif(n = N, min = 0, max = 1),
  'y' = NA
)
count <- 0
for(i in 1:m){
  if(i >= 2) count <- count + Ni[i-1]
  indexi <- count + 1:Ni[i]
  df$index[indexi] <- i
  df$y[indexi] <- get_obs(x = df$x[indexi], i = i, noisy = TRUE, s2 = 0.1)
}
table(df$index) # Ni

#
ggplot(data = df) +
  geom_point(mapping = aes(x = x, y = y, col = factor(index)))

#
my_kernel <- function(x0, x1, l){
  exp( -(2*l^2)^{-1} * (x0 - x1)^2)
}

#
compute_kernel_gp <- function(y_vector, x_vector, lambda, bandwidth){
  
  #
  x_vector_length <- length(x_vector)
  
  # compute kernel matrix "Kmat" over "x_vector"
  Kmat <- array(data = NA, dim = c(x_vector_length, x_vector_length))
  Kmat[1,1] <- 1
  for(i in 2:x_vector_length){
    Kmat[i,i] <- 1
    Kmat[i,-i] <- Kmat[-i,i] <- my_kernel(x0 = x_vector[-i], x1 = x_vector[i], l = bandwidth)
  }
  
  # compute regressor-like coefficients
  Kreg <- solve( Kmat + lambda * diag(x_vector_length) ) %*% y_vector

  #
  l <- list("x_vector" = x_vector, "y_vector" = y_vector, "kernel_matrix" = Kmat, "kernel_reg" = Kreg)
  return(l) 
  
}  

#
compute_kernel_ind_gp <- function(y_vector, x_vector, index, weights, lambda, bandwidth){
  
  #
  x_vector_length <- length(x_vector)
  
  # compute kernel matrix "Kmat" over "x_vector"
  Kmat <- array(data = NA, dim = c(x_vector_length, x_vector_length))
  Kmat[1,1] <- weights[index[1]]
  for(i in 2:x_vector_length){
    Kmat[i,i] <- weights[index[i]]
    Kmat[i,-i] <- Kmat[-i,i] <- weights[index[-i]] * my_kernel(x0 = x_vector[-i], x1 = x_vector[i], l = bandwidth)
  }
  
  # compute regressor-like coefficients
  Kreg <- solve( Kmat + lambda * diag(x_vector_length) ) %*% y_vector
  
  #
  l <- list("x_vector" = x_vector, "y_vector" = y_vector, "kernel_matrix" = Kmat, "kernel_reg" = Kreg)
  return(l) 
  
}  

#
predict_from_gp <- function(new_x_vector, gp, bandwidth){
  
  # compute distance between "new_x_vector" and training "x_vector"
  new_x_dist <- sapply(X = new_x_vector, FUN = function(x){
    my_kernel(x0 = x, x1 = gp$x_vector, l = bandwidth) 
  })
  
  # compute predictions
  new_x_pred <- t(new_x_dist) %*% gp$kernel_reg
  
  #
  l <- data.frame('x' = new_x_vector, 'y' = new_x_pred)
  l <- l[order(l$x),]
  return(l)
  
}

#
predict_from_ind_gp <- function(new_x_vector, gp, new_index, weights, bandwidth){
  
  # compute distance between "new_x_vector" and training "x_vector"
  new_x_dist <- sapply(X = 1:length(new_x_vector), FUN = function(j){
    my_kernel(x0 = new_x_vector[j], x1 = gp$x_vector, l = bandwidth) * weights[new_index[j]] 
  })
  
  # compute predictions
  new_x_pred <- t(new_x_dist) %*% gp$kernel_reg
  
  #
  l <- data.frame('x' = new_x_vector, 'y' = new_x_pred)
  l <- l[order(l$x),]
  return(l)
  
}

#
predict_from_ind_gp_2 <- function(x, x_vector, y_vector, index, weights, bandwidth){
  
  # compute set of weights
  weights <- my_kernel(x0 = x_vector, x1 = x, l = bandwidth) * weights[index]
  
  # compute prediction
  pred <- sum(weights * y_vector) / sum(weights)
  
  # return prediction
  return(pred)
  
}

# #
# LAMBDA <- 1
# BANDWIDTH_GLOBAL <- 0.1
# BANDWIDTH_IND <- 0.3
# 
# #
# gp_global <- compute_kernel_gp(y_vector = df$y, x_vector = df$x, lambda = LAMBDA, bandwidth = BANDWIDTH_GLOBAL)
# gp_global_pred <- predict_from_gp(new_x_vector = gp_global$x_vector, gp = gp_global, bandwidth = BANDWIDTH_GLOBAL)
# df_est <- gp_global_pred
# df_est$index <- 0
# df_est$categ <- 'global'
# 
# #
# gp_ind <- list()
# for(i in 1:m){
#   index_i <- (df$index == i)
#   l <- list()
#   l$gp <- compute_kernel_gp(y_vector = df$y[index_i], x_vector = df$x[index_i], lambda = LAMBDA, bandwidth = BANDWIDTH_IND)
#   l$pred <- predict_from_gp(new_x_vector = l$gp$x_vector, gp = l$gp, bandwidth = BANDWIDTH_IND)
#   gp_ind[[i]] <- l
#   l$pred$index <- i
#   l$pred$categ <- 'individual'
#   df_est <- rbind.data.frame(df_est, l$pred)
# }
# 
# #
# gp_ind_2 <- list()
# for(i in 1:m){
#   l <- list()
#   w <- rep(0, m) # /(m+1)
#   w[i] <- 1
#   l$pred <- sapply(df$x, function(x){
#     predict_from_ind_gp_2(x, x_vector = df$x, y_vector = df$y, index = df$index, weights = w, bandwidth = BANDWIDTH_IND)
#   })
#   l$pred <- data.frame('x' = df$x, 'y' = l$pred)
#   # l$gp <- compute_kernel_ind_gp(y_vector = df$y, x_vector = df$x, index = df$index, weights = w, lambda = LAMBDA, bandwidth = BANDWIDTH_IND)
#   # l$pred <- predict_from_ind_gp(new_x_vector = l$gp$x_vector, gp = l$gp, new_index = df$index, weights = w, bandwidth = BANDWIDTH_IND)
#   gp_ind_2[[i]] <- l
#   l$pred$index <- i
#   l$pred$categ <- 'mixed'
#   df_est <- rbind.data.frame(df_est, l$pred)
# }

######
BANDWIDTH_GLOBAL <- 0.1
w <- rep(1, m)
pred <- sapply(df$x, function(x){
  predict_from_ind_gp_2(x, x_vector = df$x, y_vector = df$y, index = df$index, weights = w, bandwidth = BANDWIDTH_GLOBAL)
})
pred <- data.frame('x' = df$x, 'y' = pred, 'index' = 0, 'category' = 'global')
df_est <- list()
df_est[[1]] <- pred

######
BANDWIDTH_IND <- BANDWIDTH_GLOBAL
for(i in 1:m){
  ######
  w <- rep(0, m)
  w[i] <- 1
  pred <- sapply(df$x, function(x){
    predict_from_ind_gp_2(x, x_vector = df$x, y_vector = df$y, index = df$index, weights = w, bandwidth = BANDWIDTH_IND)
  })
  pred <- data.frame('x' = df$x, 'y' = pred, 'index' = i, 'category' = 'individual')
  df_est[[length(df_est)+1]] <- pred
  ######
  w <- rep(1, m)
  w[i] <- m
  w <- w / sum(w)
  pred <- sapply(df$x, function(x){
    predict_from_ind_gp_2(x, x_vector = df$x, y_vector = df$y, index = df$index, weights = w, bandwidth = BANDWIDTH_IND)
  })
  pred <- data.frame('x' = df$x, 'y' = pred, 'index' = i, 'category' = 'mixed')
  df_est[[length(df_est)+1]] <- pred
}

#
df_est <- do.call(rbind.data.frame, df_est)
df_est$index_label <- paste0('Ind #', df_est$index)
df_est$index_label[df_est$index == 0] <- 'Base function'
df_est$index_label <- factor(df_est$index_label)
temp$category <- 'Truth'
dfff <- rbind.data.frame(df_est, temp)

#
ggplot(data = dfff[dfff$index == 0 ,]) +
  geom_line(mapping = aes(x = x, y = y, col = category)) +
  geom_rug(data = df, mapping = aes(x = x)) +
  geom_point(data = df, mapping = aes(x = x, y = y))

#
ggplot(data = dfff[dfff$index != 0 ,]) +
  geom_line(mapping = aes(x = x, y = y, col = category)) +
  geom_rug(data = df, mapping = aes(x = x)) +
  facet_wrap(~index_label)

# check 
compare <- sapply(X = 1:m, FUN = function(i){
  ref <- get_obs(x = df$x, i = i, noisy = FALSE)
  pred1 <- dfff$y[dfff$index == i & dfff$category == 'individual']
  pred2 <- dfff$y[dfff$index == i & dfff$category == 'mixed']
  pred3 <- dfff$y[dfff$index == 0 & dfff$category == 'global']
  c( mean(abs(ref-pred1)), mean(abs(ref-pred2)), mean(abs(ref-pred3)) )
})
compare <- data.frame(t(compare))
colnames(compare) <- c('individual', 'mixed', 'global')
compare
plot.ts(compare$individual, type = 'o', pch = 20, cex = 4, ylim = c(0, max(compare)))
lines(compare$mixed, type = 'o', pch = 20, cex = 4, ylim = c(0, max(compare)), col = 'red')
lines(compare$global, type = 'o', pch = 20, cex = 4, ylim = c(0, max(compare)), col = 'blue')
legend('bottomleft', legend = colnames(compare), lty = 1, pch = 20, cex = 2, col = c('black', 'red', 'blue'))

# plot(dfff[dfff$index == 9 & dfff$category == 'individual', c('x', 'y')])
# points(dfff[dfff$index == 9 & dfff$category == 'Truth', c('x', 'y')], col  ='red')
#
# plot(gp_global_pred, type = 'l', pch = 20, ylim = range(df$y))
# ind_color <- c('red', 'blue', 'green')
# for(i in 1:3){
#   lines(gp_ind[[i]]$pred, col = ind_color[i], type = 'l', pch = 20)
#   points(x = gp_ind[[i]]$gp$x_vector, y = gp_ind[[i]]$gp$y_vector, col = ind_color[i], pch = 20)
# }
# 
# #
# par(mfrow = c(3,1))
# for(i in 1:3){
#   pred_global_i <- predict_from_gp(new_x_vector = gp_ind[[i]]$gp$x_vector, gp = gp_global, bandwidth = BANDWIDTH)
#   plot(x = gp_ind[[i]]$pred$x, y = gp_ind[[i]]$pred$y - pred_global_i$y, col = ind_color[i], type = 'l', pch = 20, ylim = c(-0.5, 1))
#   lines(x = x, y = get_obs(x = x, i = i, noisy = FALSE, only.ux = TRUE),
#         # fx(x = x, i = i, noisy = FALSE) - fx(x = x, i = 0, noisy = FALSE),
#         type = 'l', col = 'red', lwd = 2, lty = 2)
# }
# 
# #
# plot(x = x, y = fx(x = x, i = i, noisy = FALSE) - fx(x = x, i = 0, noisy = FALSE),
#      type = 'l', col = 'red', lwd = 2, lty = 2)
# 
# 
# my_kernel <- function(x0, x1, type = 'exp', l){
#   exp( -(2*l^2)^{-1} * (x0 - x1)^2)
# }
# KERNEL_VALUE_EQUALS <- 1
# 
# Ntrain <- N
# kernel_bandwith <- 0.3
# Kmat <- array(data = NA, dim = c(Ntrain,Ntrain))
# Kmat[1,1] <- 1
# for(i in 2:Ntrain){
#   Kmat[i,i] <- 1
#   temp <- my_kernel(x0 = xi[-i,2], x1 = xi[i,2], l = kernel_bandwith)
#   Kmat[i,-i] <- Kmat[-i,i] <- temp
# }
# summary(c(Kmat))
# sigma2 <- 1
# normi <- solve(Kmat + diag(Ntrain) * sigma2)
# normi_times_LvecTrain <- normi %*% (yi)
# xx <- seq(0,1,length.out = 100)
# kTi <- sapply(xx, function(xx_i) my_kernel(x0 = xx_i, x1 = xi[,2], l = kernel_bandwith))
# pred0 <- as.numeric( t(kTi) %*% normi_times_LvecTrain )
# 
# xx <- xi[,2]
# kTi <- sapply(xx, function(xx_i) my_kernel(x0 = xx_i, x1 = xi[,2], l = kernel_bandwith))
# pred0 <- as.numeric( t(kTi) %*% normi_times_LvecTrain )
# 
# index1 <- (xi[,1]==1)
# plot(xi[index1,2], yi[index1])
# points(x = xi[index1,2], y = pred0[index1], col = 'red')
# plot(x = xi[index1,2], y = yi[index1] - pred0[index1], col = 'red')
# 
# Kmat1 <- Kmat[index1, index1]
# 
# # Ntrain1 <- sum(index1)
# # Kmat2 <- array(data = NA, dim = c(Ntrain1, Ntrain1))
# # Kmat2[1,1] <- 1
# # for(i in 2:Ntrain1){
# #   Kmat2[i,i] <- 1
# #   temp <- my_kernel(x0 = xi[index1,2][-i], x1 = xi[index1,2][i], l = kernel_bandwith)
# #   Kmat2[i,-i] <- Kmat2[-i,i] <- temp
# # }
# new_Kmat1 <- solve(Kmat1 + diag(sum(index1)) * sigma2) %*% ((yi - pred0)[index1])
# xx <- seq(0,1,length.out = 100)
# kTi <- sapply(xx, function(xx_i) my_kernel(x0 = xx_i, x1 = xi[index1,2], l = kernel_bandwith))
# pred1 <- as.numeric( t(kTi) %*% new_Kmat1 )
# plot.ts(pred1)
