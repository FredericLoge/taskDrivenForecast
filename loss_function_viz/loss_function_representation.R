# loss function representation
pars_grid <- expand.grid(
  'd' = with(problem_constants, seq(from = P_min, to = P_max, length.out = 50)),
  'tau' = with(problem_constants, seq(from = tau_min, to = tau_max, length.out = 50))
)
pars_grid$val <- NA
for(i in 1:nrow(pars_grid)){
  pars_grid$val[i] <- cost_foo(pars = problem_constants, D = pars_grid$d[i], tau = pars_grid$tau[i])
}
pars_grid$opt <- FALSE
levels_d <- unique(pars_grid$d)
for(i in 1:length(levels_d)){
  index <- which.min(pars_grid$val[pars_grid$d == levels_d[i]])
  pars_grid$opt[which(pars_grid$d == levels_d[i])[index]] <- TRUE
}
library(ggplot2)
ggplot(data = pars_grid) +
  geom_tile(mapping = aes(x = d, y = tau, fill = val)) +
  scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue') +
  geom_point(data = pars_grid[pars_grid$opt,], mapping = aes(x = d, y = tau, col = val), pch = 20, cex = 5) +
  xlab(expression(Customer~demand~italic(d))) +
  ylab(expression(tau)) +
  labs(fill = expression(bold(tilde(L))~group("(",list(tau, italic(d)),")")),
       col = expression(min[tau]~bold(tilde(L))~group("(",list(tau, italic(d)),")"))) + 
  ggtitle(label = expression(bold(Loss~Function~tilde(L))~group("(",list(tau, italic(d)),")"))) + # 'Loss function L~(tau,d)') + 
  theme(
    axis.text = element_text(size = 20),
    axis.title.x = element_text(size = 20, margin = margin(t = 20)),
    axis.title.y = element_text(size = 20, margin = margin(r = 20), angle = 0, vjust = 0.5),
    plot.title = element_text(size = 30, margin = margin(b = 20))
  )
ggplot(data = pars_grid[abs(pars_grid$d - 200) < 10,], mapping = aes(x = tau, y = val, group = factor(d), col = factor(round(d)))) +
  geom_point() +
  geom_line()  +
  labs(col = expression(Customer~demand~italic(d))) +
  xlab(expression(tau)) +
  ylab(expression(bold(tilde(L))~group("(",list(tau, italic(d)),")"))) + 
  ggtitle(label = expression(bold(Loss~Function~tilde(L))~group("(",list(tau, italic(d)),")"))) + # 'Loss function L~(tau,d)') + 
  theme(
    axis.text = element_text(size = 20),
    axis.title.x = element_text(size = 20, margin = margin(t = 20)),
    axis.title.y = element_text(size = 20, margin = margin(r = 20), angle = 0, vjust = 0.5),
    plot.title = element_text(size = 30, margin = margin(b = 20))
  )
