n = 100; px = 5; py = 5
x = matrix(rgamma(n * px, 1, 1), n, px)
x = x / rowSums(x)
e = matrix(rgamma(n * py, 1, 1), n, py)
e = e / rowSums(e)
m = matrix(rgamma(px * py, 1 / py, 1), px, py)
m = m / rowSums(m)
y = 0.8 * x %*% m + 0.2 * e

dd_data <- list(n = n, px = px, py = py, x = x, y = y)

library(rstan)

hmc_dd <- stan(file = "coalescence.stan", data = dd_data, chains = 1, warmup = 100, 
               iter = 500, cores = parallel::detectCores(), refresh = 0)

model_dd = stan_model('coalescence.stan')
vb_dd = vb(model_dd, data = dd_data, tol_rel_obj = 1e-4,iter=5000)

hmc_sample = hmc_dd@sim$samples[[1]]
vb_sample = vb_dd@sim$samples[[1]]

hmc_mest = matrix(0, px, py)
vb_mest = matrix(0, px, py)

for(j in 1 : px) { for(l in 1 : py) {
  # hmc_mest[j, l] = mean(hmc_sample[1+px+py + (l - 1) * px + j][[1]])
  vb_mest[j, l] = mean(vb_sample[1+px+py + (l - 1) * px + j][[1]]) 
}}
par(mfrow=c(1,2))
heatmap(m, Rowv = NA, Colv = NA)
#heatmap(hmc_mest, Rowv = NA, Colv = NA)
heatmap(vb_mest, Rowv = NA, Colv = NA)

print(hmc_dd, pars = c("a", "b", "c", "m", "lp__"), probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
