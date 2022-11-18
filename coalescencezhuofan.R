set.seed(2)
n = 100; px = 5; py = 5
x = matrix(rgamma(n * px, 1, 1), px, n)
x = t(t(x) / colSums(x))
e = matrix(rgamma(n * py, 1, 1), py, n)
e = t(t(e) / colSums(e))
m = matrix(rgamma(px * py, 1 / py, 1), py, px)
m = t(t(m) / colSums(m))
y = 0.8 * m %*% x + 0.2 * e

library(gtools)
x = t(rdirichlet(n,rep(0.01,px)))
e = t(rdirichlet(n,rep(0.1,py)))
# e = matrix(rgamma(n * py, 1, 1), py, n)
# e = t(t(e) / colSums(e))
m = matrix(rgamma(px * py, 1 / py, 1), py, px)
m = t(t(m) / colSums(m))
y = 0.8 * m %*% x + 0.2 * e
heatmap((x), Rowv = NA, Colv = NA)
# heatmap(t(oral_scale),Rowv = NA, Colv = NA)
heatmap(y, Rowv = NA, Colv = NA)

dd_data <- list(n = n, px = px, py = py, x = x, y = y)

library(rstan)

hmc_dd <- stan(file = "coalescence.stan", data = dd_data, chains = 1, 
               iter = 2000, cores = parallel::detectCores(), refresh = 0)
hmc_dd <- stan(file = "coalescencezhuofan.stan", data = dd_data, chains = 1, 
               iter = 2000, cores = parallel::detectCores(), refresh = 0)

model_dd = stan_model('coalescencezhuofan.stan')
vb_dd = vb(model_dd, data = dd_data, tol_rel_obj = 1e-5,iter=20000,seed=1)

#hmc_sample = hmc_dd@sim$samples[[1]]
vb_sample = vb_dd@sim$samples[[1]]

#hmc_mest = matrix(0, px, py)
vb_mest = matrix(0, py, px)

for(j in 1 : py) { for(l in 1 : px) {
  # hmc_mest[j, l] = mean(hmc_sample[1+px+py + (l - 1) * px + j][[1]])
  vb_mest[j, l] = mean(vb_sample[2+py*px+(l - 1) * py + j][[1]]) 
}}

heatmap(t(m), Rowv = NA, Colv = NA)
#heatmap(hmc_mest, Rowv = NA, Colv = NA)
heatmap(t(vb_mest), Rowv = NA, Colv = NA)

print(hmc_dd, pars = c("a", "b", "c", "m", "lp__"), probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
initl  = list("c" = 0.99,"m"= t(m))
initl  = list("c" = 0.99,"m"= rdirichlet(py,rep(1.0/px,px)))
initl = list(initl)
fit = sampling(model_dd,dd_data,iter=100,chains=1,cores = parallel::detectCores(),init = initl)

print(fit)
params = extract(fit)
# paramstrue = params
# rm(params)
sam_mest <- apply(params$mt,c(2,3),mean)
heatmap((m), Rowv = NA, Colv = NA)
#heatmap(hmc_mest, Rowv = NA, Colv = NA)
heatmap((sam_mest), Rowv = NA, Colv = NA)


### loglikelihood
ll = function(m,cstar,y,x){
  res = sum(log(ddirichlet(t(m),rep(1.0/px,px))))+dbeta(cstar,1,1,log = TRUE)+sum(log(ddirichlet(t(e),rep(0.1,py))))-py*log(1-cstar)
  return(res)
}
ll = function(m,cstar,rr,y,x){
  # tmp = t(y-cstar*m%*%x)/(1-cstar)
  # tmp = abs(tmp)*(tmp>=0)
  # tmp =tmp+0.0001
  # tmp = tmp/rowSums(tmp)
  res = sum(log(ddirichlet(t(m),rep(1.0/px,px))))+dbeta(cstar,1,1,log = TRUE)+sum(log(ddirichlet(rr,rep(0.1,py))))-py*log(1-cstar)
  return(res)
}
rec2 = vector()
for (i in 1:50) {
  rec2[i] = ll(params$mt[i,,],params$cstar[i],params$rr[i,,],y,x)
}
rec = vector()
for (i in 1:50) {
  rec[i] = ll(paramstrue$mt[i,,],paramstrue$cstar[i],paramstrue$rr[i,,],y,x)
}
plot(rec,type = "l")
plot(rec,type = "l",ylim = c(1716,1721))
abline(h=1720.219,lwd=1,col="blue")
truell = sum(log(ddirichlet(t(m),rep(1.0/px,px))))+dbeta(0.8,1,1,log = TRUE)+sum(log(ddirichlet(t(e),rep(0.1,py))))-py*log(0.2)


