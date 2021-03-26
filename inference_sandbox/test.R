#----Assume
# Fundamental time unit is week
# Assume discrete mosquito psuedo-generations of 1 week each immature then mature
# All mosquitoes in previous gen die 
# Human pop size >> infected humans

#============================PARAMS & FUNS for DATA
# ------------------Bird pop
min.bird = 200
max.bird = 2000
scale = 2*pi/52
noise = rep(100,5)
bird.data = function(end, min, max, scale, noise){
  dat = ((cos(seq(1,end,by=1)*scale)+1)*((max.bird-min.bird)/2))+min.bird
  noise = stats::filter(rnorm(end), filter=noise, circular=TRUE)
  dat = dat+noise
  dat[dat<0]=0
  dat
}
# ------------------Mosquito pop
target = 1e6
p1 = 1
p2 = 1
over.dis = 1
mos.data = function(end, target, p1, p2, od){
  init = c(target*p1/(p1+p2), target*p2/(p1+p2))
  ret = matrix(-1, nrow=end, ncol=2)
  ret[1,] = init
  for(i in 2:end){
    bth = rnbinom(1, mu=init[2]*1/p2, size=over.dis)
    ret[i,2] = ret[i-1,1]
    ret[i,1] = bth
  }
  ret
}
#-------------------Human cases
agg = 52
human.cases = function(end, mean, od)rnbinom(end/agg, mu=mean, size=od)
#============================CONTROL PARAMS
lam1 = 1.5 # effective weekly biting, bird
lam2 = 2 # effective weekly biting, human
pop.size = 1e6

#============================Example instance
library(rstan)
library(ggplot2)

bd = bird.data(52*10, min.bird, max.bird, scale, noise)
md = mos.data(52*10, target, p1, p2, over.dis)
hd = human.cases(52*10, 100, 0.5)

dat = list(Np=52*10,
           No=10,
           lam1=lam1, # biting birds
           lam2=lam2,  # biting humans
           pop=pop.size,
           md=md[,2],
           hd=hd)

init.fun = function(){list(zeta=rep(1e-8, 10))}

fit = stan(file = "test.stan",
           data = dat,
           chains = 1,
           cores = 4,
           iter=5000,
           init = init.fun,
           verbose=F)

zeta = get_posterior_mean(fit, pars="zeta")
plot(hd, zeta)

