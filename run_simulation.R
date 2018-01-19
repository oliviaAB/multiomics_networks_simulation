##########################################################################################################################
#                                             Test code on several runs                                                  #
##########################################################################################################################
library(parallel)
library(compiler)

rm(list = ls(all = T))
source("simulation.R")

tmax = 500
nsim = 100
dt = 0.1
simtime = seq(0, tmax, dt)

G = 10
P = 8
M = 20
MR = 25

network = rand_network(G, P, M, MR)
indiv = rand_indiv(network)

ncores = detectCores()

####  Our algorithm

system.time(simsIND <- mclapply(1:nsim, function(i){ 
  simuInd(network, indiv, tmax, dt = dt)$time_abundance
}, mc.cores = ncores))

res_simInd = vector("list", length = ncol(simsIND[[1]])-1); names(res_simInd) = colnames(simsIND[[1]])[-1]

for(sim in simsIND){
  for(m in names(res_simInd)){ res_simInd[[m]] = rbind(res_simInd[[m]], sim[,m]) }
}

mean_simInd = lapply(res_simInd, colMeans)
quant2_5_simInd = lapply(res_simInd, function(x){apply(x, 2, quantile, probs = 0.025)})
quant97_5_simInd = lapply(res_simInd, function(x){apply(x, 2, quantile, probs = 0.975)})


####  Adaptivetau 

parAT = paramAT(network, indiv)

enableJIT(1)
system.time(simsAT <- mclapply(1:nsim, function(i){
  simAT <- ssa.adaptivetau(parAT$init.values, parAT$transitions, parAT$rateFunc, parAT$params, tf = tmax)
  res = simAT[sapply(seq(0, tmax, dt), function(t){ max(which(simAT[,"time"]<=t)) }),]
  res[,1] = seq(0, tmax, dt)
  return(res)
}, mc.cores = ncores))

res_simAT = vector("list", length = ncol(simsAT[[1]])-1); names(res_simAT) = colnames(simsAT[[1]])[-1]

for(sim in simsAT){
  for(m in names(res_simAT)){ res_simAT[[m]] = rbind(res_simAT[[m]], sim[,m]) }
}

mean_simAT = lapply(res_simAT, colMeans)
quant2_5_simAT = lapply(res_simAT, function(x){apply(x, 2, quantile, probs = 0.025)})
quant97_5_simAT = lapply(res_simAT, function(x){apply(x, 2, quantile, probs = 0.975)})


### Deterministic solution

parDS = paramDeSolve(network, indiv)
system.time(simDS <- ode(parDS$y, simtime, parDS$func, parDS$parms))

for(mol in names(mean_simInd)){
  ymin = min(quant2_5_simInd[[mol]], simDS[, mol])#, quant2_5_simAT[[mol]])
  ymax = max(quant97_5_simInd[[mol]], simDS[, mol])#, quant97_5_simAT[[mol]])
  
  plot(simtime, mean_simInd[[mol]], type = 'l', col = "blue", main = mol, xlab = "time", ylab = "abundance", ylim = c(ymin, ymax))
  lines(simtime, quant2_5_simInd[[mol]], col = "blue", lty = "dotted")
  lines(simtime, quant97_5_simInd[[mol]], col = "blue", lty = "dotted")
  polygon(c(simtime, rev(simtime)), c(quant2_5_simInd[[mol]],rev(quant97_5_simInd[[mol]])), col = alpha("blue", alpha = 0.1), border = NA)
  
  # lines(simtime, mean_simAT[[mol]], col = "green")
  # lines(simtime, quant2_5_simAT[[mol]], col = "green", lty = "dotted")
  # lines(simtime, quant97_5_simAT[[mol]], col = "green", lty = "dotted")
  # polygon(c(simtime, rev(simtime)), c(quant2_5_simAT[[mol]],rev(quant97_5_simAT[[mol]])), col = alpha("green", alpha = 0.1), border = NA)
  
  lines(simDS[, "time"], simDS[, mol], col = "black", lwd = 1.5)
}

p = "P7_A"

plot(simtime[simtime<=50], simtime[simtime<=50], ylim = c(0, max(res_simInd[[p]])), type = 'n')
for(i in 1:nrow(res_simInd[[p]])) lines(simtime, res_simInd[[p]][i,])
