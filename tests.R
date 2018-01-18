rm(list = ls(all = T))
cat("\n \n \n \n \n \n \n \n")

source("simulation.R")

G = 5
P = 3
M = 3
MR = 1

network = rand_network(G, P, M, MR)
indiv = rand_indiv(network)

tmax = 500
system.time(sim <- simuInd(network, indiv, tmax))

parAT = paramAT(network, indiv)
system.time(simAT <- ssa.adaptivetau(parAT$init.values, parAT$transitions, parAT$rateFunc, parAT$params, tf = tmax))

dt = 1
for(i in 1:length(network)){assign(names(network)[i], network[[i]])}; for(i in 1:length(indiv)){assign(names(indiv)[i], indiv[[i]])}

for(m in colnames(sim$time_abundance)[-1]){
  plot(sim$time_abundance[,"time"], sim$time_abundance[,m], type = 'l', main = m, col = "blue")
  lines(simAT[,"time"], simAT[, m], col = "red")
}
