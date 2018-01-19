rm(list = ls(all = T))
cat("\n \n \n \n \n \n \n \n")

source("simulation.R")

G = 5
P = 3
M = 5
MR = 5

network = rand_network(G, P, M, MR)
indiv = rand_indiv(network)

tmax = 500
system.time(sim1 <- simuInd(network, indiv, tmax, dt = 0.1))
system.time(sim2 <- simuInd(network, indiv, tmax, dt = 0.05))

parAT = paramAT(network, indiv)
system.time(simAT <- ssa.adaptivetau(parAT$init.values, parAT$transitions, parAT$rateFunc, parAT$params, tf = tmax))

parDS = paramDeSolve(network, indiv)
system.time(simDS <- ode(parDS$y, 0:tmax, parDS$func, parDS$parms))

for(m in colnames(sim1$time_abundance)[-1]){
  # ymin = min(sim1$time_abundance[,m], simAT[,m], simDS[,m])
  # ymax = max(sim1$time_abundance[,m], simAT[,m], simDS[,m])
  ymin = min(sim1$time_abundance[,m], simDS[,m])
  ymax = max(sim1$time_abundance[,m], simDS[,m])
    plot(sim1$time_abundance[,"time"], sim1$time_abundance[,m], type = 'l', main = m, col = "blue", ylim = c(ymin, ymax))
  lines(sim2$time_abundance[,"time"], sim2$time_abundance[,m], col = 'green')
  lines(simDS[, "time"], simDS[, m], col = 'black')
  #lines(simAT[,"time"], simAT[, m], col = "red")
}


# ---------------------------------------------------------

rm(list = ls(all = T))
cat("\n \n \n \n \n \n \n \n")

source("simulation.R")

G = 5
P = 3
M = 5
MR = 5

network = rand_network(G, P, M, MR)
indiv = rand_indiv(network)

tmax = 500

system.time(sim <- simuInd(network, indiv, tmax, dt = 1))
system.time(sim1 <- simuInd(network, indiv, tmax, dt = 0.5))
system.time(sim2 <- simuInd(network, indiv, tmax, dt = 0.1))
system.time(sim3 <- simuInd(network, indiv, tmax, dt = 0.05))

parDS = paramDeSolve(network, indiv)
system.time(simDS <- ode(parDS$y, 0:tmax, parDS$func, parDS$parms))

for(m in colnames(sim$time_abundance)[-1]){
  ymin = min(sim$time_abundance[,m], sim1$time_abundance[,m], sim2$time_abundance[,m], sim3$time_abundance[,m], simDS[,m])
  ymax = max(sim$time_abundance[,m], sim1$time_abundance[,m], sim2$time_abundance[,m], sim3$time_abundance[,m], simDS[,m])
    plot(sim$time_abundance[,"time"], sim$time_abundance[,m], type = 'l', main = m, col = "red", ylim = c(ymin, ymax))
  lines(sim1$time_abundance[,"time"], sim1$time_abundance[,m], col = "orange")
  lines(sim2$time_abundance[,"time"], sim2$time_abundance[,m], col = "yellow")
  lines(sim3$time_abundance[,"time"], sim3$time_abundance[,m], col = "green")
  #lines(sim4$time_abundance[,"time"], sim4$time_abundance[,m], col = "blue")
  lines(simDS[, "time"], simDS[, m], col = 'black')
}



# ---------------------------------------------------------

rm(list = ls(all = T))
cat("\n \n \n \n \n \n \n \n")

source("simulation.R")

G = 5
P = 3
M = 5
MR = 5

network = rand_network(G, P, M, MR)
indiv = rand_indiv(network)

tmax = 500

system.time(sim <- simuInd(network, indiv, tmax, dt = 1))
system.time(simP <- simuIndPrim(network, indiv, tmax, dt = 1))

parDS = paramDeSolve(network, indiv)
system.time(simDS <- ode(parDS$y, 0:tmax, parDS$func, parDS$parms))

for(m in colnames(sim$time_abundance)[-1]){
  ymin = min(sim$time_abundance[,m], sim1$time_abundance[,m], sim2$time_abundance[,m], sim3$time_abundance[,m], simDS[,m])
  ymax = max(sim$time_abundance[,m], sim1$time_abundance[,m], sim2$time_abundance[,m], sim3$time_abundance[,m], simDS[,m])
  plot(sim$time_abundance[,"time"], sim$time_abundance[,m], type = 'l', main = m, col = "red", ylim = c(ymin, ymax))
  lines(simP$time_abundance[,"time"], simP$time_abundance[,m], col = "orange")
  lines(simDS[, "time"], simDS[, m], col = 'black')
}

