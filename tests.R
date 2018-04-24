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





# ---------------------------------------------------------
# ---------------------------------------------------------

library(adaptivetau)

x0 = c('P' = 1, 'R' = 0, 'PR' = 0, 'RNA' = 0)
trans = matrix(c(-1, -1, 1, 0,
                 1, 1, -1, 0, 
                 0, 0, 0, 1,
                 0, 1, 0, 0),
               nrow = 4)
params = c('pbind' = 0.05, 'punbind' = 0.5, 'TC' = 2, 'synthR' = 0.5)
rateFunc = function(y, params, t){
  return(c(params['pbind']*y['P']*y['R'],
           params['punbind']*y['PR'],
           params['TC']*y['PR'],
           params['synthR']))
}

rateFunc(x0, params, 1)

res1 = ssa.adaptivetau(x0, trans, rateFunc, params, 500)
plot(res1[,'time'], res1[,'RNA'], type= 'l')
lines(res1[,'time'], res1[,'R'], lty = 2)

####

x0_prim = c('R' = 0, 'RNA' = 0)
trans_prim = matrix(c(0,1, 1, 0), nrow = 2)
rateFunc_prim = function(y, params, t){
  return(c(params['TC']*y['R']/(y['R'] + params['punbind']/params['pbind']), params['synthR']))
}
rateFunc_prim(x0_prim, params, 1)

res2 = ssa.adaptivetau(x0_prim, trans_prim, rateFunc_prim, params, 500)
lines(res2[,'time'], res2[,'RNA'], col = 'blue')
lines(res2[,'time'], res2[,'R'], col = 'blue', lty = 2)


x0_prim2 = c('R' = 0, 'RNA' = 0)
trans_prim2 = matrix(c(0,1, 1, 0), nrow = 2)
rateFunc_prim2 = function(y, params, t){
  return(c(params['TC']*y['R'], params['synthR']))
}
rateFunc_prim2(x0_prim, params, 1)

res3 = ssa.adaptivetau(x0_prim2, trans_prim2, rateFunc_prim2, params, 500)
lines(res3[,'time'], res3[,'RNA'], col = 'green')
lines(res3[,'time'], res3[,'R'], col = 'green', lty = 2)


# Plot transcription rate

plot(res1[, 'time'], params['TC']*res1[,'PR'], col = 'black')
lines(res2[, 'time'], params['TC']*res2[,'R']/(res2[,'R'] + params['punbind']/params['pbind']), col = 'blue')
lines(res3[, 'time'], params['TC']*res3[,'R'], col = 'green')



