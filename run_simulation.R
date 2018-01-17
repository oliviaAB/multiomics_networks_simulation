##########################################################################################################################
#                                                        MAIN                                                            #
##########################################################################################################################

source("simulation.R")

# Random simulation with execution time
tmax = 1000  

nw1 = rand_network_null(1000,800,1)
cohort = rand_cohort(nw1,50)
system.time(sim <- simu(nw1, cohort, tmax))
visu(nw1,cohort,sim, tmax)



##########################################################################################################################
#                                             TEST NEGATIVE FEEDBACK LOOP                                                #
#                                   One gene, its protein repress its transcription                                      #
##########################################################################################################################

library(GillespieSSA)

# --------- #
# Main code #
# --------- #

tmax = 1000  

negative_feedback = rand_network_null(1,1,0)
negative_feedback$TF_sgn[1] = -1
negative_feedback$TF_th[1] = 500
negative_feedback$TF_n[1] = 2

cohort = rand_cohort(negative_feedback,1)
sim = simu(negative_feedback, cohort, tmax)
#visu(negative_feedback,cohort,sim, tmax)


# ------------ #
# GillespieSSA #
# ------------ #


# ---- Without the active/inactive state for proteins ----

x0 = c(RNA1 = cohort$rna_0, P1 = cohort$prot_tot_0)
a = c("k_TC*(1-P1^n/(P1^n+th^n))", #Transcription reaction
      "k_TL*RNA1", #Translation reaction
      "p0_DR*RNA1", #RNA decay reaction
      "p0_DP*P1") #Protein decay rate
nu = matrix(c(+1,0,0,+1,-1,0,0,-1), nrow = 2)
params = c(k_TC = negative_feedback$k_TC[[1]],
           n = negative_feedback$TF_n[[1]]*cohort$QTL_TC[[1]],
           th = negative_feedback$TF_th[[1]],
           k_TL = negative_feedback$k_TL[[1]],
           p0_DR = negative_feedback$p0_DR[[1]],
           p0_DP = negative_feedback$p0_DP[[1]])
simSSA = ssa(x0, a, nu, params, tmax, method = "ETL")


# ---- Taking into account activated and inactivated proteins ----

x0_prim = c(RNA1 = cohort$rna_0, P1_NA = cohort$prot_NA_0, P1_A = cohort$prot_A_0)
a_prim = c("k_TC*(1-P1_A^n/(P1_A^n+th^n))", #Transcription reaction
          "k_TL*RNA1", #Translation reaction
          "p0_DR*RNA1", #RNA decay reaction
          "p0_DP*P1_NA", #Non-activated protein decay rate
          "p0_DP*P1_A", #Activated protein
          "P1_NA") #Protein activation
nu_prim = matrix(c(+1,0,0, 0,+1,0, -1,0,0, 0,-1,0, 0,0,-1, 0,-1,1), nrow = 3)
params_prim = c(k_TC = negative_feedback$k_TC[[1]],
           n = negative_feedback$TF_n[[1]]*cohort$QTL_TC[[1]],
           th = negative_feedback$TF_th[[1]],
           k_TL = negative_feedback$k_TL[[1]],
           p0_DR = negative_feedback$p0_DR[[1]],
           p0_DP = negative_feedback$p0_DP[[1]])
simSSA_prim = ssa(x0_prim, a_prim, nu_prim, params_prim, tmax, method = "ETL")

plot(1:tmax, sim$time_rna$RNA1, type = "l", col = "blue") 
lines(simSSA$data[,1], simSSA$data[,"RNA1"], col = "red")
lines(simSSA_prim$data[,1], simSSA_prim$data[,"RNA1"], col = "green")

plot(1:tmax, sim$time_prot_tot$P1, type = "l", col = "blue") 
lines(simSSA$data[,1], simSSA$data[,"P1"], col = "red")
lines(simSSA_prim$data[,1], (simSSA_prim$data[,"P1_NA"]+simSSA_prim$data[,"P1_A"]), col = "green")





##########################################################################################################################
#                                           COMPARISON GILLESPIESSA VS OUR CODE                                          #
#                                                   ON RANDOM NETWORKS                                                   #
##########################################################################################################################

source("simulation.R")

# ----------------- #
#  One simulation   #
# ----------------- # ----

tmax = 1000

nw1 = rand_network(5,3,5)
cohort = rand_cohort(nw1,1)

# -- Our method --
sim = simu(nw1, cohort, tmax)

# -- GillespieSSA --
parSSA = paramSSA(nw1, cohort)

simSSA = ssa(parSSA$x0, parSSA$a, parSSA$nu, parSSA$parms, tmax, method = "OTL")

# -- Visualization --

for(g in nw1$genes){
  ymax = max(sim$time_rna[[g]], simSSA$data[,g])
  plot(0:tmax, sim$time_rna[[g]], type = "l", col = "blue", main = g, xlab = "time", ylab = "abundance", ylim = c(0, ymax))  
   lines(simSSA$data[,1], simSSA$data[,g], col = "red")
}
 
for(p in nw1$prot){
  ymax = max(sim$time_prot_tot[[p]], (simSSA$data[,paste0(p,"_NA")]+simSSA$data[,paste0(p,"_A")]))
  plot(0:tmax, sim$time_prot_tot[[p]], type = "l", col = "blue", main = p, xlab = "time", ylab = "abundance", ylim = c(0, ymax)) 
  lines(simSSA$data[,1], (simSSA$data[,paste0(p,"_NA")]+simSSA$data[,paste0(p,"_A")]), col = "red")
}

# ----

# ------------------------------ #
# Summary of several simulations #
# ------------------------------ # ---

tmax = 1000
nsim = 5

# Negative feedback on transcription ----
negative_feedback = rand_network_null(1,1,0)
negative_feedback$TF_sgn[1] = -1
negative_feedback$TF_th[1] = 500
negative_feedback$TF_n[1] = 2
nw1 = negative_feedback
cohort = rand_cohort(nw1,1)

save(nw1, cohort, file = "negative_feedback.RData")
load("negative_feedback.RData")

cohort$rna_0[1] = 100
cohort$prot_A_0[1] = 2000
cohort$prot_NA_0[1] = 0
cohort$prot_tot_0 = cohort$prot_NA_0 + cohort$prot_A_0
# ----

nw1 = rand_network_null(1,1,0)
nw1$TF_sgn[1] = -1
nw1$TF_th[1] = 1000
nw1$TF_n[1] = 2
nw1$TF_fc[1] = 1

nw1$k_TC[1] = 1/480
nw1$k_TL[1] = 1/60
nw1$p0_DR[1] = 1/300
nw1$p0_DP[1] = 1/7200

cohort = rand_cohort_null(nw1, 1)
cohort$rna_0[1] = 100
cohort$prot_NA_0[1] = 0
cohort$prot_A_0[1] = 300



# Random network ----
nw1 = rand_network(5,3,5)
cohort = rand_cohort(nw1,1)
# ----


network_plot(nw1)

G = length(nw1$genes)
P = length(nw1$prot)
M = length(nw1$met)

# run nsim times the simulation to get a confidence interval/variation interval around the abundance of each molecule ----
res_sim = vector("list", length = G+P); names(res_sim) = c(nw1$genes, nw1$prot)

for(s in 1:nsim){ 
  sim = simu(nw1, cohort, tmax)
  for(g in nw1$genes){ res_sim[[g]] = rbind(res_sim[[g]], sim$time_rna[[g]][1,]) }
  for(p in nw1$prot){ res_sim[[p]] = rbind(res_sim[[p]], sim$time_prot_tot[[p]][1,]) }
}

# Get the mean and the 2,5% and 97,5% quantiles
mean_sim = lapply(res_sim, colMeans)
quant2_5_sim = lapply(res_sim, function(x){apply(x, 2, quantile, probs = 0.025)})
quant97_5_sim = lapply(res_sim, function(x){apply(x, 2, quantile, probs = 0.975)})

# ------------------ #
#    GillespieSSA    #
# ------------------ #

parSSA = paramSSA(nw1, cohort)

# run nsim times the simulation to get a confidence interval/variation interval around the abundance of each molecule ----
res_simSSA = vector("list", length = G+P); names(res_simSSA) = c(nw1$genes, nw1$prot)

for(s in 1:nsim){ 
  simSSA = ssa(parSSA$x0, parSSA$a, parSSA$nu, parSSA$parms, tmax, method = "ETL")
  # We want to keep the state of the system for each discrete time 0,1,2 .. tmax
  # We take the state of each time step to be the state of the system at each time point closest but inferior to the time step (ie for time step 1, if we have the 
  # state of the system for the time point for 0, 0.3, 0.9 and 1.2 we choose the time point 0.9 to represent time step 1)
  for(g in nw1$genes){ res_simSSA[[g]] = rbind(res_simSSA[[g]], sapply(0:tmax, function(i){simSSA$data[max(which(simSSA$data[,1]<=i)),g]})) }
  for(p in nw1$prot){ res_simSSA[[p]] = rbind(res_simSSA[[p]], sapply(0:tmax, function(i){
      temptime =  max(which(simSSA$data[,1]<=i)) 
      return(simSSA$data[temptime,paste0(p,"_NA")] + simSSA$data[temptime,paste0(p,"_A")])})) }
}

# Get the mean and the 2,5% and 97,5% quantiles
mean_simSSA = lapply(res_simSSA, colMeans)
quant2_5_simSSA = lapply(res_simSSA, function(x){apply(x, 2, quantile, probs = 0.025)})
quant97_5_simSSA = lapply(res_simSSA, function(x){apply(x, 2, quantile, probs = 0.975)})



# ------------------ #
#      deSolve       #
# ------------------ #

parDeSolve = paramDeSolve(nw1, cohort)

simODE = ode(parDeSolve$y, 0:tmax, parDeSolve$func, parDeSolve$parms)
resODE = matrix(NA, nrow = (tmax+1), ncol = (G+P)) ; colnames(resODE) = c(nw1$genes, nw1$prot)
resODE[,nw1$genes] = simODE[,nw1$genes]
resODE[,nw1$prot] = sapply(nw1$prot, function(p){apply(simODE[,grepl(p,colnames(simODE))],1,sum)})

# ------------------ #
#   Visualization    #
# ------------------ #


for(mol in c(nw1$genes, nw1$prot)){
  # windows()
  ymin = min(quant2_5_sim[[mol]], quant2_5_simSSA[[mol]])
  ymax = max(quant97_5_sim[[mol]], quant97_5_simSSA[[mol]])
  plot(0:tmax, mean_sim[[mol]], type = 'l', col = "blue", main = mol, xlab = "time", ylab = "abundance", ylim = c(ymin,ymax))
  lines(0:tmax, mean_simSSA[[mol]], col = "red")
  lines(0:tmax, quant2_5_sim[[mol]], col = "blue", lty = "dotted")
  lines(0:tmax, quant97_5_sim[[mol]], col = "blue", lty = "dotted")
  lines(0:tmax, quant2_5_simSSA[[mol]], col = "red", lty = "dotted")
  lines(0:tmax, quant97_5_simSSA[[mol]], col = "red", lty = "dotted")
  lines(0:tmax, resODE[,mol], col = "black")
  polygon(c(0:tmax, tmax:0), c(quant2_5_sim[[mol]],rev(quant97_5_sim[[mol]])), col = alpha("blue", alpha = 0.1), border = NA)
  polygon(c(0:tmax, tmax:0), c(quant2_5_simSSA[[mol]],rev(quant97_5_simSSA[[mol]])), col = alpha("red", alpha = 0.1), border = NA)
}


# # Deterministic solution - For one protein-coding gene, no regulation
# 
# mol = "RNA1"
# ymin = min(quant2_5_sim[[mol]], quant2_5_simSSA[[mol]])
# ymax = max(quant97_5_sim[[mol]], quant97_5_simSSA[[mol]])
# plot(0:tmax, mean_sim[[mol]], type = 'l', col = "blue", main = mol, xlab = "time", ylab = "abundance", ylim = c(ymin,ymax))
# lines(0:tmax, mean_simSSA[[mol]], col = "red")
# lines(0:tmax, quant2_5_sim[[mol]], col = "blue", lty = "dotted")
# lines(0:tmax, quant97_5_sim[[mol]], col = "blue", lty = "dotted")
# lines(0:tmax, quant2_5_simSSA[[mol]], col = "red", lty = "dotted")
# lines(0:tmax, quant97_5_simSSA[[mol]], col = "red", lty = "dotted")
# pDR = nw1$p0_DR[1] * cohort$QTL_DR[1,1]; C = cohort$rna_0[1,1]-(nw1$k_TC[1]/pDR)
# lines(0:tmax, C*exp(-pDR*(0:tmax))+(nw1$k_TC/pDR), col = "black")
# 
# 
# mol = "P1"
# ymin = min(quant2_5_sim[[mol]], quant2_5_simSSA[[mol]])
# ymax = max(quant97_5_sim[[mol]], quant97_5_simSSA[[mol]])
# plot(0:tmax, mean_sim[[mol]], type = 'l', col = "blue", main = mol, xlab = "time", ylab = "abundance", ylim = c(ymin,ymax))
# lines(0:tmax, mean_simSSA[[mol]], col = "red")
# lines(0:tmax, quant2_5_sim[[mol]], col = "blue", lty = "dotted")
# lines(0:tmax, quant97_5_sim[[mol]], col = "blue", lty = "dotted")
# lines(0:tmax, quant2_5_simSSA[[mol]], col = "red", lty = "dotted")
# lines(0:tmax, quant97_5_simSSA[[mol]], col = "red", lty = "dotted")
# C2 = cohort$prot_tot_0[1,1]-cohort$rna_0[1,1]*(nw1$k_TL[1]/nw1$p0_DP[1])
# lines(0:tmax, C2*exp(-nw1$p0_DP[1]*(0:tmax))+(C*exp(-pDR*(0:tmax))+(nw1$k_TC/pDR))*(nw1$k_TL/nw1$p0_DP[1]), col = "black")



##########################################################################################################################
#                                                      TEST SIMIND                                                       #
##########################################################################################################################

source("simulation.R")

# Random network ----
nw1 = rand_network(5,3,5)
cohort = rand_cohort(nw1,1)
tmax = 1000

indiv = cohort
indiv["ind"] = NULL
for(i in names(indiv)){
  indiv[[i]] = indiv[[i]][,1]
}

indiv["indname"] = "ind1"

# --------- #
# Main code #
# --------- #

simC = simu(nw1, cohort, tmax)
simI = simuInd(nw1, indiv, tmax)

for(g in nw1$genes){
  ymin = min(simC$time_rna[[g]], simI$time_abundance[,g])
  ymax = max(simC$time_rna[[g]], simI$time_abundance[,g])
  
  plot(0:tmax, simC$time_rna[[g]], type = 'l', col = "red", ylim = c(ymin, ymax), main = g)
  lines(simI$time_abundance[,"time"], simI$time_abundance[,g], col = "blue")
}

for(p in nw1$prot){
  ymin = min(simC$time_prot_tot[[p]], simI$time_abundance[,p])
  ymax = max(simC$time_prot_tot[[p]], simI$time_abundance[,p])
  
  plot(0:tmax, simC$time_prot_tot[[p]], type = 'l', col = "red", ylim = c(ymin, ymax), main = p)
  lines(simI$time_abundance[,"time"], simI$time_abundance[,p], col = "blue")
}

# -------

nsim = 100

res_sim = vector("list", length = length(nw1$genes)+length(nw1$prot)); names(res_sim) = c(nw1$genes, nw1$prot)

system.time(for(s in 1:nsim){ 
    sim = simu(nw1, cohort, tmax)
    for(g in nw1$genes){ res_sim[[g]] = rbind(res_sim[[g]], sim$time_rna[[g]][1,]) }
    for(p in nw1$prot){ res_sim[[p]] = rbind(res_sim[[p]], sim$time_prot_tot[[p]][1,]) }
  }
)

# Get the mean and the 2,5% and 97,5% quantiles
mean_sim = lapply(res_sim, colMeans)
quant2_5_sim = lapply(res_sim, function(x){apply(x, 2, quantile, probs = 0.025)})
quant97_5_sim = lapply(res_sim, function(x){apply(x, 2, quantile, probs = 0.975)})



res_simInd = vector("list", length = length(nw1$genes)+length(nw1$prot)); names(res_simInd) = c(nw1$genes, nw1$prot)

system.time(for(s in 1:nsim){ 
  sim = simuInd(nw1, indiv, tmax)
  for(g in nw1$genes){ res_simInd[[g]] = rbind(res_simInd[[g]], sim$time_abundance[,g]) }
  for(p in nw1$prot){ res_simInd[[p]] = rbind(res_simInd[[p]], sim$time_abundance[,p]) }
}
)

# Get the mean and the 2,5% and 97,5% quantiles
mean_simInd = lapply(res_simInd, colMeans)
quant2_5_simInd = lapply(res_simInd, function(x){apply(x, 2, quantile, probs = 0.025)})
quant97_5_simInd = lapply(res_simInd, function(x){apply(x, 2, quantile, probs = 0.975)})


for(mol in c(nw1$genes, nw1$prot)){
  # windows()
  ymin = min(quant2_5_sim[[mol]], quant2_5_simInd[[mol]])
  ymax = max(quant97_5_sim[[mol]], quant97_5_simInd[[mol]])
  plot(0:tmax, mean_sim[[mol]], type = 'l', col = "blue", main = mol, xlab = "time", ylab = "abundance", ylim = c(ymin,ymax))
  lines(0:tmax, mean_simInd[[mol]], col = "red")
  lines(0:tmax, quant2_5_sim[[mol]], col = "blue", lty = "dotted")
  lines(0:tmax, quant97_5_sim[[mol]], col = "blue", lty = "dotted")
  lines(0:tmax, quant2_5_simInd[[mol]], col = "red", lty = "dotted")
  lines(0:tmax, quant97_5_simInd[[mol]], col = "red", lty = "dotted")
  #lines(0:tmax, resODE[,mol], col = "black")
  polygon(c(0:tmax, tmax:0), c(quant2_5_sim[[mol]],rev(quant97_5_sim[[mol]])), col = alpha("blue", alpha = 0.1), border = NA)
  polygon(c(0:tmax, tmax:0), c(quant2_5_simInd[[mol]],rev(quant97_5_simInd[[mol]])), col = alpha("red", alpha = 0.1), border = NA)
}



# --------------------------------------- #
#               New test                  #
# --------------------------------------- #


source("simulation.R")

nw1 = rand_network(5, 3, 0)
indiv = rand_indiv(nw1)
tmax = 1000

sim = simuInd(nw1, indiv, tmax)

for(m in c(nw1$genes, nw1$met)){
  plot(sim$time_abundance[,"time"], sim$time_abundance[,m], type = "l", main = m)
}

for(p in  nw1$prot){
  plot(sim$time_abundance[,"time"], sim$time_abundance[,p], col = "black", type = "l", main = p)
  lines(sim$time_abundance[,"time"], sim$time_abundance[,paste0(p, "_NA")], col = "blue")
  lines(sim$time_abundance[,"time"], sim$time_abundance[,paste0(p, "_A")], col = "red")
}











##########################################################################################################################
#                 USE PACKAGE ADAPTIVETAU TO COMPARE RESULTS FROM OUR ALGO TO APPROX. SSA SIMULATION                     #
##########################################################################################################################

library(adaptivetau)

source("simulation.R")

network = rand_network_null(5, 4, 2)
ind = rand_indiv(network)

G = length(network$genes)
P = length(network$prot)
M = length(network$met)

protNAA = c(sapply(network$prot, function(p){paste(p,"NA",sep = "_")}), sapply(network$prot, function(p){paste(p,"A",sep = "_")}))

tmax = 1000

parAT = paramAT(network, ind)
init.values = parAT$init.values
transitions = parAT$transitions
rateFunc = parAT$rateFunc
params = parAT$params

# -------------------------- #
#     Only 1 simulation      #
# -------------------------- #

system.time(sim <- simuInd(network, ind, tmax))

system.time(simAD <- ssa.adaptivetau(parAT$init.values, parAT$transitions, parAT$rateFunc, parAT$params, tf = tmax))

parSSA = paramSSAindiv(network, ind)
system.time(simSSA <- GillespieSSA::ssa(parSSA$x0, parSSA$a, parSSA$nu, parSSA$parms, tmax, method = "ETL") )

for(m in setdiff(colnames(sim$time_abundance)[-1], c(network$prot, network$met))){
  plot(sim$time_abundance[,"time"], sim$time_abundance[, m], type = 'l', col = 'blue', main = m)
  lines(simAD[, 'time'], simAD[,m], col = 'green')
  lines(simSSA$data[, 1], simSSA$data[,m], col = 'red')
  
}

# -------------------------- #
#       100 simulations      #
# -------------------------- #

nsim = 100

# Our algo ----
res_simInd = vector("list", length = length(network$genes)+length(network$prot)); names(res_simInd) = c(network$genes, network$prot)
system.time(for(s in 1:nsim){ 
  sim = simuInd(network, ind, tmax)
  for(g in network$genes){ res_simInd[[g]] = rbind(res_simInd[[g]], sim$time_abundance[,g]) }
  for(p in network$prot){ res_simInd[[p]] = rbind(res_simInd[[p]], sim$time_abundance[,p]) }
}
)

# Get the mean and the 2,5% and 97,5% quantiles
mean_simInd = lapply(res_simInd, colMeans)
quant2_5_simInd = lapply(res_simInd, function(x){apply(x, 2, quantile, probs = 0.025)})
quant97_5_simInd = lapply(res_simInd, function(x){apply(x, 2, quantile, probs = 0.975)})


# SSA algo ----
parSSA = paramSSAindiv(network, ind)
res_simSSA = vector("list", length = G+P); names(res_simSSA) = c(network$genes, network$prot)
system.time(for(s in 1:nsim){ 
  simSSA = ssa(parSSA$x0, parSSA$a, parSSA$nu, parSSA$parms, tmax, method = "ETL")
  # We want to keep the state of the system for each discrete time 0,1,2 .. tmax
  # We take the state of each time step to be the state of the system at each time point closest but inferior to the time step (ie for time step 1, if we have the 
  # state of the system for the time point for 0, 0.3, 0.9 and 1.2 we choose the time point 0.9 to represent time step 1)
  for(g in network$genes){ res_simSSA[[g]] = rbind(res_simSSA[[g]], sapply(0:tmax, function(i){simSSA$data[max(which(simSSA$data[,1]<=i)),g]})) }
  for(p in network$prot){ res_simSSA[[p]] = rbind(res_simSSA[[p]], sapply(0:tmax, function(i){
    temptime =  max(which(simSSA$data[,1]<=i)) 
    return(simSSA$data[temptime,paste0(p,"_NA")] + simSSA$data[temptime,paste0(p,"_A")])})) }
})

# Get the mean and the 2,5% and 97,5% quantiles
mean_simSSA = lapply(res_simSSA, colMeans)
quant2_5_simSSA = lapply(res_simSSA, function(x){apply(x, 2, quantile, probs = 0.025)})
quant97_5_simSSA = lapply(res_simSSA, function(x){apply(x, 2, quantile, probs = 0.975)})


# Adaptivetau ----
res_simAT = vector("list", length = G+P); names(res_simAT) = c(network$genes, network$prot)
enableJIT(1)
system.time(for(s in 1:nsim){
  simAT <- ssa.adaptivetau(init.values, transitions, rateFunc, params, tf = tmax)
  for(g in network$genes){ res_simAT[[g]] = rbind(res_simAT[[g]], sapply(0:tmax, function(i){simAT[max(which(simAT[,"time"]<=i)),g]})) }
  for(p in network$prot){ res_simAT[[p]] = rbind(res_simAT[[p]], sapply(0:tmax, function(i){
    temptime =  max(which(simAT[,"time"]<=i)) 
    return(simAT[temptime,paste0(p,"_NA")] + simAT[temptime,paste0(p,"_A")])})) }
  
})
# Get the mean and the 2,5% and 97,5% quantiles
mean_simAT = lapply(res_simAT, colMeans)
quant2_5_simAT = lapply(res_simAT, function(x){apply(x, 2, quantile, probs = 0.025)})
quant97_5_simAT = lapply(res_simAT, function(x){apply(x, 2, quantile, probs = 0.975)})




# Visualization ----
for(mol in c(network$genes, network$prot)){
  # windows()
  ymin = min(quant2_5_simInd[[mol]], quant2_5_simAT[[mol]])
  ymax = max(quant97_5_simInd[[mol]], quant97_5_simAT[[mol]])
  # ymin = min(quant2_5_simInd[[mol]], quant2_5_simSSA[[mol]], quant2_5_simAT[[mol]])
  # ymax = max(quant2_5_simInd[[mol]], quant2_5_simSSA[[mol]], quant2_5_simAT[[mol]])
  plot(0:tmax, mean_simInd[[mol]], type = 'l', col = "blue", main = mol, xlab = "time", ylab = "abundance", ylim = c(ymin,ymax))
  # lines(0:tmax, mean_simSSA[[mol]], col = "red")
  lines(0:tmax, mean_simAT[[mol]], col = "green")
  lines(0:tmax, quant2_5_simInd[[mol]], col = "blue", lty = "dotted")
  lines(0:tmax, quant97_5_simInd[[mol]], col = "blue", lty = "dotted")
  # lines(0:tmax, quant2_5_simSSA[[mol]], col = "red", lty = "dotted")
  # lines(0:tmax, quant97_5_simSSA[[mol]], col = "red", lty = "dotted")
  lines(0:tmax, quant2_5_simAT[[mol]], col = "green", lty = "dotted")
  lines(0:tmax, quant97_5_simAT[[mol]], col = "green", lty = "dotted")
  polygon(c(0:tmax, tmax:0), c(quant2_5_simInd[[mol]],rev(quant97_5_simInd[[mol]])), col = alpha("blue", alpha = 0.1), border = NA)
  # polygon(c(0:tmax, tmax:0), c(quant2_5_simSSA[[mol]],rev(quant97_5_simSSA[[mol]])), col = alpha("red", alpha = 0.1), border = NA)
  polygon(c(0:tmax, tmax:0), c(quant2_5_simAT[[mol]],rev(quant97_5_simAT[[mol]])), col = alpha("green", alpha = 0.1), border = NA)
}


##########################################################################################################################
#                                            PARALLELIZE CODE                                                            #
##########################################################################################################################

library(parallel)
library(adaptivetau)
library(compiler)

source("simulation.R")

ncores = detectCores()

network = rand_network(5, 4, 10, 4)
ind = rand_indiv(network)

G = length(network$genes)
P = length(network$prot)
M = length(network$met)

protNAA = c(sapply(network$prot, function(p){paste(p,"NA",sep = "_")}), sapply(network$prot, function(p){paste(p,"A",sep = "_")}))

tmax = 1000

parAT = paramAT(network, ind)
init.values = parAT$init.values
transitions = parAT$transitions
rateFunc = parAT$rateFunc
params = parAT$params


nsim = 100

# Our algo ----
system.time(sims <- mclapply(1:nsim, function(i){ 
  simuInd(network, ind, tmax)$time_abundance
}, mc.cores = ncores))

res_simInd = vector("list", length = ncol(sims[[1]])-1); names(res_simInd) = colnames(sims[[1]])[-1]

for(sim in sims){
  for(m in names(res_simInd)){ res_simInd[[m]] = rbind(res_simInd[[m]], sim[,m]) }
}

mean_simInd = lapply(res_simInd, colMeans)
quant2_5_simInd = lapply(res_simInd, function(x){apply(x, 2, quantile, probs = 0.025)})
quant97_5_simInd = lapply(res_simInd, function(x){apply(x, 2, quantile, probs = 0.975)})


# Adaptivetau ----

enableJIT(1)
system.time(sims <- mclapply(1:nsim, function(i){
  simAT <- ssa.adaptivetau(init.values, transitions, rateFunc, params, tf = tmax)
  res = simAT[sapply(0:tmax, function(t){ max(which(simAT[,"time"]<=t)) }),]
  res[,1] = 0:tmax
  res = cbind(res, sapply(network$prot, function(p){ rowSums(res[,grep(paste0('^',p,'_'), colnames(res))]) }))
  return(res)
}, mc.cores = ncores))

res_simAT = vector("list", length = ncol(sims[[1]])-1); names(res_simAT) = colnames(sims[[1]])[-1]

for(sim in sims){
  for(m in names(res_simAT)){ res_simAT[[m]] = rbind(res_simAT[[m]], sim[,m]) }
}

mean_simAT = lapply(res_simAT, colMeans)
quant2_5_simAT = lapply(res_simAT, function(x){apply(x, 2, quantile, probs = 0.025)})
quant97_5_simAT = lapply(res_simAT, function(x){apply(x, 2, quantile, probs = 0.975)})

for(mol in names(mean_simInd)){
  ymin = min(quant2_5_simInd[[mol]], quant2_5_simAT[[mol]])
  ymax = max(quant97_5_simInd[[mol]], quant97_5_simAT[[mol]])
  plot(0:tmax, mean_simInd[[mol]], type = 'l', col = "blue", main = mol, xlab = "time", ylab = "abundance", ylim = c(ymin,ymax))
  lines(0:tmax, mean_simAT[[mol]], col = "green")
  lines(0:tmax, quant2_5_simInd[[mol]], col = "blue", lty = "dotted")
  lines(0:tmax, quant97_5_simInd[[mol]], col = "blue", lty = "dotted")
  lines(0:tmax, quant2_5_simAT[[mol]], col = "green", lty = "dotted")
  lines(0:tmax, quant97_5_simAT[[mol]], col = "green", lty = "dotted")
  polygon(c(0:tmax, tmax:0), c(quant2_5_simInd[[mol]],rev(quant97_5_simInd[[mol]])), col = alpha("blue", alpha = 0.1), border = NA)
  polygon(c(0:tmax, tmax:0), c(quant2_5_simAT[[mol]],rev(quant97_5_simAT[[mol]])), col = alpha("green", alpha = 0.1), border = NA)
}






##########################################################################################################################
#                                            SIMULATION METABOLISM                                                       #
##########################################################################################################################

# for(i in 1:length(network)){assign(names(network)[i], network[[i]])}; for(i in 1:length(indiv)){assign(names(indiv)[i], indiv[[i]])}

rm(list = ls(all=T))

source("simulation.R")

tmax = 1000

network = rand_network(5, 4, 5, 10)
indiv = rand_indiv(network)


## Our algo
system.time(sim <- simuInd(network, indiv, tmax) )
sim = sim$time_abundance

## Adaptivetau
parAT = paramAT(network, indiv)
parAT$rateFunc(parAT$init.values, parAT$params, 1)
system.time(simAT <- ssa.adaptivetau(parAT$init.values, parAT$transitions, parAT$rateFunc, parAT$params, tf = tmax) )
simAT = cbind(simAT, sapply(network$prot, function(p){ rowSums(simAT[,grep(paste0(p,"_"), colnames(simAT))]) }))

## Deterministic solution
parDS = paramDeSolve(network, indiv)
system.time(simDS <- ode(parDS$y, seq(0, tmax, by = 1), parDS$func, parDS$parms))
simDS = cbind(simDS, sapply(network$prot, function(p){ rowSums(simDS[,grep(paste0(p,"_"), colnames(simDS))]) }))


## Plot
for(mol in colnames(sim)[-1]){
  ymin = min(sim[,mol], simAT[, mol], simDS[, mol])
  ymax = max(sim[,mol], simAT[, mol], simDS[, mol])
  plot(sim[,"time"], sim[,mol], col = "blue", ylim = c(ymin, ymax), type = 'l', main = mol)
  lines(simAT[,"time"], simAT[,mol], col = "red")
  lines(simDS[,"time"], simDS[,mol], col = "black")
}


