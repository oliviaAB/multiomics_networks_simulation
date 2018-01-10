##########################################################################################################################
#                                                        MAIN                                                            #
##########################################################################################################################

source("data_simulation.R")

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

source("data_simulation.R")

# ----------------- #
#  One simulation   #
# ----------------- # ----

tmax = 1000

nw1 = rand_network(5,3,1)
cohort = rand_cohort(nw1,1)

# -- Our method --
sim = simu(nw1, cohort, tmax)

# -- GillespieSSA --
parSSA = paramSSA(nw1, cohort)

simSSA = ssa(parSSA$x0, parSSA$a, parSSA$nu, parSSA$parms, tmax, method = "ETL")

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

tmax = 8000
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
nw1 = rand_network(5,3,0)
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
  simSSA = ssa(parSSA$x0, parSSA$a, parSSA$nu, parSSA$parms, tmax, method = "OTL")
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
  windows()
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
