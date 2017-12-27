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

if(!suppressWarnings(require("GillespieSSA", quietly = T))){install.packages("GillespieSSA")}
library(GillespieSSA)

# --------- #
# Main code #
# --------- #

tmax = 1000  

negative_feedback = rand_network_null(1,1,1)
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
library(GillespieSSA)

# ---------------------------------------------------------------------- #
# Transformation sampled network and cohort into GillespieSSA parameters #
# ---------------------------------------------------------------------- #

paramSSA = function(network, cohort){
  
  G = length(nw1$genes)
  P = length(nw1$prot)
  M = length(nw1$met)
  
  protNAA = c(sapply(network$prot, function(p){paste(p,"NA",sep = "_")}), sapply(network$prot, function(p){paste(p,"A",sep = "_")}))
  
  # INITIAL CONDITIONS ----
  x0 = c(cohort$rna_0[,1], cohort$prot_NA_0[,1], cohort$prot_A_0[,1], cohort$met_tot_0)
  names(x0) = c(network$genes, protNAA, network$met)
  
  
  # PROPENSITY FUNCTIONS / REACTION RATES ---- 
  regulation_law = function(id, reaction, network){
    reg = names(which(network[[paste0(reaction,"_sgn")]][id,]!=0))
    reg[grepl('^P',reg)] = paste(reg[grepl('^P',reg)],"A",sep = "_") # consider only active proteins
    parname = paste0(reaction, id)
    sup = sapply(reg, function(x){paste(c("*(1+sgn",parname,x,"*(",x,"^n",parname,x,"/(",x,"^n",parname,x,"+th",parname,x,"^n",parname,x,")))"), collapse = "")})
    return(sup)
  }
  
  regulation_law2 = function(id, reaction, network){
    reg = names(which(network[[paste0(reaction,"_sgn")]][id,]!=0))
    reg[grepl('^P',reg)] = paste(reg[grepl('^P',reg)],"A",sep = "_") # consider only active proteins
    parname = paste0(reaction, id)
    sup = sapply(reg, function(x){paste(c("*(",x,"^n",parname,x,"/(",x,"^n",parname,x,"+th",parname,x,"^n",parname,x,"))"), collapse = "")})
    return(sup)
  }
  
  a = c( # ----
         # transcription reactions
         sapply(network$genes, function(g){
           sup = regulation_law(g, 'TF', network)
           paste(c("k_TC",g,sup), collapse = "")}),
         # translation reactions
         sapply(network$prot, function(p){           
           sup = regulation_law(network$g2p[p], 'TLF', network)
           paste(c(network$g2p[p],"*k_TL",p,sup), collapse = "")}),
         # RNA decay reactions
         sapply(network$genes, function(g){
           sup = regulation_law(g, 'DR', network)
           paste(c(g,"*p0_DR",g,sup), collapse = "")}),
         # Protein decay reactions (for active proteins then for inactive proteins)
         sapply(protNAA, function(p){
           sup = regulation_law(sub("_NA|_A","",p), 'DP', network)
           paste(c(p,"*p0_DP",sub("_NA|_A","",p),sup), collapse = "")}),
         # Protein activation reactions
         sapply(protNAA[1:P], function(p){
           sup = regulation_law2(sub("_NA","",p), 'ACT', network)
           paste(c(p,sup), collapse = "")}),
         # Protein inactivation reactions
         sapply(protNAA[(P+1):(2*P)], function(p){
           sup = regulation_law2(sub("_A","",p), 'DEACT', network)
           if(length(sup) == 0){return("0")}
           else{paste(c(p,sup), collapse = "")}})
  ) # ----
  
  # STATE-CHANGE VECTOR ----
  tempTCTL = diag(1, nrow = G+P, ncol = G+P); tempTCTL = rbind(tempTCTL, matrix(0, nrow = P, ncol = G+P))
  tempDRDP = diag(-1, nrow = G+2*P, ncol = G+2*P)
  tempPos = diag(1, nrow = P, ncol = P); tempNeg = diag(-1, nrow = P, ncol = P)
  tempACT = rbind(matrix(0, nrow = G, ncol = P), tempNeg, tempPos)
  tempDEACT = rbind(matrix(0, nrow = G, ncol = P), tempPos, tempNeg)
  nu = cbind(tempTCTL, tempDRDP, tempACT, tempDEACT); nu = rbind(nu, matrix(0, nrow = M, ncol = ncol(nu))); rownames(nu) = names(x0)
  
  # REACTION RATE PARAMETERS ----
  combnames = function(par, reaction, target, reg){
    reg[grepl('^P',reg)] = paste(reg[grepl('^P',reg)],"A",sep = "_")
    comb = expand.grid(target, reg)
    if(nrow(comb)!=0){sapply(1:nrow(comb), function(i){paste0(par, reaction, comb[i,1], comb[i, 2])})}
  }
  
  parms = c( # ----
             # Transcription parameters
             network$k_TC, # Basal transcription rates
             as.vector(network$TF_sgn), # regulation direction
             as.vector(network$TF_th*cohort$QTL_TC[,1]), # regulation threshold
             as.vector(network$TF_n), # regulation power
             # Translation parameters
             network$k_TL, # Basal translation rates
             as.vector(network$TLF_sgn), # regulation direction
             as.vector(network$TLF_th*cohort$QTL_TL[,1]), # regulation threshold
             as.vector(network$TLF_n), # regulation power
             # RNA decay parameters
             network$p0_DR*cohort$QTL_TC[,1], # Basal decay rates
             as.vector(network$DR_sgn), # regulation direction
             as.vector(network$DR_th), # regulation threshold
             as.vector(network$DR_n), # regulation power
             # protein decay parameters
             network$p0_DP, # Basal decay rates
             as.vector(network$DP_sgn), # regulation direction
             as.vector(network$DP_th), # regulation threshold
             as.vector(network$DP_n), # regulation power
             # protein activation parameters
             as.vector(network$ACT_sgn), # regulation direction
             as.vector(network$ACT_th), # regulation threshold
             as.vector(network$ACT_n), # regulation power
             # protein inactivation parameters
             as.vector(network$DEACT_sgn), # regulation direction
             as.vector(network$DEACT_th), # regulation threshold
             as.vector(network$DEACT_n) # regulation power
  ) # ----
  
  names(parms) = c( # ----
                    # Transcription parameters
                    sapply(network$genes, function(x){paste0("k_TC",x)}), # Basal transcription rates
                    combnames("sgn", "TF", rownames(network$TF_sgn), colnames(network$TF_sgn)), # regulation direction
                    combnames("th", "TF", rownames(network$TF_th), colnames(network$TF_th)), # regulation threshold
                    combnames("n", "TF", rownames(network$TF_n), colnames(network$TF_n)), # regulation power
                    # Translation parameters
                    sapply(network$prot, function(x){paste0("k_TL",x)}), # Basal translation rates
                    combnames("sgn", "TLF", rownames(network$TLF_sgn), colnames(network$TLF_sgn)), # regulation direction
                    combnames("th", "TLF", rownames(network$TLF_th), colnames(network$TLF_th)), # regulation threshold
                    combnames("n", "TLF", rownames(network$TLF_n), colnames(network$TLF_n)), # regulation power
                    # RNA decay parameters
                    sapply(network$genes, function(x){paste0("p0_DR",x)}), # Basal decay rates
                    combnames("sgn", "DR", rownames(network$DR_sgn), colnames(network$DR_sgn)), # regulation direction
                    combnames("th", "DR", rownames(network$DR_th), colnames(network$DR_th)), # regulation threshold
                    combnames("n", "DR", rownames(network$DR_n), colnames(network$DR_n)), # regulation power
                    # protein decay parameters
                    sapply(network$prot, function(x){paste0("p0_DP",x)}), # Basal decay rates
                    combnames("sgn", "DP", rownames(network$DP_sgn), colnames(network$DP_sgn)), # regulation direction
                    combnames("th", "DP", rownames(network$DP_th), colnames(network$DP_th)), # regulation threshold
                    combnames("n", "DP", rownames(network$DP_n), colnames(network$DP_n)), # regulation power
                    # protein activation parameters
                    combnames("sgn", "ACT", rownames(network$ACT_sgn), colnames(network$ACT_sgn)), # regulation direction
                    combnames("th", "ACT", rownames(network$ACT_th), colnames(network$ACT_th)), # regulation threshold
                    combnames("n", "ACT", rownames(network$ACT_n), colnames(network$ACT_n)), # regulation power
                    # protein inactivation parameters
                    combnames("sgn", "DEACT", rownames(network$DEACT_sgn), colnames(network$DEACT_sgn)), # regulation direction
                    combnames("th", "DEACT", rownames(network$DEACT_th), colnames(network$DEACT_th)), # regulation threshold
                    combnames("n", "DEACT", rownames(network$DEACT_n), colnames(network$DEACT_n)) # regulation power
  ) # ----
  
  res = list("x0" = x0, "a" = a, "nu" = nu, "parms" = parms)
  return(res)
}


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

tmax = 500
nsim = 100

nw1 = rand_network_null(5,3,1)
cohort = rand_cohort(nw1,1)

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
#   Visualization    #
# ------------------ #


for(mol in c(nw1$genes, nw1$prot)){
  ymin = min(quant2_5_sim[[mol]], quant2_5_simSSA[[mol]])
  ymax = max(quant97_5_sim[[mol]], quant97_5_simSSA[[mol]])
  plot(0:tmax, mean_sim[[mol]], type = 'l', col = "blue", main = mol, xlab = "time", ylab = "abundance", ylim = c(ymin,ymax))
  lines(0:tmax, mean_simSSA[[mol]], col = "red")
  lines(0:tmax, quant2_5_sim[[mol]], col = "blue", lty = "dotted")
  lines(0:tmax, quant97_5_sim[[mol]], col = "blue", lty = "dotted")
  lines(0:tmax, quant2_5_simSSA[[mol]], col = "red", lty = "dotted")
  lines(0:tmax, quant97_5_simSSA[[mol]], col = "red", lty = "dotted")
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
