rm(list = ls())

setwd("C:/Users/oangelin/OneDrive - Massey University/Documents/multiomics_networks_simulation")

library(bioPN)
dyn.load("propensity.dll")
source("data_simulation.R")

PartitionedLeaping <- function(model, timep, delta=1 , runs=1, ect=1e-9) {
  pre <- model$pre
  post <- model$post
  h <- model$h
  M <- model$M

  storage.mode(pre) <- storage.mode(post) <- storage.mode(runs) <- "integer"
  #  storage.mode(slow) <- "integer"
  storage.mode(delta) <- storage.mode(M) <- storage.mode(ect) <- "double"
  
  if (is.null(place <- model$place)) {
    place <- paste("P",1:ncol(pre),sep="")
  }
  if (is.null(transition <- model$transition)) {
    transition <- paste("T",1:nrow(pre),sep="")
  }
  
  .Call("PartitionedLeaping", pre, post, h, M, timep, delta, runs, place, transition, ect, parent.frame())
}

# Petri Network

model = list(
  place = c("RNA", "P_NA","P_A"),
  transition = c("Transcription","Translation","RNAdecay","ProteinNAdecay","ProteinAdecay", "ProteinActivation"),
  pre = matrix(c(0,0,0,
                 1,0,0,
                 1,0,0,
                 0,1,0,
                 0,0,1,
                 0,1,0), ncol = 3, byrow = T),
  post = matrix(c(1,0,0,
                  1,1,0,
                  0,0,0,
                  0,0,0,
                  0,0,0,
                  0,0,1), ncol = 3, byrow = T),
  h = list( "Transcription" = getNativeSymbolInfo("Transcription", PACKAGE="propensity")$address,
            "Translation" = getNativeSymbolInfo("Translation", PACKAGE="propensity")$address,
            "RNAdecay" = getNativeSymbolInfo("RNAdecay", PACKAGE="propensity")$address,
            "ProteinNAdecay" = getNativeSymbolInfo("ProteinNAdecay", PACKAGE="propensity")$address,
            "ProteinAdecay" = getNativeSymbolInfo("ProteinAdecay", PACKAGE="propensity")$address,
            "ProteinActivation" = getNativeSymbolInfo("ProteinActivation", PACKAGE="propensity")$address),
  M = c("RNA" = 100, "P_NA" = 0, "P_A" = 300),
  slow = c(0 ,0 ,0)
)

model1 = model
model1$slow = c(1 ,1 ,1)

delta = 1
timep = 200

HR = HaseltineRawlings(model, timep, delta)
HR1 = HaseltineRawlings(model1, timep, delta)

GB = GibsonBruck(model, timep, delta)

# Our method

tmax = timep

nw1 = rand_network_null(1,1,0)
nw1$TF_sgn[1] = -1
nw1$TF_th[1] = 1000
nw1$TF_n[1] = 2
nw1$TF_fc[1] = 1

nw1$k_TC[1] = 5
nw1$k_TL[1] = 1
nw1$p0_DR[1] = 1/50
nw1$p0_DP[1] = 1/120

cohort = rand_cohort_null(nw1, 1)
cohort$rna_0[1] = 100
cohort$prot_NA_0[1] = 0
cohort$prot_A_0[1] = 300

sim = simu(nw1, cohort, tmax)

parSSA = paramSSA(nw1, cohort)
res_simSSA = vector("list", length = 2); names(res_simSSA) = c(nw1$genes, nw1$prot)
simSSA = ssa(parSSA$x0, parSSA$a, parSSA$nu, parSSA$parms, tmax, method = "OTL")
for(g in nw1$genes){ res_simSSA[[g]] = sapply(0:tmax, function(i){simSSA$data[max(which(simSSA$data[,1]<=i)),g]}) }
for(p in nw1$prot){ res_simSSA[[p]] = sapply(0:tmax, function(i){
  temptime =  max(which(simSSA$data[,1]<=i))
  return(simSSA$data[temptime,paste0(p,"_NA")] + simSSA$data[temptime,paste0(p,"_A")])}) }


ymin = min(HR$run[[1]]$M[[1]], HR1$run[[1]]$M[[1]], sim$time_rna[[1]], GB$run[[1]]$M[[1]])
ymax = max(HR$run[[1]]$M[[1]], HR1$run[[1]]$M[[1]], sim$time_rna[[1]], GB$run[[1]]$M[[1]])
plot(HR$dt, HR$run[[1]]$M[[1]], type = 'none', col = "red", ylim = c(ymin, ymax))
lines(HR1$dt, HR1$run[[1]]$M[[1]], col = "red")
lines(GB$dt, GB$run[[1]]$M[[1]], col = "darkred", lty = 2)
lines(0:tmax, sim$time_rna[[1]], col = "hotpink")
lines(0:tmax, res_simSSA[[1]], col = "hotpink", lty = 2)

windows()
ymin = min((HR$run[[1]]$M[[2]]+HR$run[[1]]$M[[3]]), (HR1$run[[1]]$M[[2]]+HR1$run[[1]]$M[[3]]), sim$time_prot_tot[[1]], (GB$run[[1]]$M[[2]]+GB$run[[1]]$M[[3]]))
ymax = max((HR$run[[1]]$M[[2]]+HR$run[[1]]$M[[3]]), (HR1$run[[1]]$M[[2]]+HR1$run[[1]]$M[[3]]), sim$time_prot_tot[[1]], (GB$run[[1]]$M[[2]]+GB$run[[1]]$M[[3]]))
plot(HR$dt, (HR$run[[1]]$M[[2]]+HR$run[[1]]$M[[3]]), type = 'none', col = "green", ylim = c(ymin, ymax))
lines(HR1$dt, (HR1$run[[1]]$M[[2]]+HR1$run[[1]]$M[[3]]), col = "green")
lines(GB$dt, (GB$run[[1]]$M[[2]]+GB$run[[1]]$M[[3]]), col = "green", lty = 2)
lines(0:tmax, sim$time_prot_tot[[1]], col = "blue")
lines(0:tmax, res_simSSA[[2]], col = "blue", lty = 2)



#############################################################################################################


library(ssar)
source("simulation.R")


negative_feedback = rand_network_null(1,1,0)
negative_feedback$TF_sgn[1] = -1
negative_feedback$TF_th[1] = 500
negative_feedback$TF_n[1] = 2
network = negative_feedback

cohort = rand_cohort_null(network,1)
cohort$rna_0[1] = 100
cohort$prot_A_0[1] = 100
cohort$prot_NA_0[1] = 0
cohort$prot_tot_0[1] = cohort$prot_A_0[1] + cohort$prot_NA_0[1]


# --------------------------------------- #
# Creating parameters for ssar simulation #
# --------------------------------------- #

G = length(network$genes)
P = length(network$prot)
M = length(network$met)

protNAA = c(sapply(network$prot, function(p){paste(p,"NA",sep = "_")}), sapply(network$prot, function(p){paste(p,"A",sep = "_")}))

# INITIAL CONDITIONS ----
x0 = matrix(c(cohort$rna_0[,1], cohort$prot_NA_0[,1], cohort$prot_A_0[,1], cohort$met_tot_0), nrow = 1)
colnames(x0) = c(network$genes, protNAA, network$met)

# Propensity matrix ----
tempTCTL = diag(1, nrow = G+P, ncol = G+P); tempTCTL = rbind(tempTCTL, matrix(0, nrow = P, ncol = G+P))
tempDRDP = diag(-1, nrow = G+2*P, ncol = G+2*P)
tempPos = diag(1, nrow = P, ncol = P); tempNeg = diag(-1, nrow = P, ncol = P)
tempACT = rbind(matrix(0, nrow = G, ncol = P), tempNeg, tempPos)
tempDEACT = rbind(matrix(0, nrow = G, ncol = P), tempPos, tempNeg)
v = cbind(tempTCTL, tempDRDP, tempACT, tempDEACT); v = rbind(v, matrix(0, nrow = M, ncol = ncol(v))); rownames(v) = colnames(x0)


# REACTION RATE PARAMETERS ----
combnames = function(par, reaction, target, reg){
  reg[grepl('^P',reg)] = paste(reg[grepl('^P',reg)],"A",sep = "_")
  comb = expand.grid(target, reg)
  if(nrow(comb)!=0){sapply(1:nrow(comb), function(i){paste0(par, reaction, comb[i,1], comb[i, 2])})}
}

params = c( # ----
            # Transcription parameters
            network$k_TC, # Basal transcription rates
            as.vector(network$TF_sgn), # regulation direction
            as.vector(network$TF_th*cohort$QTL_TC[,1]), # regulation threshold
            as.vector(network$TF_n), # regulation power
            as.vector(network$TF_fc), # regulation fold change
            # Translation parameters
            network$k_TL, # Basal translation rates
            as.vector(network$TLF_sgn), # regulation direction
            as.vector(network$TLF_th*cohort$QTL_TL[,1]), # regulation threshold
            as.vector(network$TLF_n), # regulation power
            as.vector(network$TLF_fc), # regulation fold change
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

names(params) = c( # ----
                   # Transcription parameters
                   sapply(network$genes, function(x){paste0("k_TC",x)}), # Basal transcription rates
                   combnames("sgn", "TF", rownames(network$TF_sgn), colnames(network$TF_sgn)), # regulation direction
                   combnames("th", "TF", rownames(network$TF_th), colnames(network$TF_th)), # regulation threshold
                   combnames("n", "TF", rownames(network$TF_n), colnames(network$TF_n)), # regulation power
                   combnames("fc", "TF", rownames(network$TF_fc), colnames(network$TF_fc)), # regulation fold change
                   # Translation parameters
                   sapply(network$prot, function(x){paste0("k_TL",x)}), # Basal translation rates
                   combnames("sgn", "TLF", rownames(network$TLF_sgn), colnames(network$TLF_sgn)), # regulation direction
                   combnames("th", "TLF", rownames(network$TLF_th), colnames(network$TLF_th)), # regulation threshold
                   combnames("n", "TLF", rownames(network$TLF_n), colnames(network$TLF_n)), # regulation power
                   combnames("fc", "TLF", rownames(network$TLF_fc), colnames(network$TLF_fc)), # regulation fold change
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

temp = c(G, P, M); names(temp) = c("G","P","M")
params = c(params, temp)

pfun = function(t, X, params){
  
  res = vector(length = (2*G + 5*P))
  
  
}







# Propensity functions ----

# Regulation law for transcription and translation reactions (takes into account the fc parameter)
regulation_law_TCTL = function(id, reaction, network){
  reg = colnames(network[[paste0(reaction,"_sgn")]])[which(network[[paste0(reaction,"_sgn")]][id,]!=0)]
  reg[grepl('^P',reg)] = paste(reg[grepl('^P',reg)],"A",sep = "_") # consider only active proteins
  parname = paste0(reaction, id)
  sup = sapply(reg, function(x){paste(c("*(1+as.numeric(params['sgn",parname,x,"'])*as.numeric(params['fc",parname,x,"'])*(X[,'",x,"']^as.numeric(params['n",parname,x,"'])/(X[,'",x,"']^as.numeric(params['n",parname,x,"'])+as.numeric(params['th",parname,x,"'])^as.numeric(params['n",parname,x,"']))))"), collapse = "")})
  return(sup)
}

# Regulation law for the decay reactions
regulation_law_decay = function(id, reaction, network){
  reg = colnames(network[[paste0(reaction,"_sgn")]])[which(network[[paste0(reaction,"_sgn")]][id,]!=0)]
  reg[grepl('^P',reg)] = paste(reg[grepl('^P',reg)],"A",sep = "_") # consider only active proteins
  parname = paste0(reaction, id)
  sup = sapply(reg, function(x){paste(c("*(1+as.numeric(params['sgn",parname,x,"'])*(X[,'",x,"']^as.numeric(params['n",parname,x,"'])/(X[,'",x,"']^as.numeric(params['n",parname,x,"'])+as.numeric(params['th",parname,x,"'])^as.numeric(params['n",parname,x,"'])))"), collapse = "")})
  return(sup)
}

# Regulation law for the activation/inactivation reactions
regulation_law_act = function(id, reaction, network){
  #reg = names(which(network[[paste0(reaction,"_sgn")]][id,]!=0))
  reg = colnames(network[[paste0(reaction,"_sgn")]])[which(network[[paste0(reaction,"_sgn")]][id,]!=0)]
  reg[grepl('^P',reg)] = paste(reg[grepl('^P',reg)],"A",sep = "_") # consider only active proteins
  parname = paste0(reaction, id)
  sup = sapply(reg, function(x){paste(c("*(X[,'",x,"']^as.numeric(params['n",parname,x,"'])/(X[,'",x,"']^as.numeric(params['n",parname,x,"'])+as.numeric(params['th",parname,x,"'])^as.numeric(params['n",parname,x,"'])))"), collapse = "")})
  return(sup)
}

a = c( # ----
       # transcription reactions
       sapply(network$genes, function(g){
         sup = regulation_law_TCTL(g, 'TF', network)
         paste(c("as.numeric(params['k_TC",g,"'])",sup), collapse = "")}),
       # translation reactions
       sapply(network$prot, function(p){           
         sup = regulation_law_TCTL(network$g2p[p], 'TLF', network)
         paste(c("X[,'",network$g2p[p],"']*as.numeric(params['k_TL",p,"'])",sup), collapse = "")}),
       # RNA decay reactions
       sapply(network$genes, function(g){
         sup = regulation_law_decay(g, 'DR', network)
         paste(c("X[,'",g,"']*as.numeric(params['p0_DR",g,"'])",sup), collapse = "")}),
       # Protein decay reactions (for active proteins then for inactive proteins)
       sapply(protNAA, function(p){
         sup = regulation_law_decay(sub("_NA|_A","",p), 'DP', network)
         paste(c("X[,'",p,"']*as.numeric(params['p0_DP",sub("_NA|_A","",p),"'])",sup), collapse = "")}),
       # Protein activation reactions
       sapply(protNAA[1:P], function(p){
         sup = regulation_law_act(sub("_NA","",p), 'ACT', network)
         paste(c("X[,'",p,"']",sup), collapse = "")}),
       # Protein inactivation reactions
       sapply(protNAA[(P+1):(2*P)], function(p){
         sup = regulation_law_act(sub("_A","",p), 'DEACT', network)
         if(length(sup) == 0){return("0")}
         else{paste(c("X[,'",p,"']",sup), collapse = "")}})
) # ----
names(a) = sapply(1:length(a), function(i){paste0("pro",i)})



params = c(a, params)


pfun(1, x0, params)

pfun = function(t, X, params){
  react = grep("pro", names(a), value = T)
  res = sapply(react, function(i){eval(parse(text = params[i]))})
  return(matrix(res, nrow = 1))
}




ssar::ssa(xinit = x0, pfun = pfun, v = v, params = params, tmax = 1000, nsim = 1)





#########################################################################################################################################################################################################################################
#########################################################################################################################################################################################################################################
#########################################################################################################################################################################################################################################


library(ssar)
source("simulation.R")

network = rand_network(5, 4, 2)
ind = rand_indiv(network)

G = length(network$genes)
P = length(network$prot)
M = length(network$met)

protNAA = c(sapply(network$prot, function(p){paste(p,"NA",sep = "_")}), sapply(network$prot, function(p){paste(p,"A",sep = "_")}))

# Initial conditions ----
x0 = matrix(c(ind$rna_0, ind$prot_NA_0, ind$prot_A_0, ind$met_tot_0), nrow = 1)
colnames(x0) = c(network$genes, protNAA, network$met)

# Propensity matrix ----
tempTCTL = diag(1, nrow = G+P, ncol = G+P); tempTCTL = rbind(tempTCTL, matrix(0, nrow = P, ncol = G+P))
tempDRDP = diag(-1, nrow = G+2*P, ncol = G+2*P)
tempPos = diag(1, nrow = P, ncol = P); tempNeg = diag(-1, nrow = P, ncol = P)
tempACT = rbind(matrix(0, nrow = G, ncol = P), tempNeg, tempPos)
tempDEACT = rbind(matrix(0, nrow = G, ncol = P), tempPos, tempNeg)
v = cbind(tempTCTL, tempDRDP, tempACT, tempDEACT); v = rbind(v, matrix(0, nrow = M, ncol = ncol(v))); rownames(v) = colnames(x0)

# Propensity scores ----

pfun = function(t, X, params){
  
  # G transcription rates, P translation rates, G RNA degradation rates, 2*P protein degradation rates, P activation rates, P inactivation rates
  res = vector(length = (2*params["G"] + 5*params["G"]))
  
  
  
  return(matrix(res, nrow = 1))
}

# Parameters ----

params = c(length(network$genes), length(network$prot), length(network$met))

names(params) = c("G", "P", "M")


#########################################################################################################################################################################################################################################

library(ssar)

params     <- c(iter = 3, a = 3, c = 2, b = 0.01, v1 = 1, v2 = 2)
X          <- matrix(c(100, 100), ncol = 2); colnames(X) = c("v1", "v2")

#Propensity function
pfun       <- function(t, X, params){ cbind(params['a']*t*X[,params['v1']] + 1, 
                                            params['b']*X[,params['v1']]*X[,params['v2']], 
                                            params['c']*X[,params['v2']]) }


pfun = function(t, X, params){
  res = vector(length = params['iter'])
  for(i in 1:params['iter']){
    res[i] = i*params["b"]*X[,params["v1"]]
  }
  return(matrix(res, nrow = 1))
}

#Propensity score
v          <- matrix(c(+1,-1,0,0,+1,-1),nrow=2,byrow=TRUE)

#Simulate
simulation <- ssar::ssa(X, pfun, v, params)




#########################################################################################################################################################################################################################################
#########################################################################################################################################################################################################################################
#########################################################################################################################################################################################################################################


library(adaptivetau)

init.values = c(prey = 1000, pred = 500)
transitions = matrix(c(+1, 0, -2, +1, 0, -1), nrow = 2, ncol = 3, dimnames = list(c("prey", "pred"), c()))
params = vector("list")
params$p = c(r = 10, beta = 0.01, delta = 10)
rateFunc = function(x, params, t){
  with(params, {  res = vector(length = 3)
  res[1] = p['r'] * x["prey"]
  res[2] = p['beta'] * x["prey"] * x["pred"] * (x["prey"] >=2)
  res[3] = p['delta'] * x["pred"]
  return(res)})
}

sim = ssa.adaptivetau(init.values, transitions, rateFunc, params, tf = 100)


# ----------------------------------


source("simulation.R")

network = rand_network(5, 4, 2)
ind = rand_indiv(network)

tmax = 1000


system.time(sim <- simuInd(network, ind, tmax))

system.time(simAD <- ssa.adaptivetau(init.values, transitions, rateFunc, params, tf = 1000))

parSSA = paramSSAindiv(network, ind)
system.time(simSSA <- GillespieSSA::ssa(parSSA$x0, parSSA$a, parSSA$nu, parSSA$parms, tmax, method = "OTL") )

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
system.time(for(s in 1:nsim){
    simAT <- ssa.adaptivetau(init.values, transitions, rateFunc, params, tf = 1000)
    for(g in network$genes){ res_simAT[[g]] = rbind(res_simAT[[g]], sapply(0:tmax, function(i){simAT[max(which(simAT[,"time"]<=i)),g]})) }
    for(p in network$prot){ res_simAT[[p]] = rbind(res_simAT[[p]], sapply(0:tmax, function(i){
      temptime =  max(which(simAT[,"time"]<=i)) 
      return(simAT[temptime,paste0(p,"_NA")] + simAT[temptime,paste0(p,"_A")])})) }
  }
)

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





#########################################################################################################################################################################################################################################
#                                                          TEST ENZYMATIC REACTION                                                             #
#########################################################################################################################################################################################################################################

library(adaptivetau)

tf = 10

##  Complete reaction ----

x0 = c(E = 100, S = 50, ES = 0, P = 0)

transitions = matrix(c(-1, -1, 1, 0, # binding reaction
                       1, 1, -1, 0, # unbinding reaction
                       1, 0, -1, 1), # transformation reaction
                     nrow = 4, dimnames = list(c("E", "S", "ES", "P"), c("binding", "unbinding", "transformation")))

# values from Lagagne, 2013
params = list(k1 = 1, # 1 molecule^(-1).s^(-1)
              kminus1 = 28300, #28300 s^(-1)
              k2 = 3.2) # k2 = kcat (median 13.7 s^(-1))

rateFunc = function(x, params, t){
  c(binding = params$k1*x["E"]*x["S"], unbinding = params$kminus1*x["ES"], transformation = params$k2*x["ES"])
}

rateFunc(x0, params, 1)


system.time(sim1 <- ssa.adaptivetau(init.values = x0, transitions = transitions, rateFunc = rateFunc, params = params, tf = tf))


## Michaelis-Menten

x02 = c(E = x0["E"], S = x0["S"], P = x0["P"]); names(x02) = c("E", "S", "P")

transitions2 = matrix(c(0, -1, 1), nrow = 3, dimnames = list(c("E", "S", "P"), c("transformation")))

params2 = list("kcat" = params$k2, "KM" = (params$kminus1 + params$k2)/params$k1)

rateFunc2 = function(x, params, t){
  c(params2$kcat * x["E"] * x["S"] / (params2$KM + x["S"]))
}
rateFunc2(x02, params2, 1)

system.time( sim2 <- ssa.adaptivetau(init.values = x02, transitions = transitions2, rateFunc = rateFunc2, params = params2, tf = tf) )

tplot = round(seq(1, nrow(sim1), by = nrow(sim1)/100))

plot(sim1[tplot, "time"], sim1[tplot, "S"], col = "red", type = 'l', ylim = c(0, max(x0)))
lines(sim2[,"time"], sim2[,"S"], col = "red", lty = 2)
lines(sim1[tplot, "time"], sim1[tplot, "P"], col = "green")
lines(sim2[,"time"], sim2[,"P"], col = "green", lty = 2)
lines(sim1[tplot, "time"], sim1[tplot, "E"], col = "orange")
lines(sim1[tplot, "time"], sim1[tplot, "ES"], col = "yellow")

for(m in names(x0)){
  plot(sim1[tplot, "time"], sim1[tplot, m], main = m)
}
