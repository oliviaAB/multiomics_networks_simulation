##########################################################################################################################
##########################################################################################################################
###                                  DATA SIMULATION - DRAFT 1 - 28/11/2017                                            ###
##########################################################################################################################
##########################################################################################################################
# Variables ----

# N individuals, n = 1 .. N
# G genes, g = 1 .. G
# P proteins, p = 1 .. P
# NC non-coding genes, NC = G - P, nc = 1 .. NC
# M metabolites, m = 1 .. M
# Q QTLs, q = 1 .. Q

# genes : vector of characters with gene ID
# prot : vector of characters with protein ID
# protcod : vector of characters with protein-coding genes ID
# noncod : vector of characters with non-coding genes ID
# met : vector of characters with metabolite ID
# ind : vactor of characters with individuals ID

# tmax = simulation time


# Graph  ----

## Topology

# TF : array for transcription regulation network, target genes (lenght of 1st dimension = G) x transcription regulators (lenght of 2nd dimension = NC + P (=G)) X edges values (lenght of 3rd dimension = 3)
# TLF : array for translation regulation network, target protein-coding genes (lenght of 1st dimension = P) x transcription regulators (lenght of 2nd dimension = NC + P (=G)) X edges values (lenght of 3rd dimension = 3)
# DR : array for RNA degradation regulation network, target genes (lenght of 1st dimension = G) X regulators non coding RNAs (lenght of 2nd dimension = NC) X edges values (lenght of 3rd dimension = 3)
# DP: array for protein degradation regulation network, target proteins (lenght of 1st dimension = P) X regulators proteins (lenght of 2nd dimension = P) X edges values (lenght of 3rd dimension = 3)
# ACT : array for protein activation network, target proteins (lenght of 1st dimension = P) x activators (lenght of 2nd dimension = P + M) X edges values (lenght of 3rd dimension = 3)
# DEACT : array for protein deactivation network, target proteins (lenght of 1st dimension = P) x deactivators (lenght of 2nd dimension = P + M) X edges values (lenght of 3rd dimension = 3)
# g2p : vector of correspondence between protein-coding genes ID and proteins ID, names = proteins ID, values = genes ID, (length P)

## Nodes parameters
# k_TC : vector of basal transcription rates for each gene, i.e. mean basal number of transcription events in per time unit (= in absence of any regulators), genes (length G)
# k_TL : vector of basal translation rates for each protein-coding gene, i.e. mean basal number of transcription events in per time unit (= in absence of any regulators), proteins (length P)
# p0_DR : vector of basal RNA degradation rate for each gene, genes (lenght G)
# p0_DP : vector of basal protein degradation rate for each gene, genes (length P)


## Genotype effects
# QTL_TC : matrix of genotype effect on transcription (TF binding) for each individual, effects (G rows) x individuals (N columns)
# QTL_TL : matrix of genotype effect on translation (TLF binding) for each individual, effects (P rows) x individuals (N columns)
# QTL_DR : matrix of genotype effect on mRNA degradation for each individual, effects (G rows) x individuals (N columns)
#   Remark : the effects can arise from non-linear combination of the effects of several QTLs, they are calculated prior to the simulation (as genotype does not vary with time)


# Dynamic quantities ----

## Rates and probabilities

# l_TC : matrix of transcription rates for each gene and each individual (i.e. mean number of transcription events in per time unit), individuals (N rows) x genes (G columns)
# l_TL : matrix of translation rates for each gene and each individual (i.e. mean number of translation events in per time unit), individuals (N rows) x protein-coding genes (P columns)
# p_DR : matrix of mRNA degradation rate for each gene and each individual, individuals (N rows) x genes (G columns)

## Molecular quantities

# rna : matrix of RNA abundance for each individual, genes (G rows) x individuals (N columns)
# prot_tot : matrix of protein abundance for each individual, proteins (P rows) x individuals (N columns)
# prot_A : matrix of active protein abundance for each individual, proteins (P rows) x individuals (N columns)
# prot_NA : matrix of inactive protein abundance for each individual, proteins (P rows) x individuals (N columns)


# ----


##########################################################################################################################
#                                                       MAIN CODE                                                        #
##########################################################################################################################


## Initialization ----

t = 1

# Vector creation
rna = matrix(NA, nrow = G, ncol = N); rownames(rna) = genes; colnames(rna) = ind
prot_tot = matrix(NA, nrow = P, ncol = N); rownames(prot_tot) = prot; colnames(prot_tot) = ind
prot_A = matrix(NA, nrow = P, ncol = N); rownames(prot_A) = prot; colnames(prot_A) = ind
prot_NA = matrix(NA, nrow = P, ncol = N); rownames(prot_NA) = prot; colnames(prot_NA) = ind
met_tot = matrix(NA, nrow = M, ncol = N); rownames(met_tot) = met; colnames(met_tot) = ind

l_TC = matrix(NA, nrow = G, ncol = N, dimnames = list(genes, ind))
l_TL = matrix(NA, nrow = P, ncol = N, dimnames = list(protcod, ind))
p_DR = matrix(NA, nrow = G, ncol = N, dimnames = list(genes, ind))
p_DP = matrix(NA, nrow = P, ncol = N, dimnames = list(prot, ind))
p_act = matrix(NA, nrow = P, ncol = N, dimnames = list(prot, ind))
p_deact = matrix(NA, nrow = P, ncol = N, dimnames = list(prot, ind))


# Initial values
# Are the initial values the same for all individuals???
rna_prev = matrix(rna_0[genes, ind], nrow = G, ncol = N); rownames(rna_prev) = genes; colnames(rna_prev) = ind
prot_A_prev = matrix(prot_A_0[prot, ind], nrow = P, ncol = N); rownames(prot_A_prev) = prot; colnames(prot_A_prev) = ind
prot_NA_prev = matrix(prot_NA_0[prot, ind], nrow = P, ncol = N); rownames(prot_NA_prev) = prot; colnames(prot_NA_prev) = ind
prot_tot_prev = matrix(prot_A_prev[prot, ind] + prot_NA_prev[prot, ind], nrow = P, ncol = N); rownames(prot_tot_prev) = prot; colnames(prot_tot_prev) = ind
met_tot_prev = matrix(met_tot_0[met, ind], nrow = M, ncol = N); rownames(met_tot_prev) = met; colnames(met_tot_prev) = ind


# For visualization

time_rna = vector("list", G); names(time_rna) = genes
time_prot_tot = vector("list", P); names(time_prot_tot) = prot
time_prot_A = vector("list", P); names(time_prot_A) = prot
time_prot_NA = vector("list", P); names(time_prot_NA) = prot
time_met_tot = vector("list", M); names(time_met_tot) = met

# Add t = 0
for(g in genes){ time_rna[[g]] = cbind(time_rna[[g]],rna_prev[g,]) }
for(p in prot){ 
  time_prot_tot[[p]] = cbind(time_prot_tot[[p]], prot_tot_prev[p,])
  time_prot_A[[p]] = cbind(time_prot_A[[p]], prot_A_prev[p,])
  time_prot_NA[[p]] = cbind(time_prot_NA[[p]], prot_NA_prev[p,])
}
for(m in met){ time_met_tot[[m]] = cbind(time_met_tot[[m]], met_tot_prev[m,]) }

# Main loop -----

while(t<tmax){
  
  mol_prev = rbind(rna_prev[noncod,ind], prot_A_prev[,ind], met_tot_prev[,ind]); rownames(mol_prev) = c(noncod,prot,met)
  
  ## 1: Update rates and probabilites ----
  
  # Transcription rate
  l_TC[genes, ind] = k_TC[genes] * t(sapply(genes, function(g){
    #reg = names(which(TF[g,,'sgn']!=0))
    reg = names(which(setNames(TF[g,,'sgn'], colnames(TF))!=0))
    if(length(reg) == 0) return(matrix(1,ncol=N,dimnames = list(NULL, ind)))
    r_A = mol_prev[reg,]
    sgn_r = TF[g,reg,'sgn']
    th_r = t(t(matrix(TF[g,reg,'th'], nrow = length(reg), ncol = N, dimnames = list(reg,ind))) * QTL_TC[g,])
    n_r = TF[g,reg,'n']
    temp = matrix(1+sgn_r*(r_A^n_r)/(th_r^n_r+r_A^n_r), ncol=N); colnames(temp) = ind 
    res = apply(temp, 2, prod)
    
    return(res)
  }))
  
  print(l_TC['G1',])
  
  # Translation rate 
  l_TL[protcod, ind] = rna_prev[protcod,] * k_TL[protcod] * t(sapply(protcod, function(g){
    reg = names(which(TLF[g,,'sgn']!=0))
    if(length(reg) == 0) return(matrix(1,ncol=N,dimnames = list(NULL, ind)))
    r_A = mol_prev[reg,]
    sgn_r = TLF[g,reg,'sgn']
    th_r = t(t(matrix(TLF[g,reg,'th'], nrow = length(reg), ncol = N, dimnames = list(reg,ind))) * QTL_TL[g,])
    n_r = TLF[g,reg,'n']
    temp = matrix(1+sgn_r*(r_A^n_r)/(th_r^n_r+r_A^n_r), ncol=N); colnames(temp) = ind 
    res = apply(temp, 2, prod)
    
    return(res)
  }))
  
  
  # RNA degradation rate
  p_DR[genes, ind] = p0_DR[genes] * QTL_DR[genes, ind] * t(sapply(genes, function(g){
    reg = names(which(setNames(DR[g,,'sgn'], colnames(DR))!=0))
    if(length(reg) == 0) return(matrix(1,ncol=N,dimnames = list(NULL, ind)))
    r_A = mol_prev[reg,]
    th_r = matrix(DR[g,reg,'th'], nrow = length(reg), ncol = N, dimnames = list(reg,ind))
    n_r = DR[g,reg,'n']
    temp = matrix(1+(r_A^n_r)/(th_r^n_r+r_A^n_r), ncol=N); colnames(temp) = ind 
    res = apply(temp, 2, prod)
    
    return(res)
  }))
  
  
  # protein degradation rate
  p_DP[prot, ind] = p0_DP[prot] * t(sapply(prot, function(g){
    reg = names(which(setNames(DP[g,,'sgn'], colnames(DP))!=0))
    if(length(reg) == 0) return(matrix(1,ncol=N,dimnames = list(NULL, ind)))
    r_A = mol_prev[reg,]
    th_r = matrix(DP[g,reg,'th'], nrow = length(reg), ncol = N, dimnames = list(reg,ind))
    n_r = DP[g,reg,'n']
    temp = matrix(1+(r_A^n_r)/(th_r^n_r+r_A^n_r), ncol=N); colnames(temp) = ind 
    res = apply(temp, 2, prod)
    
    return(res)
  }))
  
  # protein activation rate
  p_act[prot, ind] = t(sapply(prot, function(g){
    reg = names(which(setNames(ACT[g,,'sgn'], colnames(ACT))!=0))
    if(length(reg) == 0) return(matrix(1,ncol=N,dimnames = list(NULL, ind)))
    r_A = mol_prev[reg,]
    th_r = matrix(ACT[g,reg,'th'], nrow = length(reg), ncol = N, dimnames = list(reg,ind))
    n_r = ACT[g,reg,'n']
    temp = matrix(1+(r_A^n_r)/(th_r^n_r+r_A^n_r), ncol=N); colnames(temp) = ind 
    res = apply(temp, 2, prod)
    
    return(res)
  }))
  
  # protein deactivation rate
  p_deact[prot, ind] = t(sapply(prot, function(g){
    reg = names(which(setNames(DEACT[g,,'sgn'], colnames(DEACT))!=0))
    if(length(reg) == 0) return(matrix(0,ncol=N,dimnames = list(NULL, ind)))
    r_A = mol_prev[reg,]
    th_r = matrix(DEACT[g,reg,'th'], nrow = length(reg), ncol = N, dimnames = list(reg,ind))
    n_r = DEACT[g,reg,'n']
    temp = matrix(1+(r_A^n_r)/(th_r^n_r+r_A^n_r), ncol=N); colnames(temp) = ind 
    res = apply(temp, 2, prod)
    
    return(res)
  }))
  
  
  
  
  ## 2: Compute new molecule abundances ----
  
  # Update RNA abundance
  #    Gene transcription follows a Poisson law, and each RNA molecule at time t-1 has a probability p_DR of being degraded
  rna[genes, ind] = rna_prev[genes, ind] + rpois(length(l_TC), l_TC) - rbinom(N*G, rna_prev, p_DR)
  
  
  # Update protein abundance 
  
  # Activation or degradation of inactive proteins
  # Each inactive protein has 3 possible fates:
  #      - Stay inactive with a probability of (1-proba_activation)*(1-proba_degradation)
  #      - Become active with a probability of proba_activation*(1-proba_degradation)
  #      - Be degraded with a probability of proba_degradation
  probs = rbind( as.vector((1-p_act)*(1-p_DP)), as.vector(p_act*(1-p_DP)) ,as.vector(p_DP) ) # order of the columns: Prot1 for ind1, Prot2 for ind1, ... , ProtP for ind1, Prot1 for ind2, etc
  # For each protein type in each individuals these fates are sampled using a multinomial distribution
  temp = sapply(1:(N*P), function(i){rmultinom(1,prot_NA_prev[i],probs[,i])})
  Xa = temp[2,] # proteins inactive at time t-1 that are activated
  Xdna = temp[3,] # proteins inactive at time t-1 that are degraded
  
  # Each active protein has 3 possible fates:
  #      - Stay active with a probability of (1-proba_deactivation)*(1-proba_degradation)
  #      - Become inactive with a probability of proba_deactivation*(1-proba_degradation)
  #      - Be degraded with a probability of proba_degradation
  probs = rbind( as.vector((1-p_deact)*(1-p_DP)), as.vector(p_deact*(1-p_DP)) ,as.vector(p_DP) ) # order of the columns: Prot1 for ind1, Prot2 for ind1, ... , ProtP for ind1, Prot1 for ind2, etc
  # For each protein type in each individuals these fates are sampled using a multinomial distribution
  temp = sapply(1:(N*P), function(i){rmultinom(1,prot_A_prev[i],probs[,i])})
  Xna = temp[2,] # proteins active at time t-1 that are deactivated
  Xda = temp[3,] # proteins active at time t-1 that are degraded
  
  
  # Update inactive protein abundance
  prot_NA[prot, ind] = rpois(length(l_TL), l_TL) + prot_NA_prev[prot, ind] - Xa - Xdna + Xna
  
  # Update active protein abundance
  prot_A[prot, ind] = prot_A_prev[prot, ind] - Xna - Xda + Xa
  
  # Update total protein abundance
  prot_tot[prot, ind] = prot_NA[prot, ind] + prot_A[prot, ind]
  
  
  
  ## Update the '_prev' matrices ----
  rna_prev = rna
  prot_tot_prev = prot_tot
  prot_NA_prev = prot_NA
  prot_A_prev = prot_A
  
  ## Store values for visualization
  for(g in genes){ time_rna[[g]] = cbind(time_rna[[g]],rna[g,]) }
  for(p in prot){ 
    time_prot_tot[[p]] = cbind(time_prot_tot[[p]], prot_tot[p,])
    time_prot_A[[p]] = cbind(time_prot_A[[p]], prot_A[p,])
    time_prot_NA[[p]] = cbind(time_prot_NA[[p]], prot_NA[p,])
    }
  for(m in met){ time_met_tot[[m]] = cbind(time_met_tot[[m]], met_tot[m,]) }
  
  t = t + 1
}

##########################################################################################################################
#                                                       OUTPUT                                                           #
##########################################################################################################################

## Visualization ----

cols = c('blue','green', 'red'); names(cols) = ind

# plot genes
for(g in genes){
  plot(1:(tmax+1), ylim = c(0,max(unlist(time_rna[[g]]))), type = 'n', main = paste('Gene',g, sep=' '), xlab = 'time', ylab = 'abundance')
  for(i in ind){
    lines(time_rna[[g]][i,], col=cols[i])
  }
  legend('topright', col = cols, lty=1, legend = names(cols))
}

# plot genes
for(p in prot){
  plot(1:(tmax+1), ylim = c(0,max(unlist(time_prot_tot[[p]]))), type = 'n', main = paste('Protein',p, sep=' '), xlab = 'time', ylab = 'abundance')
  for(i in ind){
    lines(time_prot_tot[[p]][i,], col=cols[i])
  }
  legend('topright', col = cols, lty=1, legend = names(cols))
}