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

source("param.R")

if(!suppressWarnings(require("GillespieSSA", quietly = T))){install.packages("GillespieSSA")}
library(GillespieSSA)
if(!suppressWarnings(require("deSolve", quietly = T))){install.packages("deSolve")}
library(deSolve)
if(!suppressWarnings(require("igraph", quietly = T))){install.packages("igraph")}
library(igraph)
if(!suppressWarnings(require("scales", quietly = T))){install.packages("scales")}
library(scales)
# if(!suppressWarnings(require("devtools", quietly = T))){install.packages("devtools")}
# if(!suppressWarnings(require("ssar", quietly = T))){devtools::install_github("INSP-RH/ssar")}
# library(ssar)
if(!suppressWarnings(require("adaptivetau", quietly = T))){install.packages("adaptivetau")}
library(adaptivetau)
if(!suppressWarnings(require("magic", quietly = T))){install.packages("magic")}
library(magic)


##########################################################################################################################
#                                                  NETWORK SIMULATION                                                    #
##########################################################################################################################

rand_network = function(G, P, M, MR = NULL){
 
  # Nodes
  if(G!=0){
    genes = sapply(1:G, function(i){ return(paste0("RNA",i)) })
  } else{ genes = vector()}
  if(P!=0){
    prot = sapply(1:P, function(i){ return(paste0("P",i)) })
  } else{ prot = vector()}
  if(M!=0){
    met = sapply(1:M, function(i){ return(paste0("M",i)) })
  } else{ met = vector()}

  g2p = genes[1:P]; names(g2p) = prot
  protcod = unname(g2p); noncod = setdiff(genes, protcod); NC = G-P
  
  ## Topology
  
  # Transcription regulation ----
  # TF : array for transcription regulation network, target genes (G rows) x transcription regulators (NC + P (=G) columns)
  TF_sgn = matrix( sample(c(-1:1), size = (G*G), replace = T, prob = proba_reg_TF), nrow = G, ncol = G, dimnames = list(genes, c(noncod,prot)))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  TF_th = matrix( 0, nrow = G, ncol = G, dimnames = list(genes, c(noncod,prot))) ; TF_th[which(TF_sgn!=0)] = get(th_sampling)(length(which(TF_sgn!=0)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  TF_n = matrix( 0, nrow = G, ncol = G, dimnames = list(genes, c(noncod,prot))) ; TF_n[which(TF_sgn!=0)] = get(n_sampling)(length(which(TF_sgn!=0)))
  # The fold change of target expression induced by maximal activation !! We only consider this parameter for activators of transcription (and not repressors, for which fc = 1)
  TF_fc = matrix( 0, nrow = G, ncol = G, dimnames = list(genes, c(noncod,prot))) ; TF_fc[which(TF_sgn>0)] = get(fc_sampling)(length(which(TF_sgn>0))); TF_fc[which(TF_sgn<0)] = 1
  
  # Translation regulation ----
  # TLF : array for translation regulation network, target protein-coding genes (P rows) x transcription regulators (NC + P (=G) columns)
  TLF_sgn = matrix( sample(c(-1:1), size = (P*G), replace = T, prob = proba_reg_TLF), nrow = P, ncol = G, dimnames = list(protcod, c(noncod,prot)))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  TLF_th = matrix( 0, nrow = P, ncol = G, dimnames = list(protcod, c(noncod,prot))) ; TLF_th[which(TLF_sgn!=0)] = get(th_sampling)(length(which(TLF_sgn!=0)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  TLF_n = matrix( 0, nrow = P, ncol = G, dimnames = list(protcod, c(noncod,prot))) ; TLF_n[which(TLF_sgn!=0)] = get(n_sampling)(length(which(TLF_sgn!=0)))
  # The fold change of target translation rate induced by maximal activation !! We only consider this parameter for activators of translation (and not repressors, for which fc = 1)
  TLF_fc = matrix( 0, nrow = P, ncol = G, dimnames = list(protcod, c(noncod,prot))) ; TLF_fc[which(TLF_sgn>0)] = get(fc_sampling)(length(which(TLF_sgn>0))); TLF_fc[which(TLF_sgn<0)] = 1
  
  
  # RNA decay ----
  # DR : array for RNA degradation regulation network, target genes (lenght of 1st dimension = G) X regulators non coding RNAs (NC+Q columns)
  DR_sgn = matrix( sample(c(-1:1), size = (G*NC), replace = T, prob = proba_reg_DR), nrow = G, ncol = NC, dimnames = list(genes, noncod))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  DR_th = matrix( 0, nrow = G, ncol = NC, dimnames = list(genes, noncod)) ; DR_th[which(DR_sgn!=0)] = get(th_sampling)(length(which(DR_sgn!=0)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  DR_n = matrix( 0, nrow = G, ncol = NC, dimnames = list(genes, noncod)) ; DR_n[which(DR_sgn!=0)] = get(n_sampling)(length(which(DR_sgn!=0)))
  
  
  # Protein decay ----
  # DP : array for RNA degradation regulation network, target genes (lenght of 1st dimension = G) X regulators non coding RNAs (NC+Q columns)
  DP_sgn = matrix( sample(c(-1:1), size = (P*P), replace = T, prob = proba_reg_DP), nrow = P, ncol = P, dimnames = list(prot, prot))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  DP_th = matrix( 0, nrow = P, ncol = P, dimnames = list(prot, prot)) ; DP_th[which(DP_sgn!=0)] = get(th_sampling)(length(which(DP_sgn!=0)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  DP_n = matrix( 0, nrow = P, ncol = P, dimnames = list(prot, prot)) ; DP_n[which(DP_sgn!=0)] = get(n_sampling)(length(which(DP_sgn!=0)))
  

  # Protein activation ----
  # ACT : array for protein activation network, target proteins (P rows) x activators (P + NC + M (= G+M) columns)
  ACT_sgn = matrix( sample(c(0,1), size = (P*(G+M)), replace = T, prob = proba_reg_ACT), nrow = P, ncol = NC+P+M, dimnames = list(prot, c(prot,noncod,met)))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  ACT_th = matrix( 0, nrow = P, ncol = G+M, dimnames = list(prot, c(prot,noncod,met))) ; ACT_th[which(ACT_sgn!=0)] = get(th_sampling)(length(which(ACT_sgn!=0)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  ACT_n = matrix( 0, nrow = P, ncol = G+M, dimnames = list(prot, c(prot,noncod,met))) ; ACT_n[which(ACT_sgn!=0)] = get(n_sampling)(length(which(ACT_sgn!=0)))
  
  # Protein deactivation ----
  # DEACT : array for protein deactivation network, target proteins (P rows) x activators (P + NC + M (= G+M) columns)
  DEACT_sgn = matrix( sample(c(0,1), size = (P*(G+M)), replace = T, prob = proba_reg_ACT), nrow = P, ncol = NC+P+M, dimnames = list(prot, c(prot,noncod,met)))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  DEACT_th = matrix( 0, nrow = P, ncol = G+M, dimnames = list(prot, c(prot,noncod,met))) ; DEACT_th[which(DEACT_sgn!=0)] = get(th_sampling)(length(which(DEACT_sgn!=0)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  DEACT_n = matrix( 0, nrow = P, ncol = G+M, dimnames = list(prot, c(prot,noncod,met))) ; DEACT_n[which(DEACT_sgn!=0)] = get(n_sampling)(length(which(DEACT_sgn!=0)))
  
  # Metabolic reactions ----
  #MR: number of metabolic reactions. Can be specified by the user, otherwise MR = M-1
  if(is.null(MR)) MR = M-1
  MR = min(MR, M*(M-1)) # there is maximum M(M-1) different possible reactions (given one reaction = Mi -> Mj)
  
  if(MR>0){ # If there is at least one metabolic reaction
    
    # Name of reactions
    metreactions = sapply(1:MR, function(i){ return(paste0("metreaction",i)) })
    
    # Compute all possible reactions (there is M(M-1) different possible reactions between M metabolites) => avoids to have twice the same reaction in the stoichiometry matrix
    temp = cbind(combn(met, 2),combn(met, 2)[2:1,])
    possreactions = matrix(0, nrow = M, ncol = M*(M-1)); rownames(possreactions) = met
    for(i in 1:ncol(temp)){
      possreactions[temp[1,i],i] = -1
      possreactions[temp[2,i],i] = 1
    }
    rm(temp)
      
    S = matrix(possreactions[,sample(ncol(possreactions),MR)], nrow = M, ncol = MR) # for each reaction randomly sample one column of the possible reactions matrix
    rownames(S) = met; colnames(S) = metreactions
    
    # Enzyme activity (only one enzyme can catalyze a given reaction, but an enzyme can catalyze several reactions)
    E2R = sample(prot, size = MR, replace = T); names(E2R) =  metreactions
    
    # Enzymatic rates
    ENZ = matrix(get(enzparam_sampling)(MR), nrow = 2, ncol = MR, dimnames = list(c("k_cat","K_M"), metreactions))
    
    
    # # Enzyme activity (only one enzyme can catalyze a given reaction)
    # if(P>=MR){ E2R = sample(prot, size = MR, replace = F); names(E2R) =  metreactions; spontreactions = vector()} # if more proteins than reaction, enzyme for each reaction is randomly chosen from the proteins
    # if(P< MR){ E2R = prot; names(E2R) = metreactions[1:P]; spontreactions = setdiff(metreactions, names(E2R)) } # if more reactions than proteins, the first P reactions are catalyzed by the proteins; the other reactions are spontaneous
    # # Enzymatic rates
    # ENZ = matrix(get(enzparam_sampling)(length(E2R)), nrow = 2, ncol = MR, dimnames = list(c("k_cat","K_M"), names(E2R)))
    # # Spontaneous reaction rates
    # if(length(spontreactions)>0){
    #   SPON = matrix(get(spontmetrate_sampling)(length(spontreactions)), nrow = 1, dimnames = list(c("rates"), spontreactions))
    # }else{SPON = matrix(nrow = 0, ncol = 0)}
    
  }else{MR = 0; S = matrix(nrow = 0, ncol = 0); metreactions = vector(); E2R = vector(); ENZ = matrix(nrow = 0, ncol = 0)}


  
  ## Nodes parameters  ----
  
  # Transcription rates
  # TC rates chosen randomly between 0.01 and 0.1
  k_TC = get(basal_transcription_rate)(G) ; names(k_TC) = genes
  
  # Translation rates
  # TL rates chosen randomly between 0.5 and 5
  k_TL = get(basal_translation_rate)(P); names(k_TL) = protcod
  
  # RNA decay rates
  # RNA decay rates chosen randomly between 0.005 and 0.01
  p0_DR = get(basal_RNAdecay_rate)(G); names(p0_DR) = genes
  
  # Protein decay rates
  # protein decay rates chosen randomly between 0.01 and 0.1
  p0_DP = get(basal_proteindecay_rate)(P); names(p0_DP) = prot

      
  res = list("genes" = genes, # ----
             "protcod" = protcod,
             "noncod" = noncod,
             "prot" = prot,
             "met" = met,
             "metreactions" = metreactions,
             "g2p" = g2p,
             "TF_sgn" = TF_sgn,
             "TF_th" = TF_th,
             "TF_n" = TF_n,
             "TF_fc" = TF_fc,
             "TLF_sgn" = TLF_sgn,
             "TLF_th" = TLF_th,
             "TLF_n" = TLF_n,
             "TLF_fc" = TLF_fc,
             "DR_sgn" = DR_sgn,
             "DR_th" = DR_th,
             "DR_n" = DR_n,
             "DP_sgn" = DP_sgn,
             "DP_th" = DP_th,
             "DP_n" = DP_n,
             "ACT_sgn" = ACT_sgn,
             "ACT_th" = ACT_th,
             "ACT_n" = ACT_n,
             "DEACT_sgn" = DEACT_sgn,
             "DEACT_th" = DEACT_th,
             "DEACT_n" = DEACT_n,
             "k_TC" = k_TC,
             "k_TL" = k_TL,
             "p0_DR" = p0_DR,
             "p0_DP" = p0_DP,
             "S" = S,
             "E2R" = E2R,
             "ENZ" = ENZ
             )
    # ----
  
  return(res)
}

rand_network_null = function(G, P, M, MR = NULL){
  
  # Nodes
  if(G!=0){
    genes = sapply(1:G, function(i){ return(paste0("RNA",i)) })
  } else{ genes = vector()}
  if(P!=0){
    prot = sapply(1:P, function(i){ return(paste0("P",i)) })
  } else{ prot = vector()}
  if(M!=0){
    met = sapply(1:M, function(i){ return(paste0("M",i)) })
  } else{ met = vector()}  
  
  g2p = genes[1:P]; names(g2p) = prot
  protcod = unname(g2p); noncod = setdiff(genes, protcod); NC = G-P
  
  ## Topology
  
  # Transcription regulation ----
  # TF : array for transcription regulation network, target genes (G rows) x transcription regulators (NC + P (=G) columns)
  TF_sgn = matrix( 0, nrow = G, ncol = G, dimnames = list(genes, c(noncod,prot)))
  TF_th = matrix( 0, nrow = G, ncol = G, dimnames = list(genes, c(noncod,prot)))
  TF_n = matrix( 0, nrow = G, ncol = G, dimnames = list(genes, c(noncod,prot)))
  TF_fc = matrix( 0, nrow = G, ncol = G, dimnames = list(genes, c(noncod,prot)))
  
  # Translation regulation ----
  # TLF : array for translation regulation network, target protein-coding genes (P rows) x transcription regulators (NC + P (=G) columns)
  TLF_sgn = matrix( 0, nrow = P, ncol = G, dimnames = list(protcod, c(noncod,prot)))
  TLF_th = matrix( 0, nrow = P, ncol = G, dimnames = list(protcod, c(noncod,prot)))
  TLF_n = matrix( 0, nrow = P, ncol = G, dimnames = list(protcod, c(noncod,prot)))
  TLF_fc = matrix( 0, nrow = P, ncol = G, dimnames = list(protcod, c(noncod,prot)))
  
  # RNA decay ----
  # DR : array for RNA degradation regulation network, target genes (lenght of 1st dimension = G) X regulators non coding RNAs (NC+Q columns)
  DR_sgn = matrix( 0, nrow = G, ncol = NC, dimnames = list(genes, noncod))
  DR_th = matrix( 0, nrow = G, ncol = NC, dimnames = list(genes, noncod))
  DR_n = matrix( 0, nrow = G, ncol = NC, dimnames = list(genes, noncod))
  
  
  # Protein decay ----
  # DP : array for RNA degradation regulation network, target genes (lenght of 1st dimension = G) X regulators non coding RNAs (NC+Q columns)
  DP_sgn = matrix( 0, nrow = P, ncol = P, dimnames = list(prot, prot))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  DP_th = matrix( 0, nrow = P, ncol = P, dimnames = list(prot, prot))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  DP_n = matrix( 0, nrow = P, ncol = P, dimnames = list(prot, prot))
  
  
  # Protein activation ----
  # ACT : array for protein activation network, target proteins (P rows) x activators (P + NC + M (= G+M) columns)
  ACT_sgn = matrix( 0, nrow = P, ncol = NC+P+M, dimnames = list(prot, c(prot,noncod,met)))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  ACT_th = matrix( 0, nrow = P, ncol = G+M, dimnames = list(prot, c(prot,noncod,met)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  ACT_n = matrix( 0, nrow = P, ncol = G+M, dimnames = list(prot, c(prot,noncod,met)))
  
  # Protein deactivation ----
  # DEACT : array for protein deactivation network, target proteins (P rows) x activators (P + NC + M (= G+M) columns)
  DEACT_sgn = matrix( 0, nrow = P, ncol = NC+P+M, dimnames = list(prot, c(prot,noncod,met)))
  DEACT_th = matrix( 0, nrow = P, ncol = G+M, dimnames = list(prot, c(prot,noncod,met)))
  DEACT_n = matrix( 0, nrow = P, ncol = G+M, dimnames = list(prot, c(prot,noncod,met)))
  
  # Metabolic reactions  ----
  #MR: number of metabolic reactions. Can be specified by the user, otherwise MR = M-1
  if(is.null(MR)) MR = M-1
  MR = min(MR, M*(M-1)) # there is maximum M(M-1) different possible reactions (given one reaction = Mi -> Mj)
  
  if(MR>0){ # If there is at least one metabolic reaction
    
    # Name of reactions
    metreactions = sapply(1:MR, function(i){ return(paste0("metreaction",i)) })
    
    # Compute all possible reactions (there is M(M-1) different possible reactions between M metabolites) => avoids to have twice the same reaction in the stoichiometry matrix
    temp = cbind(combn(met, 2),combn(met, 2)[2:1,])
    possreactions = matrix(0, nrow = M, ncol = M*(M-1)); rownames(possreactions) = met
    for(i in 1:ncol(temp)){
      possreactions[temp[1,i],i] = -1
      possreactions[temp[2,i],i] = 1
    }
    rm(temp)
    
    S = matrix(possreactions[,sample(ncol(possreactions),MR)], nrow = M, ncol = MR) # for each reaction randomly sample one column of the possible reactions matrix
    rownames(S) = met; colnames(S) = metreactions
    
    # Enzyme activity (only one enzyme can catalyze a given reaction, but an enzyme can catalyze several reactions)
    E2R = sample(prot, size = MR, replace = T); names(E2R) =  metreactions
    
    # Enzymatic rates
    ENZ = matrix(get(enzparam_sampling)(MR), nrow = 2, ncol = MR, dimnames = list(c("k_cat","K_M"), metreactions))
    
  }else{MR = 0; S = matrix(nrow = 0, ncol = 0); metreactions = vector(); E2R = vector(); ENZ = matrix(nrow = 0, ncol = 0)}
  
  
  ## Nodes parameters
  
  # Transcription rates
  # TC rates chosen randomly between 0.01 and 0.1
  k_TC = get(basal_transcription_rate)(G) ; names(k_TC) = genes
  
  # Translation rates
  # TL rates chosen randomly between 0.5 and 5
  k_TL = get(basal_translation_rate)(P); names(k_TL) = protcod
  
  # RNA decay rates
  # RNA decay rates chosen randomly between 0.005 and 0.01
  p0_DR = get(basal_RNAdecay_rate)(G); names(p0_DR) = genes
  
  # Protein decay rates
  # protein decay rates chosen randomly between 0.01 and 0.1
  p0_DP = get(basal_proteindecay_rate)(P); names(p0_DP) = prot
  
  res = list("genes" = genes, # ----
             "protcod" = protcod,
             "noncod" = noncod,
             "prot" = prot,
             "met" = met,
             "metreactions" = metreactions,
             "g2p" = g2p,
             "TF_sgn" = TF_sgn,
             "TF_th" = TF_th,
             "TF_n" = TF_n,
             "TF_fc" = TF_fc,
             "TLF_sgn" = TLF_sgn,
             "TLF_th" = TLF_th,
             "TLF_n" = TLF_n,
             "TLF_fc" = TLF_fc,
             "DR_sgn" = DR_sgn,
             "DR_th" = DR_th,
             "DR_n" = DR_n,
             "DP_sgn" = DP_sgn,
             "DP_th" = DP_th,
             "DP_n" = DP_n,
             "ACT_sgn" = ACT_sgn,
             "ACT_th" = ACT_th,
             "ACT_n" = ACT_n,
             "DEACT_sgn" = DEACT_sgn,
             "DEACT_th" = DEACT_th,
             "DEACT_n" = DEACT_n,
             "k_TC" = k_TC,
             "k_TL" = k_TL,
             "p0_DR" = p0_DR,
             "p0_DP" = p0_DP,
             "S" = S,
             "E2R" = E2R,
             "ENZ" = ENZ
  )
  # ----
  
  return(res)
}


##########################################################################################################################
#                                                INDIVIDUALS SIMULATION                                                  #
##########################################################################################################################

rand_cohort = function(network, N){
# 
#   G = length(network$genes)
#   P = length(network$prot)
#   M = length(network$met)
#   
#   # Individuals
#   ind = sapply(1:N, function(i){ return(paste0("ind",i)) })
# 
#   # Genotype effects
#   
#   # QTL_TC : matrix of genotype effect on transcription (TF binding) for each individual, effects (G rows) x individuals (N columns)
#   QTL_TC = matrix(get(qtl_effect_transcription)(G*N), nrow = G, ncol = N); rownames(QTL_TC) = network$genes; colnames(QTL_TC) = ind
# 
#   # QTL_TL : matrix of genotype effect on translation (TLF binding) for each individual, effects (P rows) x individuals (N columns)
#   QTL_TL = matrix(get(qtl_effect_translation)(P*N), nrow = P, ncol = N); rownames(QTL_TL) = network$protcod; colnames(QTL_TL) = ind
#   
#   # QTL_DR : matrix of genotype effect on mRNA degradation for each individual, effects (G rows) x individuals (N columns)
#   QTL_DR = matrix(get(qtl_effect_RNAdecay)(G*N), nrow = G, ncol = N); rownames(QTL_DR) = network$genes; colnames(QTL_DR) = ind
#   
#   if(G!=0){
#     rna_0 = matrix(get(initial_abundanceCohort)(G,N), nrow = G, ncol = N); rownames(rna_0) = network$genes; colnames(rna_0) = ind
#   }else{rna_0 = vector()}
#   if(P!=0){
#     prot_A_0 = matrix(get(initial_abundanceCohort)(P,N), nrow = P, ncol = N); rownames(prot_A_0) = network$prot; colnames(prot_A_0) = ind
#     prot_NA_0 = matrix(get(initial_abundanceCohort)(P,N), nrow = P, ncol = N); rownames(prot_NA_0) = network$prot; colnames(prot_NA_0) = ind
#     prot_tot_0 = prot_NA_0 + prot_A_0; rownames(prot_tot_0) = network$prot; colnames(prot_tot_0) = ind
#   }else{prot_A_0 = vector(); prot_NA_0 = vector(); prot_tot_0 = vector()}
#   if(M!=0){
#     met_tot_0 = matrix(get(initial_abundanceCohort)(M,N), nrow = M, ncol = N); rownames(met_tot_0) = network$met; colnames(met_tot_0) = ind
#   }else{met_tot_0 = vector()}
#   
#   res = list("ind" = ind,
#              "QTL_TC" = QTL_TC,
#              "QTL_TL" = QTL_TL,
#              "QTL_DR" = QTL_DR,
#              "rna_0" = rna_0,
#              "prot_tot_0" = prot_tot_0,
#              "prot_A_0" = prot_A_0,
#              "prot_NA_0" = prot_NA_0,
#              "met_tot_0" = met_tot_0)
#   return(res)
} 

rand_cohort_null = function(network, N){
#   
#   G = length(network$genes)
#   P = length(network$prot)
#   M = length(network$met)
#   
#   # Individuals
#   ind = sapply(1:N, function(i){ return(paste0("ind",i)) })
#   
#   # Genotype effects
#   
#   # QTL_TC : matrix of genotype effect on transcription (TF binding) for each individual, effects (G rows) x individuals (N columns)
#   QTL_TC = matrix(1, nrow = G, ncol = N); rownames(QTL_TC) = network$genes; colnames(QTL_TC) = ind
#   
#   # QTL_TL : matrix of genotype effect on translation (TLF binding) for each individual, effects (P rows) x individuals (N columns)
#   QTL_TL = matrix(1, nrow = P, ncol = N); rownames(QTL_TL) = network$protcod; colnames(QTL_TL) = ind
#   
#   # QTL_DR : matrix of genotype effect on mRNA degradation for each individual, effects (G rows) x individuals (N columns)
#   QTL_DR = matrix(1, nrow = G, ncol = N); rownames(QTL_DR) = network$genes; colnames(QTL_DR) = ind
#   
#   if(G!=0){
#     rna_0 = matrix(0, nrow = G, ncol = N); rownames(rna_0) = network$genes; colnames(rna_0) = ind
#   }else{rna_0 = vector()}
#   if(P!=0){
#     prot_A_0 = matrix(0, nrow = P, ncol = N); rownames(prot_A_0) = network$prot; colnames(prot_A_0) = ind
#     prot_NA_0 = matrix(0, nrow = P, ncol = N); rownames(prot_NA_0) = network$prot; colnames(prot_NA_0) = ind
#     prot_tot_0 = prot_NA_0 + prot_A_0; rownames(prot_tot_0) = network$prot; colnames(prot_tot_0) = ind
#   }else{prot_A_0 = vector(); prot_NA_0 = vector(); prot_tot_0 = vector()}
#   if(M!=0){
#     met_tot_0 = matrix(0, nrow = M, ncol = N); rownames(met_tot_0) = network$met; colnames(met_tot_0) = ind
#   }else{met_tot_0 = vector()}
#   
#   res = list("ind" = ind,
#              "QTL_TC" = QTL_TC,
#              "QTL_TL" = QTL_TL,
#              "QTL_DR" = QTL_DR,
#              "rna_0" = rna_0,
#              "prot_tot_0" = prot_tot_0,
#              "prot_A_0" = prot_A_0,
#              "prot_NA_0" = prot_NA_0,
#              "met_tot_0" = met_tot_0)
#   return(res)
}

rand_indiv = function(network, num = 1){
  
  G = length(network$genes)
  P = length(network$prot)
  M = length(network$met)
  
  # Individual name
  indname = paste0("ind",num)
  
  # Genotype effects
  
  # QTL_TC : vector of genotype effect on transcription (TF binding)
  QTL_TC = get(qtl_effect_transcription)(G); names(QTL_TC) = network$genes

  # QTL_TL : vector of genotype effect on translation (TLF binding)
  QTL_TL = get(qtl_effect_translation)(P); names(QTL_TL) = network$protcod
  
  # QTL_DR : vector of genotype effect on mRNA degradation
  QTL_DR = get(qtl_effect_RNAdecay)(G); names(QTL_DR) = network$genes
  
  # Initial abundance
  rna_0 = get(initial_abundance)(G); names(rna_0) = network$genes
  prot_NA_0 = get(initial_abundance)(P); names(prot_NA_0) = network$prot
  prot_A_0 = get(initial_abundance)(P); names(prot_A_0) = network$prot
  prot_tot_0 = prot_NA_0 + prot_A_0; names(prot_tot_0) = network$prot
  met_tot_0 = get(initial_abundance_met)(M); names(met_tot_0) = network$met
  
  res = list("indname" = indname,
             "QTL_TC" = QTL_TC,
             "QTL_TL" = QTL_TL,
             "QTL_DR" = QTL_DR,
             "rna_0" = rna_0,
             "prot_tot_0" = prot_tot_0,
             "prot_A_0" = prot_A_0,
             "prot_NA_0" = prot_NA_0,
             "met_tot_0" = met_tot_0)
  return(res)
}

rand_indiv_null = function(network, num = 1){
  
  G = length(network$genes)
  P = length(network$prot)
  M = length(network$met)
  
  # Individual name
  indname = paste0("ind",num)
  
  # Genotype effects
  
  # QTL_TC : vector of genotype effect on transcription (TF binding)
  QTL_TC = rep(1, G); names(QTL_TC) = network$genes
  
  # QTL_TL : vector of genotype effect on translation (TLF binding)
  QTL_TL = rep(1, P); names(QTL_TL) = network$protcod
  
  # QTL_DR : vector of genotype effect on mRNA degradation
  QTL_DR = rep(1, G); names(QTL_DR) = network$genes
  
  # Initial abundance
  rna_0 = get(initial_abundance)(G); names(rna_0) = network$genes
  prot_NA_0 = get(initial_abundance)(P); names(prot_NA_0) = network$prot
  prot_A_0 = get(initial_abundance)(P); names(prot_A_0) = network$prot
  prot_tot_0 = prot_NA_0 + prot_A_0; names(prot_tot_0) = network$prot
  met_tot_0 = get(initial_abundance)(M); names(met_tot_0) = network$met
  
  res = list("indname" = indname,
             "QTL_TC" = QTL_TC,
             "QTL_TL" = QTL_TL,
             "QTL_DR" = QTL_DR,
             "rna_0" = rna_0,
             "prot_tot_0" = prot_tot_0,
             "prot_A_0" = prot_A_0,
             "prot_NA_0" = prot_NA_0,
             "met_tot_0" = met_tot_0)
  return(res)
}


##########################################################################################################################
#                                                       MAIN CODE                                                        #
##########################################################################################################################

simu = function(network, cohort, tmax){ 
#   with(as.list(c(network, cohort)),{
#   
#   # # ----
#   # ind = network$ind
#   # genes = network$genes
#   # protcod = network$protcod
#   # noncod = network$noncod
#   # prot = network$prot
#   # met = network$met
#   # g2p = network$g2p
#   # TF_sgn = network$TF_sgn
#   # TF_th = network$TF_th
#   # TF_n = network$TF_n
#   # TF_fc = network$TF_fc
#   # TLF_sgn = network$TLF_sgn
#   # TLF_th = network$TLF_th
#   # TLF_n = network$TLF_n
#   # TLF_fc = network$TLF_fc
#   # DR_sgn = network$DR_sgn
#   # DR_th = network$DR_th
#   # DR_n = network$DR_n
#   # DP_sgn = network$DP_sgn
#   # DP_th = network$DP_th
#   # DP_n = network$DP_n
#   # ACT_sgn = network$ACT_sgn
#   # ACT_th = network$ACT_th
#   # ACT_n = network$ACT_n
#   # DEACT_sgn = network$DEACT_sgn
#   # DEACT_th = network$DEACT_th
#   # DEACT_n = network$DEACT_n
#   # k_TC = network$k_TC
#   # k_TL = network$k_TL
#   # p0_DR = network$p0_DR
#   # p0_DP = network$p0_DP
#   # 
#   # ind = cohort$ind
#   # QTL_TC = cohort$QTL_TC
#   # QTL_TL = cohort$QTL_TL
#   # QTL_DR = cohort$QTL_DR
#   # rna_0 = cohort$rna_0
#   # prot_NA_0 = cohort$prot_NA_0
#   # prot_A_0 = cohort$prot_A_0
#   # met_tot_0 = cohort$met_tot_0
#   # # ----
#   
#   ## Initialization ----
#   
#   t = 1
#   
#   N= length(ind)
#   G = length(genes)
#   P = length(prot)
#   NC = length(noncod)
#   M = length(met)
#   
#   # Vector creation
#   rna = matrix(NA, nrow = G, ncol = N); rownames(rna) = genes; colnames(rna) = ind
#   prot_tot = matrix(NA, nrow = P, ncol = N); rownames(prot_tot) = prot; colnames(prot_tot) = ind
#   prot_A = matrix(NA, nrow = P, ncol = N); rownames(prot_A) = prot; colnames(prot_A) = ind
#   prot_NA = matrix(NA, nrow = P, ncol = N); rownames(prot_NA) = prot; colnames(prot_NA) = ind
#   met_tot = matrix(NA, nrow = M, ncol = N); rownames(met_tot) = met; colnames(met_tot) = ind
#   
#   l_TC = matrix(NA, nrow = G, ncol = N, dimnames = list(genes, ind))
#   l_TL = matrix(NA, nrow = P, ncol = N, dimnames = list(protcod, ind))
#   p_DR = matrix(NA, nrow = G, ncol = N, dimnames = list(genes, ind))
#   p_DP = matrix(NA, nrow = P, ncol = N, dimnames = list(prot, ind))
#   p_act = matrix(NA, nrow = P, ncol = N, dimnames = list(prot, ind))
#   p_deact = matrix(NA, nrow = P, ncol = N, dimnames = list(prot, ind))
#   
#   
#   # Initial values
#   # Are the initial values the same for all individuals???
#   rna_prev = matrix(rna_0, nrow = G, ncol = N); rownames(rna_prev) = genes; colnames(rna_prev) = ind
#   prot_A_prev = matrix(prot_A_0, nrow = P, ncol = N); rownames(prot_A_prev) = prot; colnames(prot_A_prev) = ind
#   prot_NA_prev = matrix(prot_NA_0, nrow = P, ncol = N); rownames(prot_NA_prev) = prot; colnames(prot_NA_prev) = ind
#   prot_tot_prev = matrix(prot_A_prev[prot, ind] + prot_NA_prev[prot, ind], nrow = P, ncol = N); rownames(prot_tot_prev) = prot; colnames(prot_tot_prev) = ind
#   met_tot_prev = matrix(met_tot_0, nrow = M, ncol = N); rownames(met_tot_prev) = met; colnames(met_tot_prev) = ind
#   
#   
#   # For visualization
#   
#   time_rna = vector("list", G); names(time_rna) = genes
#   time_prot_tot = vector("list", P); names(time_prot_tot) = prot
#   time_prot_A = vector("list", P); names(time_prot_A) = prot
#   time_prot_NA = vector("list", P); names(time_prot_NA) = prot
#   time_met_tot = vector("list", M); names(time_met_tot) = met
#   
#   time_TC_rate = vector("list", G); names(time_TC_rate) = genes
#   
#   # Add t = 0
#   for(g in genes){ time_rna[[g]] = cbind(time_rna[[g]],rna_prev[g,]) }
#   for(p in prot){ 
#     time_prot_tot[[p]] = cbind(time_prot_tot[[p]], prot_tot_prev[p,])
#     time_prot_A[[p]] = cbind(time_prot_A[[p]], prot_A_prev[p,])
#     time_prot_NA[[p]] = cbind(time_prot_NA[[p]], prot_NA_prev[p,])
#   }
#   for(m in met){ time_met_tot[[m]] = cbind(time_met_tot[[m]], met_tot_prev[m,]) }
#   
#   # Main loop -----
#   
#   while(t<=tmax){
#     
#     mol_prev = rbind(as.matrix(rna_prev[noncod, ind], ncol = length(ind)), as.matrix(prot_A_prev[, ind], ncol = length(ind)), as.matrix(met_tot_prev[, ind], ncol = length(ind))); rownames(mol_prev) = c(noncod,prot,met); colnames(mol_prev) = ind
#     
#     ## 1: Update rates and probabilites ----
#     
#     # Transcription rate
#     l_TC[genes, ind] = k_TC[genes] * t(sapply(genes, function(g){
#       reg = mol_prev[colnames(TF_n),]^TF_n[g,]
#       theta = t(QTL_TC[g,] * matrix(TF_th[g,], nrow = N, ncol = ncol(TF_th), byrow = T)) ^ TF_n[g,]
#       temp = 1 + TF_sgn[g,] * TF_fc[g,] * (reg/(reg + theta)); rownames(temp) = colnames(TF_sgn); colnames(temp) = ind
#       return(apply(temp, 2, prod))
#     }))
#     
#     
#     # Translation rate 
#       l_TL[protcod, ind] = rna_prev[protcod,] * k_TL[protcod] * t(sapply(protcod, function(g){
#         reg = mol_prev[colnames(TLF_n),]^TLF_n[g,]
#         theta = t(QTL_TL[g,] * matrix(TLF_th[g,], nrow = N, ncol = ncol(TLF_th), byrow = T)) ^ TLF_n[g,]
#         temp = 1 + TLF_sgn[g,] * TLF_fc[g,] * (reg/(reg + theta)); rownames(temp) = colnames(TLF_sgn); colnames(temp) = ind
#         return(apply(temp, 2, prod))
#       }))
#     
#     
#     # RNA degradation rate
#       p_DR[genes, ind] = p0_DR[genes] * QTL_DR[genes, ind] * t(sapply(genes, function(g){
#         reg = mol_prev[colnames(DR_n),]^DR_n[g,]
#         theta = matrix(DR_th[g,], nrow = ncol(DR_th), ncol = N) ^ DR_n[g,]
#         temp = 1 + DR_sgn[g,] * (reg/(reg + theta)); rownames(temp) = colnames(DR_sgn); colnames(temp) = ind
#         return(apply(temp, 2, prod))
#       }))  
#       
#     
#     
#     # protein degradation rate
#       p_DP[prot, ind] = p0_DP[prot] * t(sapply(prot, function(g){
#         reg = mol_prev[colnames(DP_n),]^DP_n[g,]
#         theta = matrix(DP_th[g,], nrow = ncol(DP_th), ncol = N) ^ DP_n[g,]
#         temp = 1 + DP_sgn[g,] * (reg/(reg + theta)); rownames(temp) = colnames(DP_sgn); colnames(temp) = ind
#         return(apply(temp, 2, prod))
#       }))
#       
#     
#     # protein activation rate
#     p_act[prot, ind] = t(sapply(prot, function(g){
#       # if(sum(ACT_sgn[g,]) == 0){return(matrix(1, nrow = 1, ncol = N, dimnames = list(g, ind)))}
#       reg = mol_prev[colnames(ACT_n),]^ACT_n[g,]
#       theta = matrix(ACT_th[g,], nrow = ncol(ACT_th), ncol = N) ^ ACT_n[g,]
#       temp = ACT_sgn[g,] * (reg/(reg + theta)) + 1 - ACT_sgn[g,]; rownames(temp) = colnames(ACT_sgn); colnames(temp) = ind # + 1 - ACT_sgn[g,] allows to set non-regulators to 1 = neutral in the product of all regulator contributions
#       return(apply(temp, 2, prod))
#     }))
# 
#     
#     # protein deactivation rate
#     p_deact[prot, ind] = t(sapply(prot, function(g){
#       if(sum(DEACT_sgn[g,]) == 0){return(matrix(0, nrow = 1, ncol = N, dimnames = list(g, ind)))}
#       reg = mol_prev[colnames(DEACT_n),]^DEACT_n[g,]
#       theta = matrix(DEACT_th[g,], nrow = ncol(DEACT_th), ncol = N) ^ DEACT_n[g,]
#       temp = DEACT_sgn[g,] * (reg/(reg + theta)) + 1 - DEACT_sgn[g,]; rownames(temp) = colnames(DEACT_sgn); colnames(temp) = ind
#       return(apply(temp, 2, prod))
#     }))
#     
#     
#     
#     
#     ## 2: Compute new molecule abundances ----
#     
#     # Update RNA abundance
#     #    Gene transcription follows a Poisson law, and each RNA molecule at time t-1 has a probability p_DR of being degraded
#     rna[genes, ind] = rna_prev[genes, ind] + rpois(length(l_TC), l_TC) - rbinom(N*G, rna_prev, p_DR)
# 
#     # Update protein abundance 
#     
#     # Activation or degradation of inactive proteins
#     # Each inactive protein has 3 possible fates:
#     #      - Stay inactive with a probability of (1-proba_activation)*(1-proba_degradation)
#     #      - Become active with a probability of proba_activation*(1-proba_degradation)
#     #      - Be degraded with a probability of proba_degradation
#     probs = rbind( as.vector((1-p_act)*(1-p_DP)), as.vector(p_act*(1-p_DP)) ,as.vector(p_DP) ) # order of the columns: Prot1 for ind1, Prot2 for ind1, ... , ProtP for ind1, Prot1 for ind2, etc
#     # For each protein type in each individuals these fates are sampled using a multinomial distribution
#     temp = sapply(1:(N*P), function(i){rmultinom(1,prot_NA_prev[i],probs[,i])})
#     Xa = temp[2,] # proteins inactive at time t-1 that are activated
#     Xdna = temp[3,] # proteins inactive at time t-1 that are degraded
#     
#     # Each active protein has 3 possible fates:
#     #      - Stay active with a probability of (1-proba_deactivation)*(1-proba_degradation)
#     #      - Become inactive with a probability of proba_deactivation*(1-proba_degradation)
#     #      - Be degraded with a probability of proba_degradation
#     probs = rbind( as.vector((1-p_deact)*(1-p_DP)), as.vector(p_deact*(1-p_DP)) ,as.vector(p_DP) ) # order of the columns: Prot1 for ind1, Prot2 for ind1, ... , ProtP for ind1, Prot1 for ind2, etc
#     # For each protein type in each individuals these fates are sampled using a multinomial distribution
#     temp = sapply(1:(N*P), function(i){rmultinom(1,prot_A_prev[i],probs[,i])})
#     Xna = temp[2,] # proteins active at time t-1 that are deactivated
#     Xda = temp[3,] # proteins active at time t-1 that are degraded
#     
#     
#     # Update inactive protein abundance
#     prot_NA[prot, ind] = rpois(length(l_TL), l_TL) + prot_NA_prev[prot, ind] - Xa - Xdna + Xna
#     
#     # Update active protein abundance
#     prot_A[prot, ind] = prot_A_prev[prot, ind] - Xna - Xda + Xa
#     
#     # Update total protein abundance
#     prot_tot[prot, ind] = prot_NA[prot, ind] + prot_A[prot, ind]
#     
#     
#     
#     ## Update the '_prev' matrices ----
#     rna_prev = rna
#     prot_tot_prev = prot_tot
#     prot_NA_prev = prot_NA
#     prot_A_prev = prot_A
#     
#     ## Store values for visualization
#     for(g in genes){ 
#       time_rna[[g]] = cbind(time_rna[[g]],rna[g,])
#       time_TC_rate[[g]] = cbind(time_TC_rate[[g]],l_TC[g,])
#       }
#     for(p in prot){ 
#       time_prot_tot[[p]] = cbind(time_prot_tot[[p]], prot_tot[p,])
#       time_prot_A[[p]] = cbind(time_prot_A[[p]], prot_A[p,])
#       time_prot_NA[[p]] = cbind(time_prot_NA[[p]], prot_NA[p,])
#       }
#     for(m in met){ time_met_tot[[m]] = cbind(time_met_tot[[m]], met_tot[m,]) }
#     
#     t = t + 1
#   }
#   
#   for(g in genes){ 
#     rownames(time_rna[[g]]) = ind
#     rownames(time_TC_rate[[g]]) = ind
#   }
#   for(p in prot){ 
#     rownames(time_prot_tot[[p]]) = ind
#     rownames(time_prot_A[[p]]) = ind
#     rownames(time_prot_NA[[p]]) = ind
#   }
#   for(m in met){ rownames(time_met_tot[[m]]) = ind }
#   
#   res = list("time_rna" = time_rna,
#              "time_prot_tot" = time_prot_tot,
#              "time_prot_A" = time_prot_A,
#              "time_prot_NA" = time_prot_NA,
#              "time_met_tot" = time_met_tot,
#              "time_TC_rate" = time_TC_rate)
# 
#   return(res)
#   })
}

simuInd = function(network, indiv, tmax, dt = 1){
  with(as.list(c(network, indiv)),{
    
    ## ONLY FOR TESTS INSIDE THE FUNCTION (create the variables of the environment network, indiv)
    # for(i in 1:length(network)){assign(names(network)[i], network[[i]])}
    # for(i in 1:length(indiv)){assign(names(indiv)[i], indiv[[i]])}
    
    ## Initialization ----
    
    t = 0
    
    G = length(genes)
    P = length(prot)
    NC = length(noncod)
    M = length(met)
    MR = length(metreactions)
    
    # Vector creation
    
    rna = vector(length = G); names(rna) = genes
    prot_NA = vector(length = P); names(prot_NA) = prot
    prot_A = vector(length = P); names(prot_A) = prot
    prot_tot = vector(length = P); names(prot_tot) = prot
    met_tot = vector(length = M); names(met_tot) = met
    
    l_TC = vector(length = G); names(l_TC) = genes
    l_TL = vector(length = P); names(l_TL) = protcod
    p_DR = vector(length = G); names(p_DR) = genes
    p_DP = vector(length = P); names(p_DP) = prot
    p_act = vector(length = P); names(p_act) = prot
    p_deact = vector(length = P); names(p_deact) = prot
        
    # SSA parameters creation (for metabolic reactions simulation)
    
    if(MR>0){
      rateFunc = function(x, params, t){
        sapply(params$metreactions, function(r){  x[params$E2R[r]] * params$ENZ["k_cat",r] * x[rownames(params$S)[params$S[,r] == -1]] / ( params$ENZ["K_M",r] + x[rownames(params$S)[params$S[,r] == -1]] ) })
      }
      params = as.list(c(network, indiv))
      Sprim = rbind(S, matrix(0, nrow = P, ncol = MR, dimnames = list(prot, metreactions)))
    }
 
    # Initial values
    
    rna_prev = rna_0; names(rna_prev) = genes
    prot_NA_prev = prot_NA_0; names(prot_NA_prev) = prot
    prot_A_prev = prot_A_0; names(prot_A_prev) = prot
    prot_tot_prev = prot_NA_0 + prot_A_0; names(prot_tot_prev) = prot
    met_tot_prev = met_tot_0; names(met_tot_prev) = met

      
    # For visualization
    
    time_abundance = matrix(c(t, rna_prev, prot_NA_prev, prot_A_prev, prot_tot_prev, met_tot_prev), nrow = 1, dimnames = list(c(),c("time",genes, paste0(prot,"_NA"), paste0(prot,"_A"), prot, met)))
    time_rate = matrix(nrow = 0, ncol = length(c(t, l_TC, l_TL, p_DR, p_DP, p_act, p_deact)))
    colnames(time_rate) = c( "time",
      sapply(genes, function(g){paste("l_TC",g, sep = "_")}),
      sapply(protcod, function(g){paste("l_TL",g, sep = "_")}),
      sapply(genes, function(g){paste("p_DR",g, sep = "_")}),
      sapply(prot, function(p){paste("p_DP",p, sep = "_")}),
      sapply(prot, function(p){paste("p_act",p, sep = "_")}),
      sapply(prot, function(p){paste("p_deact",p, sep = "_")})
      )
    
    # Main loop -----
    
    while(t<tmax){
      
      t = t + dt
      
      ## 1: Update reaction rates ----
      
      regprev = c(rna_prev[noncod], prot_A_prev[prot], met_tot_prev[met])
      
      regTC = colnames(TF_sgn)
      l_TC = k_TC * apply((1 + TF_sgn * TF_fc * t(regprev[regTC]^t(TF_n))/( (QTL_TC * TF_th)^TF_n +  t(regprev[regTC]^t(TF_n)))), 1, prod)
      
      regTL = colnames(TLF_sgn)
      l_TL = k_TL * apply((1 + TLF_sgn * TLF_fc * t(regprev[regTL]^t(TLF_n))/( (QTL_TL * TLF_th)^TLF_n +  t(regprev[regTL]^t(TLF_n)))), 1, prod)
      
      regDR = colnames(DR_sgn)
      p_DR = QTL_DR * p0_DR * apply((1 + DR_sgn * t(regprev[regDR]^t(DR_n))/( DR_th^DR_n +  t(regprev[regDR]^t(DR_n)))), 1, prod)
      
      regDP = colnames(DP_sgn)
      p_DP = p0_DP * apply((1 + DP_sgn * t(regprev[regDP]^t(DP_n))/( DP_th^DP_n +  t(regprev[regDP]^t(DP_n)))), 1, prod)
      
      regACT = colnames(ACT_sgn)
      p_act = apply((ACT_sgn * t(regprev[regACT]^t(ACT_n))/( ACT_th^ACT_n +  t(regprev[regACT]^t(ACT_n))) + 1 -ACT_sgn), 1, prod)
 
      regDEACT = colnames(DEACT_sgn)
      p_deact = (rowSums(DEACT_sgn) != 0) * apply((DEACT_sgn * t(regprev[regDEACT]^t(DEACT_n))/( DEACT_th^DEACT_n +  t(regprev[regDEACT]^t(DEACT_n))) + 1 -DEACT_sgn), 1, prod)
      
           
     
      
      ## 2: Compute new molecule abundances ----
      
      # Update RNA abundance
      #    Gene transcription follows a Poisson law, and each RNA molecule at time t-1 has a probability p_DR of being degraded
      rna[genes] = rna_prev[genes] + rpois(length(l_TC), l_TC*dt) - rbinom(G, rna_prev, p_DR*dt)
      
      # Update protein abundance 
      
      # Activation or degradation of inactive proteins
      # Each inactive protein has 3 possible fates:
      #      - Stay inactive with a probability of (1-proba_activation)*(1-proba_degradation)
      #      - Become active with a probability of proba_activation*(1-proba_degradation)
      #      - Be degraded with a probability of proba_degradation
      probs = rbind( (1-p_act)*(1-p_DP), p_act*(1-p_DP) ,p_DP ) # order of the columns: Prot1 for ind1, Prot2 for ind1, ... , ProtP for ind1, Prot1 for ind2, etc
      # For each protein type in each individuals these fates are sampled using a multinomial distribution
      temp = sapply(1:P, function(i){rmultinom(1,prot_NA_prev[i],probs[,i])})
      Xa = temp[2,] # proteins inactive at time t-1 that are activated
      Xdna = temp[3,] # proteins inactive at time t-1 that are degraded
      
      # Each active protein has 3 possible fates:
      #      - Stay active with a probability of (1-proba_deactivation)*(1-proba_degradation) (row 1)
      #      - Become inactive with a probability of proba_deactivation*(1-proba_degradation) (row 2)
      #      - Be degraded with a probability of proba_degradation (row 3)
      probs = rbind( (1-p_deact)*(1-p_DP), p_deact*(1-p_DP) ,p_DP ) # columns = proteins
      # For each protein type in each individuals these fates are sampled using a multinomial distribution
      temp = sapply(1:P, function(i){rmultinom(1,prot_A_prev[i],probs[,i])})
      Xna = temp[2,] # proteins active at time t-1 that are deactivated
      Xda = temp[3,] # proteins active at time t-1 that are degraded
      
      
      # Update inactive protein abundance
      prot_NA[prot] = rpois(length(l_TL), rna_prev[g2p[prot]] * l_TL) + prot_NA_prev[prot] - Xa - Xdna + Xna
      
      # Update active protein abundance
      prot_A[prot] = prot_A_prev[prot] - Xna - Xda + Xa
      
      # Update total protein abundance
      prot_tot[prot] = prot_NA[prot] + prot_A[prot]
      
      # Update metabolite abundance
      if(MR>0){ 
        SSAmetabo = ssa.adaptivetau(init.values = c(met_tot_prev[met], prot_A_prev[prot]), transitions = Sprim, rateFunc = rateFunc, params = params, tf = 1)
        met_tot[met] = SSAmetabo[max(which(SSAmetabo[, "time"] <= 1)), met] # we take the iteration that is closest to time = 1 (we don't take into account reactions occuring after one second)
      }else{met_tot[met] = met_tot_prev[met]}
      
      ## Update the '_prev' matrices ----
      rna_prev = rna
      prot_tot_prev = prot_tot
      prot_NA_prev = prot_NA
      prot_A_prev = prot_A
      met_tot_prev = met_tot
      
      
      ## Store values for visualization ----
      time_abundance = rbind(time_abundance, c(t, rna, prot_NA, prot_A, prot_tot, met_tot))
      time_rate = rbind(time_rate, c(t, l_TC, l_TL, p_DR, p_DP, p_act, p_deact))
      
    }
      
    res = list("time_abundance" = time_abundance,
               "time_rate" = time_rate)
    
    return(res)
  })
}

 
  
##########################################################################################################################
#                                                       OUTPUT                                                           #
##########################################################################################################################

visu_abundance = function(sim, network){
  
  # Plot RNA abundance for each gene
  for(g in network$genes){
    plot(sim$time_abundance[,"time"], sim$time_abundance[,g], col = "black", type = "l", main = paste(g, "abundance", sep = " "), xlab = "time", ylab = "Number of molecules")
  }
  
  # Plot total, active and inactive abundance for each protein
  for(p in  network$prot){
    plot(sim$time_abundance[,"time"], sim$time_abundance[,p], col = "black", type = "l", main = paste(p, "abundance", sep = " "), xlab = "time", ylab = "Number of molecules")
    lines(sim$time_abundance[,"time"], sim$time_abundance[,paste0(p, "_NA")], col = "blue")
    lines(sim$time_abundance[,"time"], sim$time_abundance[,paste0(p, "_A")], col = "red")
    legend("topleft", c("Total", "Active", "Inactive"), col = c("black", "red", "blue"), lty = c(1,1,1))
  }

  # Plot abundance for each metabolite
  for(m in network$met){
    plot(sim$time_abundance[,"time"], sim$time_abundance[,m], col = "black", type = "l", main = paste(m, "abundance", sep = " "), xlab = "time", ylab = "Number of molecules")
  }
  
}  

visu_rate = function(sim, network){
  
  # Transcription rates
  TC_rates = grep("l_TC", colnames(sim$time_rate), value = T)
  col = rainbow(length(TC_rates)); names(col) = TC_rates
  ymin = min(sim$time_rate[,TC_rates])
  ymax = max(sim$time_rate[,TC_rates])
  plot(sim$time_rate[, "time"], sim$time_rate[, "time"], type = 'n', ylim = c(ymin, ymax), main = "Transcription rates", xlab = "time", ylab = "rates")
  for(t in TC_rates){
    lines(sim$time_rate[, "time"], sim$time_rate[, t], col = col[t])
  }
  legend('topleft', gsub("l_TC_", "", TC_rates), col = col, lty = 1, horiz = T)
  
  # Translation rates
  TL_rates = grep("l_TL", colnames(sim$time_rate), value = T)
  col = rainbow(length(TL_rates)); names(col) = TL_rates
  ymin = min(sim$time_rate[,TL_rates])
  ymax = max(sim$time_rate[,TL_rates])
  plot(sim$time_rate[, "time"], sim$time_rate[, "time"], type = 'n', ylim = c(ymin, ymax), main = "Translation rates", xlab = "time", ylab = "rates")
  for(t in TL_rates){
    lines(sim$time_rate[, "time"], sim$time_rate[, t], col = col[t])
  }
  legend('topleft', gsub("l_TL_", "", TL_rates), col = col, lty = 1, horiz = T)
  
  # RNA degradation rates
  DR_rates = grep("p_DR", colnames(sim$time_rate), value = T)
  col = rainbow(length(DR_rates)); names(col) = DR_rates
  ymin = min(sim$time_rate[,DR_rates])
  ymax = max(sim$time_rate[,DR_rates])
  plot(sim$time_rate[, "time"], sim$time_rate[, "time"], type = 'n', ylim = c(ymin, ymax), main = "RNA degradation rates", xlab = "time", ylab = "rates")
  for(t in DR_rates){
    lines(sim$time_rate[, "time"], sim$time_rate[, t], col = col[t])
  }
  legend('topleft', gsub("p_DR_", "", DR_rates), col = col, lty = 1, horiz = T)
  
  # Protein degradation rates
  DP_rates = grep("p_DP", colnames(sim$time_rate), value = T)
  col = rainbow(length(DP_rates)); names(col) = DP_rates
  ymin = min(sim$time_rate[,DP_rates])
  ymax = max(sim$time_rate[,DP_rates])
  plot(sim$time_rate[, "time"], sim$time_rate[, "time"], type = 'n', ylim = c(ymin, ymax), main = "Protein degradation rates", xlab = "time", ylab = "rates")
  for(t in DP_rates){
    lines(sim$time_rate[, "time"], sim$time_rate[, t], col = col[t])
  }
  legend('topleft', gsub("p_DP_", "", DP_rates), col = col, lty = 1, horiz = T)
  
  # Protein activation rates
  act_rates = grep("p_act", colnames(sim$time_rate), value = T)
  col = rainbow(length(act_rates)); names(col) = act_rates
  ymin = min(sim$time_rate[,act_rates])
  ymax = max(sim$time_rate[,act_rates])
  plot(sim$time_rate[, "time"], sim$time_rate[, "time"], type = 'n', ylim = c(ymin, ymax), main = "Protein activation rates", xlab = "time", ylab = "rates")
  for(t in act_rates){
    lines(sim$time_rate[, "time"], sim$time_rate[, t], col = col[t])
  }
  legend('topleft', gsub("p_act_", "", act_rates), col = col, lty = 1, horiz = T)
  
  # Protein inactivation rates
  deact_rates = grep("p_deact", colnames(sim$time_rate), value = T)
  col = rainbow(length(deact_rates)); names(col) = deact_rates
  ymin = min(sim$time_rate[,deact_rates])
  ymax = max(sim$time_rate[,deact_rates])
  plot(sim$time_rate[, "time"], sim$time_rate[, "time"], type = 'n', ylim = c(ymin, ymax), main = "Protein inactivation rates", xlab = "time", ylab = "rates")
  for(t in deact_rates){
    lines(sim$time_rate[, "time"], sim$time_rate[, t], col = col[t])
  }
  legend('topleft', gsub("p_deact_", "", deact_rates), col = col, lty = 1, horiz = T)
  
}


##########################################################################################################################
#                                              NETWORK VISUALIZATION                                                     #
##########################################################################################################################

GRN_plot = function(nw){
  G = length(nw$genes)
  P = length(nw$prot)
  M = length(nw$met)
  NC = length(nw$noncod)
  nwnodes = data.frame("ID" = c(nw$genes, nw$prot, nw$met), "mol_type" = rep(c("RNA","Protein", "Metabolite"),c(G,P,M)))
  nwedges = data.frame("source" = vector(), "target" = vector(), "reg_type" = vector(), "sign" = vector())
  
  # Adding synthesis edges
  nwedges = rbind(nwedges, data.frame("source" = unname(nw$g2p), "target" = names(nw$g2p), "reg_type" = rep("synthesis",P), "sign" = rep(1,P)))
  
  # Adding transcription regulation
  if(sum(dim(nw$TF_sgn)==0) == 0){
    nwedges = rbind(nwedges, data.frame("source" = unname(expand.grid(rownames(nw$TF_sgn), colnames(nw$TF_sgn))[,2]), "target" = unname(expand.grid(rownames(nw$TF_sgn), colnames(nw$TF_sgn))[,1]), "reg_type" = rep("transcription", length(nw$TF_sgn)), "sign" = as.vector(nw$TF_sgn)))
  }
  
  # Adding translation regulation
  if(sum(dim(nw$TLF_sgn)==0) == 0){
    nwedges = rbind(nwedges, data.frame("source" = unname(expand.grid(rownames(nw$TLF_sgn), colnames(nw$TLF_sgn))[,2]), "target" = unname(expand.grid(rownames(nw$TLF_sgn), colnames(nw$TLF_sgn))[,1]), "reg_type" = rep("translation", length(nw$TLF_sgn)), "sign" = as.vector(nw$TLF_sgn)))
  }
  
  # Adding RNA decay regulation
  if(sum(dim(nw$DR_sgn)==0) == 0){
    nwedges = rbind(nwedges, data.frame("source" = unname(expand.grid(rownames(nw$DR_sgn), colnames(nw$DR_sgn))[,2]), "target" = unname(expand.grid(rownames(nw$DR_sgn), colnames(nw$DR_sgn))[,1]), "reg_type" = rep("decay", length(nw$DR_sgn)), "sign" = as.vector(nw$DR_sgn)))
  }
  
  # Adding protein decay regulation
  if(sum(dim(nw$DP_sgn)==0) == 0){
    nwedges = rbind(nwedges, data.frame("source" = unname(expand.grid(rownames(nw$DP_sgn), colnames(nw$DP_sgn))[,2]), "target" = unname(expand.grid(rownames(nw$DP_sgn), colnames(nw$DP_sgn))[,1]), "reg_type" = rep("decay", length(nw$DP_sgn)), "sign" = as.vector(nw$DP_sgn)))
  }
  
  # Adding protein activation regulation
  if(sum(dim(nw$ACT_sgn)==0) == 0){
    nwedges = rbind(nwedges, data.frame("source" = unname(expand.grid(rownames(nw$ACT_sgn), colnames(nw$ACT_sgn))[,2]), "target" = unname(expand.grid(rownames(nw$ACT_sgn), colnames(nw$ACT_sgn))[,1]), "reg_type" = rep("protein_activation", length(nw$ACT_sgn)), "sign" = as.vector(nw$ACT_sgn)))
  }
  
  # Adding protein inactivation regulation
  if(sum(dim(nw$DEACT_sgn)==0) == 0){
    #nwedges = rbind(nwedges, data.frame("source" = unname(expand.grid(rownames(nw$DEACT_sgn), colnames(nw$DEACT_sgn))[,2]), "target" = unname(expand.grid(rownames(nw$DEACT_sgn), colnames(nw$DEACT_sgn))[,1]), "reg_type" = rep("protein_inactivation", length(nw$DEACT_sgn)), "sign" = as.vector(nw$DEACT_sgn)))
    nwedges = rbind(nwedges, data.frame("source" = unname(expand.grid(rownames(nw$DEACT_sgn), colnames(nw$DEACT_sgn))[,2]), "target" = unname(expand.grid(rownames(nw$DEACT_sgn), colnames(nw$DEACT_sgn))[,1]), "reg_type" = rep("protein_activation", length(nw$DEACT_sgn)), "sign" = (-1*as.vector(nw$DEACT_sgn))))
  }
  
  #Removing null edges
  nwedges = nwedges[nwedges$sign!=0,]
  
  # Plot network
  
  col_nodes = c("orange", "limegreen", "gold"); names(col_nodes) = c("RNA", "Protein", "Metabolite")
  
  col_edges = c("red", "blue"); names(col_edges) = c("1","-1")
  width_edges = c(2, 2, 2, 1, 1); names(width_edges) = c("synthesis", "transcription", "translation", "decay", "protein_activation")
  lty_edges = c(1, 2, 2, 1, 3); names(lty_edges) = c("synthesis", "transcription", "translation", "decay", "protein_activation")

  nwgraph = graph_from_data_frame(d=nwedges, vertices=nwnodes, directed=T) 
  
  # Delete from vertices metabolites that do not act as regulators ( = degree=0)
  nwgraph = delete_vertices(nwgraph,V(nwgraph)[(V(nwgraph)$name %in% nw$met & degree(nwgraph) == 0)])

  V(nwgraph)$color = unname(col_nodes[V(nwgraph)$mol_type])
  E(nwgraph)$color = unname(col_edges[as.character(E(nwgraph)$sign)])
  E(nwgraph)[E(nwgraph)$reg_type == "synthesis"]$color = "black"
  E(nwgraph)$width = unname(width_edges[E(nwgraph)$reg_type])
  E(nwgraph)$lty = unname(lty_edges[E(nwgraph)$reg_type])
  
  l <- layout_nicely(nwgraph)
  
  plot(nwgraph, edge.arrow.size = .4, vertex.label.cex = .8, vertex.size = 20, edge.curved =T, layout = l)
  
  #  legend(x = 0, y = -1.1, legend = names(lty_edges), lty = lty_edges, lwd = width_edges, ncol = 3, xjust = 0.5, text.width = 0.5, cex = 1, bty = "n", x.intersp = 0.3, y.intersp = 0.5)
  legend(x = -2.5, y = 0, legend = c(names(lty_edges),"", "positive regulation", "negative regulation"), lty = c(lty_edges, 1, 1, 1), lwd = c(width_edges, 0, 2, 2), col = c(rep("black",length(lty_edges)),"white","red","blue"), ncol = 1, yjust = 0.5, text.width = 0.5, cex = 1, bty = "n", x.intersp = 0.3, y.intersp = 0.5)
  
}

Metabolism_plot = function(nw){
  P = length(nw$prot)
  M = length(nw$met)
  MR = length(nw$metreactions)
  
  nwnodes = data.frame("ID" = c(nw$met, unname(nw$E2R)), "mol_type" = rep(c("Metabolite", "Protein"), c(M, length(nw$E2R))))
  nwedges = data.frame("source" = vector(), "target" = vector(), "type" = vector()) # type can be: toenz (from a metabolite to an enzyme = no arrow),
                                                                                    #              fromenz (from a metabolite to an enzyme = arrow)
                                                                                    #           or direct (from a metabolite to a metabolite = arrow)
  
  # Adding reactions catalyzed by an enzyme 
  nwedges = rbind(nwedges, data.frame("source" = sapply(names(nw$E2R), function(r){rownames(nw$S)[which(nw$S[,r] == -1)]}), "target" = unname(nw$E2R), "type" = rep("toenz", length(nw$E2R))))
  nwedges = rbind(nwedges, data.frame("source" = unname(nw$E2R), "target" = sapply(names(nw$E2R), function(r){rownames(nw$S)[which(nw$S[,r] == 1)]}), "type" = rep("fromenz", length(nw$E2R))))
  
  # Adding spontaneous reactions
  spontreact = setdiff(nw$metreactions, names(nw$E2R))
  nwedges = rbind(nwedges, data.frame("source" = sapply(spontreact, function(r){rownames(nw$S)[which(nw$S[,r] == -1)]}), "target" = sapply(spontreact, function(r){rownames(nw$S)[which(nw$S[,r] == 1)]}), "type" = rep("direct", length(spontreact))))
  
  # Plot network
  
  col_nodes = c("limegreen", "gold"); names(col_nodes) = c("Protein", "Metabolite")
  arrow_type = c(0, 2, 2); names(arrow_type) = c("toenz", "fromenz", "direct")
  
  nwgraph = graph_from_data_frame(d=nwedges, vertices=nwnodes, directed=T) 
  
  V(nwgraph)$color = unname(col_nodes[V(nwgraph)$mol_type])
  E(nwgraph)$arrow.mode = unname(arrow_type[E(nwgraph)$type])
  
  l <- layout_nicely(nwgraph)
  
  plot(nwgraph, edge.color = "black", edge.arrow.size = .6, vertex.label.cex = .8, vertex.size = 20, layout = l)
   
}

global_plot = function(nw){
  par(mfrow=c(1,2))
  GRN_plot(nw)
  Metabolism_plot(nw)
  par(mfrow=c(1,1))
}

##########################################################################################################################
#                      TRANSFORMATION OF SAMPLED NETWORK AND COHORT INTO GILLESPIESSA PARAMETERS                         #
##########################################################################################################################

paramSSA = function(network, cohort){
  
  # G = length(network$genes)
  # P = length(network$prot)
  # M = length(network$met)
  # 
  # protNAA = c(sapply(network$prot, function(p){paste(p,"NA",sep = "_")}), sapply(network$prot, function(p){paste(p,"A",sep = "_")}))
  # 
  # # INITIAL CONDITIONS ----
  # x0 = c(cohort$rna_0[,1], cohort$prot_NA_0[,1], cohort$prot_A_0[,1], cohort$met_tot_0)
  # names(x0) = c(network$genes, protNAA, network$met)
  # 
  # 
  # # PROPENSITY FUNCTIONS / REACTION RATES ---- 
  # 
  # # Regulation law for transcription and translation reactions (takes into account the fc parameter)
  # regulation_law_TCTL = function(id, reaction, network){
  #   reg = colnames(network[[paste0(reaction,"_sgn")]])[which(network[[paste0(reaction,"_sgn")]][id,]!=0)]
  #   reg[grepl('^P',reg)] = paste(reg[grepl('^P',reg)],"A",sep = "_") # consider only active proteins
  #   parname = paste0(reaction, id)
  #   sup = sapply(reg, function(x){paste(c("*(1+sgn",parname,x,"*fc",parname,x,"*(",x,"^n",parname,x,"/(",x,"^n",parname,x,"+th",parname,x,"^n",parname,x,")))"), collapse = "")})
  #   return(sup)
  # }
  # 
  # # Regulation law for the decay reactions
  # regulation_law_decay = function(id, reaction, network){
  #   reg = colnames(network[[paste0(reaction,"_sgn")]])[which(network[[paste0(reaction,"_sgn")]][id,]!=0)]
  #   reg[grepl('^P',reg)] = paste(reg[grepl('^P',reg)],"A",sep = "_") # consider only active proteins
  #   parname = paste0(reaction, id)
  #   sup = sapply(reg, function(x){paste(c("*(1+sgn",parname,x,"*(",x,"^n",parname,x,"/(",x,"^n",parname,x,"+th",parname,x,"^n",parname,x,")))"), collapse = "")})
  #   return(sup)
  # }
  # 
  # # Regulation law for the activation/inactivation reactions
  # regulation_law_act = function(id, reaction, network){
  #   #reg = names(which(network[[paste0(reaction,"_sgn")]][id,]!=0))
  #   reg = colnames(network[[paste0(reaction,"_sgn")]])[which(network[[paste0(reaction,"_sgn")]][id,]!=0)]
  #   reg[grepl('^P',reg)] = paste(reg[grepl('^P',reg)],"A",sep = "_") # consider only active proteins
  #   parname = paste0(reaction, id)
  #   sup = sapply(reg, function(x){paste(c("*(",x,"^n",parname,x,"/(",x,"^n",parname,x,"+th",parname,x,"^n",parname,x,"))"), collapse = "")})
  #   return(sup)
  # }
  # 
  # a = c( # ----
  #        # transcription reactions
  #        sapply(network$genes, function(g){
  #          sup = regulation_law_TCTL(g, 'TF', network)
  #          paste(c("k_TC",g,sup), collapse = "")}),
  #        # translation reactions
  #        sapply(network$prot, function(p){           
  #          sup = regulation_law_TCTL(network$g2p[p], 'TLF', network)
  #          paste(c(network$g2p[p],"*k_TL",p,sup), collapse = "")}),
  #        # RNA decay reactions
  #        sapply(network$genes, function(g){
  #          sup = regulation_law_decay(g, 'DR', network)
  #          paste(c(g,"*p0_DR",g,sup), collapse = "")}),
  #        # Protein decay reactions (for active proteins then for inactive proteins)
  #        sapply(protNAA, function(p){
  #          sup = regulation_law_decay(sub("_NA|_A","",p), 'DP', network)
  #          paste(c(p,"*p0_DP",sub("_NA|_A","",p),sup), collapse = "")}),
  #        # Protein activation reactions
  #        sapply(protNAA[1:P], function(p){
  #          sup = regulation_law_act(sub("_NA","",p), 'ACT', network)
  #          paste(c(p,sup), collapse = "")}),
  #        # Protein inactivation reactions
  #        sapply(protNAA[(P+1):(2*P)], function(p){
  #          sup = regulation_law_act(sub("_A","",p), 'DEACT', network)
  #          if(length(sup) == 0){return("0")}
  #          else{paste(c(p,sup), collapse = "")}})
  # ) # ----
  # 
  # names(a) = c(
  #   # transcription reactions
  #   sapply(network$genes, function(g){paste0("TRANSCRIPTION",g)}),
  #   # translation reactions
  #   sapply(network$prot, function(p){paste0("TRANSLATION",p)}),
  #   # RNA decay reactions
  #   sapply(network$genes, function(g){paste0("DECAY",g)}),
  #   # Protein decay reactions
  #   sapply(protNAA, function(p){paste0("DECAY",p)}),
  #   # protein activation reactions
  #   sapply(protNAA[1:P], function(p){paste0("ACTIVATION",p)}),
  #   # protein inactivation reactions
  #   sapply(protNAA[(P+1):(2*P)], function(p){paste0("DEACTIVATION",p)})
  # )
  # 
  # # STATE-CHANGE VECTOR ----
  # tempTCTL = diag(1, nrow = G+P, ncol = G+P); tempTCTL = rbind(tempTCTL, matrix(0, nrow = P, ncol = G+P))
  # tempDRDP = diag(-1, nrow = G+2*P, ncol = G+2*P)
  # tempPos = diag(1, nrow = P, ncol = P); tempNeg = diag(-1, nrow = P, ncol = P)
  # tempACT = rbind(matrix(0, nrow = G, ncol = P), tempNeg, tempPos)
  # tempDEACT = rbind(matrix(0, nrow = G, ncol = P), tempPos, tempNeg)
  # nu = cbind(tempTCTL, tempDRDP, tempACT, tempDEACT); nu = rbind(nu, matrix(0, nrow = M, ncol = ncol(nu))); rownames(nu) = names(x0)
  # 
  # # REACTION RATE PARAMETERS ----
  # combnames = function(par, reaction, target, reg){
  #   reg[grepl('^P',reg)] = paste(reg[grepl('^P',reg)],"A",sep = "_")
  #   comb = expand.grid(target, reg)
  #   if(nrow(comb)!=0){sapply(1:nrow(comb), function(i){paste0(par, reaction, comb[i,1], comb[i, 2])})}
  # }
  # 
  # parms = c( # ----
  #            # Transcription parameters
  #            network$k_TC, # Basal transcription rates
  #            as.vector(network$TF_sgn), # regulation direction
  #            as.vector(network$TF_th*cohort$QTL_TC[,1]), # regulation threshold
  #            as.vector(network$TF_n), # regulation power
  #            as.vector(network$TF_fc), # regulation fold change
  #            # Translation parameters
  #            network$k_TL, # Basal translation rates
  #            as.vector(network$TLF_sgn), # regulation direction
  #            as.vector(network$TLF_th*cohort$QTL_TL[,1]), # regulation threshold
  #            as.vector(network$TLF_n), # regulation power
  #            as.vector(network$TLF_fc), # regulation fold change
  #            # RNA decay parameters
  #            network$p0_DR*cohort$QTL_TC[,1], # Basal decay rates
  #            as.vector(network$DR_sgn), # regulation direction
  #            as.vector(network$DR_th), # regulation threshold
  #            as.vector(network$DR_n), # regulation power
  #            # protein decay parameters
  #            network$p0_DP, # Basal decay rates
  #            as.vector(network$DP_sgn), # regulation direction
  #            as.vector(network$DP_th), # regulation threshold
  #            as.vector(network$DP_n), # regulation power
  #            # protein activation parameters
  #            as.vector(network$ACT_sgn), # regulation direction
  #            as.vector(network$ACT_th), # regulation threshold
  #            as.vector(network$ACT_n), # regulation power
  #            # protein inactivation parameters
  #            as.vector(network$DEACT_sgn), # regulation direction
  #            as.vector(network$DEACT_th), # regulation threshold
  #            as.vector(network$DEACT_n) # regulation power
  # ) # ----
  # 
  # names(parms) = c( # ----
  #                   # Transcription parameters
  #                   sapply(network$genes, function(x){paste0("k_TC",x)}), # Basal transcription rates
  #                   combnames("sgn", "TF", rownames(network$TF_sgn), colnames(network$TF_sgn)), # regulation direction
  #                   combnames("th", "TF", rownames(network$TF_th), colnames(network$TF_th)), # regulation threshold
  #                   combnames("n", "TF", rownames(network$TF_n), colnames(network$TF_n)), # regulation power
  #                   combnames("fc", "TF", rownames(network$TF_fc), colnames(network$TF_fc)), # regulation fold change
  #                   # Translation parameters
  #                   sapply(network$prot, function(x){paste0("k_TL",x)}), # Basal translation rates
  #                   combnames("sgn", "TLF", rownames(network$TLF_sgn), colnames(network$TLF_sgn)), # regulation direction
  #                   combnames("th", "TLF", rownames(network$TLF_th), colnames(network$TLF_th)), # regulation threshold
  #                   combnames("n", "TLF", rownames(network$TLF_n), colnames(network$TLF_n)), # regulation power
  #                   combnames("fc", "TLF", rownames(network$TLF_fc), colnames(network$TLF_fc)), # regulation fold change
  #                   # RNA decay parameters
  #                   sapply(network$genes, function(x){paste0("p0_DR",x)}), # Basal decay rates
  #                   combnames("sgn", "DR", rownames(network$DR_sgn), colnames(network$DR_sgn)), # regulation direction
  #                   combnames("th", "DR", rownames(network$DR_th), colnames(network$DR_th)), # regulation threshold
  #                   combnames("n", "DR", rownames(network$DR_n), colnames(network$DR_n)), # regulation power
  #                   # protein decay parameters
  #                   sapply(network$prot, function(x){paste0("p0_DP",x)}), # Basal decay rates
  #                   combnames("sgn", "DP", rownames(network$DP_sgn), colnames(network$DP_sgn)), # regulation direction
  #                   combnames("th", "DP", rownames(network$DP_th), colnames(network$DP_th)), # regulation threshold
  #                   combnames("n", "DP", rownames(network$DP_n), colnames(network$DP_n)), # regulation power
  #                   # protein activation parameters
  #                   combnames("sgn", "ACT", rownames(network$ACT_sgn), colnames(network$ACT_sgn)), # regulation direction
  #                   combnames("th", "ACT", rownames(network$ACT_th), colnames(network$ACT_th)), # regulation threshold
  #                   combnames("n", "ACT", rownames(network$ACT_n), colnames(network$ACT_n)), # regulation power
  #                   # protein inactivation parameters
  #                   combnames("sgn", "DEACT", rownames(network$DEACT_sgn), colnames(network$DEACT_sgn)), # regulation direction
  #                   combnames("th", "DEACT", rownames(network$DEACT_th), colnames(network$DEACT_th)), # regulation threshold
  #                   combnames("n", "DEACT", rownames(network$DEACT_n), colnames(network$DEACT_n)) # regulation power
  # ) # ----
  # 
  # res = list("x0" = x0, "a" = a, "nu" = nu, "parms" = parms)
  # return(res)
} 

paramSSAindiv = function(network, indiv){
  
  G = length(network$genes)
  P = length(network$prot)
  M = length(network$met)
  
  protNAA = c(sapply(network$prot, function(p){paste(p,"NA",sep = "_")}), sapply(network$prot, function(p){paste(p,"A",sep = "_")}))
  
  # INITIAL CONDITIONS ----
  x0 = c(indiv$rna_0, indiv$prot_NA_0, indiv$prot_A_0, indiv$met_tot_0)
  names(x0) = c(network$genes, protNAA, network$met)
  
  
  # PROPENSITY FUNCTIONS / REACTION RATES ---- 
  
  # Regulation law for transcription and translation reactions (takes into account the fc parameter)
  regulation_law_TCTL = function(id, reaction, network){
    reg = colnames(network[[paste0(reaction,"_sgn")]])[which(network[[paste0(reaction,"_sgn")]][id,]!=0)]
    reg[grepl('^P',reg)] = paste(reg[grepl('^P',reg)],"A",sep = "_") # consider only active proteins
    parname = paste0(reaction, id)
    sup = sapply(reg, function(x){paste(c("*(1+sgn",parname,x,"*fc",parname,x,"*(",x,"^n",parname,x,"/(",x,"^n",parname,x,"+th",parname,x,"^n",parname,x,")))"), collapse = "")})
    return(sup)
  }
  
  # Regulation law for the decay reactions
  regulation_law_decay = function(id, reaction, network){
    reg = colnames(network[[paste0(reaction,"_sgn")]])[which(network[[paste0(reaction,"_sgn")]][id,]!=0)]
    reg[grepl('^P',reg)] = paste(reg[grepl('^P',reg)],"A",sep = "_") # consider only active proteins
    parname = paste0(reaction, id)
    sup = sapply(reg, function(x){paste(c("*(1+sgn",parname,x,"*(",x,"^n",parname,x,"/(",x,"^n",parname,x,"+th",parname,x,"^n",parname,x,")))"), collapse = "")})
    return(sup)
  }
  
  # Regulation law for the activation/inactivation reactions
  regulation_law_act = function(id, reaction, network){
    #reg = names(which(network[[paste0(reaction,"_sgn")]][id,]!=0))
    reg = colnames(network[[paste0(reaction,"_sgn")]])[which(network[[paste0(reaction,"_sgn")]][id,]!=0)]
    reg[grepl('^P',reg)] = paste(reg[grepl('^P',reg)],"A",sep = "_") # consider only active proteins
    parname = paste0(reaction, id)
    sup = sapply(reg, function(x){paste(c("*(",x,"^n",parname,x,"/(",x,"^n",parname,x,"+th",parname,x,"^n",parname,x,"))"), collapse = "")})
    return(sup)
  }
  
  a = c( # ----
         # transcription reactions
         sapply(network$genes, function(g){
           sup = regulation_law_TCTL(g, 'TF', network)
           paste(c("k_TC",g,sup), collapse = "")}),
         # translation reactions
         sapply(network$prot, function(p){           
           sup = regulation_law_TCTL(network$g2p[p], 'TLF', network)
           paste(c(network$g2p[p],"*k_TL",p,sup), collapse = "")}),
         # RNA decay reactions
         sapply(network$genes, function(g){
           sup = regulation_law_decay(g, 'DR', network)
           paste(c(g,"*p0_DR",g,sup), collapse = "")}),
         # Protein decay reactions (for active proteins then for inactive proteins)
         sapply(protNAA, function(p){
           sup = regulation_law_decay(sub("_NA|_A","",p), 'DP', network)
           paste(c(p,"*p0_DP",sub("_NA|_A","",p),sup), collapse = "")}),
         # Protein activation reactions
         sapply(protNAA[1:P], function(p){
           sup = regulation_law_act(sub("_NA","",p), 'ACT', network)
           paste(c(p,sup), collapse = "")}),
         # Protein inactivation reactions
         sapply(protNAA[(P+1):(2*P)], function(p){
           sup = regulation_law_act(sub("_A","",p), 'DEACT', network)
           if(length(sup) == 0){return("0")}
           else{paste(c(p,sup), collapse = "")}})
  ) # ----
  
  names(a) = c(
    # transcription reactions
    sapply(network$genes, function(g){paste0("TRANSCRIPTION",g)}),
    # translation reactions
    sapply(network$prot, function(p){paste0("TRANSLATION",p)}),
    # RNA decay reactions
    sapply(network$genes, function(g){paste0("DECAY",g)}),
    # Protein decay reactions
    sapply(protNAA, function(p){paste0("DECAY",p)}),
    # protein activation reactions
    sapply(protNAA[1:P], function(p){paste0("ACTIVATION",p)}),
    # protein inactivation reactions
    sapply(protNAA[(P+1):(2*P)], function(p){paste0("DEACTIVATION",p)})
  )
  
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
             as.vector(network$TF_th*indiv$QTL_TC), # regulation threshold
             as.vector(network$TF_n), # regulation power
             as.vector(network$TF_fc), # regulation fold change
             # Translation parameters
             network$k_TL, # Basal translation rates
             as.vector(network$TLF_sgn), # regulation direction
             as.vector(network$TLF_th*indiv$QTL_TL), # regulation threshold
             as.vector(network$TLF_n), # regulation power
             as.vector(network$TLF_fc), # regulation fold change
             # RNA decay parameters
             network$p0_DR*indiv$QTL_TC, # Basal decay rates
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
  
  res = list("x0" = x0, "a" = a, "nu" = nu, "parms" = parms)
  return(res)
}

##########################################################################################################################
#                         TRANSFORMATION OF SAMPLED NETWORK AND COHORT INTO DESOLVE PARAMETERS                           #
##########################################################################################################################

paramDeSolve = function(network, indiv){
  G = length(network$genes)
  P = length(network$prot)
  M = length(network$met)
  
  protNAA = c(sapply(network$prot, function(p){paste(p,"NA",sep = "_")}), sapply(network$prot, function(p){paste(p,"A",sep = "_")}))
  
  # Initial conditions ----
  x0 = c(indiv$rna_0, indiv$prot_NA_0, indiv$prot_A_0, indiv$met_tot_0)
  names(x0) = c(network$genes, protNAA, network$met)
  
  # Function giving the derivatives for each molecule ----
  
  func = function(t, y ,parms){
    with(parms, {
      protNAA = c(sapply(prot, function(p){paste(p,"NA",sep = "_")}), sapply(prot, function(p){paste(p,"A",sep = "_")}))
      names(protNAA) = rep(prot,2)
      
      res = vector(length = (length(genes)+2*length(prot)+length(met))); names(res) = c(genes, protNAA, met)
      
      ## Reaction rates
      regprev = y[c(noncod, paste0(prot, "_A"), met)]
      
      regTC = colnames(TF_sgn); regTC[grepl('^P',regTC)] = paste(regTC[grepl('^P',regTC)],"A",sep = "_")
      l_TC = k_TC * apply((1 + TF_sgn * TF_fc * t(regprev[regTC]^t(TF_n))/( (QTL_TC * TF_th)^TF_n +  t(regprev[regTC]^t(TF_n)))), 1, prod)
      
      regTL = colnames(TLF_sgn); regTL[grepl('^P',regTL)] = paste(regTL[grepl('^P',regTL)],"A",sep = "_")
      l_TL = k_TL * apply((1 + TLF_sgn * TLF_fc * t(regprev[regTL]^t(TLF_n))/( (QTL_TL * TLF_th)^TLF_n +  t(regprev[regTL]^t(TLF_n)))), 1, prod)
      
      regDR = colnames(DR_sgn); regDR[grepl('^P',regDR)] = paste(regDR[grepl('^P',regDR)],"A",sep = "_")
      p_DR = QTL_DR * p0_DR * apply((1 + DR_sgn * t(regprev[regDR]^t(DR_n))/( DR_th^DR_n +  t(regprev[regDR]^t(DR_n)))), 1, prod)
      
      regDP = colnames(DP_sgn); regDP[grepl('^P',regDP)] = paste(regDP[grepl('^P',regDP)],"A",sep = "_")
      p_DP = p0_DP * apply((1 + DP_sgn * t(regprev[regDP]^t(DP_n))/( DP_th^DP_n +  t(regprev[regDP]^t(DP_n)))), 1, prod)
      
      regACT = colnames(ACT_sgn); regACT[grepl('^P',regACT)] = paste(regACT[grepl('^P',regACT)],"A",sep = "_")
      p_act = apply((ACT_sgn * t(regprev[regACT]^t(ACT_n))/( ACT_th^ACT_n +  t(regprev[regACT]^t(ACT_n))) + 1 -ACT_sgn), 1, prod)
      
      regDEACT = colnames(DEACT_sgn); regDEACT[grepl('^P',regDEACT)] = paste(regDEACT[grepl('^P',regDEACT)],"A",sep = "_")
      p_deact = (rowSums(DEACT_sgn) != 0) * apply((DEACT_sgn * t(regprev[regDEACT]^t(DEACT_n))/( DEACT_th^DEACT_n +  t(regprev[regDEACT]^t(DEACT_n))) + 1 -DEACT_sgn), 1, prod)
      
      metrates = sapply(metreactions, function(r){  y[E2R[r]] * ENZ["k_cat",r] * y[rownames(S)[S[,r] == -1]] / ( ENZ["K_M",r] + y[rownames(S)[S[,r] == -1]] ) })
      names(metrates) = metreactions

      P = length(prot)
      
      res[genes] = l_TC[genes] - y[genes]*p_DR[genes]
      res[protNAA[1:P]] = y[g2p[prot]]*l_TL[g2p[prot]] - y[protNAA[1:P]]*p_act[prot]*(1-p_DP[prot]) + y[protNAA[(P+1):(2*P)]]*p_deact[prot]*(1-p_DP[prot]) - y[protNAA[1:P]]*p_DP[prot]
      res[protNAA[(P+1):(2*P)]] = y[protNAA[1:P]]*p_act[prot]*(1-p_DP[prot]) - y[protNAA[(P+1):(2*P)]]*p_deact[prot]*(1-p_DP[prot]) - y[protNAA[(P+1):(2*P)]]*p_DP[prot]
      res[met] = sapply(met, function(i){ sum(metrates[colnames(S)[S[i, ] == 1]]) - sum(metrates[colnames(S)[S[i, ] == -1]])})
      
      return(list(res))
    })
  }
  

  # Parameters ----
  parms = c(network, indiv)
  parms$E2R[] = paste0(parms$E2R,"_A")
  
  return(list("y" = x0, "func" = func, "parms" = parms))
}

##########################################################################################################################
#                       TRANSFORMATION OF SAMPLED NETWORK AND COHORT INTO ADAPTIVETAU PARAMETERS                        #
##########################################################################################################################

paramAT = function(network, ind){
  
  G = length(network$genes)
  P = length(network$prot)
  M = length(network$met)
  MR = length(network$metreactions)
  
  protNAA = c(sapply(network$prot, function(p){paste(p,"NA",sep = "_")}), sapply(network$prot, function(p){paste(p,"A",sep = "_")}))
  
  # Initial values ----
  init.values = c(ind$rna_0, ind$prot_NA_0, ind$prot_A_0, ind$met_tot_0)
  names(init.values) = c(network$genes, protNAA, network$met)
  
  # Propensity matrix ----
  tempTCTL = diag(1, nrow = G+P, ncol = G+P); tempTCTL = rbind(tempTCTL, matrix(0, nrow = P, ncol = G+P))
  tempDRDP = diag(-1, nrow = G+2*P, ncol = G+2*P)
  tempPos = diag(1, nrow = P, ncol = P); tempNeg = diag(-1, nrow = P, ncol = P)
  tempACT = rbind(matrix(0, nrow = G, ncol = P), tempNeg, tempPos)
  tempDEACT = rbind(matrix(0, nrow = G, ncol = P), tempPos, tempNeg)
  if(MR == 0){ transitions = cbind(tempTCTL, tempDRDP, tempACT, tempDEACT); transitions = rbind(transitions, matrix(0, nrow = M, ncol = ncol(transitions))) }
  if(MR>0){ transitions = cbind(tempTCTL, tempDRDP, tempACT, tempDEACT); transitions = adiag(transitions, network$S) }
  rownames(transitions) = names(init.values)
  
  # Parameters ----
  params = as.list(c(network, ind))
  params$E2R[] = paste0(params$E2R,"_A")
  
  # Rate function ----
  rateFunc = function(x, params, t){
    with(params, {
      
      regprev = x[c(noncod, paste0(prot, "_A"), met)]
      
      regTC = colnames(TF_sgn); regTC[grepl('^P',regTC)] = paste(regTC[grepl('^P',regTC)],"A",sep = "_")
      l_TC = k_TC * apply((1 + TF_sgn * TF_fc * t(regprev[regTC]^t(TF_n))/( (QTL_TC * TF_th)^TF_n +  t(regprev[regTC]^t(TF_n)))), 1, prod)
      
      regTL = colnames(TLF_sgn); regTL[grepl('^P',regTL)] = paste(regTL[grepl('^P',regTL)],"A",sep = "_")
      l_TL = k_TL * apply((1 + TLF_sgn * TLF_fc * t(regprev[regTL]^t(TLF_n))/( (QTL_TL * TLF_th)^TLF_n +  t(regprev[regTL]^t(TLF_n)))), 1, prod)
      
      regDR = colnames(DR_sgn); regDR[grepl('^P',regDR)] = paste(regDR[grepl('^P',regDR)],"A",sep = "_")
      p_DR = QTL_DR * p0_DR * apply((1 + DR_sgn * t(regprev[regDR]^t(DR_n))/( DR_th^DR_n +  t(regprev[regDR]^t(DR_n)))), 1, prod)
      
      regDP = colnames(DP_sgn); regDP[grepl('^P',regDP)] = paste(regDP[grepl('^P',regDP)],"A",sep = "_")
      p_DP = p0_DP * apply((1 + DP_sgn * t(regprev[regDP]^t(DP_n))/( DP_th^DP_n +  t(regprev[regDP]^t(DP_n)))), 1, prod)
      
      regACT = colnames(ACT_sgn); regACT[grepl('^P',regACT)] = paste(regACT[grepl('^P',regACT)],"A",sep = "_")
      p_act = apply((ACT_sgn * t(regprev[regACT]^t(ACT_n))/( ACT_th^ACT_n +  t(regprev[regACT]^t(ACT_n))) + 1 -ACT_sgn), 1, prod)
      
      regDEACT = colnames(DEACT_sgn); regDEACT[grepl('^P',regDEACT)] = paste(regDEACT[grepl('^P',regDEACT)],"A",sep = "_")
      p_deact = (rowSums(DEACT_sgn) != 0) * apply((DEACT_sgn * t(regprev[regDEACT]^t(DEACT_n))/( DEACT_th^DEACT_n +  t(regprev[regDEACT]^t(DEACT_n))) + 1 -DEACT_sgn), 1, prod)
      
      metrates = sapply(metreactions, function(r){  x[E2R[r]] * ENZ["k_cat",r] * x[rownames(S)[S[,r] == -1]] / ( ENZ["K_M",r] + x[rownames(S)[S[,r] == -1]] ) })
      
      return(c(l_TC,
               x[protcod] * l_TL,
               x[genes] * p_DR,
               x[paste0(prot, "_NA")] * p_DP,
               x[paste0(prot, "_A")] * p_DP,
               x[paste0(prot, "_NA")] * p_act,
               x[paste0(prot, "_A")] * p_deact,
               unlist(metrates)))
    })
  }
  
  return(list("init.values" = init.values,
              "transitions" = transitions,
              "params" = params,
              rateFunc = rateFunc))
  
}

