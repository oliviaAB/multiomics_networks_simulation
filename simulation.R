##########################################################################################################################
##########################################################################################################################
###                                  DATA SIMULATION - DRAFT 2 - 18/01/2018                                            ###
##########################################################################################################################
##########################################################################################################################



source("param.R")

if(!suppressWarnings(require("truncnorm", quietly = T))){install.packages("truncnorm")}
library(truncnorm)
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


rand_network = function(G, P, M , RM = NULL){
  
  
  ## Creating molecules ----

  if(G!=0){
    g.id = sapply(1:G, function(i){ return(paste0("RNA",i)) })
  } else{ g.id = vector()}
  if(P!=0){
    protspec.id = sapply(1:P, function(i){ return(paste0("P",i)) })
  } else{ p.id = vector()}
  if(M!=0){
    m.id = sapply(1:M, function(i){ return(paste0("M",i)) })
  } else{ m.id = vector()}
  
  if(G<P) stop("There must be at least as many genes as proteins")
  
  # Assigning a coding gene to each protein
  g2p = g.id[1:P]; names(g2p) = protspec.id
  protcod = unname(g2p); noncod = setdiff(g.id, protcod); NC = G-P
  
  # Creating active and inactive forms of proteins
  #   A protein can have active and inactive form (but not always the case)
  #   If it is the case the protein is synthesized in the inactive form an later can be activated/deactivated
  
  if(param.actprot>1 | param.actprot<0) stop("In param.R param.actprot value must be between 0 and 1")
  pNN.id = protspec.id[1:round(P*(1-param.actprot))]
  pNA.id = sapply(setdiff(protspec.id, pNN.id), function(i){paste0(i, "_NA")})
  pA.id = sapply(setdiff(protspec.id, pNN.id), function(i){paste0(i, "_A")})
  names(pA.id) = names(pNA.id) = setdiff(protspec.id, pNN.id)
    
  # p.id = vector of all possible protein states (considered as different species in the simulation)
  p.id = c(pNN.id, pNA.id, pA.id)
  
  ## Creating regulation network --
  
  # Transcription regulation ----
  # TF : array for transcription regulation network, target genes (G rows) x transcription regulators (NC + P (=G) columns)
  TF_sgn = matrix( sample(c(-1:1), size = (G*G), replace = T, prob = proba_reg_TF), nrow = G, ncol = G, dimnames = list(g.id, c(noncod, pNN.id, pA.id)))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  TF_th = matrix( 0, nrow = G, ncol = G, dimnames = list(g.id, c(noncod, pNN.id, pA.id))) ; TF_th[which(TF_sgn!=0)] = get(th_sampling)(length(which(TF_sgn!=0)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  TF_n = matrix( 0, nrow = G, ncol = G, dimnames = list(g.id, c(noncod, pNN.id, pA.id))) ; TF_n[which(TF_sgn!=0)] = get(n_sampling)(length(which(TF_sgn!=0)))
  # The fold change of target expression induced by maximal activation !! We only consider this parameter for activators of transcription (and not repressors, for which fc = 1)
  TF_fc = matrix( 0, nrow = G, ncol = G, dimnames = list(g.id, c(noncod, pNN.id, pA.id))) ; TF_fc[which(TF_sgn>0)] = get(fc_sampling)(length(which(TF_sgn>0))); TF_fc[which(TF_sgn<0)] = 1
  
  # Translation regulation ----
  # TLF : array for translation regulation network, target protein-coding genes (P rows) x transcription regulators (NC + P (=G) columns)
  TLF_sgn = matrix( sample(c(-1:1), size = (P*G), replace = T, prob = proba_reg_TLF), nrow = P, ncol = G, dimnames = list(protcod, c(noncod, pNN.id, pA.id)))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  TLF_th = matrix( 0, nrow = P, ncol = G, dimnames = list(protcod, c(noncod, pNN.id, pA.id))) ; TLF_th[which(TLF_sgn!=0)] = get(th_sampling)(length(which(TLF_sgn!=0)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  TLF_n = matrix( 0, nrow = P, ncol = G, dimnames = list(protcod, c(noncod, pNN.id, pA.id))) ; TLF_n[which(TLF_sgn!=0)] = get(n_sampling)(length(which(TLF_sgn!=0)))
  # The fold change of target translation rate induced by maximal activation !! We only consider this parameter for activators of translation (and not repressors, for which fc = 1)
  TLF_fc = matrix( 0, nrow = P, ncol = G, dimnames = list(protcod, c(noncod, pNN.id, pA.id))) ; TLF_fc[which(TLF_sgn>0)] = get(fc_sampling)(length(which(TLF_sgn>0))); TLF_fc[which(TLF_sgn<0)] = 1
 
  # RNA decay ----
  # DR : array for RNA degradation regulation network, target genes (lenght of 1st dimension = G) X regulators non coding RNAs (NC+Q columns)
  DR_sgn = matrix( sample(c(-1:1), size = (G*NC), replace = T, prob = proba_reg_DR), nrow = G, ncol = NC, dimnames = list(g.id, noncod))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  DR_th = matrix( 0, nrow = G, ncol = NC, dimnames = list(g.id, noncod)) ; DR_th[which(DR_sgn!=0)] = get(th_sampling)(length(which(DR_sgn!=0)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  DR_n = matrix( 0, nrow = G, ncol = NC, dimnames = list(g.id, noncod)) ; DR_n[which(DR_sgn!=0)] = get(n_sampling)(length(which(DR_sgn!=0)))
  
  # Protein decay ----
  # DP : array for RNA degradation regulation network, target genes (lenght of 1st dimension = G) X regulators non coding RNAs (NC+Q columns)
  DP_sgn = matrix( sample(c(-1:1), size = (P*P), replace = T, prob = proba_reg_DP), nrow = P, ncol = P, dimnames = list(protspec.id, c(pNN.id, pA.id)))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  DP_th = matrix( 0, nrow = P, ncol = P, dimnames = list(protspec.id, c(pNN.id, pA.id))) ; DP_th[which(DP_sgn!=0)] = get(th_sampling)(length(which(DP_sgn!=0)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  DP_n = matrix( 0, nrow = P, ncol = P, dimnames = list(protspec.id, c(pNN.id, pA.id))) ; DP_n[which(DP_sgn!=0)] = get(n_sampling)(length(which(DP_sgn!=0)))
  # Duplicate regulation for active and inactive states of proteins
  DP_sgn = rbind(DP_sgn[c(pNN.id, names(pA.id)),], DP_sgn[names(pA.id),]); rownames(DP_sgn) = c(pNN.id, pNA.id, pA.id)
  DP_th = rbind(DP_th[c(pNN.id, names(pA.id)),], DP_th[names(pA.id),]); rownames(DP_th) = c(pNN.id, pNA.id, pA.id)
  DP_n = rbind(DP_n[c(pNN.id, names(pA.id)),], DP_n[names(pA.id),]); rownames(DP_n) = c(pNN.id, pNA.id, pA.id)
  
  # Protein activation ----
  # ACT : array for protein activation network, target proteins (P rows) x activators (P + NC + M (= G+M) columns)
  ACT_sgn = matrix( sample(c(0,1), size = (length(pNA.id) * (P+NC+M)), replace = T, prob = proba_reg_ACT), nrow = length(pNA.id), ncol = NC+P+M, dimnames = list(pNA.id, c(pNN.id, pA.id, noncod,m.id)))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  ACT_th = matrix( 0, nrow = length(pNA.id), ncol = G+M, dimnames = list(pNA.id, c(pNN.id, pA.id, noncod,m.id))) ; ACT_th[which(ACT_sgn!=0)] = get(th_sampling)(length(which(ACT_sgn!=0)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  ACT_n = matrix( 0, nrow = length(pNA.id), ncol = G+M, dimnames = list(pNA.id, c(pNN.id, pA.id, noncod,m.id))) ; ACT_n[which(ACT_sgn!=0)] = get(n_sampling)(length(which(ACT_sgn!=0)))
  
  # Protein deactivation ----
  # DEACT : array for protein deactivation network, target proteins (P rows) x activators (P + NC + M (= G+M) columns)
  DEACT_sgn = matrix( sample(c(0,1), size = (length(pA.id) * (P+NC+M)), replace = T, prob = proba_reg_DEACT), nrow = length(pA.id), ncol = NC+P+M, dimnames = list(pA.id, c(pNN.id, pA.id, noncod,m.id)))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  DEACT_th = matrix( 0, nrow = length(pA.id), ncol = G+M, dimnames = list(pA.id, c(pNN.id, pA.id, noncod,m.id))) ; DEACT_th[which(DEACT_sgn!=0)] = get(th_sampling)(length(which(DEACT_sgn!=0)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  DEACT_n = matrix( 0, nrow = length(pA.id), ncol = G+M, dimnames = list(pA.id, c(pNN.id, pA.id, noncod,m.id))) ; DEACT_n[which(DEACT_sgn!=0)] = get(n_sampling)(length(which(DEACT_sgn!=0)))
  
  # Metabolic reactions ----
  #MR: number of metabolic reactions. Can be specified by the user, otherwise MR = M-1
  if(is.null(MR)) MR = M-1
  MR = min(MR, M*(M-1)) # there is maximum M(M-1) different possible reactions (given one reaction = Mi -> Mj)
  
  if(MR>0){ # If there is at least one metabolic reaction
    
    # Name of reactions
    mr.id = sapply(1:MR, function(i){ return(paste0("metreaction",i)) })
    
    # Compute all possible reactions (there is M(M-1) different possible reactions between M metabolites) => avoids to have twice the same reaction in the stoichiometry matrix
    temp = cbind(combn(m.id, 2),combn(m.id, 2)[2:1,])
    possreactions = matrix(0, nrow = M, ncol = M*(M-1)); rownames(possreactions) = m.id
    for(i in 1:ncol(temp)){
      possreactions[temp[1,i],i] = -1
      possreactions[temp[2,i],i] = 1
    }
    rm(temp)
    
    S = matrix(possreactions[,sample(ncol(possreactions),MR)], nrow = M, ncol = MR) # for each reaction randomly sample one column of the possible reactions matrix
    rownames(S) = m.id; colnames(S) = mr.id
    
    # Enzyme activity (only one enzyme can catalyze a given reaction, but an enzyme can catalyze several reactions)
    E2R = sample(c(pNN.id, pA.id), size = MR, replace = T); names(E2R) =  mr.id
    
    # Enzymatic rates
    ENZ = matrix(get(enzparam_sampling)(MR), nrow = 2, ncol = MR, dimnames = list(c("k_cat","K_M"), mr.id))
    
  }else{MR = 0; S = matrix(nrow = 0, ncol = 0); mr.id = vector(); E2R = vector(); ENZ = matrix(nrow = 0, ncol = 0)}
 
  
  ## Nodes parameters  ----
  
  # Transcription rates
  # TC rates chosen randomly between 0.01 and 0.1
  k_TC = get(basal_transcription_rate)(G) ; names(k_TC) = g.id
  
  # Translation rates
  # TL rates chosen randomly between 0.5 and 5
  k_TL = get(basal_translation_rate)(P); names(k_TL) = protcod
  
  # RNA decay rates
  #p0_DR = get(basal_RNAdecay_rate)(G); names(p0_DR) = g.id
  LT0_R = get(basal_RNAlifetime)(G); names(LT0_R) = g.id
  
  # Protein decay rates
  #p0_DP = get(basal_proteindecay_rate)(P); names(p0_DP) = protspec.id
  #p0_DP = c(p0_DP[c(pNN.id, names(pNA.id), names(pA.id))])
  #names(p0_DP) = c(pNN.id, pNA.id, pA.id)
  
  LT0_P = get(basal_protlifetime)(P); names(LT0_P) = protspec.id
  LT0_P = c(LT0_P[c(pNN.id, names(pNA.id), names(pA.id))])
  names(LT0_P) = c(pNN.id, pNA.id, pA.id)
  
  
  res = list("g.id" = g.id,
             "p.id" = p.id,
             "pNN.id" = pNN.id,
             "pNA.id" = pNA.id,
             "pA.id" = pA.id,
             "protspec.id" = protspec.id,
             "protcod" = protcod,
             "noncod" = noncod,
             "m.id" = m.id,
             "mr.id" = mr.id,
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
             # "p0_DR" = p0_DR,
             # "p0_DP" = p0_DP,
             "LT0_R" = LT0_R,
             "LT0_P" = LT0_P,
             "S" = S,
             "E2R" = E2R,
             "ENZ" = ENZ)
   
}

# Generates a network without regulation
rand_network_null = function(G, P, M , RM = NULL){
  
  
  ## Creating molecules ----
  
  if(G!=0){
    g.id = sapply(1:G, function(i){ return(paste0("RNA",i)) })
  } else{ g.id = vector()}
  if(P!=0){
    protspec.id = sapply(1:P, function(i){ return(paste0("P",i)) })
  } else{ p.id = vector()}
  if(M!=0){
    m.id = sapply(1:M, function(i){ return(paste0("M",i)) })
  } else{ m.id = vector()}
  
  if(G<P) stop("There must be at least as many genes as proteins")
  
  # Assigning a coding gene to each protein
  g2p = g.id[1:P]; names(g2p) = protspec.id
  protcod = unname(g2p); noncod = setdiff(g.id, protcod); NC = G-P
  
  # Creating active and inactive forms of proteins
  #   A protein can have active and inactive form (but not always the case)
  #   If it is the case the protein is synthesized in the inactive form an later can be activated/deactivated
  
  pNN.id = protspec.id
  pNA.id = vector()
  pA.id = vector()

  # p.id = vector of all possible protein states (considered as different species in the simulation)
  p.id = c(pNN.id, pNA.id, pA.id)
  
  ## Creating regulation network --
  
  # Transcription regulation ----
  # TF : array for transcription regulation network, target genes (G rows) x transcription regulators (NC + P (=G) columns)
  TF_sgn = matrix( 0, nrow = G, ncol = G, dimnames = list(g.id, c(noncod, pNN.id, pA.id)))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  TF_th = matrix( 0, nrow = G, ncol = G, dimnames = list(g.id, c(noncod, pNN.id, pA.id)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  TF_n = matrix( 0, nrow = G, ncol = G, dimnames = list(g.id, c(noncod, pNN.id, pA.id)))
  # The fold change of target expression induced by maximal activation !! We only consider this parameter for activators of transcription (and not repressors, for which fc = 1)
  TF_fc = matrix( 0, nrow = G, ncol = G, dimnames = list(g.id, c(noncod, pNN.id, pA.id)))
  
  # Translation regulation ----
  # TLF : array for translation regulation network, target protein-coding genes (P rows) x transcription regulators (NC + P (=G) columns)
  TLF_sgn = matrix( 0, nrow = P, ncol = G, dimnames = list(protcod, c(noncod, pNN.id, pA.id)))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  TLF_th = matrix( 0, nrow = P, ncol = G, dimnames = list(protcod, c(noncod, pNN.id, pA.id)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  TLF_n = matrix( 0, nrow = P, ncol = G, dimnames = list(protcod, c(noncod, pNN.id, pA.id)))
  # The fold change of target translation rate induced by maximal activation !! We only consider this parameter for activators of translation (and not repressors, for which fc = 1)
  TLF_fc = matrix( 0, nrow = P, ncol = G, dimnames = list(protcod, c(noncod, pNN.id, pA.id)))
  
  # RNA decay ----
  # DR : array for RNA degradation regulation network, target genes (lenght of 1st dimension = G) X regulators non coding RNAs (NC+Q columns)
  DR_sgn = matrix( 0, nrow = G, ncol = NC, dimnames = list(g.id, noncod))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  DR_th = matrix( 0, nrow = G, ncol = NC, dimnames = list(g.id, noncod))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  DR_n = matrix( 0, nrow = G, ncol = NC, dimnames = list(g.id, noncod))
  
  # Protein decay ----
  # DP : array for RNA degradation regulation network, target genes (lenght of 1st dimension = G) X regulators non coding RNAs (NC+Q columns)
  DP_sgn = matrix( 0, nrow = P, ncol = P, dimnames = list(protspec.id, c(pNN.id, pA.id)))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  DP_th = matrix( 0, nrow = P, ncol = P, dimnames = list(protspec.id, c(pNN.id, pA.id)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  DP_n = matrix( 0, nrow = P, ncol = P, dimnames = list(protspec.id, c(pNN.id, pA.id)))

  # Protein activation ----
  # ACT : array for protein activation network, target proteins (P rows) x activators (P + NC + M (= G+M) columns)
  ACT_sgn = matrix(nrow = 0, ncol = 0)
  ACT_th = matrix(nrow = 0, ncol = 0)
  ACT_n = matrix(nrow = 0, ncol = 0)
  
  # Protein deactivation ----
  # DEACT : array for protein deactivation network, target proteins (P rows) x activators (P + NC + M (= G+M) columns)
  DEACT_sgn = matrix(nrow = 0, ncol = 0)
  DEACT_th = matrix(nrow = 0, ncol = 0)
  DEACT_n = matrix(nrow = 0, ncol = 0)  
  # Metabolic reactions ----
  #MR: number of metabolic reactions. Can be specified by the user, otherwise MR = M-1
  if(is.null(MR)) MR = M-1
  MR = min(MR, M*(M-1)) # there is maximum M(M-1) different possible reactions (given one reaction = Mi -> Mj)
  
  if(MR>0){ # If there is at least one metabolic reaction
    
    # Name of reactions
    mr.id = sapply(1:MR, function(i){ return(paste0("metreaction",i)) })
    
    # Compute all possible reactions (there is M(M-1) different possible reactions between M metabolites) => avoids to have twice the same reaction in the stoichiometry matrix
    temp = cbind(combn(m.id, 2),combn(m.id, 2)[2:1,])
    possreactions = matrix(0, nrow = M, ncol = M*(M-1)); rownames(possreactions) = m.id
    for(i in 1:ncol(temp)){
      possreactions[temp[1,i],i] = -1
      possreactions[temp[2,i],i] = 1
    }
    rm(temp)
    
    S = matrix(possreactions[,sample(ncol(possreactions),MR)], nrow = M, ncol = MR) # for each reaction randomly sample one column of the possible reactions matrix
    rownames(S) = m.id; colnames(S) = mr.id
    
    # Enzyme activity (only one enzyme can catalyze a given reaction, but an enzyme can catalyze several reactions)
    E2R = sample(c(pNN.id, pA.id), size = MR, replace = T); names(E2R) =  mr.id
    
    # Enzymatic rates
    ENZ = matrix(get(enzparam_sampling)(MR), nrow = 2, ncol = MR, dimnames = list(c("k_cat","K_M"), mr.id))
    
  }else{MR = 0; S = matrix(nrow = 0, ncol = 0); mr.id = vector(); E2R = vector(); ENZ = matrix(nrow = 0, ncol = 0)}
  
  
  ## Nodes parameters  ----
  
  # Transcription rates
  # TC rates chosen randomly between 0.01 and 0.1
  k_TC = get(basal_transcription_rate)(G) ; names(k_TC) = g.id
  
  # Translation rates
  # TL rates chosen randomly between 0.5 and 5
  k_TL = get(basal_translation_rate)(P); names(k_TL) = protcod
  
  # RNA decay rates
  #p0_DR = get(basal_RNAdecay_rate)(G); names(p0_DR) = g.id
  LT0_R = get(basal_RNAlifetime)(G); names(LT0_R) = g.id
  
  # Protein decay rates
  #p0_DP = get(basal_proteindecay_rate)(P); names(p0_DP) = protspec.id
  #p0_DP = c(p0_DP[c(pNN.id, names(pNA.id), names(pA.id))])
  #names(p0_DP) = c(pNN.id, pNA.id, pA.id)
  
  LT0_P = get(basal_protlifetime)(P); names(LT0_P) = protspec.id
  LT0_P = c(LT0_P[c(pNN.id, names(pNA.id), names(pA.id))])
  names(LT0_P) = c(pNN.id, pNA.id, pA.id)
  
  res = list("g.id" = g.id,
             "p.id" = p.id,
             "pNN.id" = pNN.id,
             "pNA.id" = pNA.id,
             "pA.id" = pA.id,
             "protspec.id" = protspec.id,
             "protcod" = protcod,
             "noncod" = noncod,
             "m.id" = m.id,
             "mr.id" = mr.id,
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
             # "p0_DR" = p0_DR,
             # "p0_DP" = p0_DP,
             "LT0_R" = LT0_R,
             "LT0_P" = LR0_P,
             "S" = S,
             "E2R" = E2R,
             "ENZ" = ENZ)
  
}


##########################################################################################################################
#                                                INDIVIDUALS SIMULATION                                                  #
##########################################################################################################################

rand_indiv = function(network, num = 1){
  
  G = length(network$g.id)
  P = length(network$p.id)
  M = length(network$m.id)
  MR = length(network$mr.id)
  
  # Individual name
  indname = paste0("ind",num)
  
  # Genotype effects
  
  # QTL_TC : vector of genotype effect on transcription (TF binding)
  QTL_TC = get(qtl_effect_transcription)(G); names(QTL_TC) = network$g.id
  
  # QTL_TL : vector of genotype effect on translation (TLF binding)
  QTL_TL = get(qtl_effect_translation)(length(network$protspec.id)); names(QTL_TL) = network$protcod
  
  # QTL_DR : vector of genotype effect on mRNA degradation
  QTL_DR = get(qtl_effect_RNAdecay)(G); names(QTL_DR) = network$g.id
  
  # QTL_ENZ: vector of genotype effect on enzymatic activity
  QTL_ENZ = get(qtl_effect_enzact)(MR); names(QTL_ENZ) = network$mr.id
  
  # Initial abundance
  rna_0 = get(initial_rna_abundance)(G); names(rna_0) = network$g.id
  prot_0 = get(initial_prot_abundance)(P); names(prot_0) = network$p.id
  met_0 = get(initial_met_abundance)(M); names(met_0) = network$m.id
  
  res = list("indname" = indname,
             "QTL_TC" = QTL_TC,
             "QTL_TL" = QTL_TL,
             "QTL_DR" = QTL_DR,
             "QTL_ENZ" = QTL_ENZ,
             "rna_0" = rna_0,
             "prot_0" = prot_0,
             "met_0" = met_0)
  return(res)
}

# Generates an individual with no QTL effect
rand_indiv_null = function(network, num = 1){
  
  G = length(network$g.id)
  P = length(network$p.id)
  M = length(network$m.id)
  MR = length(network$mr.id)
  
  # Individual name
  indname = paste0("ind",num)
  
  # Genotype effects
  
  # QTL_TC : vector of genotype effect on transcription (TF binding)
  QTL_TC = rep(0, G); names(QTL_TC) = network$g.id
  
  # QTL_TL : vector of genotype effect on translation (TLF binding)
  QTL_TL = rep(0, length(network$protspec.id)); names(QTL_TL) = network$protcod
  
  # QTL_DR : vector of genotype effect on mRNA degradation
  QTL_DR = rep(0, G); names(QTL_DR) = network$g.id
  
  # QTL_ENZ: vector of genotype effect on enzymatic activity
  QTL_ENZ = rep(0, MR); names(QTL_ENZ) = network$mr.id
  
  # Initial abundance
  rna_0 = get(initial_rna_abundance)(G); names(rna_0) = network$g.id
  prot_0 = get(initial_prot_abundance)(P); names(prot_0) = network$p.id
  met_0 = get(initial_met_abundance)(M); names(met_0) = network$m.id
  
  res = list("indname" = indname,
             "QTL_TC" = QTL_TC,
             "QTL_TL" = QTL_TL,
             "QTL_DR" = QTL_DR,
             "QTL_ENZ" = QTL_ENZ,
             "rna_0" = rna_0,
             "prot_0" = prot_0,
             "met_0" = met_0)
  return(res)
}


##########################################################################################################################
#                                                       MAIN CODE                                                        #
##########################################################################################################################


simuInd = function(network, indiv, tmax, dt = 1){
  with(as.list(c(network, indiv)),{
    
    ## ONLY FOR TESTS INSIDE THE FUNCTION (create the variables of the environment network, indiv)
    # for(i in 1:length(network)){assign(names(network)[i], network[[i]])}; for(i in 1:length(indiv)){assign(names(indiv)[i], indiv[[i]])}
    
    ## Initialization ----
    
    t = 0
    
    G = length(g.id)
    P = length(p.id)
    NC = length(noncod)
    M = length(m.id)
    MR = length(mr.id)
    
    # Vector creation
    
    rna = vector(length = G); names(rna) = g.id
    prot = vector(length = P); names(prot) = p.id
    met = vector(length = M); names(met) = m.id
    
    l_TC = vector(length = G); names(l_TC) = g.id
    l_TL = vector(length = G-NC); names(l_TL) = protcod
    p_DR = vector(length = G); names(p_DR) = g.id
    p_DP = vector(length = P); names(p_DP) = p.id
    p_act = vector(length = length(pNA.id)); names(p_act) = pNA.id
    p_deact = vector(length = length(pA.id)); names(p_deact) = pA.id
    
    # SSA parameters creation (for metabolic reactions simulation)
    
    if(MR>0){
      rateFunc = function(x, params, t){
        sapply(params$metreactions, function(r){  x[params$E2R[r]] * params$ENZ["k_cat",r] * x[rownames(params$S)[params$S[,r] == -1]] / ( params$ENZ["K_M",r] + x[rownames(params$S)[params$S[,r] == -1]] ) })
      }
      params = as.list(E2R, ENZ, S)
      Sprim = rbind(S, matrix(0, nrow = P, ncol = MR, dimnames = list(p.id, mr.id)))
    }
    
    # Initial values
    
    rna_prev = rna_0; names(rna_prev) = g.id
    prot_prev = prot_0; names(prot_prev) = p.id
    met_prev = met_0; names(met_prev) = m.id
    
    
    # For visualization
    
    time_abundance = matrix(c(t, rna_prev, prot_prev, met_prev), nrow = 1, dimnames = list(c(),c("time",g.id, p.id, m.id)))
    time_rate = matrix(nrow = 0, ncol = length(c(t, l_TC, l_TL, p_DR, p_DP, p_act, p_deact)))
    colnames(time_rate) = c( "time",
                             sapply(g.id, function(g){paste("l_TC",g, sep = "_")}),
                             sapply(protcod, function(g){paste("l_TL",g, sep = "_")}),
                             sapply(g.id, function(g){paste("p_DR",g, sep = "_")}),
                             sapply(p.id, function(p){paste("p_DP",p, sep = "_")}),
                             sapply(pNA.id, function(p){paste("p_act",p, sep = "_")}),
                             sapply(pA.id, function(p){paste("p_deact",p, sep = "_")})
    )
    
    # Main loop -----
    
    while(t<tmax){
      
      t = t + dt
      
      ## 1: Update reaction rates ----
      
      regprev = c(rna_prev[noncod], prot_prev[c(pNN.id, pA.id)], met_prev[m.id])
      
      regTC = colnames(TF_sgn)
      l_TC = k_TC * (1 + rowSums((TF_sgn > 0)*TF_fc*t(regprev[regTC]^t(TF_n))/( (QTL_TC * TF_th)^TF_n +  t(regprev[regTC]^t(TF_n))))) * apply( (1-(TF_sgn<0)*(t(regprev[regTC]^t(TF_n))/( (QTL_TC * TF_th)^TF_n +  t(regprev[regTC]^t(TF_n))))), 1, prod)

      regTL = colnames(TLF_sgn)
      l_TL = k_TL * (1 + rowSums((TLF_sgn > 0)*TLF_fc*t(regprev[regTL]^t(TLF_n))/( (QTL_TL * TLF_th)^TLF_n +  t(regprev[regTL]^t(TLF_n))))) * apply( (1-(TLF_sgn<0)*(t(regprev[regTL]^t(TLF_n))/( (QTL_TL * TLF_th)^TLF_n +  t(regprev[regTL]^t(TLF_n))))), 1, prod)
      
      regDR = colnames(DR_sgn)
      p_DR = (1/(LT0_R*QTL_DR)) * (1 + ((LT0_R*QTL_DR) - 1)*(1 - apply( 1 - (DR_sgn>0)*(t(regprev[regDR]^t(DR_n))/( DR_th^DR_n +  t(regprev[regDR]^t(DR_n)))) , 1, prod))) * apply((1 + (DR_sgn<0) * t(regprev[regDR]^t(DR_n))/( DR_th^DR_n +  t(regprev[regDR]^t(DR_n)))), 1, prod)
      
      
      regDP = colnames(DP_sgn)
      p_DP = (1/LT0_P) * (1 + (LT0_P - 1)*(1 - apply( 1 - (DP_sgn>0)*(t(regprev[regDP]^t(DP_n))/( DP_th^DP_n +  t(regprev[regDP]^t(DP_n)))) , 1, prod))) * apply((1 + (DP_sgn<0) * t(regprev[regDP]^t(DP_n))/( DP_th^DP_n +  t(regprev[regDP]^t(DP_n)))), 1, prod)
      
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


