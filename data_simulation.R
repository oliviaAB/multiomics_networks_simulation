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

##########################################################################################################################
#                                                  NETWORK SIMULATION                                                    #
##########################################################################################################################


rand_network = function(G, P, M){
 
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
  
  # Transcription regulation
  # TF : array for transcription regulation network, target genes (G rows) x transcription regulators (NC + P (=G) columns)
  TF_sgn = matrix( sample(c(-1:1), size = (G*G), replace = T, prob = proba_reg_TF), nrow = G, ncol = G, dimnames = list(genes, c(noncod,prot)))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  TF_th = matrix( 0, nrow = G, ncol = G, dimnames = list(genes, c(noncod,prot))) ; TF_th[which(TF_sgn!=0)] = get(th_sampling)(length(which(TF_sgn!=0)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  TF_n = matrix( 0, nrow = G, ncol = G, dimnames = list(genes, c(noncod,prot))) ; TF_n[which(TF_sgn!=0)] = get(n_sampling)(length(which(TF_sgn!=0)))
  
  # Translation regulation
  # TLF : array for translation regulation network, target protein-coding genes (P rows) x transcription regulators (NC + P (=G) columns)
  TLF_sgn = matrix( sample(c(-1:1), size = (P*G), replace = T, prob = proba_reg_TLF), nrow = P, ncol = G, dimnames = list(protcod, c(noncod,prot)))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  TLF_th = matrix( 0, nrow = P, ncol = G, dimnames = list(protcod, c(noncod,prot))) ; TLF_th[which(TLF_sgn!=0)] = get(th_sampling)(length(which(TLF_sgn!=0)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  TLF_n = matrix( 0, nrow = P, ncol = G, dimnames = list(protcod, c(noncod,prot))) ; TLF_n[which(TLF_sgn!=0)] = get(n_sampling)(length(which(TLF_sgn!=0)))
  
  # RNA decay
  # DR : array for RNA degradation regulation network, target genes (lenght of 1st dimension = G) X regulators non coding RNAs (NC+Q columns)
  DR_sgn = matrix( sample(c(-1:1), size = (G*NC), replace = T, prob = proba_reg_DR), nrow = G, ncol = NC, dimnames = list(genes, noncod))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  DR_th = matrix( 0, nrow = G, ncol = NC, dimnames = list(genes, noncod)) ; DR_th[which(DR_sgn!=0)] = get(th_sampling)(length(which(DR_sgn!=0)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  DR_n = matrix( 0, nrow = G, ncol = NC, dimnames = list(genes, noncod)) ; DR_n[which(DR_sgn!=0)] = get(n_sampling)(length(which(DR_sgn!=0)))
  
  
  # Protein decay
  # DP : array for RNA degradation regulation network, target genes (lenght of 1st dimension = G) X regulators non coding RNAs (NC+Q columns)
  DP_sgn = matrix( sample(c(-1:1), size = (P*P), replace = T, prob = proba_reg_DP), nrow = P, ncol = P, dimnames = list(prot, prot))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  DP_th = matrix( 0, nrow = P, ncol = P, dimnames = list(prot, prot)) ; DP_th[which(DP_sgn!=0)] = get(th_sampling)(length(which(DP_sgn!=0)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  DP_n = matrix( 0, nrow = P, ncol = P, dimnames = list(prot, prot)) ; DP_n[which(DP_sgn!=0)] = get(n_sampling)(length(which(DP_sgn!=0)))
  

  # Protein activation
  # ACT : array for protein activation network, target proteins (P rows) x activators (P + NC + M (= G+M) columns)
  ACT_sgn = matrix( sample(c(0,1), size = (P*(G+M)), replace = T, prob = proba_reg_ACT), nrow = P, ncol = NC+P+M, dimnames = list(prot, c(prot,noncod,met)))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  ACT_th = matrix( 0, nrow = P, ncol = G+M, dimnames = list(prot, c(prot,noncod,met))) ; ACT_th[which(ACT_sgn!=0)] = get(th_sampling)(length(which(ACT_sgn!=0)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  ACT_n = matrix( 0, nrow = P, ncol = G+M, dimnames = list(prot, c(prot,noncod,met))) ; ACT_n[which(ACT_sgn!=0)] = get(n_sampling)(length(which(ACT_sgn!=0)))
  
  # Protein deactivation
  # DEACT : array for protein deactivation network, target proteins (P rows) x activators (P + NC + M (= G+M) columns)
  DEACT_sgn = matrix( sample(c(0,1), size = (P*(G+M)), replace = T, prob = proba_reg_ACT), nrow = P, ncol = NC+P+M, dimnames = list(prot, c(prot,noncod,met)))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  DEACT_th = matrix( 0, nrow = P, ncol = G+M, dimnames = list(prot, c(prot,noncod,met))) ; DEACT_th[which(DEACT_sgn!=0)] = get(th_sampling)(length(which(DEACT_sgn!=0)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  DEACT_n = matrix( 0, nrow = P, ncol = G+M, dimnames = list(prot, c(prot,noncod,met))) ; DEACT_n[which(DEACT_sgn!=0)] = get(n_sampling)(length(which(DEACT_sgn!=0)))
  
  
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
             "g2p" = g2p,
             "TF_sgn" = TF_sgn,
             "TF_th" = TF_th,
             "TF_n" = TF_n,
             "TLF_sgn" = TLF_sgn,
             "TLF_th" = TLF_th,
             "TLF_n" = TLF_n,
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
             "p0_DP" = p0_DP
             )
    # ----
  
  return(res)
}

rand_network_null = function(G, P, M){
  
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
  
  # Transcription regulation
  # TF : array for transcription regulation network, target genes (G rows) x transcription regulators (NC + P (=G) columns)
  TF_sgn = matrix( 0, nrow = G, ncol = G, dimnames = list(genes, c(noncod,prot)))
  TF_th = matrix( 0, nrow = G, ncol = G, dimnames = list(genes, c(noncod,prot)))
  TF_n = matrix( 0, nrow = G, ncol = G, dimnames = list(genes, c(noncod,prot)))
  
  # Translation regulation
  # TLF : array for translation regulation network, target protein-coding genes (P rows) x transcription regulators (NC + P (=G) columns)
  TLF_sgn = matrix( 0, nrow = P, ncol = G, dimnames = list(protcod, c(noncod,prot)))
  TLF_th = matrix( 0, nrow = P, ncol = G, dimnames = list(protcod, c(noncod,prot)))
  TLF_n = matrix( 0, nrow = P, ncol = G, dimnames = list(protcod, c(noncod,prot)))
  
  # RNA decay
  # DR : array for RNA degradation regulation network, target genes (lenght of 1st dimension = G) X regulators non coding RNAs (NC+Q columns)
  DR_sgn = matrix( 0, nrow = G, ncol = NC, dimnames = list(genes, noncod))
  DR_th = matrix( 0, nrow = G, ncol = NC, dimnames = list(genes, noncod))
  DR_n = matrix( 0, nrow = G, ncol = NC, dimnames = list(genes, noncod))
  
  
  # Protein decay
  # DP : array for RNA degradation regulation network, target genes (lenght of 1st dimension = G) X regulators non coding RNAs (NC+Q columns)
  DP_sgn = matrix( 0, nrow = P, ncol = P, dimnames = list(prot, prot))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  DP_th = matrix( 0, nrow = P, ncol = P, dimnames = list(prot, prot))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  DP_n = matrix( 0, nrow = P, ncol = P, dimnames = list(prot, prot))
  
  
  # Protein activation
  # ACT : array for protein activation network, target proteins (P rows) x activators (P + NC + M (= G+M) columns)
  ACT_sgn = matrix( 0, nrow = P, ncol = NC+P+M, dimnames = list(prot, c(prot,noncod,met)))
  # The th parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 50 and 100
  ACT_th = matrix( 0, nrow = P, ncol = G+M, dimnames = list(prot, c(prot,noncod,met)))
  # The n parameter for the Hill function of each regulator is sampled from a uniform distribution on discrete numbers between 1 and 4
  ACT_n = matrix( 0, nrow = P, ncol = G+M, dimnames = list(prot, c(prot,noncod,met)))
  
  # Protein deactivation
  # DEACT : array for protein deactivation network, target proteins (P rows) x activators (P + NC + M (= G+M) columns)
  DEACT_sgn = matrix( 0, nrow = P, ncol = NC+P+M, dimnames = list(prot, c(prot,noncod,met)))
  DEACT_th = matrix( 0, nrow = P, ncol = G+M, dimnames = list(prot, c(prot,noncod,met)))
  DEACT_n = matrix( 0, nrow = P, ncol = G+M, dimnames = list(prot, c(prot,noncod,met)))
  
  
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
             "g2p" = g2p,
             "TF_sgn" = TF_sgn,
             "TF_th" = TF_th,
             "TF_n" = TF_n,
             "TLF_sgn" = TLF_sgn,
             "TLF_th" = TLF_th,
             "TLF_n" = TLF_n,
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
             "p0_DP" = p0_DP
  )
  # ----
  
  return(res)
}

##########################################################################################################################
#                                                INDIVIDUALS SIMULATION                                                  #
##########################################################################################################################


rand_cohort = function(network, N){

  G = length(network$genes)
  P = length(network$prot)
  M = length(network$met)
  
  # Individuals
  ind = sapply(1:N, function(i){ return(paste0("ind",i)) })

  # Genotype effects
  
  # QTL_TC : matrix of genotype effect on transcription (TF binding) for each individual, effects (G rows) x individuals (N columns)
  QTL_TC = matrix(get(qtl_effect_transcription)(G*N), nrow = G, ncol = N); rownames(QTL_TC) = network$genes; colnames(QTL_TC) = ind

  # QTL_TL : matrix of genotype effect on translation (TLF binding) for each individual, effects (P rows) x individuals (N columns)
  QTL_TL = matrix(get(qtl_effect_translation)(P*N), nrow = P, ncol = N); rownames(QTL_TL) = network$protcod; colnames(QTL_TL) = ind
  
  # QTL_DR : matrix of genotype effect on mRNA degradation for each individual, effects (G rows) x individuals (N columns)
  QTL_DR = matrix(get(qtl_effect_RNAdecay)(G*N), nrow = G, ncol = N); rownames(QTL_DR) = network$genes; colnames(QTL_DR) = ind
  
  if(G!=0){
    rna_0 = matrix(get(initial_abundance)(G,N), nrow = G, ncol = N); rownames(rna_0) = network$genes; colnames(rna_0) = ind
  }else{rna_0 = vector()}
  if(P!=0){
    prot_A_0 = matrix(get(initial_abundance)(P,N), nrow = P, ncol = N); rownames(prot_A_0) = network$prot; colnames(prot_A_0) = ind
    prot_NA_0 = matrix(get(initial_abundance)(P,N), nrow = P, ncol = N); rownames(prot_NA_0) = network$prot; colnames(prot_NA_0) = ind
    prot_tot_0 = prot_NA_0 + prot_A_0; rownames(prot_tot_0) = network$prot; colnames(prot_tot_0) = ind
  }else{prot_A_0 = vector(); prot_NA_0 = vector(); prot_tot_0 = vector()}
  if(M!=0){
    met_tot_0 = matrix(get(initial_abundance)(M,N), nrow = M, ncol = N); rownames(met_tot_0) = network$met; colnames(met_tot_0) = ind
  }else{met_tot_0 = vector()}
  
  res = list("ind" = ind,
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
  with(as.list(c(network, cohort)),{
  # # ----
  # ind = network$ind
  # genes = network$genes
  # protcod = network$protcod
  # noncod = network$noncod
  # prot = network$prot
  # met = network$met
  # g2p = network$g2p
  # TF_sgn = network$TF_sgn
  # TF_th = network$TF_th
  # TF_n = network$TF_n
  # TLF_sgn = network$TLF_sgn
  # TLF_th = network$TLF_th
  # TLF_n = network$TLF_n
  # DR_sgn = network$DR_sgn
  # DR_th = network$DR_th
  # DR_n = network$DR_n
  # DP_sgn = network$DP_sgn
  # DP_th = network$DP_th
  # DP_n = network$DP_n
  # ACT_sgn = network$ACT_sgn
  # ACT_th = network$ACT_th
  # ACT_n = network$ACT_n
  # DEACT_sgn = network$DEACT_sgn
  # DEACT_th = network$DEACT_th
  # DEACT_n = network$DEACT_n
  # k_TC = network$k_TC
  # k_TL = network$k_TL
  # p0_DR = network$p0_DR
  # p0_DP = network$p0_DP
  # 
  # ind = cohort$ind
  # QTL_TC = cohort$QTL_TC
  # QTL_TL = cohort$QTL_TL
  # QTL_DR = cohort$QTL_DR
  # # ----
  
  ## Initialization ----
  
  t = 1
  
  N= length(ind)
  G = length(genes)
  P = length(prot)
  NC = length(noncod)
  M = length(met)
  
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
  rna_prev = matrix(rna_0, nrow = G, ncol = N); rownames(rna_prev) = genes; colnames(rna_prev) = ind
  prot_A_prev = matrix(prot_A_0, nrow = P, ncol = N); rownames(prot_A_prev) = prot; colnames(prot_A_prev) = ind
  prot_NA_prev = matrix(prot_NA_0, nrow = P, ncol = N); rownames(prot_NA_prev) = prot; colnames(prot_NA_prev) = ind
  prot_tot_prev = matrix(prot_A_prev[prot, ind] + prot_NA_prev[prot, ind], nrow = P, ncol = N); rownames(prot_tot_prev) = prot; colnames(prot_tot_prev) = ind
  met_tot_prev = matrix(met_tot_0, nrow = M, ncol = N); rownames(met_tot_prev) = met; colnames(met_tot_prev) = ind
  
  
  # For visualization
  
  time_rna = vector("list", G); names(time_rna) = genes
  time_prot_tot = vector("list", P); names(time_prot_tot) = prot
  time_prot_A = vector("list", P); names(time_prot_A) = prot
  time_prot_NA = vector("list", P); names(time_prot_NA) = prot
  time_met_tot = vector("list", M); names(time_met_tot) = met
  
  time_TC_rate = vector("list", G); names(time_TC_rate) = genes
  
  # Add t = 0
  for(g in genes){ time_rna[[g]] = cbind(time_rna[[g]],rna_prev[g,]) }
  for(p in prot){ 
    time_prot_tot[[p]] = cbind(time_prot_tot[[p]], prot_tot_prev[p,])
    time_prot_A[[p]] = cbind(time_prot_A[[p]], prot_A_prev[p,])
    time_prot_NA[[p]] = cbind(time_prot_NA[[p]], prot_NA_prev[p,])
  }
  for(m in met){ time_met_tot[[m]] = cbind(time_met_tot[[m]], met_tot_prev[m,]) }
  
  # Main loop -----
  
  while(t<=tmax){
    
    mol_prev = rbind(as.matrix(rna_prev[noncod, ind], ncol = length(ind)), as.matrix(prot_A_prev[, ind], ncol = length(ind)), as.matrix(met_tot_prev[, ind], ncol = length(ind))); rownames(mol_prev) = c(noncod,prot,met); colnames(mol_prev) = ind
    
    ## 1: Update rates and probabilites ----
    
    # Transcription rate
    l_TC[genes, ind] = k_TC[genes] * t(sapply(genes, function(g){
      reg = mol_prev[colnames(TF_n),]^TF_n[g,]
      theta = t(QTL_TC[g,] * matrix(TF_th[g,], nrow = N, ncol = ncol(TF_th), byrow = T)) ^ TF_n[g,]
      temp = 1 + TF_sgn[g,] * (reg/(reg + theta)); rownames(temp) = colnames(TF_sgn); colnames(temp) = ind
      return(apply(temp, 2, prod))
    }))
    
    
    # Translation rate 
      l_TL[protcod, ind] = rna_prev[protcod,] * k_TL[protcod] * t(sapply(protcod, function(g){
        reg = mol_prev[colnames(TLF_n),]^TLF_n[g,]
        theta = t(QTL_TL[g,] * matrix(TLF_th[g,], nrow = N, ncol = ncol(TLF_th), byrow = T)) ^ TLF_n[g,]
        temp = 1 + TLF_sgn[g,] * (reg/(reg + theta)); rownames(temp) = colnames(TLF_sgn); colnames(temp) = ind
        return(apply(temp, 2, prod))
      }))
    
    
    # RNA degradation rate
      p_DR[genes, ind] = p0_DR[genes] * QTL_DR[genes, ind] * t(sapply(genes, function(g){
        reg = mol_prev[colnames(DR_n),]^DR_n[g,]
        theta = matrix(DR_th[g,], nrow = ncol(DR_th), ncol = N) ^ DR_n[g,]
        temp = 1 + DR_sgn[g,] * (reg/(reg + theta)); rownames(temp) = colnames(DR_sgn); colnames(temp) = ind
        return(apply(temp, 2, prod))
      }))  
      
    
    
    # protein degradation rate
      p_DP[prot, ind] = p0_DP[prot] * t(sapply(prot, function(g){
        reg = mol_prev[colnames(DP_n),]^DP_n[g,]
        theta = matrix(DP_th[g,], nrow = ncol(DP_th), ncol = N) ^ DP_n[g,]
        temp = 1 + DP_sgn[g,] * (reg/(reg + theta)); rownames(temp) = colnames(DP_sgn); colnames(temp) = ind
        return(apply(temp, 2, prod))
      }))
      
    
    # protein activation rate
    p_act[prot, ind] = t(sapply(prot, function(g){
      # if(sum(ACT_sgn[g,]) == 0){return(matrix(1, nrow = 1, ncol = N, dimnames = list(g, ind)))}
      reg = mol_prev[colnames(ACT_n),]^ACT_n[g,]
      theta = matrix(ACT_th[g,], nrow = ncol(ACT_th), ncol = N) ^ ACT_n[g,]
      temp = ACT_sgn[g,] * (reg/(reg + theta)) + 1 - ACT_sgn[g,]; rownames(temp) = colnames(ACT_sgn); colnames(temp) = ind # + 1 - ACT_sgn[g,] allows to set non-regulators to 1 = neutral in the product of all regulator contributions
      return(apply(temp, 2, prod))
    }))

    
    # protein deactivation rate
    p_deact[prot, ind] = t(sapply(prot, function(g){
      if(sum(DEACT_sgn[g,]) == 0){return(matrix(0, nrow = 1, ncol = N, dimnames = list(g, ind)))}
      reg = mol_prev[colnames(DEACT_n),]^DEACT_n[g,]
      theta = matrix(DEACT_th[g,], nrow = ncol(DEACT_th), ncol = N) ^ DEACT_n[g,]
      temp = DEACT_sgn[g,] * (reg/(reg + theta)) + 1 - DEACT_sgn[g,]; rownames(temp) = colnames(DEACT_sgn); colnames(temp) = ind
      return(apply(temp, 2, prod))
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
    for(g in genes){ 
      time_rna[[g]] = cbind(time_rna[[g]],rna[g,])
      time_TC_rate[[g]] = cbind(time_TC_rate[[g]],l_TC[g,])
      }
    for(p in prot){ 
      time_prot_tot[[p]] = cbind(time_prot_tot[[p]], prot_tot[p,])
      time_prot_A[[p]] = cbind(time_prot_A[[p]], prot_A[p,])
      time_prot_NA[[p]] = cbind(time_prot_NA[[p]], prot_NA[p,])
      }
    for(m in met){ time_met_tot[[m]] = cbind(time_met_tot[[m]], met_tot[m,]) }
    
    t = t + 1
  }
  
  for(g in genes){ 
    rownames(time_rna[[g]]) = ind
    rownames(time_TC_rate[[g]]) = ind
  }
  for(p in prot){ 
    rownames(time_prot_tot[[p]]) = ind
    rownames(time_prot_A[[p]]) = ind
    rownames(time_prot_NA[[p]]) = ind
  }
  for(m in met){ rownames(time_met_tot[[m]]) = ind }
  
  res = list("time_rna" = time_rna,
             "time_prot_tot" = time_prot_tot,
             "time_prot_A" = time_prot_A,
             "time_prot_NA" = time_prot_NA,
             "time_met_tot" = time_met_tot,
             "time_TC_rate" = time_TC_rate)

  return(res)
  })
}
  
  
##########################################################################################################################
#                                                       OUTPUT                                                           #
##########################################################################################################################

visu = function(network, cohort, sim, tmax){
  
  #cols = c('blue','green', 'red'); names(cols) = cohort$ind
  cols = rainbow(length(cohort$ind)); names(cols) = cohort$ind
  
  # plot genes
  for(g in network$genes){
    plot(1:(tmax+1), ylim = c(0,max(unlist(sim$time_rna[[g]]))), type = 'n', main = paste('Gene',g, sep=' '), xlab = 'time', ylab = 'abundance')
    for(i in cohort$ind){
      lines(sim$time_rna[[g]][i,], col=cols[i])
    }
    legend('topright', col = cols, lty=1, legend = names(cols))
  }
  
  # plot proteins
  for(p in network$prot){
    plot(1:(tmax+1), ylim = c(0,max(unlist(sim$time_prot_tot[[p]]))), type = 'n', main = paste('Protein',p, sep=' '), xlab = 'time', ylab = 'abundance')
    for(i in cohort$ind){
      lines(sim$time_prot_tot[[p]][i,], col=cols[i])
    }
    legend('topright', col = cols, lty=1, legend = names(cols))
  }
  
   # plot transcription rates
  for(g in network$genes){
    plot(1:(tmax+1), ylim = c(min(unlist(sim$time_TC_rate[[g]])),max(unlist(sim$time_TC_rate[[g]]))), type = 'n', main = paste('TC rate - Gene',g, sep=' '), xlab = 'time', ylab = 'Transcription rate')
    for(i in cohort$ind){
      lines(sim$time_TC_rate[[g]][i,], col=cols[i])
    }
    legend('topright', col = cols, lty=1, legend = names(cols))
    }
  
}  

