# N = 4
# G = 5
# P = 3
# NC = G - P
# M = 6
# Q = 10


# Variables ----

ind = c("ind1", "ind2", "ind3"); N = length(ind)

genes = c('G1'); G = length(genes)
prot = c('P1'); P = length(prot)
met = c('M1'); M = length(met)

g2p = genes[1:P]; names(g2p) = prot
protcod = unname(g2p); noncod = setdiff(genes, protcod); NC = G-P

# Graph  ----

## Topology 

# TF : array for transcription regulation network, target genes (lenght of 1st dimension = G) x transcription regulators (lenght of 2nd dimension = NC + P (=G)) X edges values (lenght of 3rd dimension = 3)
nar = list (genes, c(noncod, prot), c("sgn", "th", "n"))
names(nar) = c("targets", "regulators", "values")
TF = array(data = 0, dim = c(G, G, 3), dimnames = nar)
# TF[,,'sgn'] = sample(c(-1:1),G*G,replace=T, prob=c(0.2,0.6,0.2))
# tofill=which(TF[, , "sgn"] != 0, arr.ind = T)
# TF[,,'th'][tofill] = sample(c(5:10),nrow(tofill),replace=T)
# TF[,,'n'][tofill] = sample(c(1:2),nrow(tofill),replace=T)

TF[1,1,'sgn'] = - 1
TF[1,1,'th'] = 50
TF[1,1,'n'] = 2

# TLF : array for translation regulation network, target protein-coding genes (lenght of 1st dimension = P) x transcription regulators (lenght of 2nd dimension = NC + P (=G)) X edges values (lenght of 3rd dimension = 3)
nar = list (protcod, c(noncod, prot), c("sgn", "th", "n"))
names(nar) = c("targets", "regulators", "values")
TLF = array(data = 0, dim = c(P, G, 3), dimnames = nar)
# TLF[,,'sgn'] = sample(c(-1:1),P*G,replace=T, prob=c(0.2,0.6,0.2))
# tofill=which(TLF[, , "sgn"] != 0, arr.ind = T)
# TLF[,,'th'][tofill] = sample(c(5:10),nrow(tofill),replace=T)
# TLF[,,'n'][tofill] = sample(c(0:2),nrow(tofill),replace=T)

# DR : array for RNA degradation regulation network, target genes (lenght of 1st dimension = G) X regulators non coding RNAs (lenght of 2nd dimension = NC) X edges values (lenght of 3rd dimension = 3)
nar = list (genes, c(noncod), c("sgn", "th", "n"))
names(nar) = c("targets", "regulators", "values")
DR = array(data = 0, dim = c(G, NC, 3), dimnames = nar)
# DR[,,'sgn'] = rbinom(G*NC,1,prob = 0.3)
# tofill=which(DR[, , "sgn"] != 0, arr.ind = T)
# DR[,,'th'][tofill] = sample(c(5:10),nrow(tofill),replace=T)
# DR[,,'n'][tofill] = sample(c(0:2),nrow(tofill),replace=T)

# DP: array for protein degradation regulation network, target proteins (lenght of 1st dimension = P) X regulators proteins (lenght of 2nd dimension = P) X edges values (lenght of 3rd dimension = 3)
nar = list (prot, c(prot), c("sgn", "th", "n"))
names(nar)=c("targets", "regulators", "values")
DP = array(data = 0, dim = c(P, P, 3), dimnames = nar)
# DP[,,'sgn']=rbinom(P*P,1,prob = 0.3)
# tofill=which(DP[, , "sgn"] != 0, arr.ind = T)
# DP[,,'th'][tofill] = sample(c(5:10),nrow(tofill),replace=T)
# DP[,,'n'][tofill] = sample(c(0:2),nrow(tofill),replace=T)

# ACT : array for protein activation network, target proteins (lenght of 1st dimension = P) x activators (lenght of 2nd dimension = P + M) X edges values (lenght of 3rd dimension = 3)
nar = list (prot, c(prot, met), c("sgn", "th", "n"))
names(nar)=c("targets", "regulators", "values")
ACT = array(data = 0, dim = c(P, P+M, 3), dimnames = nar)
# ACT[,,'sgn']=rbinom(P*(P+M),1,prob = 0.3)
# tofill=which(ACT[, , "sgn"] != 0, arr.ind = T)
# ACT[,,'th'][tofill] = sample(c(5:10),nrow(tofill),replace=T)
# ACT[,,'n'][tofill] = sample(c(0:2),nrow(tofill),replace=T)

# DEACT : array for protein deactivation network, target proteins (lenght of 1st dimension = P) x deactivators (lenght of 2nd dimension = P + M) X edges values (lenght of 3rd dimension = 3)
nar = list (prot, c(prot, met), c("sgn", "th", "n"))
names(nar)=c("targets", "regulators", "values")
DEACT = array(data = 0, dim = c(P, P+M, 3), dimnames = nar)
# DEACT[,,'sgn']=rbinom(P*(P+M),1,prob = 0.3)
# tofill=which(DEACT[, , "sgn"] != 0, arr.ind = T)
# DEACT[,,'th'][tofill] = sample(c(5:10),nrow(tofill),replace=T)
# DEACT[,,'n'][tofill] = sample(c(0:2),nrow(tofill),replace=T)



## Nodes parameters

# k_TC : vector of basal transcription rates for each gene, i.e. mean basal number of transcription events in per time unit (= in absence of any regulators), genes (length G)
# k_TC=runif(G,0,1); names(k_TC) = genes
k_TC=c(5); names(k_TC) = genes

# k_TL : vector of basal translation rates for each protein-coding gene, i.e. mean basal number of transcription events in per time unit (= in absence of any regulators), proteins (length P)
#k_TL=runif(P,0,1); names(k_TL) = protcod
k_TL=c(5); names(k_TL) = genes

# p0_DR : vector of basal RNA degradation rate for each gene, genes (lenght G)
# p0_DR=runif(G,0,1); names(p0_DR) = genes
p0_DR=c(0.1); names(p0_DR) = genes

# p0_DP : vector of basal protein degradation rate for each gene, genes (length P)
# p0_DP=runif(P,0,1); names(p0_DP) = prot
p0_DP=c(0.3); names(p0_DP) = prot



## Genotype effects

# QTL_TC : matrix of genotype effect on transcription (TF binding) for each individual, effects (G rows) x individuals (N columns)
# QTL_TC = matrix(runif(G * N ,0.7,1.3), nrow = G, ncol = N); rownames(QTL_TC) = genes; colnames(QTL_TC) = ind
QTL_TC = matrix(1, nrow = G, ncol = N); rownames(QTL_TC) = genes; colnames(QTL_TC) = ind
# QTL_TC[c(1,3,5)] = c(0.2,1,2)

# QTL_TL : matrix of genotype effect on translation (TLF binding) for each individual, effects (P rows) x individuals (N columns)
# QTL_TL = matrix(runif(P * N ,0.7,1.3), nrow = P, ncol = N); rownames(QTL_TL) = protcod; colnames(QTL_TL) = ind
QTL_TL = matrix(1, nrow = P, ncol = N); rownames(QTL_TL) = protcod; colnames(QTL_TL) = ind

# QTL_DR : matrix of genotype effect on mRNA degradation for each individual, effects (G rows) x individuals (N columns)
# QTL_DR = matrix(runif(G * N ,0.7,1.3), nrow = G, ncol = N); rownames(QTL_DR) = genes; colnames(QTL_DR) = ind
QTL_DR = matrix(1, nrow = G, ncol = N); rownames(QTL_DR) = genes; colnames(QTL_DR) = ind

## Initial values
# rna_0 = matrix(round(runif(G*N,0,10),0),nrow=G,ncol=N); rownames(rna_0) = genes; colnames(rna_0) = ind
# prot_A_0 = matrix(round(runif(P*N,0,10),0),nrow=P,ncol=N); rownames(prot_A_0) = prot; colnames(prot_A_0) = ind
# prot_NA_0 = matrix(round(runif(P*N,0,10),0),nrow=P,ncol=N); rownames(prot_NA_0) = prot; colnames(prot_NA_0) = ind
# met_tot_0 = matrix(round(runif(M*N,0,10),0),nrow=M,ncol=N); rownames(met_tot_0) = met; colnames(met_tot_0) = ind

rna_0 = matrix(c(5,5,5),nrow=G,ncol=N, byrow = T); rownames(rna_0) = genes; colnames(rna_0) = ind
prot_A_0 = matrix(c(5,5,5),nrow=P,ncol=N, byrow = T); rownames(prot_A_0) = prot; colnames(prot_A_0) = ind
prot_NA_0 = matrix(0,nrow=P,ncol=N); rownames(prot_NA_0) = prot; colnames(prot_NA_0) = ind
met_tot_0 = matrix(round(runif(M*N,0,10),0),nrow=M,ncol=N); rownames(met_tot_0) = met; colnames(met_tot_0) = ind



tmax = 100
