##########################################################################################################################
##########################################################################################################################
###                                       NETWORK GENERATION  - using Julia                                            ###
##########################################################################################################################
##########################################################################################################################

if(!suppressWarnings(require("XRJulia", quietly = T))){install.packages("XRJulia")}
library(XRJulia)
if(!suppressWarnings(require("igraph", quietly = T))){install.packages("igraph")}
library(igraph)
if(!suppressWarnings(require("ggplot2", quietly = T))){install.packages("ggplot2")}
library(ggplot2)
if(!suppressWarnings(require("truncnorm", quietly = T))){install.packages("truncnorm")}
library(truncnorm)

## Test if Julia is installed on the computer
if(!findJulia(test = T)) stop("Julia is not installed on the computer or not accessible by R. Check that Julia is correcly installed and/or in the PATH variable\n")


##########################################################################################################################
###                                               PARAMETERS                                                           ###
##########################################################################################################################


## Temporary
#setwd("~/winData/multiomics_networks_simulation")



##########################################################################################################################
###                                        DEFINE JULIA AND R FUNCTIONS                                                ###
##########################################################################################################################

## Get Julia functions from source code
juliaSource(paste0(getwd(),"/julia_functions.jl"))


# ------------------------------------------------------------------------------------------------------------ #
#                                            TEMPORARY FUNCTIONS                                               #
# ------------------------------------------------------------------------------------------------------------ # 


howmanyautoreg = function(edg){
  cat("There are  ")
  cat(sum(edg[,1] == edg[,2]))
  cat("  self-regulatory edges in the graph.\n")
}

howmanyloops = function(edg){
  cat("There are  ")
  res = edg[(edg[,1]!=edg[,2]), c("from", "to")] ## remove self-loops
  cat(sum(duplicated.matrix(data.frame("1" = c(res$from, res$to), "2" = c(res$to, res$from))))/2)
  cat("  2-nodes loops in the graph.\n")
}

ggbarplot = function(x, title, labx, laby){
  #data = data.frame("values" = factor(x, level = factor(0:max(x))))
  #g = ggplot(data, aes(x = values))+geom_bar() +ggtitle(title) + xlab(labx) + ylab(laby) + scale_x_discrete(drop=FALSE)
  
  data = data.frame("values" = x)
  g = ggplot(data, aes(x = values)) + geom_histogram() +ggtitle(title) + xlab(labx) + ylab(laby)
  return(g)
  
}

shownw = function(nw){
  
  reaction = edge_attr(nw, "TargetReaction", index = 1) # which reaction is regulated in the network (used to retrieve the regulators)
  
  ## Plot network
  colrs = c("target"="gold","PCreg"="yellowgreen", "NCreg"="violetred") # target genes are gold, protein regulators are yellowgreen, noncoding regulators are violetred
  V(nw)$color = colrs[V(nw)$nodetype]
  ## Color differently activating and inhibiting regulations
  colrs = c("1" = "brown1", "-1" = "dodgerblue")
  E(nw)$color = colrs[E(nw)$RegSign]
  
  plot(nw, edge.arrow.size=.4, layout = layout_with_fr, vertex.size = 7, vertex.label = NA)
  
  ## Histogram of the in- and out-degree
  #hist(degree(nw, which(V(nw)$nodetype != "target"), mode = "out"), main = "Out-degree of regulators", xlab = "Number of targets per regulator", ylab = "Regulators")
  print(ggbarplot(degree(nw, which(V(nw)$nodetype != "target"), mode = "out"), "Out-degree of regulators", "Number of targets per regulator", "Regulators"))
  
  #hist(degree(nw, mode = "in"), main = "In-degree of genes", xlab = "Number of regulators", ylab = "Targets")
  print(ggbarplot(degree(nw, mode = "in"), "In-degree of genes", "Number of regulators", "Targets"))
  
  ## Separate regulation by protein regulators and noncoding reglators
  nwPC = delete.edges(nw, E(nw)[E(nw)$RegBy == "NCreg"])
  #hist(degree(nwPC, mode = "in"), main = "In-degree of genes (only protein regulators)", xlab = "Number of protein regulators per target", ylab = "Targets")
  print(ggbarplot(degree(nwPC, mode = "in"), "In-degree of genes (only protein regulators)", "Number of protein regulators per target", "Targets"))
  #hist(degree(nwPC, which(V(nwPC)$nodetype == "PCreg"), mode = "out"), main = "Out-degree of protein regulators", xlab = "Number of targets per protein regulator", ylab = "Protein regulators")
  print(ggbarplot(degree(nwPC, which(V(nwPC)$nodetype == "PCreg"), mode = "out"), "Out-degree of protein regulators", "Number of targets per protein regulator", "Protein regulators"))
  
  nwNC = delete.edges(nw, E(nw)[E(nw)$RegBy == "PCreg"])
  #hist(degree(nwNC, mode = "in"), main = "In-degree of genes (only noncoding regulators)", xlab = "Number of noncoding regulators per target", ylab = "Targets")
  print(ggbarplot(degree(nwNC, mode = "in"), "In-degree of genes (only noncoding regulators)", "Number of noncoding regulators per target", "Targets"))
  #hist(degree(nwNC, which(V(nwNC)$nodetype == "NCreg"), mode = "out"), main = "Out-degree of noncoding regulators", xlab = "Number of targets per noncoding regulator", ylab = "Noncoding regulator")
  data = data.frame("values" = as.factor(degree(nwNC, which(V(nwNC)$nodetype == "NCreg"), mode = "out")))
  g = ggplot(data, aes(x = values))+geom_bar()+ggtitle("Out-degree of noncoding regulators") + xlab("Number of targets per noncoding regulator") + ylab("Noncoding regulator")
  print(ggbarplot(degree(nwNC, which(V(nwNC)$nodetype == "NCreg"), mode = "out"), "Out-degree of noncoding regulators", "Number of targets per noncoding regulator", "Noncoding regulator"))
  
  data = data.frame("protreg" = degree(nwPC, mode = "in"), "ncreg" = degree(nwNC, mode = "in"))  
  g = ggplot(data, aes(x = protreg, y = ncreg)) + geom_jitter() + ggtitle("Respective number of protein and noncoding regulators for each gene") + xlab("Number of protein regulators") + ylab("Number of noncoding regulators")
  print(g)
}

systemvisualize = function(mysystem){
  
  sys.nw = igraph::graph_from_data_frame(d = mysystem$edg, directed = T, vertices = mysystem$nod)
  hist(degree(sys.nw, mode = "out"), main = "Out-degree distribution of the global regulatory network", xlab = "Number of targets")
  hist(degree(sys.nw, mode = "in"), main = "In-degree distribution of the global regulatory network", xlab = "Number of regulators")
  
  for(n in c("TCRN", "TLRN", "RDRN", "PDRN")){
    hist(degree(mysystem[["RN.nw"]][[paste0(n,".nw")]], mode = "out"), main = paste("Out-degree distribution of the", n, "network", sep = " "), xlab = "Number of targets")
    hist(degree(mysystem[["RN.nw"]][[paste0(n,".nw")]], mode = "in"), main = paste("In-degree distribution of the", n, "network", sep = " "), xlab = "Number of regulators")
  }
  
}

# ------------------------------------------------------------------------------------------------------------ #
#                                     PARAMETERS FOR NETWORK GENERATION                                        #
# ------------------------------------------------------------------------------------------------------------ # 

## Constructor function for the insilicosystemargs class
## An object insilicosystemsargs contains a list of all arguments necessary for the generation of an in silico system
insilicosystemargs = function( ## ----
                       G = 50, ## G : number of genes in the system
                       ## protein-coding ratios
                       PC.p = 0.7, ## PC.p : ratio of protein-coding genes in the system
                       PC.TC.p = 0.4, ## PC.TC.p : ratio of regulators of transcription among the protein-coding genes
                       PC.TL.p = 0.3, ## PC.TL.p : ratio of regulators of translation among the protein-coding genes
                       PC.RD.p = 0.1, ## PC.RD.p : ratio of regulators of RNA decay among the protein-coding genes
                       PC.PD.p = 0.1, ## PC.PD.p : ratio of regulators of protein decay among the protein-coding genes
                       PC.PTM.p = 0.05, ## PC.PTM.p : ratio of regulators of protein post-translational modification among the protein-coding genes
                       PC.MR.p = 0.05, ## PC.MR.p : ratio of metabolic enzymes among the protein-coding genes
                       PC.PTM.form.p = 0.6, ## PC.PTM.form.p: for protein coding genes, ratio of protein having a PTM form
                       ## noncoding ratios
                       NC.TC.p = 0.3, ## NC.TC.p : ratio of regulators of transcription among the noncoding genes
                       NC.TL.p = 0.35, ## NC.TL.p : ratio of regulators of translation among the noncoding genes
                       NC.RD.p = 0.3, ## NC.RD.p : ratio of regulators of RNA decay among the noncoding genes
                       NC.PD.p = 0, ## NC.PD.p : ratio of regulators of protein decay among the noncoding genes
                       NC.PTM.p = 0.05, ## NC.PTM.p : ratio of regulators of protein post-translational modification among the noncoding genes
                       ## Ratio of postive/negative regulation
                       TC.PC.pos.p = 0.5, ## TC.PC.pos.p: probability that the transcription is positively regulated by protein regulators
                       TC.NC.pos.p = 0.5, ## TC.NC.pos.p: probability that the transcription is positively regulated by noncoding regulators 
                       TL.PC.pos.p = 0.5, ## TL.PC.pos.p: probability that the translation is positively regulated by protein regulators 
                       TL.NC.pos.p = 0.05, ## TL.NC.pos.p: probability that the translation is positively regulated by noncoding regulators 
                       # PTM.PC.pos.p = 0.5, ## PTM.PC.pos.p: probability that the protein decay is is activated by the post-translational modification by protein regulators
                       # PTM.NC.pos.p = 0.5, ## PTM.NC.pos.p: probability that the protein decay is is activated by the post-translational modification by noncoding regulators
                       #### Distribution of the different kinetic parameters of genes 
                       basal_transcription_rate_samplingfct = function(x){ runif(x, 0.01, 0.1) }, ## Function from which the transcription rates of genes are sampled (input x is the required sample size)
                       basal_translation_rate_samplingfct = function(x){ runif(x, 0.5, 5) }, ## Function from which the translation rates of genes are sampled (input x is the required sample size)
                       basal_RNAlifetime_samplingfct = function(x){ sample(60:3600, x, replace = T) }, ## Function from which the transcript lifetimes are sampled (input x is the required sample size)
                       basal_protlifetime_samplingfct = function(x){ sample(5400:14400, x, replace = T) }, ## Function from which the protein lifetime are sampled (input x is the required sample size)
                       #### TC reg. network properties
                       TC.PC.outdeg.distr = "powerlaw", ## Form of the distribution of the number of targets of transcription factors (can be either "powerlaw" or "exponential")
                       TC.NC.outdeg.distr = "powerlaw", ## Form of the the distribution of the number of targets of noncoding RNAs regulating transcription (can be either "powerlaw" or "exponential")
                       TC.PC.outdeg.exp = 3, ## exponent of the distribution for the out-degree of the protein regulators in the transcription graph
                       TC.NC.outdeg.exp = 1,   ## exponent of the distribution for the out-degree of the noncoding regulators in the transcription graph
                       TC.PC.indeg.distr = "powerlaw", ## Type of preferential attachment for the targets of protein regulators in the transcription graph
                       TC.NC.indeg.distr = "powerlaw", ## Type of preferential attachment for the targets of ncRNAs in the transcription graph
                       TC.PC.autoregproba = 0.2, ## Probability of protein regulators to perform autoregulation
                       TC.NC.autoregproba = 0, ## Probability of ncRNAs to perform autoregulation
                       TC.PC.twonodesloop = FALSE, ## Are 2-nodes loops authorised in the transcription network with protein regulators?
                       TC.NC.twonodesloop = FALSE, ## Are 2-nodes loops authorised in the transcription network with noncoding regulators?
                       TCbindingrate_samplingfct = function(x){ runif(x, 0.001, 0.01) }, ## Function from which the binding rate of regulators on target are sampled (input x is the required sample size)
                       TCunbindingrate_samplingfct = function(x){ runif(x, 0.001, 0.01) }, ## Function from which the unbinding rate of regulators from target are sampled (input x is the required sample size)
                       TCfoldchange_samplingfct = function(x){ sample(2:30, x, replace = T)  }, ## Function from which the fold change induced by a bound regulator are sampled (input x is the required sample size)
                       #### TL reg. network properties
                       TL.PC.outdeg.distr = "powerlaw", ## Form of the distribution of the number of targets of translation factors (can be either "powerlaw" or "exponential")
                       TL.NC.outdeg.distr = "powerlaw", ## Form of the the distribution of the number of targets of noncoding RNAs regulating translation (can be either "powerlaw" or "exponential")
                       TL.PC.outdeg.exp = 3, ## exponent of the distribution for the out-degree of the protein regulators in the translation graph
                       TL.NC.outdeg.exp = 1,   ## exponent of the distribution for the out-degree of the noncoding regulators in the translation graph
                       TL.PC.indeg.distr = "powerlaw", ## Type of preferential attachment for the targets of protein regulators in the translation graph
                       TL.NC.indeg.distr = "powerlaw", ## Type of preferential attachment for the targets of ncRNAs in the translation graph
                       TL.PC.autoregproba = 0.2, ## Probability of protein regulators to perform autoregulation
                       TL.NC.autoregproba = 0, ## Probability of ncRNAs to perform autoregulation
                       TL.PC.twonodesloop = FALSE, ## Are 2-nodes loops authorised in the translation network with protein regulators?
                       TL.NC.twonodesloop = FALSE, ## Are 2-nodes loops authorised in the translation network with noncoding regulators?
                       TLbindingrate_samplingfct = function(x){ runif(x, 0.001, 0.01) }, ## Function from which the binding rate of regulators on target are sampled (input x is the required sample size)
                       TLunbindingrate_samplingfct = function(x){ runif(x, 0.001, 0.01) }, ## Function from which the unbinding rate of regulators from target are sampled (input x is the required sample size)
                       TLfoldchange_samplingfct = function(x){ sample(2:30, x, replace = T)  }, ## Function from which the fold change induced by a bound regulator are sampled (input x is the required sample size)
                       #### RD reg. network properties
                       RD.PC.outdeg.distr = "powerlaw", ## Form of the distribution of the number of targets of RNA decay factors (can be either "powerlaw" or "exponential")
                       RD.NC.outdeg.distr = "powerlaw", ## Form of the the distribution of the number of targets of noncoding RNAs regulating RNA decay (can be either "powerlaw" or "exponential")
                       RD.PC.outdeg.exp = 3, ## exponent of the distribution for the out-degree of the protein regulators in the RNA decay graph
                       RD.NC.outdeg.exp = 1,   ## exponent of the distribution for the out-degree of the noncoding regulators in the RNA decay graph
                       RD.PC.indeg.distr = "powerlaw", ## Type of preferential attachment for the targets of protein regulators in the RNA decay graph
                       RD.NC.indeg.distr = "powerlaw", ## Type of preferential attachment for the targets of ncRNAs in the RNA decay graph
                       RD.PC.autoregproba = 0.2, ## Probability of protein regulators to perform autoregulation
                       RD.NC.autoregproba = 0, ## Probability of ncRNAs to perform autoregulation
                       RD.PC.twonodesloop = FALSE, ## Are 2-nodes loops authorised in the RNA decay network with protein regulators?
                       RD.NC.twonodesloop = FALSE, ## Are 2-nodes loops authorised in the RNA decay network with noncoding regulators?
                       RDbindingrate_samplingfct = function(x){ runif(x, 0.001, 0.01) }, ## Function from which the binding rate of regulators on target are sampled (input x is the required sample size)
                       #### PD reg. network properties
                       PD.PC.outdeg.distr = "powerlaw", ## Form of the distribution of the number of targets of protein decay factors (can be either "powerlaw" or "exponential")
                       PD.NC.outdeg.distr = "powerlaw", ## Form of the the distribution of the number of targets of noncoding RNAs regulating protein decay (can be either "powerlaw" or "exponential")
                       PD.PC.outdeg.exp = 3, ## exponent of the distribution for the out-degree of the protein regulators in the protein decay graph
                       PD.NC.outdeg.exp = 1,   ## exponent of the distribution for the out-degree of the noncoding regulators in the protein decay graph
                       PD.PC.indeg.distr = "powerlaw", ## Type of preferential attachment for the targets of protein regulators in the protein decay graph
                       PD.NC.indeg.distr = "powerlaw", ## Type of preferential attachment for the targets of ncRNAs in the protein decay graph
                       PD.PC.autoregproba = 0.2, ## Probability of protein regulators to perform autoregulation
                       PD.NC.autoregproba = 0, ## Probability of ncRNAs to perform autoregulation
                       PD.PC.twonodesloop = FALSE, ## Are 2-nodes loops authorised in the protein decay network with protein regulators?
                       PD.NC.twonodesloop = FALSE, ## Are 2-nodes loops authorised in the protein decay network with noncoding regulators?
                       PDbindingrate_samplingfct = function(x){ runif(x, 0.001, 0.01) }, ## Function from which the binding rate of regulators on target are sampled (input x is the required sample size)
                       #### PTM reg. network properties
                       PTM.PC.outdeg.distr = "powerlaw", ## Form of the distribution of the number of targets of protein post-translational modification factors (can be either "powerlaw" or "exponential")
                       PTM.NC.outdeg.distr = "powerlaw", ## Form of the the distribution of the number of targets of noncoding RNAs regulating protein post-translational modification (can be either "powerlaw" or "exponential")
                       PTM.PC.outdeg.exp = 3, ## exponent of the distribution for the out-degree of the protein regulators in the protein post-translational modification graph
                       PTM.NC.outdeg.exp = 1,   ## exponent of the distribution for the out-degree of the noncoding regulators in the protein post-translational modification graph
                       PTM.PC.indeg.distr = "powerlaw", ## Type of preferential attachment for the targets of protein regulators in the protein post-translational modification graph
                       PTM.NC.indeg.distr = "powerlaw", ## Type of preferential attachment for the targets of ncRNAs in the protein post-translational modification graph
                       PTM.PC.autoregproba = 0.2, ## Probability of protein regulators to perform autoregulation
                       PTM.NC.autoregproba = 0, ## Probability of ncRNAs to perform autoregulation
                       PTM.PC.twonodesloop = FALSE, ## Are 2-nodes loops authorised in the protein post-translational modification network with protein regulators?
                       PTM.NC.twonodesloop = FALSE, ## Are 2-nodes loops authorised in the protein post-translational modification network with noncoding regulators?
                       PTMbindingrate_samplingfct = function(x){ runif(x, 0.001, 0.01) }, ## Function from which the binding rate of regulators on target are sampled (input x is the required sample size)
                       ## Regulatory complexes
                       regcomplexes = 'prot', ## Can the regulators controlling a common target form regulatory complexes? can be 'none', 'prot' (only protein can form regulatory complexes) or 'both' (both RNAs and proteins can form regulatory complexes)
                       regcomplexes.p = 1, ## probability that regulators controlling a common target form regulatory complexes; ignore if regcomplexes = 'none'
                       regcomplexes.size = 2, ## number of components of a regulatory complex; ignore if regcomplexes = 'none'
                       complexesformationrate_samplingfct = function(x){ runif(x, 0.001, 0.01) }, ## Function from which the formation rate of regulatory complexes are sampled (input x is the required sample size)
                       complexesdissociationrate_samplingfct = function(x){ runif(x, 0.001, 0.01) } ## Function from which the dissociation rate of regulatory complexes are sampled (input x is the required sample size)
){
  
  NC.p = 1 - PC.p
  
  temp = sum(PC.TC.p + PC.TL.p + PC.RD.p + PC.PTM.p + PC.MR.p)
  PC.TC.p = PC.TC.p/temp
  PC.TL.p = PC.TL.p/temp
  PC.RD.p = PC.RD.p/temp
  PC.PD.p = PC.PD.p/temp
  PC.PTM.p = PC.PTM.p/temp
  PC.MR.p = PC.MR.p/temp
  
  temp = sum(NC.TC.p + NC.TL.p + NC.RD.p + NC.PTM.p)
  NC.TC.p = NC.TC.p/temp
  NC.TL.p = NC.TL.p/temp
  NC.RD.p = NC.RD.p/temp
  NC.PD.p = NC.PD.p/temp
  NC.PTM.p = NC.PTM.p/temp
  
  ## There cannot be negative regulation for transcript and protein decay
  RD.PC.pos.p = 1 ## RD.PC.pos.p: probability that the RNA decay is positively regulated by protein regulators (faster decay)
  RD.NC.pos.p = 1 ## RD.NC.pos.p: probability that the RNA decay is positively regulated by noncoding regulators (faster decay)
  PD.PC.pos.p = 1 ## PD.PC.pos.p: probability that the protein decay is positively regulated by protein regulators (faster decay)
  PD.NC.pos.p = 1 ## PD.NC.pos.p: probability that the protein decay is positively regulated by noncoding regulators (faster decay)
  
  value = list(  "G" = G,
                 "PC.p" = PC.p,
                 "PC.TC.p" = PC.TC.p,
                 "PC.TL.p" = PC.TL.p,
                 "PC.RD.p" = PC.RD.p,
                 "PC.PD.p" = PC.PD.p,
                 "PC.PTM.p" = PC.PTM.p,
                 "PC.MR.p" = PC.MR.p,
                 "PC.PTM.form.p" = PC.PTM.form.p,
                 "NC.p" = NC.p,
                 "NC.TC.p" = NC.TC.p,
                 "NC.TL.p" = NC.TL.p,
                 "NC.RD.p" = NC.RD.p,
                 "NC.PD.p" = NC.PD.p,
                 "NC.PTM.p" = NC.PTM.p,
                 "TC.PC.pos.p" = TC.PC.pos.p,
                 "TC.NC.pos.p" = TC.NC.pos.p,
                 "TL.PC.pos.p" = TL.PC.pos.p,
                 "TL.NC.pos.p" = TL.NC.pos.p,
                 "RD.PC.pos.p" = RD.PC.pos.p,
                 "RD.NC.pos.p" = RD.NC.pos.p,
                 "PD.PC.pos.p" = PD.PC.pos.p,
                 "PD.NC.pos.p" = PD.NC.pos.p,
                 # "PTM.PC.pos.p" = PTM.PC.pos.p,
                 # "PTM.NC.pos.p" = PTM.NC.pos.p,
                 "basal_transcription_rate_samplingfct" = basal_transcription_rate_samplingfct,
                 "basal_translation_rate_samplingfct" = basal_translation_rate_samplingfct,
                 "basal_RNAlifetime_samplingfct" = basal_RNAlifetime_samplingfct,
                 "basal_protlifetime_samplingfct" = basal_protlifetime_samplingfct,
                 "TC.PC.outdeg.distr" = TC.PC.outdeg.distr,
                 "TC.NC.outdeg.distr" = TC.NC.outdeg.distr,
                 "TC.PC.outdeg.exp" = TC.PC.outdeg.exp, 
                 "TC.NC.outdeg.exp" = TC.NC.outdeg.exp,
                 "TC.PC.indeg.distr" = TC.PC.indeg.distr,
                 "TC.NC.indeg.distr" = TC.NC.indeg.distr,
                 "TC.PC.autoregproba" = TC.PC.autoregproba,
                 "TC.NC.autoregproba" = TC.NC.autoregproba,
                 "TC.PC.twonodesloop" = TC.PC.twonodesloop,
                 "TC.NC.twonodesloop" = TC.NC.twonodesloop,
                 "TCbindingrate_samplingfct" = TCbindingrate_samplingfct,
                 "TCunbindingrate_samplingfct" = TCunbindingrate_samplingfct,
                 "TCfoldchange_samplingfct" = TCfoldchange_samplingfct,
                 "TL.PC.outdeg.distr" = TL.PC.outdeg.distr,
                 "TL.NC.outdeg.distr" = TL.NC.outdeg.distr,
                 "TL.PC.outdeg.exp" = TL.PC.outdeg.exp,
                 "TL.NC.outdeg.exp" = TL.NC.outdeg.exp,
                 "TL.PC.indeg.distr" = TL.PC.indeg.distr,
                 "TL.NC.indeg.distr" = TL.NC.indeg.distr,
                 "TL.PC.autoregproba" = TL.PC.autoregproba,
                 "TL.NC.autoregproba" = TL.NC.autoregproba,
                 "TL.PC.twonodesloop" = TL.PC.twonodesloop,
                 "TL.NC.twonodesloop" = TL.NC.twonodesloop,
                 "TLbindingrate_samplingfct" = TLbindingrate_samplingfct,
                 "TLunbindingrate_samplingfct" = TLunbindingrate_samplingfct,
                 "TLfoldchange_samplingfct" = TLfoldchange_samplingfct,
                 "RD.PC.outdeg.distr" = RD.PC.outdeg.distr,
                 "RD.NC.outdeg.distr" = RD.NC.outdeg.distr,
                 "RD.PC.outdeg.exp" = RD.PC.outdeg.exp,
                 "RD.NC.outdeg.exp" = RD.NC.outdeg.exp, 
                 "RD.PC.indeg.distr" = RD.PC.indeg.distr,
                 "RD.NC.indeg.distr" = RD.NC.indeg.distr,
                 "RD.PC.autoregproba" = RD.PC.autoregproba,
                 "RD.NC.autoregproba" = RD.NC.autoregproba,
                 "RD.PC.twonodesloop" = RD.PC.twonodesloop,
                 "RD.NC.twonodesloop" = RD.NC.twonodesloop,
                 "RDbindingrate_samplingfct" = RDbindingrate_samplingfct,
                 "PD.PC.outdeg.distr" = PD.PC.outdeg.distr,
                 "PD.NC.outdeg.distr" = PD.NC.outdeg.distr,
                 "PD.PC.outdeg.exp" = PD.PC.outdeg.exp,
                 "PD.NC.outdeg.exp" = PD.NC.outdeg.exp,
                 "PD.PC.indeg.distr" = PD.PC.indeg.distr,
                 "PD.NC.indeg.distr" = PD.NC.indeg.distr,
                 "PD.PC.autoregproba" = PD.PC.autoregproba,
                 "PD.NC.autoregproba" = PD.NC.autoregproba,
                 "PD.PC.twonodesloop" = PD.PC.twonodesloop,
                 "PD.NC.twonodesloop" = PD.NC.twonodesloop,
                 "PDbindingrate_samplingfct" = PDbindingrate_samplingfct,
                 "PTM.PC.outdeg.distr" = PTM.PC.outdeg.distr,
                 "PTM.NC.outdeg.distr" = PTM.NC.outdeg.distr,
                 "PTM.PC.outdeg.exp" = PTM.PC.outdeg.exp,
                 "PTM.NC.outdeg.exp" = PTM.NC.outdeg.exp,
                 "PTM.PC.indeg.distr" = PTM.PC.indeg.distr,
                 "PTM.NC.indeg.distr" = PTM.NC.indeg.distr,
                 "PTM.PC.autoregproba" = PTM.PC.autoregproba,
                 "PTM.NC.autoregproba" = PTM.NC.autoregproba,
                 "PTM.PC.twonodesloop" = PTM.PC.twonodesloop,
                 "PTM.NC.twonodesloop" = PTM.NC.twonodesloop,
                 "PTMbindingrate_samplingfct" = PTMbindingrate_samplingfct,
                 "regcomplexes" = regcomplexes,
                 "regcomplexes.p" = regcomplexes.p,
                 "regcomplexes.size" = regcomplexes.size ,
                 "complexesformationrate_samplingfct" = complexesformationrate_samplingfct, 
                 "complexesdissociationrate_samplingfct" = complexesdissociationrate_samplingfct)
  
  
  attr(value, "class") = "insilicosystemargs"
  
  return(value)
  
}  ##----


# ------------------------------------------------------------------------------------------------------------ #
#                                   PARAMETERS FOR INDIVIDUAL GENERATION                                       #
# ------------------------------------------------------------------------------------------------------------ # 

## Constructor function for the insilicosystemargs class
## An object insilicosystemsargs contains a list of all arguments necessary for the generation of an in silico system
insilicoindividualargs = function( ## ----
                            ploidy = 4, ## ploidy = number of alleles for each gene
                            gcnname = "GCN", ## gcnname = name to give to each allele version
                            ngenevariants = 5, ## number of alleles existing for each gene
                            qtleffect_samplingfct = function(x){rtruncnorm(x, a = 0, b = Inf, mean = 1, sd = 0.1)} ## function from which is sampled the effect of a QTL (input x is the required sample size)
){
  
  gcnList = sapply(1:ploidy, function(x){paste0(gcnname, x)})
 
  value = list("ploidy" = ploidy,
               "gcnname" = gcnname,
               "gcnList" = gcnList,
               "ngenevariants" = ngenevariants,
               "qtleffect_samplingfct" = qtleffect_samplingfct)
  
  
  attr(value, "class") = "insilicoindividualargs"
  
  return(value)
  
}  ##----


# ------------------------------------------------------------------------------------------------------------ #
#                                      GENERATE VARIANTS FOR THE GENES                                         #
# ------------------------------------------------------------------------------------------------------------ # 

## Function that generates the QTL effects of each gene variant 
## Input
##    - nod: the data.frame
createVariants = function(genes, indargs){
  
  G = nrow(genes)
  
  variants = vector("list", G)
  names(variants) = genes$id
  
  ## names of the qtl effects
  ## The first 5 are the qtl affecting all genes, the last 4 only affect protein coding genes
  qtlnames = c("qtlTCrate", "qtlRDrate", "qtlTCregbind", "qtlRDbindreg", "qtlactivity", "qtlTLrate", "qtlPDrate", "qtlTLregbind", "qtlPDregbind")
    
  if(indargs$ngenevariants > 1){  
    for(i in genes$id){
      potentialqtls = 1:(5 + 4*(genes[i, "coding"] == "PC")) ## protein-coding genes can have mutations affecting more parameters compared to noncoding genes
      nchanges = sample(potentialqtls, indargs$ngenevariants - 1, replace = T) ## Select the number of "mutations" from the original allele for each variant (minimum 1 otherwise we would have several copies of the original allele)
                                                                               ## the 1st variant is the original allele - no mutation
      qtlchanges = unlist(sapply(1:(indargs$ngenevariants-1), function(x){ length(qtlnames)*x + sample(potentialqtls, nchanges[x], replace = F)})) ## sample which qtl are affected by mutations for each variant, and convert it into matrix coordinates
      ## For each gene, the variants are stocked in the form of a matrix, rows being the different potential qtls and columns being the different gene variaents. Element i,j correspond to the effect of mutation at QTL i for variant j
      ## The matrix is filled with 1 (no effect, allele identical to the "original allele", and 0 for QTLs that do not exist in the gene (e.g. QTL affecting the translation rate for noncoding genes))
      temp = matrix(rep(c(1, 0), indargs$ngenevariants * c(length(potentialqtls), length(qtlnames) - length(potentialqtls))), byrow = T, nrow = length(qtlnames), ncol = indargs$ngenevariants, dimnames = list(qtlnames, 1:indargs$ngenevariants))
      temp[qtlchanges] = indargs$qtleffect_samplingfct(sum(nchanges))
      variants[[i]] = temp
    }
  }else{
    for(i in genes$id){
      variants[[i]] = matrix(rep(1, length(qtlnames)), byrow = T, nrow = length(qtlnames), ncol = 1, dimnames = list(qtlnames, 1))
    }
  }
  
  return(variants)
  
}

createIndividual = function(variantsList, indargs){

  G = length(variantsList)
  QTLeffects = vector("list", indargs$ploidy)
  names(QTLeffects) = indargs$gcnList
  ## individualvariants: data frame where rows are genes and columns are "copy number" ids i.e. each columns represent a homolog chromosom
  ## Element i, j in the data frame corresponds to the variant of gene i present in the homolog chromosom j
  individualvariants = as.data.frame(matrix(sample(1:indargs$ngenevariants, G*indargs$ploidy, replace = T), nrow = G, ncol = indargs$ploidy))
  names(individualvariants) = indargs$gcnList
  qtlnames = c("qtlTCrate", "qtlRDrate", "qtlTCregbind", "qtlRDbindreg", "qtlactivity", "qtlTLrate", "qtlPDrate", "qtlTLregbind", "qtlPDregbind")
  
  ## Work for each gene copy (here in the sense each homolog chromosome)
  for(gcn in indargs$gcnList){
    QTLeffects[[gcn]] = vector("list", length(qtlnames))
    names(QTLeffects[[gcn]]) = qtlnames
    for(q in qtlnames){
      for(g in 1:G){
        QTLeffects[[gcn]][[q]][g] = variantsList[[g]][q, individualvariants[g, gcn]]
      }
    }
  }
  
  value = list("QTLeffects" = QTLeffects, "haplotype" = individualvariants)
  attr(value, "class") = "insilicoindividual"
  
  return(value)
}


createPopulation = function(nind, insilicosystem, indargs){
  
  genvariants = createVariants(insilicosystem$genes, indargs)
  indnames = sapply(1:nind, function(x){paste0("Ind", x)})
  individualsList = vector("list", nind)
  names(individualsList) = indnames
  
  for(i in indnames){
    individualsList[[i]] = createIndividual(genvariants, indargs)
  }
  
  value = list("GenesVariants" = genvariants, "individualsList" = individualsList, "indargs" = indargs)
  
}

# ------------------------------------------------------------------------------------------------------------ #
#                                     GENERATE GENES FOR IN SILICO SYSTEM                                      #
# ------------------------------------------------------------------------------------------------------------ # 

## Generate the genes in the system and their attributes, according to the user parameters
## Inputs:
##  - nod: data frame created by the function createGenes
##  - sysargs: an object of class insiliosystemargs, i.e. list of all parameters for in silico network generation
## Outputs:
##    - nod: a data frame of genes (rows) and their attributes
createGenes = function(sysargs){
  
  ## G.name : genes name 
  ##G.nameid = sapply(1:G, function(i){paste0("G_",i)})
  
  ## nod is the vertices data frame (1st column "id" = an integer value (faster for computation), 2nd column "name", 3rd column "coding" (values "NC" or "PC"),  
  ##  4rd column "TargetReaction" (values "TC", "TL", "RD", "PTM", "PD", "MR"), next columns: kinetic parameters (transcription rate, translation rate, RNA decay rate, protein decay rate))
  nod = data.frame("id" = 1:sysargs[["G"]], "coding" = rep("", sysargs[["G"]]), "TargetReaction" = rep("", sysargs[["G"]]),  "PTMform" = rep("", sysargs[["G"]]), "ActiveForm" = rep("", sysargs[["G"]]),
                   "TCrate" = rep(0,sysargs[["G"]]), "TLrate" = rep(0,sysargs[["G"]]), "RDrate" = rep(0,sysargs[["G"]]), "PDrate" = rep(0,sysargs[["G"]]), stringsAsFactors = F)
  rownames(nod) = nod$id
  
  ## Deciding gene status
  nod$coding = sample(c("PC", "NC"), sysargs[["G"]], prob = c(sysargs[["PC.p"]], sysargs[["NC.p"]]), replace = T)
  
  ## Deciding gene function (reaction to be regulated)
  nod$TargetReaction[nod$coding == "PC"] = sample(c("TC", "TL", "RD", "PD", "PTM", "RD"), sum(nod$coding == "PC"), prob = c(sysargs[["PC.TC.p"]], sysargs[["PC.TL.p"]], sysargs[["PC.RD.p"]], sysargs[["PC.PD.p"]], sysargs[["PC.PTM.p"]], sysargs[["PC.MR.p"]]), replace = T)
  nod$TargetReaction[nod$coding == "NC"] = sample(c("TC", "TL", "RD", "PD", "PTM"), sum(nod$coding == "NC"), prob = c(sysargs[["NC.TC.p"]], sysargs[["NC.TL.p"]], sysargs[["NC.RD.p"]], sysargs[["NC.PD.p"]], sysargs[["NC.PTM.p"]]), replace = T)
  
  ## Choose which proteins (from protein-coding genes) have a PTM form
  nod$PTMform[nod$coding == "PC"] = sample(c("1","0"), sum(nod$coding == "PC"), prob = c(sysargs[["PC.PTM.form.p"]], 1-sysargs[["PC.PTM.form.p"]]), replace = T)
  nod$PTMform[nod$coding == "NC"] = "0"
  
  ## In nod, state what is the active form of each gene, i.e. which form (i.e. RNA, protein, activated protein) is performing the regulation
  nod$ActiveForm[nod$coding == "NC"] = "R" ## noncoding genes act through their RNA
  nod$ActiveForm[nod$coding == "PC" & nod$PTMform == "0"] = "P" ## protein-coding genes act through their protein
  nod$ActiveForm[nod$coding == "PC" & nod$PTMform == "1"] = "Pm" ## For proteins undergoing a post-translational modification, only the PTM form is active
  nod$ActiveForm = sapply(1:nrow(nod), function(x){paste0(nod$ActiveForm[x],".",nod$id[x])})
  
  ## Sample the kinetic parameters of the genes
  
  ## Transcription rate: applicable to all genes
  nod$TCrate = sysargs[["basal_transcription_rate_samplingfct"]](sysargs[["G"]])
  
  ## RNA decay rate: applicable to all genes
  ## Sample the lifetime, the decay rate is defined as 1/lifetime
  nod$RDrate = 1/sysargs[["basal_RNAlifetime_samplingfct"]](sysargs[["G"]])
  
  ## Translation rate: applicable to protein-coding genes
  nod$TLrate[nod$coding == "PC"] = sysargs[["basal_translation_rate_samplingfct"]](sum(nod$coding == "PC"))
  
  ## Protein coding rate: applicable to protein-coding genes
  nod$PDrate[nod$coding == "PC"] = 1/sysargs[["basal_protlifetime_samplingfct"]](sum(nod$coding == "PC"))
  
  return(nod)
  
}


# ------------------------------------------------------------------------------------------------------------ #
#                                       GENERATE A REGULATORY NETWORK                                          #
# ------------------------------------------------------------------------------------------------------------ # 

## Inputs:
##  - regsList: a list of vector with regulator ids. 2 elements: 1st vector is for protein-coding regulators, the 2nd for noncoding regulators (as each can have a different set of target + have a different out-degree distribution)
##  - tarsList: a list of vector with target ids. 2 elements: 1st vector is for targets of protein-coding regulators, the 2nd for targets of noncoding regulators (as each type of regulator can have a different set of target + have a different out-degree distribution)
##  - which reaction is regulated? (used to retrieve automatically the variables)
##  - the data frame of nodes in the system
##  - sysargs: arguments of the system
## Outputs:
##  - nod: data frame of nodes (and ther attributes) in the network
##  - edg: data frame of edges (and their attributes)
createRegulatoryNetwork = function(regsList, tarsList, reaction, nod, sysargs){
  
  ## Construct the regulatory network nodes and edges data.frame
  nwnod = nod[nod$id %in% c(unlist(regsList), unlist(tarsList)),]
  nwnod = data.frame(nwnod, "nodetype" = rep("target", nrow(nwnod)), stringsAsFactors = F) # add a node attribute specifying if the node is a target, a protein regulator or a noncoding regulator (regulators can be also targets but will be labeled as regulators)
  nwnod[nwnod$id %in% regsList[["PC"]], "nodetype"] = "PCreg"
  nwnod[nwnod$id %in% regsList[["NC"]], "nodetype"] = "NCreg"
  
  ## Call the julia function nwgeneration to generate the regulatory network where the protein regulators are the regulatory nodes
  edgPC = juliaGet(juliaCall("nwgeneration", regsList[["PC"]], tarsList[["PC"]], sysargs[[paste(reaction, "PC", "indeg.distr", sep = ".")]], sysargs[[paste(reaction, "PC", "outdeg.distr", sep = ".")]], sysargs[[paste(reaction, "PC", "outdeg.exp", sep = ".")]], sysargs[[paste(reaction, "PC", "autoregproba", sep = ".")]], sysargs[[paste(reaction, "PC", "twonodesloop", sep = ".")]]))
  
  ## Call the julia function nwgeneration to generate the regulatory network where the noncoding regulators are the regulatory nodes
  edgNC = juliaGet(juliaCall("nwgeneration", regsList[["NC"]], tarsList[["NC"]], sysargs[[paste(reaction, "NC", "indeg.distr", sep = ".")]], sysargs[[paste(reaction, "NC", "outdeg.distr", sep = ".")]], sysargs[[paste(reaction, "NC", "outdeg.exp", sep = ".")]], sysargs[[paste(reaction, "NC", "autoregproba", sep = ".")]], sysargs[[paste(reaction, "NC", "twonodesloop", sep = ".")]]))
  
  ## create the edge dataframe
  nwedg = data.frame("from" = c(edgPC[,1], edgNC[,1]), "to" = c(edgPC[,2], edgNC[,2]), "TargetReaction" = rep(reaction,nrow(edgPC)+nrow(edgNC)), "RegSign" = rep("",nrow(edgPC)+nrow(edgNC)), "RegBy" = rep(c("PC", "NC"), c(nrow(edgPC), nrow(edgNC))), stringsAsFactors = F)
  
  ## Choose the sign (activation or repression) of each regulation (=edge)
  ## First for the regulatory interactions exerted by protein regulators
  if(nrow(edgPC) > 0) nwedg$RegSign[1:nrow(edgPC)] = sample(c("1","-1"), nrow(edgPC), prob = c(sysargs[[paste(reaction, "PC", "pos.p", sep = ".")]], 1 - sysargs[[paste(reaction, "PC", "pos.p", sep = ".")]]), replace = T)
  ## Then for the regulatory interactions exerted by noncoding regulators
  if(nrow(edgNC) > 0) nwedg$RegSign[(nrow(edgPC)+1):(nrow(edgPC) + nrow(edgNC))] = sample(c("1","-1"), nrow(edgNC), prob = c(sysargs[[paste(reaction, "NC", "pos.p", sep = ".")]], 1 - sysargs[[paste(reaction, "NC", "pos.p", sep = ".")]]), replace = T) 
  
  ## Create corresponding igraph object (to save the regulatory interactions as is, before the creation of the combinatorial regulation)
  nw = igraph::graph_from_data_frame(d = nwedg, directed = T, vertices = nwnod)
  
  ## Creation of combinatorial regulation
  ## if regcomplexes != 'none', if several regulators control a common target they can form regulatory complexes
  ## The composition of each complex is stored in complexes
  ## Complexes can be composed only of proteins if recomplexes = "prot" or protein and noncoding regulators if regcomplexes = "both"
  
  if(sysargs[["regcomplexes"]] == "none"){ ## If regulators controlling a same target are not allowed to form a regulatory complex, simply reformat the edg and nod dataframes
    
    nwedgcomp = nwedg
    nwedgcomp$from = sapply(nwedgcomp$from, toString) ## transform the id of regulators from integer to string (not the id of target because we need it to be integer for computational speed later)
    complexes = list()
    
  }else if(sysargs[["regcomplexes"]] == "prot"){ ## If the regulatory complexes can only be protein complexes
    
    temp = nwedg[nwedg$RegBy =="PC", ]
    tempregcom = juliaGet(juliaCall("combreg", temp$from, temp$to, temp$RegSign, sysargs[["regcomplexes.p"]], sysargs[["regcomplexes.size"]], reaction))
    ## only keep the noncoding regulators (the regulation from protein-coding regulators is given by the Julia function combreg)
    nwedgcomp = nwedg[nwedg$RegBy =="NC", c("from", "to", "TargetReaction", "RegSign")]
    nwedgcomp = rbind(nwedgcomp, data.frame("from" = unlist(tempregcom$newedg[,1]), "to" = unlist(tempregcom$newedg[,2]), "TargetReaction" = rep(reaction, nrow(tempregcom$newedg)), "RegSign" = unlist(tempregcom$newedg[,3]), stringsAsFactors = F))
    complexes = lapply(tempregcom$Complexes, unlist)
    
  }else if(sysargs[["regcomplexes"]] == "both"){ ## If the regulatory complexes can be protein/noncoding complexes
    
    tempregcom = juliaGet(juliaCall("combreg", nwedg$from, nwedg$to, nwedg$RegSign, regcomplexes.p, regcomplexes.size, reaction))
    nwedgcomp = data.frame("from" = unlist(tempregcom$newedg[,1]), "to" = unlist(tempregcom$newedg[,2]), "TargetReaction" = rep(reaction, nrow(tempregcom$newedg)), "RegSign" = unlist(tempregcom$newedg[,3]), stringsAsFactors = F)
    complexes = lapply(tempregcom$Complexes, unlist)
    
  }
  
  nwedg = nwedg[order(nwedg$to),]
  nwedgcomp = nwedgcomp[order(nwedgcomp$to),]
  rownames(nwedg) = NULL
  rownames(nwedgcomp) = NULL
  return(list("edg" = nwedg, "edgcomp" = nwedgcomp, "complexes" = complexes, "igraph" = nw))
  
}


# ------------------------------------------------------------------------------------------------------------ #
#                          GENERATE THE MULTI-OMIC NETWORK (LIST OF REGULATORY NETWORK)                        #
# ------------------------------------------------------------------------------------------------------------ # 

## This function generates the different regulatory networks, one for each gene expression step potentially targeted for regulation
## Inputs:
##  - nod: data frame created by the function createGenes
##  - sysargs: an object of class insiliosystemargs, i.e. list of all parameters for in silico network generation
## Outputs:
##    -
createMultiOmicNetwork = function(nod, sysargs){
  

  ##edg is the edges data frame (1st column "from", 2nd column "to", 3rd column "TargetReaction" (values "TC", "TL", "RD", "PTM", "PD", "MR"), 4th column "RegSign" (value +1 or -1))
  edg = data.frame("from" = NULL, "to" = NULL, "TargetReaction" = NULL, "RegSign" = NULL, stringsAsFactors = F)
  
  ## complexes is a list where each element gives the components of a regulatory complex, whose id is the name of the element of the list
  complexes = list()
  
  
  ####  Define Transcriptional regulatory network (TCRN)

  ## Identify TFs in the system
  PCreg.id = nod$id[nod$coding == "PC" & nod$TargetReaction == "TC"]
  # Identify non-coding genes regulating transcription in the system
  NCreg.id = nod$id[nod$coding == "NC" & nod$TargetReaction == "TC"]
  
  ## which genes are targeted by TFs? here any gene (including the TFs because a TF can autoregulate himself)
  PCtarget.id = nod$id
  ## which genes are targeted by noncoding RNAs? here any protein-coding gene
  NCtarget.id = nod$id[nod$coding == "PC"]
  
  ## Construct the regulatory network
  TCRN = createRegulatoryNetwork(regsList = list("PC" = PCreg.id, "NC" = NCreg.id), tarsList = list("PC" = PCtarget.id, "NC" = NCtarget.id), reaction = "TC", nod = nod, sysargs)
  TCRN.edg = TCRN[["edgcomp"]]

  complexes = c(complexes, TCRN[["complexes"]])
  
  ## Sample the kinetic parameters of each regulatory interaction
  ##    Kinetic parameters for transcription regulation include the binding and unbinding rate of regulators to gene promoter, and the fold change induced on transcription rate by a regulator bound to the promoter
  TCRN.edg = data.frame(TCRN.edg, "TCbindingrate" = sysargs[["TCbindingrate_samplingfct"]](nrow(TCRN.edg)), "TCunbindingrate" = sysargs[["TCunbindingrate_samplingfct"]](nrow(TCRN.edg)), "TCfoldchange" = rep(0, nrow(TCRN.edg)), stringsAsFactors = F)
  TCRN.edg$TCfoldchange[TCRN.edg$RegSign == "1"] = sysargs[["TCfoldchange_samplingfct"]](sum(TCRN.edg$RegSign == "1"))  ## Repressors induce a fold change of 0
  

  #### Define Translational regulatory network (TLRN) 

  ## Identify TLFs in the system
  PCreg.id = nod$id[nod$coding == "PC" & nod$TargetReaction == "TL"]
  # Identify non-coding genes regulating translation in the system
  NCreg.id = nod$id[nod$coding == "NC" & nod$TargetReaction == "TL"]
  
  ## which genes are targeted by TLFs? here any protein-coding gene
  PCtarget.id = nod$id[nod$coding == "PC"]
  ## which genes are targeted by noncoding RNAs? here any protein-coding gene 
  NCtarget.id = nod$id[nod$coding == "PC"]
  
  ## Construct the regulatory network
  TLRN = createRegulatoryNetwork(regsList = list("PC" = PCreg.id, "NC" = NCreg.id), tarsList = list("PC" = PCtarget.id, "NC" = NCtarget.id), reaction = "TL", nod = nod, sysargs = sysargs)
  TLRN.edg = TLRN[["edgcomp"]]

  complexes = c(complexes, TLRN[["complexes"]])
  
  
  ## Sample the kinetic parameters of each regulatory interaction
  ##    Kinetic parameters for translation regulation include the binding and unbinding rate of regulators to mRNA binding sequence, and the fold change induced on transcription rate by a regulator bound to the mRNA
  TLRN.edg = data.frame(TLRN.edg, "TLbindingrate" = sysargs[["TLbindingrate_samplingfct"]](nrow(TLRN.edg)), "TLunbindingrate" = sysargs[["TLunbindingrate_samplingfct"]](nrow(TLRN.edg)), "TLfoldchange" = rep(0, nrow(TLRN.edg)), stringsAsFactors = F)
  TLRN.edg$TLfoldchange[TLRN.edg$RegSign == "1"] = sysargs[["TLfoldchange_samplingfct"]](sum(TLRN.edg$RegSign == "1")) ## Repressors induce a fold change of 0
  
  
  #### Define RNA decay regulatory network (RDRN)

  ## Identify proteins regulating RNA decay in the system
  PCreg.id = nod$id[nod$coding == "PC" & nod$TargetReaction == "RD"]
  ## Identify noncoding RNAs regulating RNA decay (miRNAs or siRNAs for ex)
  NCreg.id = nod$id[nod$coding == "NC" & nod$TargetReaction == "RD"]
  
  ## which genes are targeted by coding regulators? here any gene (temporary)
  PCtarget.id = nod$id
  ## which genes are targeted by noncoding regulators? here any gene (temporary)
  NCtarget.id = nod$id
  
  ## Construct the regulatory network
  RDRN = createRegulatoryNetwork(regsList = list("PC" = PCreg.id, "NC" = NCreg.id), tarsList = list("PC" = PCtarget.id, "NC" = NCtarget.id), reaction = "RD", nod = nod, sysargs = sysargs)
  RDRN.edg = RDRN[["edgcomp"]]

  complexes = c(complexes, RDRN[["complexes"]])
  
  ## Sample the kinetic parameters of each regulatory interaction
  ##    Kinetic parameters for RNA decay includes the binding (and unbinding for repressors of decay) rate of the regulator on the RNA 
  RDRN.edg = data.frame(RDRN.edg, "RDbindingrate" = sysargs[["RDbindingrate_samplingfct"]](nrow(RDRN.edg)), stringsAsFactors = F)

  
  #### Define protein decay regulatory network (PDRN)

  ## Identify proteins regulating RNA decay in the system
  PCreg.id = nod$id[nod$coding == "PC" & nod$TargetReaction == "PD"]
  ## Identify noncoding RNAs regulating RNA decay (miRNAs or siRNAs for ex)
  NCreg.id = nod$id[nod$coding == "NC" & nod$TargetReaction == "PD"]
  
  ## which genes are targeted by coding regulators? here any gene (temporary)
  PCtarget.id = nod$id[nod$coding == "PC"]
  ## which genes are targeted by noncoding regulators? here any gene (temporary)
  NCtarget.id = nod$id[nod$coding == "PC"]
  
  ## Construct the regulatory network
  PDRN = createRegulatoryNetwork(regsList = list("PC" = PCreg.id, "NC" = NCreg.id), tarsList = list("PC" = PCtarget.id, "NC" = NCtarget.id), reaction = "PD", nod = nod, sysargs = sysargs)
  PDRN.edg = PDRN[["edgcomp"]]

  complexes = c(complexes, PDRN[["complexes"]])  

  ## Sample the kinetic parameters of each regulatory interaction
  ##    Kinetic parameters for protein decay includes the binding (and unbinding for repressors of decay) rate of the regulator on the protein 
  PDRN.edg = data.frame(PDRN.edg, "PDbindingrate" = sysargs[["PDbindingrate_samplingfct"]](nrow(PDRN.edg)), stringsAsFactors = F)

  
  #### Define protein post-translational modification regulatory network (PTMRN) ----

  PTMRN.edg = data.frame("from" = character(), "to" = integer(), "TargetReaction" = character(), "RegSign" = character(), "PTMbindingrate" = numeric(), stringsAsFactors = F)
  PTMRN.nw = NULL
  PTMRN = list("edg" = data.frame("from" = character(), "to" = integer(), "TargetReaction" = character(), "RegSign" = character(), "PTMbindingrate" = numeric(), "RegBy" = character(), stringsAsFactors = F))
  
  
  #### Define regulatory complexes kinetic parameters

  complexeskinetics = list()
  if(length(complexes)>0){
    formrates = sysargs[["complexesformationrate_samplingfct"]](length(complexes))
    dissrates = sysargs[["complexesdissociationrate_samplingfct"]](length(complexes))
    for(c in 1:length(complexes)){
      complexeskinetics[[names(complexes)[c]]] = list("formationrate" = formrates[c], "dissociationrate" = dissrates[c])
    }}
  
  # -----------------------------------------------------------------
  ####                          RETURN                           ----
  # -----------------------------------------------------------------
  
  
  
  ## save all the interactions in edg
  temp = list(TCRN, TLRN, RDRN, PDRN, PTMRN)
  for(t in temp){
    edg = rbind(edg, t[["edg"]][,c("from", "to", "TargetReaction", "RegSign", "RegBy")])
  }
  
  ## Create the lists to be sent to Julia
  # temp = list("TCRN.edg", "TLRN.edg", "RDRN.edg", "PDRN.edg")
  # for(t in temp){
  #   new = list()
  #   for(cols in colnames(get(t))){
  #     new[[cols]] = get(t)[,cols]
  #   }
  #   assign(t, new)
  # }
  
  ## Return
  res = list("edg" = edg,
             "complexes" = complexes,
             "complexeskinetics" = complexeskinetics,
             "TCRN.edg" = TCRN.edg,
             "TLRN.edg" = TLRN.edg,
             "RDRN.edg" = RDRN.edg,
             "RDRN.edg" = RDRN.edg,
             "PDRN.edg" = PDRN.edg,
             "PTMRN.edg" = PTMRN.edg,
             "TCRN.nw" = TCRN[["nw"]],
             "TLRN.nw" = TLRN[["nw"]],
             "RDRN.nw" = RDRN[["nw"]],
             "RDRN.nw" = RDRN[["nw"]],
             "PDRN.nw" = PDRN[["nw"]],
             "PTMRN.nw" = PTMRN[["nw"]])
  
  return(res)
}


createInSilicoSystem = function(sysargs){
  
  genes = createGenes(sysargs)
  
  mosystem = createMultiOmicNetwork(genes, sysargs)
  
  value = list("sysargs" = sysargs, "genes" = genes, "mosystem" = mosystem)
  attr(value, "class") = "insilicosystem"
  
  return(value)
}


# ------------------------------------------------------------------------------------------------------------ #
#                   GENERATE THE LIST OF SPECIES AND REACTIONS FOR THE STOCHASTIC SIMULATION                   #
# ------------------------------------------------------------------------------------------------------------ # 

createStochSystem = function(insilicosystem, indargs){
  
  ## Create the network and regulatory complexes lists to be sent to Julia (converted to dictionaries in Julia)
  temp = list("TCRN.edg", "TLRN.edg", "RDRN.edg", "PDRN.edg", "PTMRN.edg")
  for(t in temp){
   new = list()
   for(cols in colnames(insilicosystem[["mosystem"]][[t]])){
     new[[cols]] = insilicosystem[["mosystem"]][[t]][,cols]
   }
  assign(t, new)
  }
  

  ## Create the gene list to be sent to Julia (converted to dictionaries in Julia)
  new = list()
  for(cols in colnames(insilicosystem["genes"])){
    new[[cols]] = insilicosystem["genes"][,cols]
  }
  assign("nod", new)
  
  stochsystem = juliaCall("generateReactionList", nod, TCRN.edg, TLRN.edg, RDRN.edg, PDRN.edg, PTMRN.edg, insilicosystem$mosystem$complexes, insilicosystem$mosystem$complexeskinetics, as.integer(insilicosystem$sysargs$regcomplexes.size), indargs$gcnList)
  
  species = unlist(juliaGet(juliaCall("getDictfromKey", stochsystem, "species")))
  reactions = unlist(juliaGet(juliaCall("getDictfromKey", stochsystem, "reactions")))
  reactionsnames = unlist(juliaGet(juliaCall("getDictfromKey", stochsystem, "reactionsnames")))
  
  return(list("JuliaObject" = stochsystem, "species" = species, "reactions" = reactions, "reactionsnames" = reactionsnames))
  
}

# howmanyautoreg(edg)
# howmanyloops(edg)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
# 
# ############################################################################################################################
# ############################################################################################################################
# 
# 
# ## Color differently TFs and other genes
# colrs = c("gold","tomato", "blue")
# V(TCRN.nw)$color = colrs[1 + (V(TCRN.nw)$name %in% PCreg.id) + (V(TCRN.nw)$name %in% NCreg.id)*2]
# 
# ## Color differently activating and inhibiting regulations
# # colrs = c("1" = "red", "-1" = "blue")
# # E(TCRN.nw)$color = colrs[E(TCRN.nw)$RegSign]
# 
# ## Visualisation
# plot(TCRN.nw, edge.arrow.size=.4, layout = layout_with_fr, vertex.size = 7, vertex.label = NA)
# 
# #fitPowerLaw(degree(TCRN.nw, mode = "out")[degree(TCRN.nw, mode = "out")!=0])
# hist(degree(TCRN.nw, mode = "out"))
# 
# #fitPowerLaw(degree(TCRN.nw, mode = "in"))
# hist(degree(TCRN.nw, mode = "in"))
# 
# 
# 
# ## Color differently TFs and other genes
# colrs = c("gold","green", "blue")
# V(TLRN.nw)$color = colrs[1 + (V(TLRN.nw)$name %in% PCreg.id) + (V(TLRN.nw)$name %in% NCreg.id)*2]
# 
# ## Color differently activating and inhibiting regulations
# # colrs = c("1" = "red", "-1" = "blue")
# # E(TCRN.nw)$color = colrs[E(TCRN.nw)$RegSign]
# 
# ## Visualisation
# plot(TLRN.nw, edge.arrow.size=.4, layout = layout_with_fr, vertex.size = 7, vertex.label = NA)
# 
# #fitPowerLaw(degree(TCRN.nw, mode = "out")[degree(TCRN.nw, mode = "out")!=0])
# hist(degree(TLRN.nw, mode = "out"))
# 
# #fitPowerLaw(degree(TCRN.nw, mode = "in"))
# hist(degree(TLRN.nw, mode = "in"))
# 
# 
# ##############################################################
# 
# reg1 = 1:50
# reg2 = 51:100
# target = 101:1100
# 
# edg = juliaGet(juliaCall("nwgeneration", reg1, target, "exponential", "powerlaw", 1.5, 0, FALSE))
# 
# ## Create corresponding igraph object
# nw = igraph::graph_from_data_frame(d = edg, directed = T, vertices =c(reg1, target)) 
# 
# plot(nw, edge.arrow.size=.4, layout = layout_with_fr, vertex.size = 7, vertex.label = NA)
# 
# hist(degree(nw, mode = "out"))
# hist(degree(nw, mode = "in"))
# 
# edg = juliaGet(juliaCall("nwgeneration", reg2, target, "powerlaw", "powerlaw", 1, 0, FALSE, edg))
# nw2 = igraph::graph_from_data_frame(d = edg, directed = T, vertices = c(reg1, reg2, target))
# 
# plot(nw2, edge.arrow.size=.4, layout = layout_with_fr, vertex.size = 7, vertex.label = NA)
# 
# hist(degree(nw2, mode = "out"))
# hist(degree(nw2, mode = "in"))
# 
# 
# nw3 = induced_subgraph(nw2, vids = c(reg2, target))
# hist(degree(nw3, mode = "out"))
# hist(degree(nw3, mode = "in"))
# 
# plot(degree(nw, v = (target-50), mode = "in"), degree(nw3, v = (target-50), mode = "in"))
# 
# ##############################################################
# 
# iter = 100
# 
# Edist = vector(length = iter)
# tarav = vector(length = iter)
# for(i in 1:iter){
#   out = juliaGet(juliaCall("samplepowerlaw", 124, 1.4, 4410))
#   Edist[i] = sum(out)
#   tarav[i] = mean(out)
# }
# 
# 
# 
# hist(Edist)
# hist(tarav)
# hist(2*Edist/(157+4410))
