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

## Test if Julia is installed on the computer
if(!findJulia(test = T)) stop("Julia is not installed on the computer or not accessible by R. Check that Julia is correcly installed and/or in the PATH variable\n")


##########################################################################################################################
###                                               PARAMETERS                                                           ###
##########################################################################################################################


## Temporary - load parameters from param_nw.R
setwd("~/winData/multiomics_networks_simulation")
#source("param_nw.R")

## Test the validity of the input parameters
## TO DO


##########################################################################################################################
###                                        DEFINE JULIA AND R FUNCTIONS                                                ###
##########################################################################################################################

## Get Julia functions from source code
juliaSource(paste0(getwd(),"/julia_functions.jl"))

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


getregnw = function(regsList, tarsList, reaction, nod){
  ## Inputs:
  ##  - regsList: a list of vector with regulator ids. 2 elements: 1st vector is for protein-coding regulators, the 2nd for noncoding regulators (as each can have a different set of target + have a different out-degree distribution)
  ##  - tarsList: a list of vector with target ids. 2 elements: 1st vector is for protein-coding regulators, the 2nd for noncoding regulators (as each can have a different set of target + have a different out-degree distribution)
  ##  - which reaction is regulated? (used to retrieve automatically the variables)
  ##  - the list of nodes from the system
  
  ## Construct the network nodes and edges data.frame
  nwnod = nod[nod$id %in% c(unlist(regsList), unlist(tarsList)),]
  nwnod = data.frame(nwnod, "nodetype" = rep("target", nrow(nwnod)), stringsAsFactors = F) # add a node attribute specifying if the node is a target, a protein regulator or a noncoding regulator (regulators can be also targets but will be labeled as regulators)
  nwnod[nwnod$id %in% regsList[["PC"]], "nodetype"] = "PCreg"
  nwnod[nwnod$id %in% regsList[["NC"]], "nodetype"] = "NCreg"
  
  ## Call the julia function nwgeneration to generate the regulatory network where the protein regulators are the regulatory nodes
  edgPC = juliaGet(juliaCall("nwgeneration", regsList[["PC"]], tarsList[["PC"]], get(paste(reaction, "PC", "indeg.distr", sep = ".")), get(paste(reaction, "PC", "outdeg.distr", sep = ".")), get(paste(reaction, "PC", "outdeg.exp", sep = ".")), get(paste(reaction, "PC", "autoregproba", sep = ".")), get(paste(reaction, "PC", "twonodesloop", sep = "."))))

  ## Call the julia function nwgeneration to generate the regulatory network where the noncoding regulators are the regulatory nodes
  edgNC = juliaGet(juliaCall("nwgeneration", regsList[["NC"]], tarsList[["NC"]], get(paste(reaction, "NC", "indeg.distr", sep = ".")), get(paste(reaction, "NC", "outdeg.distr", sep = ".")), get(paste(reaction, "NC", "outdeg.exp", sep = ".")), get(paste(reaction, "NC", "autoregproba", sep = ".")), get(paste(reaction, "NC", "twonodesloop", sep = "."))))
  
  ## create the edge dataframe
  nwedg = data.frame("from" = c(edgPC[,1], edgNC[,1]), "to" = c(edgPC[,2], edgNC[,2]), "TargetReaction" = rep(reaction,nrow(edgPC)+nrow(edgNC)), "RegSign" = rep("",nrow(edgPC)+nrow(edgNC)), "RegBy" = rep(c("PCreg", "NCreg"), c(nrow(edgPC), nrow(edgNC))), stringsAsFactors = F)
  
  ## Choose the sign (activation or repression) of each regulation (=edge)
  ## First for the regulatory interactions exerted by protein regulators
  if(nrow(edgPC) > 0) nwedg$RegSign[1:nrow(edgPC)] = sample(c("1","-1"), nrow(edgPC), prob = c(get(paste(reaction, "PC", "pos.p", sep = ".")), 1 - get(paste(reaction, "PC", "pos.p", sep = "."))), replace = T)
  ## Then for the regulatory interactions exerted by noncoding regulators
  if(nrow(edgNC) > 0) nwedg$RegSign[(nrow(edgPC)+1):(nrow(edgPC) + nrow(edgNC))] = sample(c("1","-1"), nrow(edgNC), prob = c(get(paste(reaction, "NC", "pos.p", sep = ".")), 1 - get(paste(reaction, "NC", "pos.p", sep = "."))), replace = T) 
  
  return(list("nod" = nwnod, "edg" = nwedg))
  
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



##########################################################################################################################
###                                                      R CODE                                                        ###
##########################################################################################################################


## THis function generates 
creationsystem = function(){
  
  source("param_nw.R")
  
  #### 0: CREATION OF THE DIFFERENT VARIABLES ----
  
  ## G.name : genes name 
  G.nameid = sapply(1:G, function(i){paste0("G_",i)})
  
  ## nod is the vertices data frame (1st column "id" = an integer value (faster for computation), 2nd column "name", 3rd column "coding" (values "NC" or "PC"),  
  ##  4rd column "TargetReaction" (values "TC", "TL", "RD", "PTM", "PD", "MR"), next columns: kinetic parameters (transcription rate, translation rate, RNA decay rate, protein decay rate))
  nod = data.frame("id" = 1:G, "nameid" = G.nameid, "coding" = rep("", G), "TargetReaction" = rep("", G),  "PTMform" = rep("", G), "ActiveForm" = rep("", G),
                   "TCrate" = rep(0,G), "TLrate" = rep(0,G), "RDrate" = rep(0,G), "PDrate" = rep(0,G), stringsAsFactors = F)
  rownames(nod) = nod$id
  
  ##edg is the edges data frame (1st column "from", 2nd column "to", 3rd column "TargetReaction" (values "TC", "TL", "RD", "PTM", "PD", "MR"), 4th column "RegSign" (value +1 or -1))
  edg = data.frame("from" = NULL, "to" = NULL, "TargetReaction" = NULL, "RegSign" = NULL, stringsAsFactors = F)
  
  ## complexes is a list where each element gives the components of a regulatory complex, whose id is the name of the element of the list
  complexes = list()
  
  
  
  # ------------------------------------------------
  #### STEP 1: decide the function of each gene ----
  # ------------------------------------------------
  
  ## Deciding gene function
  func = sample(c("NC.TC", "NC.TL", "NC.RD", "NC.PTM", "PC.TC", "PC.TL", "PC.RD", "PC.PD", "PC.PTM", "PC.MR"), G, replace = T,
                prob = c(NC.p*c(NC.TC.p, NC.TL.p, NC.RD.p, NC.PTM.p), PC.p*c(PC.TC.p, PC.TL.p, PC.RD.p, PC.PD.p, PC.PTM.p, PC.MR.p)))
  
  ## Extract from func the coding class of each gene (remove everything after the .)
  nod$coding = sub("\\.[[:alpha:]]+$", "", func)
  ## Extract from func the targeted reaction for each gene (remove everything after the .)
  nod$TargetReaction = sub("^[[:alpha:]]+\\.", "", func)
  
  ## Choose which proteins (from protein-coding genes) have a PTM form
  nod$PTMform[nod$coding == "PC"] = sample(c("1","0"), sum(nod$coding == "PC"), prob = c(PC.PTM.form.p, 1-PC.PTM.form.p), replace = T)
  nod$PTMform[nod$coding == "NC"] = "0"
  
  ## In nod, state which form (i.e. RNA, protein, activated protein) is performing the regulation
  nod$ActiveForm[nod$coding == "NC"] = "R" ## noncoding genes act through their RNA
  nod$ActiveForm[nod$coding == "PC" & nod$PTMform == "0"] = "P" ## protein-coding genes act through their protein
  nod$ActiveForm[nod$coding == "PC" & nod$PTMform == "1"] = "Pm" ## For proteins undergoing a post-translational modification, only the PTM form is active
  nod$ActiveForm = sapply(1:nrow(nod), function(x){paste0(nod$ActiveForm[x],".",nod$id[x])})
  
  
  # ----------------------------------------------------------
  #### STEP 2: sample the kinetic parameters of each gene ----
  # ----------------------------------------------------------
  
  ## Transcription rate: applicable to all genes
  nod$TCrate = get(basal_transcription_rate)(G)
  
  ## RNA decay rate: applicable to all genes
  ## Sample the lifetime, the decay rate is defined as 1/lifetime
  nod$RDrate = 1/get(basal_RNAlifetime)(G)
  
  ## Translation rate: applicable to protein-coding genes
  nod$TLrate[nod$coding == "PC"] = get(basal_translation_rate)(sum(nod$coding == "PC"))
  
  ## Protein coding rate: applicable to protein-coding genes
  nod$PDrate[nod$coding == "PC"] = 1/get(basal_protlifetime)(sum(nod$coding == "PC"))
  
  
  # -----------------------------------------------------------------
  #### STEP 3 : Define Transcriptional regulatory network (TCRN) ----
  # -----------------------------------------------------------------
  
  ## Identify TFs in the system
  PCreg.id = nod$id[nod$coding == "PC" & nod$TargetReaction == "TC"]
  # Identify non-coding genes regulating transcription in the system
  NCreg.id = nod$id[nod$coding == "NC" & nod$TargetReaction == "TC"]
  
  ## which genes are targeted by TFs? here any gene (including the TFs because a TF can autoregulate himself)
  PCtarget.id = nod$id
  ## which genes are targeted by noncoding RNAs? here any protein-coding gene 
  NCtarget.id = nod$id[nod$coding == "PC"]
  
  ## Construct the regulatory network
  TCRN = getregnw(regsList = list("PC" = PCreg.id, "NC" = NCreg.id), tarsList = list("PC" = PCtarget.id, "NC" = NCtarget.id), reaction = "TC", nod = nod)
  TCRN.nod = TCRN[["nod"]]
  TCRN.edg = TCRN[["edg"]]
  
  ## Create corresponding igraph object (to save the regulatory interactions as is, before the creation of the combinatorial regulation)
  TCRN.nw = igraph::graph_from_data_frame(d = TCRN.edg, directed = T, vertices = TCRN.nod)
  
  ## Creation of combinatiorial regulation
  ## if regcomplexes != 'none', if several regulators control a common target they can form regulatory complexes
  ## The composition of each complex is stored in complexes
  ## Complexes can be composed only of proteins if recomplexes = "prot" or protein and noncoding regulators if regcomplexes = "both"
  
  if(regcomplexes == "none"){ ## If regulators controlling a same target are not allowed to form a regulatory complex, simply reformat the edg and nod dataframes
    
    TCRN.edg$from = sapply(TCRN.edg$from, toString) ## transform the id of regulators from integer to string (not the id of target because we need it to be integer for computational speed later)

  }else if(regcomplexes == "prot"){ ## If the regulatory complexes can only be protein complexes
    
      temp = TCRN.edg[TCRN.edg$RegBy =="PCreg", ]
      tempregcom = juliaGet(juliaCall("combreg", PCtarget.id, temp$from, temp$to, temp$RegSign, regcomplexes.p, regcomplexes.size, "TC"))
      ## only keep the noncoding regulators (the regulation from protein-coding regulators is given by the Julia function combreg)
      TCRN.edg = TCRN.edg[TCRN.edg$RegBy =="NCreg", c("from", "to", "TargetReaction", "RegSign")]
      TCRN.edg = rbind(TCRN.edg, data.frame("from" = unlist(tempregcom$newedg[,1]), "to" = unlist(tempregcom$newedg[,2]), "TargetReaction" = rep("TC", nrow(tempregcom$newedg)), "RegSign" = unlist(tempregcom$newedg[,3])))
      rownames(TCRN.edg) = NULL
      complexes = c(complexes, lapply(tempregcom$Complexes, unlist))
      
  }else if(regcomplexes == "both"){ ## If the regulatory complexes can be protein/noncoding complexes
    
    tempregcom = juliaGet(juliaCall("combreg", unique(c(PCtarget.id, NCtarget.id)), TCRN.edg$from, TCRN.edg$to, TCRN.edg$RegSign, regcomplexes.p, regcomplexes.size, "TC"))
    TCRN.edg = data.frame("from" = unlist(tempregcom$newedg[,1]), "to" = unlist(tempregcom$newedg[,2]), "TargetReaction" = rep("TC", nrow(tempregcom$newedg)), "RegSign" = unlist(tempregcom$newedg[,3]))
    rownames(TCRN.edg) = NULL
    complexes = c(complexes, lapply(tempregcom$Complexes, unlist))
    
  }
    
  
  
  ## Sample the kinetic parameters of each regulatory interaction
  ##    Kinetic parameters for transcription regulation include the binding and unbinding rate of regulators to gene promoter, and the fold change induced on transcription rate by a regulator bound to the promoter
  TCRN.edg = data.frame(TCRN.edg, "TCbindingrate" = get(TCbindingrate)(nrow(TCRN.edg)), "TCunbindingrate" = get(TCunbindingrate)(nrow(TCRN.edg)), "TCfoldchange" = rep(0, nrow(TCRN.edg)), stringsAsFactors = F)
  TCRN.edg$TCfoldchange[TCRN.edg$RegSign == "1"] = get(TCfoldchange)(sum(TCRN.edg$RegSign == "1"))  ## Repressors induce a fold change of 0
  

  # ## Construct the network nodes and edges data.frame
  # TCRN.nod = nod[nod$id %in% c(TF.id, ncRNA.id, targetTF.id, targetncRNA.id),]
  # 
  # # Call the julia function nwgeneration to construct the TF-mediated regulatory network
  # edg = juliaGet(juliaCall("nwgeneration", TF.id, targetTF.id, TC.TF.indeg.distr, TC.TF.outdeg.distr, TC.TF.outdeg.exp, TC.TF.autoregproba, TC.twonodesloop))
  # temp = nrow(edg) #keep in memory the number of edges coming from protein-coding genes (to be used when sampling the sign of the regulation for each edge)
  # #edg = juliaGet(juliaCall("nwgenerationSR", list(TF.id, ncRNA.id), target.id, TC.indeg.distr, list(TC.TF.outdeg.distr, TC.ncRNA.outdeg.distr), list(TC.TF.outdeg.exp, TC.ncRNA.outdeg.exp)))
  # edg = juliaGet(juliaCall("nwgeneration", ncRNA.id, targetncRNA.id , TC.ncRNA.indeg.distr, TC.ncRNA.outdeg.distr, TC.ncRNA.outdeg.exp, TC.ncRNA.autoregproba, TC.twonodesloop, edg))
  # 
  # TCRN.edg = data.frame("from" = edg[,1], "to" = edg[,2], "TargetReaction" = rep("TC",nrow(edg)), "RegSign" = rep("",nrow(edg)), 
  #                       "bindingrate" = rep(0,nrow(edg)), "unbindingrate" = rep(0,nrow(edg)), "FoldChange" = rep(0,nrow(edg)), stringsAsFactors = F)
  # 
  # ## Choose the sign (activation or repression) of each regulation (=edge)
  # TCRN.edg$RegSign[1:temp] = sample(c("1","-1"), temp, prob = c(TC.PC.pos.p, 1 - TC.PC.pos.p), replace = T) ## the first rows are edges coming from protein-coding genes
  # TCRN.edg$RegSign[(temp+1):nrow(TCRN.edg)] = sample(c("1","-1"), nrow(TCRN.edg)-temp, prob = c(TC.NC.pos.p, 1 - TC.PC.pos.p), replace = T) ## the rest of the edges are from noncoding genes
  # 
  # ## Sample the kinetic parameters for each regulatory reaction
  # TCRN.edg$bindingrate = get(TCbindingrate)(nrow(TCRN.edg))
  # TCRN.edg$unbindingrate = get(TCunbindingrate)(nrow(TCRN.edg))
  # TCRN.edg$FoldChange[TCRN.edg$RegSign == "1"] = get(TCunbindingrate)(sum(TCRN.edg$RegSign == "1"))
  
  
  # -----------------------------------------------------------------
  #### STEP 4 : Define Translational regulatory network (TLRN) ----
  # -----------------------------------------------------------------
  
  ## Identify TLFs in the system
  PCreg.id = nod$id[nod$coding == "PC" & nod$TargetReaction == "TL"]
  # Identify non-coding genes regulating translation in the system
  NCreg.id = nod$id[nod$coding == "NC" & nod$TargetReaction == "TL"]
  
  ## which genes are targeted by TLFs? here any protein-coding gene
  PCtarget.id = nod$id[nod$coding == "PC"]
  ## which genes are targeted by noncoding RNAs? here any protein-coding gene 
  NCtarget.id = nod$id[nod$coding == "PC"]
  
  ## Construct the regulatory network
  TLRN = getregnw(regsList = list("PC" = PCreg.id, "NC" = NCreg.id), tarsList = list("PC" = PCtarget.id, "NC" = NCtarget.id), reaction = "TL", nod = nod)
  TLRN.nod = TLRN[["nod"]]
  TLRN.edg = TLRN[["edg"]]
  
  ## Create corresponding igraph object
  TLRN.nw = igraph::graph_from_data_frame(d = TLRN.edg, directed = T, vertices = TLRN.nod)
  
  ## Creation of combinatorial regulation
  ## if regcomplexes != 'none', if several regulators control a common target they can form regulatory complexes
  ## The composition of each complex is stored in complexes
  ## Complexes can be composed only of proteins if recomplexes = "prot" or protein and noncoding regulators if regcomplexes = "both"
  
  if(regcomplexes == "none"){ ## If regulators controlling a same target are not allowed to form a regulatory complex, simply reformat the edg and nod dataframes
    
    TLRN.edg$from = sapply(TLRN.edg$from, toString) ## transform the id of regulators from integer to string (not the id of target because we need it to be integer for computational speed later)

  }else if(regcomplexes == "prot"){ ## If the regulatory complexes can only be protein complexes
    
    temp = TLRN.edg[TLRN.edg$RegBy =="PCreg", ]
    tempregcom = juliaGet(juliaCall("combreg", PCtarget.id, temp$from, temp$to, temp$RegSign, regcomplexes.p, regcomplexes.size, "TL"))
    ## only keep the noncoding regulators (the regulation from protein-coding regulators is given by the Julia function combreg)
    TLRN.edg = TLRN.edg[TLRN.edg$RegBy =="NCreg", c("from", "to", "TargetReaction", "RegSign")]
    TLRN.edg = rbind(TLRN.edg, data.frame("from" = unlist(tempregcom$newedg[,1]), "to" = unlist(tempregcom$newedg[,2]), "TargetReaction" = rep("TL", nrow(tempregcom$newedg)), "RegSign" = unlist(tempregcom$newedg[,3])))
    rownames(TLRN.edg) = NULL
    complexes = c(complexes, lapply(tempregcom$Complexes, unlist))
    
  }else if(regcomplexes == "both"){ ## If the regulatory complexes can be protein/noncoding complexes
    
    tempregcom = juliaGet(juliaCall("combreg", unique(c(PCtarget.id, NCtarget.id)), TLRN.edg$from, TLRN.edg$to, TLRN.edg$RegSign, regcomplexes.p, regcomplexes.size, "TL"))
    TLRN.edg = data.frame("from" = unlist(tempregcom$newedg[,1]), "to" = unlist(tempregcom$newedg[,2]), "TargetReaction" = rep("TL", nrow(tempregcom$newedg)), "RegSign" = unlist(tempregcom$newedg[,3]))
    rownames(TLRN.edg) = NULL
    complexes = c(complexes, lapply(tempregcom$Complexes, unlist))
    
  }
  
  
  ## Sample the kinetic parameters of each regulatory interaction
  ##    Kinetic parameters for translation regulation include the binding and unbinding rate of regulators to mRNA binding sequence, and the fold change induced on transcription rate by a regulator bound to the mRNA
  TLRN.edg = data.frame(TLRN.edg, "TLbindingrate" = get(TLbindingrate)(nrow(TLRN.edg)), "TLunbindingrate" = get(TLunbindingrate)(nrow(TLRN.edg)), "TLfoldchange" = rep(0, nrow(TLRN.edg)), stringsAsFactors = F)
  TLRN.edg$TLfoldchange[TLRN.edg$RegSign == "1"] = get(TLfoldchange)(sum(TLRN.edg$RegSign == "1")) ## Repressors induce a fold change of 0
  
  
  # -----------------------------------------------------------------
  #### STEP 5 : Define RNA decay regulatory network (RDRN) ----
  # -----------------------------------------------------------------
  
  ## Identify proteins regulating RNA decay in the system
  PCreg.id = nod$id[nod$coding == "PC" & nod$TargetReaction == "RD"]
  ## Identify noncoding RNAs regulating RNA decay (miRNAs or siRNAs for ex)
  NCreg.id = nod$id[nod$coding == "NC" & nod$TargetReaction == "RD"]
  
  ## which genes are targeted by coding regulators? here any gene (temporary)
  PCtarget.id = nod$id
  ## which genes are targeted by noncoding regulators? here any gene (temporary)
  NCtarget.id = nod$id
  
  ## Construct the regulatory network
  RDRN = getregnw(regsList = list("PC" = PCreg.id, "NC" = NCreg.id), tarsList = list("PC" = PCtarget.id, "NC" = NCtarget.id), reaction = "RD", nod = nod)
  RDRN.nod = RDRN[["nod"]]
  RDRN.edg = RDRN[["edg"]]
  
  ## Create corresponding igraph object
  RDRN.nw = igraph::graph_from_data_frame(d = RDRN.edg, directed = T, vertices = RDRN.nod)
  
  ## Creation of combinatorial regulation
  ## if regcomplexes != 'none', if several regulators control a common target they can form regulatory complexes
  ## The composition of each complex is stored in complexes
  ## Complexes can be composed only of proteins if recomplexes = "prot" or protein and noncoding regulators if regcomplexes = "both"
  
  if(regcomplexes == "none"){ ## If regulators controlling a same target are not allowed to form a regulatory complex, simply reformat the edg and nod dataframes
    
    RDRN.edg$from = sapply(RDRN.edg$from, toString) ## transform the id of regulators from integer to string (not the id of target because we need it to be integer for computational speed later)
    
  }else if(regcomplexes == "prot"){ ## If the regulatory complexes can only be protein complexes
    
    temp = RDRN.edg[RDRN.edg$RegBy =="PCreg", ]
    tempregcom = juliaGet(juliaCall("combreg", PCtarget.id, temp$from, temp$to, temp$RegSign, regcomplexes.p, regcomplexes.size, "RD"))
    ## only keep the noncoding regulators (the regulation from protein-coding regulators is given by the Julia function combreg)
    RDRN.edg = RDRN.edg[RDRN.edg$RegBy =="NCreg", c("from", "to", "TargetReaction", "RegSign")]
    RDRN.edg = rbind(RDRN.edg, data.frame("from" = unlist(tempregcom$newedg[,1]), "to" = unlist(tempregcom$newedg[,2]), "TargetReaction" = rep("RD", nrow(tempregcom$newedg)), "RegSign" = unlist(tempregcom$newedg[,3])))
    rownames(RDRN.edg) = NULL
    complexes = c(complexes, lapply(tempregcom$Complexes, unlist))
    
  }else if(regcomplexes == "both"){ ## If the regulatory complexes can be protein/noncoding complexes
    
    tempregcom = juliaGet(juliaCall("combreg", unique(c(PCtarget.id, NCtarget.id)), RDRN.edg$from, RDRN.edg$to, RDRN.edg$RegSign, regcomplexes.p, regcomplexes.size, "RD"))
    RDRN.edg = data.frame("from" = unlist(tempregcom$newedg[,1]), "to" = unlist(tempregcom$newedg[,2]), "TargetReaction" = rep("RD", nrow(tempregcom$newedg)), "RegSign" = unlist(tempregcom$newedg[,3]))
    rownames(RDRN.edg) = NULL
    complexes = c(complexes, lapply(tempregcom$Complexes, unlist))
    
  }
  
  ## Sample the kinetic parameters of each regulatory interaction
  ##    Kinetic parameters for RNA decay includes the binding (and unbinding for repressors of decay) rate of the regulator on the RNA 
  RDRN.edg = data.frame(RDRN.edg, "RDbindingrate" = get(RDbindingrate)(nrow(RDRN.edg)), "RDunbindingrate" = rep(0, nrow(RDRN.edg)))
  RDRN.edg$RDunbindingrate[RDRN.edg$RegSign == "-1"] = get(RDunbindingrate)(sum(RDRN.edg$RegSign == "-1")) ## only the repressors (i.e. protecting the target from decay) can unbind the target
  
  
  
  # -----------------------------------------------------------------
  #### STEP 6 : Define protein decay regulatory network (PDRN) ----
  # -----------------------------------------------------------------
  
  ## Identify proteins regulating RNA decay in the system
  PCreg.id = nod$id[nod$coding == "PC" & nod$TargetReaction == "PD"]
  ## Identify noncoding RNAs regulating RNA decay (miRNAs or siRNAs for ex)
  NCreg.id = nod$id[nod$coding == "NC" & nod$TargetReaction == "PD"]
  
  ## which genes are targeted by coding regulators? here any gene (temporary)
  PCtarget.id = nod$id[nod$coding == "PC"]
  ## which genes are targeted by noncoding regulators? here any gene (temporary)
  NCtarget.id = nod$id[nod$coding == "PC"]
  
  ## Construct the regulatory network
  PDRN = getregnw(regsList = list("PC" = PCreg.id, "NC" = NCreg.id), tarsList = list("PC" = PCtarget.id, "NC" = NCtarget.id), reaction = "PD", nod = nod)
  PDRN.nod = PDRN[["nod"]]
  PDRN.edg = PDRN[["edg"]]
  
  ## Create corresponding igraph object
  PDRN.nw = igraph::graph_from_data_frame(d = PDRN.edg, directed = T, vertices = PDRN.nod)
  
  
  ## Creation of combinatorial regulation
  ## if regcomplexes != 'none', if several regulators control a common target they can form regulatory complexes
  ## The composition of each complex is stored in complexes
  ## Complexes can be composed only of proteins if recomplexes = "prot" or protein and noncoding regulators if regcomplexes = "both"
  
  if(regcomplexes != "none"){ ## If regulators controlling a same target are not allowed to form a regulatory complex, simply reformat the edg and nod dataframes
    
    PDRN.edg$from = sapply(PDRN.edg$from, toString) ## transform the id of regulators from integer to string (not the id of target because we need it to be integer for computational speed later)
    
  }else if(regcomplexes == "prot"){ ## If the regulatory complexes can only be protein complexes
    
    temp = PDRN.edg[PDRN.edg$RegBy =="PCreg", ]
    tempregcom = juliaGet(juliaCall("combreg", PCtarget.id, temp$from, temp$to, temp$RegSign, regcomplexes.p, regcomplexes.size, "PD"))
    ## only keep the noncoding regulators (the regulation from protein-coding regulators is given by the Julia function combreg)
    PDRN.edg = PDRN.edg[PDRN.edg$RegBy =="NCreg", c("from", "to", "TargetReaction", "RegSign")]
    PDRN.edg = rbind(PDRN.edg, data.frame("from" = unlist(tempregcom$newedg[,1]), "to" = unlist(tempregcom$newedg[,2]), "TargetReaction" = rep("PD", nrow(tempregcom$newedg)), "RegSign" = unlist(tempregcom$newedg[,3])))
    rownames(PDRN.edg) = NULL
    complexes = c(complexes, lapply(tempregcom$Complexes, unlist))
    
  }else if(regcomplexes == "both"){ ## If the regulatory complexes can be protein/noncoding complexes
    
    tempregcom = juliaGet(juliaCall("combreg", unique(c(PCtarget.id, NCtarget.id)), PDRN.edg$from, PDRN.edg$to, PDRN.edg$RegSign, regcomplexes.p, regcomplexes.size, "PD"))
    PDRN.edg = data.frame("from" = unlist(tempregcom$newedg[,1]), "to" = unlist(tempregcom$newedg[,2]), "TargetReaction" = rep("PD", nrow(tempregcom$newedg)), "RegSign" = unlist(tempregcom$newedg[,3]))
    rownames(PDRN.edg) = NULL
    complexes = c(complexes, lapply(tempregcom$Complexes, unlist))
    
  }
  
  ## Sample the kinetic parameters of each regulatory interaction
  ##    Kinetic parameters for protein decay includes the binding (and unbinding for repressors of decay) rate of the regulator on the protein 
  PDRN.edg = data.frame(PDRN.edg, "PDbindingrate" = get(PDbindingrate)(nrow(PDRN.edg)), "PDunbindingrate" = rep(0, nrow(PDRN.edg)))
  PDRN.edg$PDunbindingrate[PDRN.edg$RegSign == "-1"] = get(PDunbindingrate)(sum(PDRN.edg$RegSign == "-1")) ## only the repressors (i.e. protecting the target from decay) can unbind the target
  
  
  # ------------------------------------------------------------------------------------------
  #### STEP 7 : Define protein post-translational modification regulatory network (PTMRN) ----
  # ------------------------------------------------------------------------------------------
  
  
  # ---------------------------------------------------------------
  #### STEP 8 : Define regulatory complexes kinetic parameters ----
  # ---------------------------------------------------------------
  
  complexeskinetics = list()
  formrates = get(complexesformationrate)(length(complexes))
  dissrates = get(complexesdissociationrate)(length(complexes))
  for(c in 1:length(complexes)){
    complexeskinetics[[names(complexes)[c]]] = list("formationrate" = formrates[c], "dissociationrate" = dissrates[c])
  }
  
  # -----------------------------------------------------------------
  ####                          RETURN                           ----
  # -----------------------------------------------------------------
  
  
  
  ## save all the interactions in edg
  temp = list(TCRN.edg, TLRN.edg, RDRN.edg, PDRN.edg)
  for(t in temp){
    edg = rbind(edg, t[,c("from", "to", "TargetReaction", "RegSign")])
  }
  
  ## Create the lists to be sent to Julia
  nodDF = nod
  temp = list("nod", "TCRN.edg", "TLRN.edg", "RDRN.edg", "PDRN.edg")
  for(t in temp){
    new = list()
    for(cols in colnames(get(t))){
      new[[cols]] = get(t)[,cols]
    }
    assign(t, new)
  }
  
  ## Return
  res = list("nod" = nodDF, 
             "edg" = edg,
             "complexes" = complexes,
             "RN.edg" = list("TCRN.edg" = TCRN.edg,
                            "TLRN.edg" = TLRN.edg,
                            "RDRN.edg" = RDRN.edg,
                            "RDRN.edg" = RDRN.edg,
                            "PDRN.edg" = PDRN.edg),
             "RN.nw" = list("TCRN.nw" = TCRN.nw,
                            "TLRN.nw" = TLRN.nw,
                            "RDRN.nw" = RDRN.nw,
                            "RDRN.nw" = RDRN.nw,
                            "PDRN.nw" = PDRN.nw))
  
  return(res)
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
