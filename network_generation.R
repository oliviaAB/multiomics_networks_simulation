##########################################################################################################################
##########################################################################################################################
###                                       NETWORK GENERATION  - using Julia                                            ###
##########################################################################################################################
##########################################################################################################################

if(!suppressWarnings(require("XRJulia", quietly = T))){install.packages("XRJulia")}
library(XRJulia)

if(!suppressWarnings(require("igraph", quietly = T))){install.packages("igraph")}
library(igraph)

## Test if Julia is installed on the computer
if(!findJulia(test = T)) stop("Julia is not installed on the computer or not accessible by R. Check that Julia is correcly installed and/or in the PATH variable\n")


##########################################################################################################################
###                                               PARAMETERS                                                           ###
##########################################################################################################################


## Temporary - load parameters from param_nw.R
setwd("~/winData/multiomics_networks_simulation")
source("param_nw.R")

## Test the validity of the input parameters
## TO DO


##########################################################################################################################
###                                        DEFINE JULIA AND R FUNCTIONS                                                ###
##########################################################################################################################

## Get Julia functions from source code
juliaSource("julia_functions.jl")

howmanyautoreg = function(edg){
  cat("There are  ")
  cat(sum(edg[,1] == edg[,2]))
  cat("  self-regulatory edges in the graph.\n")
}

howmanyloops = function(edg){
  cat("There are  ")
  res = edg[(edg[,1]!=edg[,2]), ]
  cat(sum(duplicated.matrix(rbind(res, res[,2:1])))/2)
  cat("  2-nodes loops in the graph.\n")
}

##########################################################################################################################
###                                                      R CODE                                                        ###
##########################################################################################################################

#### 0: CREATION OF THE DIFFERENT VARIABLES ----

## G.name : genes name 
G.nameid = sapply(1:G, function(i){paste0("G_",i)})

## nod is the vertices data frame (1st column "id" = an integer value (faster for computation), 2nd column "name", 3rd column "coding" (values "NC" or "PC"),  
##  4rd column "TargetReaction" (values "TC", "TL", "RD", "PTM", "PD", "MR"), next columns: kinetic parameters (transcription rate, translation rate, RNA decay rate, protein decay rate))
nod = data.frame("id" = 1:G, "nameid" = G.nameid, "coding" = rep("", G), "TargetReaction" = rep("", G),  "PTMform" = rep("", G),
                 "TCrate" = rep(0,G), "TLrate" = rep(0,G), "RDrate" = rep(0,G), "PDrate" = rep(0,G), stringsAsFactors = F)
rownames(nod) = nod$id

##edg is the edges data frame (1st column "from", 2nd column "to", 3rd column "TargetReaction" (values "TC", "TL", "RD", "PTM", "PD", "MR"), 4th column "RegSign" (value +1 or -1))
edg = data.frame("from" = NULL, "to" = NULL, "TargetReaction" = NULL, "RegSign" = NULL, stringsAsFactors = F)





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
TF.id = nod$id[nod$coding == "PC" & nod$TargetReaction == "TC"]

## which genes are targeted by TFs? here any gene (including the TFs because a TF can autoregulate himself)
target.id = nod$id

## Construct the network nodes and edges data.frame
TCRN.nod = nod[nod$id %in% c(TF.id, target.id),]

edg = juliaGet(juliaCall("nwgeneration", TF.id, target.id, "exponential", "powerlaw", TC.outdeg.power))
TCRN.edg = data.frame("from" = edg[,1], "to" = edg[,2], "TargetReaction" = rep("TC",nrow(edg)), "RegSign" = rep("",nrow(edg)), 
                      "bindingrate" = rep(0,nrow(edg)), "unbindingrate" = rep(0,nrow(edg)), "FoldChange" = rep(0,nrow(edg)), stringsAsFactors = F)

## Choose the sign (activation or repression) of each regulation (=edge)
TCRN.edg$RegSign = sample(c("1","-1"), nrow(TCRN.edg), prob = c(TC.pos.p, 1 - TC.pos.p), replace = T)

## Sample the kinetic parameters for each regulatory reaction
TCRN.edg$bindingrate = get(TCbindingrate)(nrow(TCRN.edg))
TCRN.edg$unbindingrate = get(TCunbindingrate)(nrow(TCRN.edg))
TCRN.edg$FoldChange[TCRN.edg$RegSign == "1"] = get(TCunbindingrate)(sum(TCRN.edg$RegSign == "1"))

## Create corresponding igraph object
TCRN.nw = igraph::graph_from_data_frame(d = TCRN.edg, directed = T, vertices = TCRN.nod)



# -----------------------------------------------------------------
#### STEP 4 : Define Translational regulatory network (TLRN) ----
# -----------------------------------------------------------------

## Identify TLFs in the system
TLF.id = nod$id[nod$coding == "PC" & nod$TargetReaction == "TL"]

## which genes are targeted by TFs? here any gene (including the TFs because a TF can autoregulate himself)
target.id = nod$id

## Construct the network nodes and edges data.frame
TCRN.nod = nod[nod$id %in% c(TF.id, target.id),]

edg = juliaGet(juliaCall("nwgeneration", TF.id, target.id, "exponential", "powerlaw", 0.8))
TCRN.edg = data.frame("from" = edg[,1], "to" = edg[,2], "TargetReaction" = rep("TC",nrow(edg)), "RegSign" = rep("",nrow(edg)), 
                      "bindingrate" = rep(0,nrow(edg)), "unbindingrate" = rep(0,nrow(edg)), "FoldChange" = rep(0,nrow(edg)), stringsAsFactors = F)

## Choose the sign (activation or repression) of each regulation (=edge)
TCRN.edg$RegSign = sample(c("1","-1"), nrow(TCRN.edg), prob = c(TC.pos.p, 1 - TC.pos.p), replace = T)

## Sample the kinetic parameters for each regulatory reaction
TCRN.edg$bindingrate = get(TCbindingrate)(nrow(TCRN.edg))
TCRN.edg$unbindingrate = get(TCunbindingrate)(nrow(TCRN.edg))
TCRN.edg$FoldChange[TCRN.edg$RegSign == "1"] = get(TCunbindingrate)(sum(TCRN.edg$RegSign == "1"))

## Create corresponding igraph object
TCRN.nw = igraph::graph_from_data_frame(d = TCRN.edg, directed = T, vertices = TCRN.nod)



howmanyautoreg(edg)
howmanyloops(edg)

############################################################################################################################
############################################################################################################################


## Color differently TFs and other genes
colrs = c("gold", "tomato")
V(TCRN.nw)$color = colrs[V(TCRN.nw)$name %in% TF.id +1]

# colrs = rainbow(Mod)
# V(TCRN.nw)$color = colrs[as.factor(V(TCRN.nw)$Module)]

## Color differently activating and inhibiting regulations
colrs = c("1" = "red", "-1" = "blue")
E(TCRN.nw)$color = colrs[E(TCRN.nw)$RegSign]


## Visualisation
plot(TCRN.nw, edge.arrow.size=.4, layout = layout_with_fr, vertex.size = 7, vertex.label = NA)

fitPowerLaw(degree(TCRN.nw, mode = "out")[degree(TCRN.nw, mode = "out")!=0])
hist(degree(TCRN.nw, mode = "out"))

fitPowerLaw(degree(TCRN.nw, mode = "in"))
hist(degree(TCRN.nw, mode = "in"))

##############################################################

reg = sapply(1:157, function(x){paste0("G",x)})
target = sapply(1:4410, function(x){paste0("G",x)})

edg = juliaGet(juliaCall("nwgeneration", reg, target, "exponential", "powerlaw", 0.8))
edg = data.frame("from" = edg[,1], "to" = edg[,2], stringsAsFactors = F)

## Create corresponding igraph object
nw = igraph::graph_from_data_frame(d = edg, directed = T)


##############################################################

iter = 100

Edist = vector(length = iter)
tarav = vector(length = iter)
for(i in 1:iter){
  out = juliaGet(juliaCall("samplepowerlaw", 124, 1.4, 4410))
  Edist[i] = sum(out)
  tarav[i] = mean(out)
}



hist(Edist)
hist(tarav)
hist(2*Edist/(157+4410))
