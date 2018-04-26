##########################################################################################################################
##########################################################################################################################
###                                       NETWORK GENERATION  - using Julia                                            ###
##########################################################################################################################
##########################################################################################################################

if(!suppressWarnings(require("XRJulia", quietly = T))){install.packages("XRJulia")}
library(XRJulia)

if(!suppressWarnings(require("igraph", quietly = T))){install.packages("igraph")}
library(igraph)

## Temporary - load parameters from param_nw.R
setwd("~/winData/multiomics_networks_simulation")
source("param_nw.R")

## Test if Julia is installed on the computer
if(!findJulia(test = T)) stop("Julia is not installed on the computer or not accessible by R. Check that Julia is correcly installed and/or in the PATH variable\n")


## Test the validity of the input parameters
## TO DO


##########################################################################################################################
###                                        DEFINE JULIA AND R FUNCTIONS                                                ###
##########################################################################################################################

## Get Julia functions from source code
juliaSource("julia_functions.jl")

## Fit a power law to the frequency distribution of the values given in obs using nls
fitPowerLaw = function(obs, cstart = 1, gstart= 1, cmin= 0, gmin= 0.1){
  cat("\n Fitting the distribution freq(obs = x) ~ c*x^(-g)\n\n")
  x = 1:max(obs)
  y = sapply(x, function(i){sum(obs == i)})
  y = y/sum(y)
  
  m = nls(y~c*x^(-g), algorithm = "port", start = list(c = cstart, g = gstart), lower = c(cmin, gmin))
  print(m)
  cat("\n Correlation between observed and predicted values:   ")
  cat(cor(y, predict(m)))
}

## Fit an exponential to the frequency distribution of the values given in obs using nls
fitExponential = function(obs, cstart = 1, lstart= 1, cmin= 0, lmin= 0.1){
  cat("\n Fitting the distribution freq(obs = x) ~ (c/l)*exp(-x/l)\n\n")
  x = 1:max(obs)
  y = sapply(x, function(i){sum(obs == i)})
  y = y/sum(y)
  
  m = nls(y~(c/l)*exp(-x/l), algorithm = "port", start = list(c = cstart, l = lstart), lower = c(cmin,lmin))
  print(m)
  cat("\n Correlation between observed and predicted values: \t ")
  cat(cor(y, predict(m)))
}

##########################################################################################################################
###                                                      R CODE                                                        ###
##########################################################################################################################

#### 0: CREATION OF THE DIFFERENT VARIABLES ----

## G.id : ID of genes 
G.id = sapply(1:G, function(i){paste0("G_",i)})

##Mod.id : ID of modules
Mod.id = sapply(1:Mod, function(i){paste0("Module_",i)})

##nod is the vertices data frame (1st column "id", 2nd column "coding" (values "NC" or "PC"), 3rd column "TargetReaction" (values "TC", "TL", "RD", "PTM", "PD", "MR"), 4th column "Module" (which module the gene belongs to))
nod = data.frame("id" = G.id, "coding" = rep("", G), "TargetReaction" = rep("", G),  "PTMform" = rep("", G), "Module" = rep("", G), stringsAsFactors = F)
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

## Choose the module each gene belongs to
nod$Module = sample(Mod.id, G, prob = rep(1.0/Mod, Mod), replace = T)




# -----------------------------------------------------------------
#### STEP 2 : Define Transcriptional regulatory network (TCRN) ----
# -----------------------------------------------------------------

## Identify TFs in the system
TF.id = nod$id[nod$coding == "PC" & nod$TargetReaction == "TC"]
TF = length(TF.id)

## which genes are targeted by TFs? here any gene (including the TFs because a TF can autoregulate himself)
TF.target = nod$id
## out.max: maximum out-degree, equal to the number of possible targets
out.max = length(TF.target)

## Construct the network nodes and edges data.frame
TCRN.nod = nod[nod$id %in% c(TF.id, TF.target),]

edg = juliaGet(juliaCall("nwgeneration", TF.id, TF.target, "powerlaw", "powerlaw", 0.8))
TCRN.edg = data.frame("from" = edg[,1], "to" = edg[,2], "TargetReaction" = rep("TC",nrow(edg)), "RegSign" = rep("",nrow(edg)), stringsAsFactors = F)

## Choose the sign (activation or repression) of each regulation (=edge)
TCRN.edg$RegSign = sample(c("1","-1"), nrow(TCRN.edg), prob = c(TC.pos.p, 1 - TC.pos.p), replace = T)

## Create corresponding igraph object
TCRN.nw = igraph::graph_from_data_frame(d = TCRN.edg, directed = T, vertices = TCRN.nod)


##############################################################

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

