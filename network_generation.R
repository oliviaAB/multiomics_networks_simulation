##########################################################################################################################
##########################################################################################################################
###                                       NETWORK GENERATION  - using Julia                                            ###
##########################################################################################################################
##########################################################################################################################

if(!suppressWarnings(require("XRJulia", quietly = T))){install.packages("XRJulia")}
library(XRJulia)

## Temporary - load parameters from param_nw.R
setwd("~/winData/multiomics_networks_simulation")
source("param_nw.R")

## Test if Julia is installed on the computer
if(!findJulia(test = T)) stop("Julia is not installed on the computer or not accessible by R. Check that Julia is correcly installed and/or in the PATH variable\n")


## Test the validity of the input parameters
## TO DO


##########################################################################################################################
###                                              DEFINE JULIA FUNCTIONS                                                ###
##########################################################################################################################

j.sample_expon = juliaEval("
    function sample_expon(n, lambda, max)
      prob = (1/lambda)*exp.(-(1:max)/lambda)
      cumprob = cumsum(prob / sum(prob))
      rnb = rand(Int(n))
      res = Int64[]
      for nb in rnb
        push!(res, findfirst(x -> x >= nb, cumprob))
      end 
      return res
    end
                  ")


j.sample_powerlaw = juliaEval("
    function sample_powerlaw(n, gamma, max)
      prob = (1:max).^(-gamma)
      cumprob = cumsum(prob / sum(prob))
      rnb = rand(Int(n))
      res = Int64[]
      for nb in rnb
        push!(res, findfirst(x -> x >= nb, cumprob))
      end 
      return res
    end
                  ")




j.sample_expon_fct = JuliaFunction(j.sample_expon)
j.sample_powerlaw_fct = JuliaFunction(j.sample_powerlaw)

myres = juliaGet(j.sample_powerlaw_fct(500,0.8,60))

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
TCRN.edg = data.frame("from" = NULL, "to" = NULL, "TargetReaction" = NULL, "RegSign" = NULL, stringsAsFactors = F)


