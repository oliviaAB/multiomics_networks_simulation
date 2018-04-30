# ------------------------------------------------------------------------------------------------------------ #
#                                     PARAMETERS FOR NETWORK GENERATION                                        #
# ------------------------------------------------------------------------------------------------------------ # 


## G : number of genes in the system
G = 100

## M : number of metabolites in the system
M = 5

## Mod : number of modules in the system
Mod = 15

#### Ratio of the different gene functions in the system ----

## NC : non coding genes (ratio NC.p)
## NC.TC (targeting transcription) (ratio among noncoding genes NC.TC.p)
## NC.TL (targeting translation) (ratio among noncoding genes NC.TL.p)
## NC.RD (targeting RNA decay) (ratio among noncoding genes NC.RD.p)
## NC.PTM (targeting post-translational modification) (ratio among noncoding genes NC.PTM.p)
NC.p = 0.3

NC.TC.p = 0.3
NC.TL.p = 0.35
NC.RD.p = 0.3
NC.PTM.p = 0.05

(NC.TC.p + NC.TL.p + NC.RD.p + NC.PTM.p) == 1

## PC: protein coding genes (ratio PC.p)
## PC.TC/TF (transcription factor targeting transcription) (ratio among noncoding genes PC.TC.p)
## PC.TL/TLF (targeting translation) (ratio among noncoding genes PC.TL.p)
## PC.RD (targeting RNA decay) (ratio among noncoding genes PC.RD.p)
## PC.PD (targeting protein decay) (ratio among noncoding genes PC.PD.p)
## PC.PTM (targeting post-translational modification) (ratio among noncoding genes PC.PTM.p)
## PC.MR (targeting metabolic reaction) (ratio among noncoding genes PC.MR.p)

PC.p = 0.7

PC.TC.p = 0.4
PC.TL.p = 0.3
PC.RD.p = 0.1
PC.PD.p = 0.1
PC.PTM.p = 0.05
PC.MR.p = 0.05

(NC.p + PC.p) == 1
(PC.TC.p + PC.TL.p + PC.RD.p + PC.PD.p + PC.PTM.p + PC.MR.p) == 1


## For protein coding genes, ratio of protein having a PTM form

PC.PTM.form.p = 0.6

#### Ratio of the different type (activation/repression) of each reaction ----

##TC.pos.p: probability that the transcription is positively regulated
##TL.pos.p: probability that the translation is positively regulated
##RD.pos.p: probability that the RNA decay is positively regulated (faster decay)
##PD.pos.p: probability that the protein decay is positively regulated (faster decay)
##PTM.pos.p: probability that the protein is activated by the post-translational modification

TC.pos.p = 0.5 
TL.pos.p = 0.5
RD.pos.p = 0.6
PD.pos.p = 0.6
PTM.pos.p = 0.5


#### Distribution of the different kinetic parameters of genes ----

##     Given as sampling functions to allow the user to control the distribution from which are sampled the parameters
##     Input parameter of each function is the sample size required

basal_transcription_rate_default = function(x){ runif(x, 0.01, 0.1) } ## default: basal TC rate chosen from a uniform distribution ranging from 0.01 to 0.1
basal_transcription_rate = "basal_transcription_rate_default"

basal_translation_rate_default = function(x){ runif(x, 0.5, 5) } ## default: basal TL rate chosen from a uniform distribution ranging from 0.5 to 5
basal_translation_rate = "basal_translation_rate_default"

# basal_RNAdecay_rate_default = function(x){ runif(x, 0.005, 0.01) } ## default: basal RNA decay rate chosen from a uniform distribution ranging from 0.005 to 0.1
# basal_RNAdecay_rate = "basal_RNAdecay_rate_default"
# 
# basal_proteindecay_rate_default = function(x){ runif(x, 0.01, 0.1) } ## default: basal protein decay rate chosen from a uniform distribution ranging from 0.01 to 0.1
# basal_proteindecay_rate = "basal_proteindecay_rate_default"


basal_RNAlifetime_default = function(x){ sample(60:3600, x, replace = T) } ## default: basal RNA lifetime chosen from a discrete uniform distribution ranging from 1 minutes to 1 hour (in seconds)
basal_RNAlifetime = "basal_RNAlifetime_default"

basal_protlifetime_default = function(x){ sample(5400:14400, x, replace = T) } ## default: basal protein lifetime chosen from a discrete uniform distribution ranging from 90 minutes to 4 hours (in seconds)
basal_protlifetime = "basal_protlifetime_default"



#### TC reg. network properties ----

## TC.outdeg.power : power of the power-law distribution for the out-degree of the transcription graph
TC.outdeg.power = 2.2
## TC.indeg.exp : power of the exponent distribution for the in-degree of the transcription graph
# TC.indeg.exp = 0.6


#### Distribution of the different kinetic parameters of transcription regulatory reactions ----
TCbindingrate_default = function(x){ runif(x, 0.001, 0.01) } ## default: transcription factor binding rate chosen from a uniform distribution ranging from ??????
TCbindingrate = "TCbindingrate_default"

TCunbindingrate_default = function(x){ runif(x, 0.001, 0.01) } ## default: transcription factor unbinding rate chosen from a uniform distribution ranging from ??????
TCunbindingrate = "TCunbindingrate_default"

TCFoldChange_default = function(x){ sample(2:30, x, replace = T)  } ## default: transcription factor fold change from a discrete uniform distribution ranging from ??????
TCFoldChange = "TCFoldChange_default"


# ----
