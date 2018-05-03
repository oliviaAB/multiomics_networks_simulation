# ------------------------------------------------------------------------------------------------------------ #
#                                     PARAMETERS FOR NETWORK GENERATION                                        #
# ------------------------------------------------------------------------------------------------------------ # 


## G : number of genes in the system
G = 50

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

##TC.PC.pos.p: probability that the transcription is positively regulated by protein regulators
##TC.NC.pos.p: probability that the transcription is positively regulated by noncoding regulators
##TL.PC.pos.p: probability that the translation is positively regulated by protein regulators
##TL.NC.pos.p: probability that the translation is positively regulated by noncoding regulators
##RD.PC.pos.p: probability that the RNA decay is positively regulated by protein regulators (faster decay)
##RD.NC.pos.p: probability that the RNA decay is positively regulated by noncoding regulators (faster decay)
##PD.PC.pos.p: probability that the protein decay is positively regulated by protein regulators (faster decay)
##PD.NC.pos.p: probability that the protein decay is positively regulated by noncoding regulators (faster decay)
##PTM.PC.pos.p: probability that the protein decay is is activated by the post-translational modification by protein regulators
##PTM.NC.pos.p: probability that the protein decay is is activated by the post-translational modification by noncoding regulators


TC.PC.pos.p = 0.5 
TC.NC.pos.p = 0.5 
TL.PC.pos.p = 0.5 
TL.NC.pos.p = 0.05 
RD.PC.pos.p = 0.7 
RD.NC.pos.p = 0.9 
PD.PC.pos.p = 0.5
PD.NC.pos.p = 0.5 
PTM.PC.pos.p = 0.5
PTM.NC.pos.p = 0.5


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


## Form of the distribution of the number of targets of transcription factors (can be either "powerlaw" or "exponential")
TC.PC.outdeg.distr = "powerlaw"
## Form of the the distribution of the number of targets of noncoding RNAs regulating transcription (can be either "powerlaw" or "exponential")
TC.NC.outdeg.distr = "powerlaw"

## exponent of the distribution for the out-degree of the TFs in the transcription graph
TC.PC.outdeg.exp = 2.2
## exponent of the distribution for the out-degree of the noncoding RNAs in the transcription graph
TC.NC.outdeg.exp = 1
## Type of preferential attachment for the targets of TFs in the transcription graph
TC.PC.indeg.distr = "exponential"
## Type of preferential attachment for the targets of ncRNAs in the transcription graph
TC.NC.indeg.distr = "powerlaw"

## Probability of TFs to perform autoregulation
TC.PC.autoregproba = 0.2
## Probability of ncRNAs to perform autoregulation
TC.NC.autoregproba = 0

## Are 2-nodes loops authorised in the transcription network with protein regulators?
TC.PC.twonodesloop = FALSE
## Are 2-nodes loops authorised in the transcription network with noncoding regulators?
TC.NC.twonodesloop = FALSE


#### Distribution of the different kinetic parameters of transcription regulatory reactions ----
TCbindingrate_default = function(x){ runif(x, 0.001, 0.01) } ## default: transcription factor binding rate chosen from a uniform distribution ranging from ??????
TCbindingrate = "TCbindingrate_default"

TCunbindingrate_default = function(x){ runif(x, 0.001, 0.01) } ## default: transcription factor unbinding rate chosen from a uniform distribution ranging from ??????
TCunbindingrate = "TCunbindingrate_default"

TCFoldChange_default = function(x){ sample(2:30, x, replace = T)  } ## default: transcription factor fold change from a discrete uniform distribution ranging from ??????
TCFoldChange = "TCFoldChange_default"





#### TL reg. network properties ----


## Form of the distribution of the number of targets of translation factors (can be either "powerlaw" or "exponential")
TL.PC.outdeg.distr = "powerlaw"
## Form of the the distribution of the number of targets of noncoding RNAs regulating translation (can be either "powerlaw" or "exponential")
TL.NC.outdeg.distr = "powerlaw"

## exponent of the distribution for the out-degree of the TLFs in the translation graph
TL.PC.outdeg.exp = 3
## exponent of the distribution for the out-degree of the noncoding RNAs in the translation graph
TL.NC.outdeg.exp = 1
## Type of preferential attachment for the targets of TLFs in the translation graph
TL.PC.indeg.distr = "exponential"
## Type of preferential attachment for the targets of ncRNAs in the translation graph
TL.NC.indeg.distr = "powerlaw"

## Probability of TLFs to perform autoregulation
TL.PC.autoregproba = 0.2
## Probability of ncRNAs to perform autoregulation
TL.NC.autoregproba = 0

## Are 2-nodes loops authorised in the translation network with protein regulators?
TL.PC.twonodesloop = FALSE
## Are 2-nodes loops authorised in the translation network with noncoding regulators?
TL.NC.twonodesloop = FALSE


#### Distribution of the different kinetic parameters of translation regulatory reactions ----
TLbindingrate_default = function(x){ runif(x, 0.001, 0.01) } ## default: translation factor binding rate chosen from a uniform distribution ranging from ??????
TLbindingrate = "TLbindingrate_default"

TLunbindingrate_default = function(x){ runif(x, 0.001, 0.01) } ## default: translation factor unbinding rate chosen from a uniform distribution ranging from ??????
TLunbindingrate = "TLunbindingrate_default"

TLFoldChange_default = function(x){ sample(2:30, x, replace = T)  } ## default: translation factor fold change from a discrete uniform distribution ranging from ??????
TLFoldChange = "TLFoldChange_default"


#### RD reg. network properties ----


## Form of the distribution of the number of targets of protein regulating RNA decay (can be either "powerlaw" or "exponential")
RD.PC.outdeg.distr = "exponential"
## Form of the the distribution of the number of targets of noncoding RNAs regulating RNA decay (can be either "powerlaw" or "exponential")
RD.NC.outdeg.distr = "powerlaw"

## exponent of the distribution for the out-degree of the regulatory proteins in the RNA decay graph
RD.PC.outdeg.exp = 2
## exponent of the distribution for the out-degree of the noncoding RNAs in the RNA decay graph
RD.NC.outdeg.exp = 2
## Type of preferential attachment for the targets of regulatory proteins in the RNA decay graph
RD.PC.indeg.distr = "exponential"
## Type of preferential attachment for the targets of ncRNAs in the RNA decay graph
RD.NC.indeg.distr = "powerlaw"

## Probability of regulatory proteins to perform autoregulation
RD.PC.autoregproba = 0.2
## Probability of ncRNAs to perform autoregulation
RD.NC.autoregproba = 0

## Are 2-nodes loops authorised in the RNA decay network with protein regulators?
RD.PC.twonodesloop = FALSE
## Are 2-nodes loops authorised in the RNA decay network with noncoding regulators?
RD.NC.twonodesloop = FALSE


#### Distribution of the different kinetic parameters of RNA decay regulatory reactions ----
RDbindingrate_default = function(x){ runif(x, 0.001, 0.01) } ## default: regulator binding rate chosen from a uniform distribution ranging from ??????
RDbindingrate = "RDbindingrate_default"

RDunbindingrate_default = function(x){ runif(x, 0.001, 0.01) } ## default: regulator unbinding rate chosen from a uniform distribution ranging from ??????
RDunbindingrate = "RDunbindingrate_default"



#### PD reg. network properties ----


## Form of the distribution of the number of targets of protein regulating protein decay (can be either "powerlaw" or "exponential")
PD.PC.outdeg.distr = "exponential"
## Form of the the distribution of the number of targets of noncoding RNAs regulating protein decay (can be either "powerlaw" or "exponential")
PD.NC.outdeg.distr = "powerlaw"

## exponent of the distribution for the out-degree of the regulatory proteins in the protein decay graph
PD.PC.outdeg.exp = 2
## exponent of the distribution for the out-degree of the noncoding RNAs in the protein decay graph
PD.NC.outdeg.exp = 2
## Type of preferential attachment for the targets of regulatory proteins in the protein decay graph
PD.PC.indeg.distr = "exponential"
## Type of preferential attachment for the targets of ncRNAs in the protein decay graph
PD.NC.indeg.distr = "powerlaw"

## Probability of regulatory proteins to perform autoregulation
PD.PC.autoregproba = 0.2
## Probability of ncRNAs to perform autoregulation
PD.NC.autoregproba = 0

## Are 2-nodes loops authorised in the protein decay network with protein regulators?
PD.PC.twonodesloop = FALSE
## Are 2-nodes loops authorised in the protein decay network with noncoding regulators?
PD.NC.twonodesloop = FALSE


#### Distribution of the different kinetic parameters of protein decay regulatory reactions ----
PDbindingrate_default = function(x){ runif(x, 0.001, 0.01) } ## default: regulator binding rate chosen from a uniform distribution ranging from ??????
PDbindingrate = "PDbindingrate_default"

PDunbindingrate_default = function(x){ runif(x, 0.001, 0.01) } ## default: regulator unbinding rate chosen from a uniform distribution ranging from ??????
PDunbindingrate = "PDunbindingrate_default"


# ----
