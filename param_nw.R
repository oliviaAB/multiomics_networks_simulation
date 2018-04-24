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

#### TC reg. network properties ----

## TC.outdeg.power : power of the power-law distribution for the out-degree of the transcription graph
TC.outdeg.power = 2.2
## TC.indeg.exp : power of the exponent distribution for the in-degree of the transcription graph
TC.indeg.exp = 0.6

# ----
