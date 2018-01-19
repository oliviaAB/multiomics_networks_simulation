##########################################################################################################################
##########################################################################################################################
###                                              PARAMETER FILE                                                        ###
##########################################################################################################################
##########################################################################################################################

# ---------------------------------------------------------------------------------------------------------------------- #
##                                         Network sampling parameters                                                  ##
# ---------------------------------------------------------------------------------------------------------------------- #

# Volume of a cell
Volcel = 10^(-15)


## Percentage of proteins that will have an active and an inactive state
param.actprot = 0.8 # Between 0 and 1


# Regulation probabilities ----
# proba_reg_TF = c(0.05,0.85,0.1) # probability of having a regulation of the sign resp. -1, 0, 1 in the transcription regulation matrix
# proba_reg_TLF = c(0.05,0.85,0.1) # probability of having a regulation of the sign resp. -1, 0, 1 in the translation regulation matrix
# proba_reg_DR = c(0.05,0.85,0.1) # probability of having a regulation of the sign resp. -1, 0, 1 in the RNA decay regulation matrix
# proba_reg_DP = c(0.05,0.85,0.1) # probability of having a regulation of the sign resp. -1, 0, 1 in the protein decay regulation matrix
# proba_reg_ACT = c(0.7,0.3) # probability of having a regulation of the sign resp. 0 and 1 in the protein activation regulation matrix
# proba_reg_DEACT = c(0.7,0.3) # probability of having a regulation of the sign resp. 0 and 1 in the protein deactivation regulation matrix

proba_reg_TF = c(0.05,0.75,0.2) # probability of having a regulation of the sign resp. -1, 0, 1 in the transcription regulation matrix
proba_reg_TLF = c(0.05,0.75,0.2) # probability of having a regulation of the sign resp. -1, 0, 1 in the translation regulation matrix
proba_reg_DR = c(0.05,0.75,0.2) # probability of having a regulation of the sign resp. -1, 0, 1 in the RNA decay regulation matrix
proba_reg_DP = c(0.05,0.75,0.2) # probability of having a regulation of the sign resp. -1, 0, 1 in the protein decay regulation matrix
proba_reg_ACT = c(0.7,0.3) # probability of having a regulation of the sign resp. 0 and 1 in the protein activation regulation matrix
proba_reg_DEACT = c(0.7,0.3) # probability of having a regulation of the sign resp. 0 and 1 in the protein deactivation regulation matrix


# Hill function parameters ----
#     Given as sampling functions to allow the user to control the distribution from which are sampled the parameters
#     Input parameter of each function is the sample size required

th_sampling_default = function(x){ sample(50:100, x, replace = T) } # default: th param sampled from a integer uniform function bewteen 50 and 100
th_sampling = "th_sampling_default"

n_sampling_default = function(x){ sample(1:4, x, replace = T) } # default: n param sampled from a integer uniform function bewteen 1 and 4
n_sampling = "n_sampling_default"

fc_sampling_default = function(x){ sample(1:10, x, replace = T) } # default: fc param sampled from a integer uniform function bewteen 1 and 10
fc_sampling = "fc_sampling_default"


# Enzymatic kinetic parameters ----
#     Given as sampling functions to allow the user to control the distribution from which are sampled the parameters
#     Input parameter of each function is the sample size required

enzparam_sampling_default = function(x){ # values for mean and sd used for the sampling function are from "Cell Biology by the numbers" (Milo and Phillips, 2016)
  logkcat = rnorm(x, 1, 1.1) # STEP 1: sample kcat
  kcat = 10^logkcat
  logratio = rtruncnorm(x, a = 0, b = 9, mean = 5, sd = 1.6) # STEP 2: sample kcat/KM
  KM = kcat/(10^logratio) # STEP 3: compute KM using the sampled values of kcat and ratio
  KM = round(KM * 6 * 10^23 * Volcel)
  return(rbind(kcat,KM))
}
enzparam_sampling = "enzparam_sampling_default"


# Basal rates ----
#     Given as sampling functions to allow the user to control the distribution from which are sampled the parameters
#     Input parameter of each function is the sample size required

basal_transcription_rate_default = function(x){ runif(x, 1, 5) } # default: basal TC rate chosen from a uniform distribution ranging from 0.01 to 0.1
basal_transcription_rate = "basal_transcription_rate_default"

basal_translation_rate_default = function(x){ runif(x, 0.5, 5) } # default: basal TL rate chosen from a uniform distribution ranging from 0.5 to 5
basal_translation_rate = "basal_translation_rate_default"

# basal_RNAdecay_rate_default = function(x){ runif(x, 0.005, 0.01) } # default: basal RNA decay rate chosen from a uniform distribution ranging from 0.005 to 0.1
# basal_RNAdecay_rate = "basal_RNAdecay_rate_default"
# 
# basal_proteindecay_rate_default = function(x){ runif(x, 0.01, 0.1) } # default: basal protein decay rate chosen from a uniform distribution ranging from 0.01 to 0.1
# basal_proteindecay_rate = "basal_proteindecay_rate_default"


basal_RNAlifetime_default = function(x){ sample(60:3600, x, replace = T) } # default: basal RNA lifetime chosen from a discrete uniform distribution ranging from 1 minutes to 1 hour (in seconds)
basal_RNAlifetime = "basal_RNAlifetime_default"

basal_protlifetime_default = function(x){ sample(5400:14400, x, replace = T) } # default: basal protein lifetime chosen from a discrete uniform distribution ranging from 90 minutes to 4 hours (in seconds)
basal_protlifetime = "basal_protlifetime_default"


# --------------------------------- #
##  Individual sampling parameters ##
# --------------------------------- #

# Genotype effects ----
#     Given as sampling functions to allow the user to control the distribution from which are sampled the parameters
#     Input parameter of each function is the sample size required

qtl_effect_transcription_default = function(x){ rtruncnorm(x , a = 0, b = Inf , mean = 1, sd = 0.05) } # default: QTL effect on transcription rate sampled from a normal distribution with mean 1 and sd 0.05   
# rtruncnorm ensures that the sampled values are positive
qtl_effect_transcription = "qtl_effect_transcription_default"

qtl_effect_translation_default = function(x){ rtruncnorm(x , a = 0, b = Inf , mean = 1, sd = 0.05) } # default: QTL effect on translation rate sampled from a normal distribution with mean 1 and sd 0.05   
# rtruncnorm ensures that the sampled values are positive
qtl_effect_translation = "qtl_effect_translation_default"

qtl_effect_RNAdecay_default = function(x){ rtruncnorm(x , a = 0, b = Inf , mean = 1, sd = 0.05) } # default: QTL effect on RNA decay rate sampled from a normal distribution with mean 1 and sd 0.05   
# rtruncnorm ensures that the sampled values are positive
qtl_effect_RNAdecay = "qtl_effect_RNAdecay_default"

qtl_effect_enzact_default = function(x){ rtruncnorm(x , a = 0, b = Inf , mean = 1, sd = 0.05) } # default: QTL effect on RNA decay rate sampled from a normal distribution with mean 1 and sd 0.05   
# rtruncnorm ensures that the sampled values are positive
qtl_effect_enzact = "qtl_effect_enzact_default"


# Initial abundance ----
#     Given as sampling functions to allow the user to control the distribution from which are sampled the parameters
#     Input parameter of each function is the number of different molecules x

initial_rna_abundance_default = function(x){ sample(0:100, size = x, replace = T)  }
initial_rna_abundance = "initial_rna_abundance_default"

initial_prot_abundance_default = function(x){ sample(10:500, size = x, replace = T)  }
initial_prot_abundance = "initial_prot_abundance_default"

initial_met_abundance_default = function(x){ sample(1000:5000, size = x, replace = T)  }
initial_met_abundance = "initial_met_abundance_default"

