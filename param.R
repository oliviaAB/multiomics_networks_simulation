##########################################################################################################################
#                                           PARAMETERS FOR NETWORK SIMULATION                                            #
##########################################################################################################################

require("truncnorm")
library(truncnorm)

# ----------------------------- #
## Network sampling parameters ##
# ----------------------------- #

# Regulation probabilities ----
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


# Basal rates ----
#     Given as sampling functions to allow the user to control the distribution from which are sampled the parameters
#     Input parameter of each function is the sample size required

basal_transcription_rate_default = function(x){ runif(x, 1, 5) } # default: basal TC rate chosen from a uniform distribution ranging from 0.01 to 0.1
basal_transcription_rate = "basal_transcription_rate_default"

basal_translation_rate_default = function(x){ runif(x, 0.5, 5) } # default: basal TL rate chosen from a uniform distribution ranging from 0.5 to 5
basal_translation_rate = "basal_translation_rate_default"

basal_RNAdecay_rate_default = function(x){ runif(x, 0.005, 0.01) } # default: basal RNA decay rate chosen from a uniform distribution ranging from 0.005 to 0.1
basal_RNAdecay_rate = "basal_RNAdecay_rate_default"

basal_proteindecay_rate_default = function(x){ runif(x, 0.01, 0.1) } # default: basal protein decay rate chosen from a uniform distribution ranging from 0.01 to 0.1
basal_proteindecay_rate = "basal_proteindecay_rate_default"

# ----



# ----------------------------- #
##  Cohort sampling parameters ##
# ----------------------------- #

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


# Initial abundance
#     Given as sampling functions to allow the user to control the distribution from which are sampled the parameters
#     Input parameter of each function is the number of different molecules m and the number of individuals i
initial_abundance_default = function(m, i){
  mol_mean = round(sample(0:500, m, replace = T))                                       # each molecule has an mean initial abundance across the individuals sampled from a discrete uniform distribution ranging from 0 to 500 
  init = round(rtruncnorm(m*i, a = 0, b = Inf, mean = mol_mean, sd = round(mol_mean/10)))    # for each individual the initial abundance of a given molecule is sampled from a discrete positive normal distribution with corresponding mean and sd = mean/10
                                                                                        # so non abundant molecules have a smaller initial variation around the mean while abundant molecules display a larger variation across individuals
  return(init)
}
initial_abundance = "initial_abundance_default"


# ----