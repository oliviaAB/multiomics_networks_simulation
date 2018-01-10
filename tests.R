rm(list = ls())

setwd("C:/Users/oangelin/OneDrive - Massey University/Documents/multiomics_networks_simulation")

library(bioPN)
dyn.load("propensity.dll")
source("data_simulation.R")

PartitionedLeaping <- function(model, timep, delta=1 , runs=1, ect=1e-9) {
  pre <- model$pre
  post <- model$post
  h <- model$h
  M <- model$M

  storage.mode(pre) <- storage.mode(post) <- storage.mode(runs) <- "integer"
  #  storage.mode(slow) <- "integer"
  storage.mode(delta) <- storage.mode(M) <- storage.mode(ect) <- "double"
  
  if (is.null(place <- model$place)) {
    place <- paste("P",1:ncol(pre),sep="")
  }
  if (is.null(transition <- model$transition)) {
    transition <- paste("T",1:nrow(pre),sep="")
  }
  
  .Call("PartitionedLeaping", pre, post, h, M, timep, delta, runs, place, transition, ect, parent.frame())
}

# Petri Network

model = list(
  place = c("RNA", "P_NA","P_A"),
  transition = c("Transcription","Translation","RNAdecay","ProteinNAdecay","ProteinAdecay", "ProteinActivation"),
  pre = matrix(c(0,0,0,
                 1,0,0,
                 1,0,0,
                 0,1,0,
                 0,0,1,
                 0,1,0), ncol = 3, byrow = T),
  post = matrix(c(1,0,0,
                  1,1,0,
                  0,0,0,
                  0,0,0,
                  0,0,0,
                  0,0,1), ncol = 3, byrow = T),
  h = list( "Transcription" = getNativeSymbolInfo("Transcription", PACKAGE="propensity")$address,
            "Translation" = getNativeSymbolInfo("Translation", PACKAGE="propensity")$address,
            "RNAdecay" = getNativeSymbolInfo("RNAdecay", PACKAGE="propensity")$address,
            "ProteinNAdecay" = getNativeSymbolInfo("ProteinNAdecay", PACKAGE="propensity")$address,
            "ProteinAdecay" = getNativeSymbolInfo("ProteinAdecay", PACKAGE="propensity")$address,
            "ProteinActivation" = getNativeSymbolInfo("ProteinActivation", PACKAGE="propensity")$address),
  M = c("RNA" = 100, "P_NA" = 0, "P_A" = 300),
  slow = c(0 ,0 ,0)
)

model1 = model
model1$slow = c(1 ,1 ,1)

delta = 1
timep = 200

HR = HaseltineRawlings(model, timep, delta)
HR1 = HaseltineRawlings(model1, timep, delta)

GB = GibsonBruck(model, timep, delta)

# Our method

tmax = timep

nw1 = rand_network_null(1,1,0)
nw1$TF_sgn[1] = -1
nw1$TF_th[1] = 1000
nw1$TF_n[1] = 2
nw1$TF_fc[1] = 1

nw1$k_TC[1] = 5
nw1$k_TL[1] = 1
nw1$p0_DR[1] = 1/50
nw1$p0_DP[1] = 1/120

cohort = rand_cohort_null(nw1, 1)
cohort$rna_0[1] = 100
cohort$prot_NA_0[1] = 0
cohort$prot_A_0[1] = 300

sim = simu(nw1, cohort, tmax)

parSSA = paramSSA(nw1, cohort)
res_simSSA = vector("list", length = 2); names(res_simSSA) = c(nw1$genes, nw1$prot)
simSSA = ssa(parSSA$x0, parSSA$a, parSSA$nu, parSSA$parms, tmax, method = "OTL")
for(g in nw1$genes){ res_simSSA[[g]] = sapply(0:tmax, function(i){simSSA$data[max(which(simSSA$data[,1]<=i)),g]}) }
for(p in nw1$prot){ res_simSSA[[p]] = sapply(0:tmax, function(i){
  temptime =  max(which(simSSA$data[,1]<=i))
  return(simSSA$data[temptime,paste0(p,"_NA")] + simSSA$data[temptime,paste0(p,"_A")])}) }


ymin = min(HR$run[[1]]$M[[1]], HR1$run[[1]]$M[[1]], sim$time_rna[[1]], GB$run[[1]]$M[[1]])
ymax = max(HR$run[[1]]$M[[1]], HR1$run[[1]]$M[[1]], sim$time_rna[[1]], GB$run[[1]]$M[[1]])
plot(HR$dt, HR$run[[1]]$M[[1]], type = 'none', col = "red", ylim = c(ymin, ymax))
lines(HR1$dt, HR1$run[[1]]$M[[1]], col = "red")
lines(GB$dt, GB$run[[1]]$M[[1]], col = "darkred", lty = 2)
lines(0:tmax, sim$time_rna[[1]], col = "hotpink")
lines(0:tmax, res_simSSA[[1]], col = "hotpink", lty = 2)

windows()
ymin = min((HR$run[[1]]$M[[2]]+HR$run[[1]]$M[[3]]), (HR1$run[[1]]$M[[2]]+HR1$run[[1]]$M[[3]]), sim$time_prot_tot[[1]], (GB$run[[1]]$M[[2]]+GB$run[[1]]$M[[3]]))
ymax = max((HR$run[[1]]$M[[2]]+HR$run[[1]]$M[[3]]), (HR1$run[[1]]$M[[2]]+HR1$run[[1]]$M[[3]]), sim$time_prot_tot[[1]], (GB$run[[1]]$M[[2]]+GB$run[[1]]$M[[3]]))
plot(HR$dt, (HR$run[[1]]$M[[2]]+HR$run[[1]]$M[[3]]), type = 'none', col = "green", ylim = c(ymin, ymax))
lines(HR1$dt, (HR1$run[[1]]$M[[2]]+HR1$run[[1]]$M[[3]]), col = "green")
lines(GB$dt, (GB$run[[1]]$M[[2]]+GB$run[[1]]$M[[3]]), col = "green", lty = 2)
lines(0:tmax, sim$time_prot_tot[[1]], col = "blue")
lines(0:tmax, res_simSSA[[2]], col = "blue", lty = 2)



#############################################################################################################




