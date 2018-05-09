library(tictoc)

setwd("~/winData/multiomics_networks_simulation")

source("network_generation.R")

tic(); mysystem = creationsystem(); toc()

tic()
for(i in 1:10) mysystem = creationsystem()
toc()
