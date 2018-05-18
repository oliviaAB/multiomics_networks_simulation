library(tictoc)

setwd("~/winData/multiomics_networks_simulation")

source("network_generation.R")

mysystemargs = insilicosystemargs()
mysystem = creationsystem(mysystemargs)

tic()
for(i in 1:10) mysystem = creationsystem()
toc()


####
sysargs = insilicosystemargs()
for(i in names(sysargs)){
  assign(i, sysargs[[i]])
}

