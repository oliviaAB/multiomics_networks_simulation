library(tictoc)

setwd("~/winData/multiomics_networks_simulation")
setwd("~/GitHub/multiomics_networks_simulation")
source("network_generation.R")

mysystemargs = insilicosystemargs(G = 100)

mysystemgenes = createGenes(mysystemargs)

mysystem = createMultiOmicNetwork(mysystemgenes, mysystemargs)

tic()
for(i in 1:100){
  mysystemargs = insilicosystemargs(G = 100)
  mysystemgenes = createGenes(mysystemargs)
  mysystem = createMultiOmicNetwork(mysystemgenes, mysystemargs)
  }
toc()


####
sysargs = insilicosystemargs()
for(i in names(sysargs)){
  assign(i, sysargs[[i]])
}

