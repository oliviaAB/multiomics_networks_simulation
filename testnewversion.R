
####################################################################################################################################################3
####################################################################################################################################################

setwd("~/winData/multiomics_networks_simulation")
#setwd("~/GitHub/multiomics_networks_simulation")

source("network_generation.R")

#anotherev = newJuliaEvaluator()
mysystemargs = insilicosystemargs(G = 3, RD.NC.outdeg.exp = 3, PC.PTM.p = 0.5)
insilicosystem = createInSilicoSystem(mysystemargs, empty = T)


# plotGlobalSystem(insilicosystem, show = T)

# plotRegulationSystem(insilicosystem, show = T)

myindivargs = insilicoindividualargs()
insilicopopulation = createPopulation(15, insilicosystem, myindivargs)

tic()
res = simulateSystemStochastic(insilicosystem, insilicopopulation, simtime = 3600, nepochs = 20, ntrialsPerInd = 1, simalgorithm = "ODM", returnStochModel = F)
toc()

tic()
res2 = simulateSystemStochasticParallel(insilicosystem, insilicopopulation, simtime = 1, nepochs = 20, ntrialsPerInd = 1, simalgorithm = "ODM", returnStochModel = F)
toc()


sapply(1:length(insilicopopulation$individualsList), function(i){identical(res$resTable[[i]][1,], res2[[i]][1,])})

# resTable = res$resTable

plotExpressionProfiles(insilicosystem, insilicopopulation, res2)


myfunc = function(i){
  
  # myev = newJuliaEvaluator()
  myev = evList[[i]]
  mysystemargs = insilicosystemargs(G = 3, RD.NC.outdeg.exp = 3, PC.PTM.p = 0.5)
  insilicosystem = createInSilicoSystem(mysystemargs, ev = myev)
  
  myindivargs = insilicoindividualargs()
  insilicopopulation = createPopulation(2, insilicosystem, myindivargs)
  
  res = simulateSystemStochastic(insilicosystem, insilicopopulation, simtime = 1, nepochs = 20, ntrialsPerInd = 1, simalgorithm = "ODM", returnStochModel = F, ev = myev)
  #removeJuliaEvaluator(myev)
  return(list("simID" = i, "insilicosystem" = insilicosystem, "insilicopopulation" = insilicopopulation, "res" = res))
}

nsim = 100
evList = sapply(1:50, function(x){
  myev = newJuliaEvaluator()
  print(myev)
  print(showConnections())
  return(myev)})
test = mclapply(1:nsim, myfunc, mc.cores = (detectCores()-1))
sapply(evList, removeJuliaEvaluator)
# ------------

setwd("~/winData/multiomics_networks_simulation")

juliaCommand("
if !haskey(Pkg.installed(), \"ClobberingReload\") 
	Pkg.clone(\"git://github.com/cstjean/ClobberingReload.jl.git\")
end
")

juliaCommand("addprocs(1)")
juliaCommand("@everywhere sinclude(\"julia_functions.jl\")")
# ------------

setwd("~/winData/multiomics_networks_simulation")
source("network_generation.R")

load("/home/oangelin/Documents/noerror.RData")
stochmodel = createStochSystem(insilicosystem, insilicopopulation$indargs, returnList = F)

evaluator = XR::getInterface(getClass("JuliaInterface"))
expr = gettextf("%s(%s)","stochasticsimulation", evaluator$ServerArglist(stochmodel$JuliaObject, insilicopopulation$individualsList[[1]]$QTLeffects, insilicopopulation$individualsList[[1]]$InitVar, df2list(insilicosystem$genes), 0.00001, modelname = "Ind1", ntrials = 1, nepochs = 1, simalgorithm = "SSA"))
#expr = "jat()"
key = RJulia()$ProxyName()
cmd = jsonlite::toJSON(c("eval", expr, key, T))
writeLines(cmd, evaluator$connection)
#evaluator$ServerQuit()
for(try in 1:10) {
  value <- readLines(evaluator$connection, 1)
  #print(value)
  if(length(value) == 0)  # But shouldn't happen?
    Sys.sleep(1)
  else
    break
}
test = XR::valueFromServer(value, key, T, evaluator)


# ------------

mysystemargs = insilicosystemargs(G = 10, PC.PTM.form.p = 0)
myinsilicosystemEmpty = createInSilicoSystem(mysystemargs, empty = T)


plotGlobalSystem(myinsilicosystemEmpty, show = T)

plotRegulationSystem(myinsilicosystem, c("TC"))

myindivargs = insilicoindividualargs()
mypopulation = createPopulation(20, myinsilicosystemEmpty, myindivargs)

test = createStochSystem(myinsilicosystemEmpty, myindivargs, returnList = T)



##############################################
library(tictoc)

setwd("~/winData/multiomics_networks_simulation")
#setwd("~/GitHub/multiomics_networks_simulation")
source("network_generation.R")

mysystemargs = insilicosystemargs(G = 500)
tic(); myinsilicosystem = createInSilicoSystem(mysystemargs); toc()

## plotmosystem

test = plotGlobalSystem(myinsilicosystem, show = T)

ggsave("/media/sf_data/globalPanel1.png", plot = test$globalPanel1, width = 33.9, height = 19.1, units = "cm")
ggsave("/media/sf_data/globalPanel2.png", plot = test$globalPanel2, width = 33.9, height = 19.1, units = "cm")
ggsave("/media/sf_data/globalPanel3.png", plot = test$globalPanel3, width = 33.9, height = 19.1, units = "cm")
ggsave("/media/sf_data/globalPanel4.png", plot = test$globalPanel4, width = 33.9, height = 19.1, units = "cm")


test2 = plotRegulationSystem(myinsilicosystem, c("TC"))
ggsave("/media/sf_data/TCPanel1.png", plot = test2$TCPanel1, width = 33.9, height = 19.1, units = "cm")
ggsave("/media/sf_data/TCPanel2.png", plot = test2$TCPanel2, width = 33.9, height = 19.1, units = "cm")
ggsave("/media/sf_data/TCPanel3.png", plot = test2$TCPanel3, width = 33.9, height = 19.1, units = "cm")




##############################################
library(tictoc)

setwd("~/winData/multiomics_networks_simulation")
#setwd("~/GitHub/multiomics_networks_simulation")
source("network_generation.R")

nsim = 100
sizesim = c(50, 100, 200, 500)
size_telapsed = vector("list", length(sizesim))
size_simnod = vector("list", length(sizesim))
size_simedg = vector("list", length(sizesim))

for(s in sizesim){
  print(s)
    mysystemargs = insilicosystemargs(G = s)
    telapsed = vector("numeric", nsim)
    simnod = vector("list", nsim)
    simedg = vector("list", nsim)
    
    for(i in 1:nsim){
      tic(); myinsilicosystem = createInSilicoSystem(mysystemargs); itel = toc(quiet = T)
      telapsed[i] = itel$toc - itel$tic
      simnod[[i]] = myinsilicosystem$genes
      simedg[[i]] = myinsilicosystem$mosystem$edg
    }
    
    size_telapsed[[s]] = telapsed
    size_simnod[[s]] = simnod
    size_simedg[[s]] = simedg
    
}  

#plot(rep(sizesim, rep(nsim, length(sizesim))), unlist(size_telapsed), xlab = "System size (G: number of nodes)", ylab = c("Running time (sec)"))
runningtime = data.frame("G" = rep(sizesim, rep(nsim, length(sizesim))), "runningtime" = unlist(size_telapsed))

runtime = ggplot(runningtime, aes(x = G, y = runningtime)) + geom_point() + xlab("System size (number of genes)") + ylab("Running time (s)")
ggsave("/media/sf_data/runtime.png", plot = runtime, width = 33.9, height = 19.1, units = "cm")



ressimdf = data.frame("simid" = numeric(), "G" = numeric(), "ratioPC" = numeric(), "E" = numeric())

for(s in sizesim){
  ratioPC = sapply(size_simnod[[s]], function(x){sum(x$coding =="PC")/nrow(x)})
  E = sapply(size_simedg[[s]], nrow)
  ressimdf = rbind(ressimdf, data.frame("simid" = 1:nsim, "G" = rep(s, nsim), "ratioPC" = ratioPC, "E" = E))
}

g1 = ggplot(ressimdf, aes(x = ratioPC)) + geom_histogram() + facet_grid(G~.) + xlab("Ratio of protein coding genes in the system")
g2 = ggplot(ressimdf, aes(x = E)) + geom_histogram() + facet_grid(G~.) + annotate("text", x = 5000, y = 70, label = sapply(sizesim, function(x){paste("mean =", mean(ressimdf[ressimdf$G == x, "E"]), sep = " ")})) + xlab("Number of regulatory interactions")

sevsimplot2 = ggarrange(g1, g2, ncol = 2)
ggsave("/media/sf_data/sevsimplot2.png", plot = sevsimplot2, width = 33.9, height = 19.1, units = "cm")


#################################################################################################################
#################################################################################################################

setwd("~/winData/multiomics_networks_simulation")
#setwd("~/GitHub/multiomics_networks_simulation")

source("network_generation.R")


mybaseport = RJulia()$port

cat("Starting simulations at", format(Sys.time(), usetz = T), "\n")

tic()
test = mclapply(1:20, function(i){
  myev = newJuliaEvaluator(port = mybaseport + i) ## create a new Julia evaluator with a port number equal to mybaseport + i (id of the simulation)
  sleeprand = sample(1:50, 1)
  #print(paste0("Port ", mybaseport+i," sleeping ", sleeprand, " seconds\n"))
  juliaCommand("println(\"port %s sleeping for %s seconds\")", mybaseport + i, sleeprand, evaluator = myev)
  juliaCommand("sleep(%s)", sleeprand, evaluator = myev)
  juliaCommand("println(\"port %s done\")", mybaseport + i, evaluator = myev)
  #print(paste0("Port ", mybaseport+i," done!"))
  removeJuliaEvaluator(myev)
  return(sleeprand)
}, mc.cores = detectCores()-1)
toc()


