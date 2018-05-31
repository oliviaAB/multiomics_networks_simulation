library(tictoc)

setwd("~/winData/multiomics_networks_simulation")
#setwd("~/GitHub/multiomics_networks_simulation")



source("network_generation.R")


mysystemargs = insilicosystemargs(G = 10, RD.NC.outdeg.exp = 3)
insilicosystem = createInSilicoSystem(mysystemargs)


#plotGlobalSystem(myinsilicosystem, show = T)

#plotRegulationSystem(myinsilicosystem, c("TC"))

myindivargs = insilicoindividualargs()
insilicopopulation = createPopulation(5, insilicosystem, myindivargs)

simtime = 100
nepochs = 200
ntrialsPerInd = 1
simalgorithm = "SSA"
returnStochModel = F
ind = names(insilicopopulation$individualsList)[1]
library(viridis)

resTable = simulateSystemStochastic(insilicosystem, insilicopopulation, simtime, nepochs, ntrialsPerInd, simalgorithm , returnStochModel)

mycolsvariants = viridis(insilicopopulation$indargs$ngenevariants)
names(mycolsvariants) = as.character(1:insilicopopulation$indargs$ngenevariants)

for(ind in names(insilicopopulation$individualsList)){

    profiles = resTable[[ind]]
    for(g in insilicosystem$genes$id){
      RNAspecies = grep(paste0("^R",g,"GCN"), names(profiles), value = T)
      alleleVariants = sapply(RNAspecies, function(x){
        id = rev(strsplit(as.character(x), "GCN")[[1]])[1]
        allelevar = insilicopopulation$individualsList[[ind]]$haplotype[g, paste0("GCN", id)]
        return(as.character(allelevar))
      })
      
      colsAllelevariants = mycolsvariants[alleleVariants]
      names(colsAllelevariants) = names(alleleVariants)
      RNAprofile = melt(profiles, id.vars = c("time", "trial"), measure.vars = RNAspecies)
      
      plotRNA = ggplot(RNAprofile, aes(x = time, y = value, color = variable)) +
        geom_line() + 
        scale_color_manual(values = colsAllelevariants, breaks = names(colsAllelevariants), drop = F, name = "Allele variant") + 
        ggtitle(paste(c("Expression profile of gene", as.character(g), "- RNA"), collapse = " "))

      if(insilicosystem$genes[g, "coding"] == "PC"){
        colsAllelevariantsP = colsAllelevariants
        names(colsAllelevariantsP) = sub("^R", "P", names(colsAllelevariants))
        protspecies = grep(paste0("^P",g,"GCN"), names(profiles), value = T)
        protprofile = melt(profiles, id.vars = c("time", "trial"), measure.vars = protspecies)
        plotprot = ggplot(protprofile, aes(x = time, y = value, color = variable)) +
          geom_line() + 
          scale_color_manual(values = colsAllelevariantsP, breaks = names(colsAllelevariantsP), drop = F, name = "Allele variant") + 
          ggtitle(paste(c("Expression profile of gene", as.character(g), "- Protein"), collapse = " "))

        ggarrange(plotRNA, plotprot, ncol = 2, draw = T)
      }else{ggarrange(plotRNA, draw = T)}
  }
}  



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
