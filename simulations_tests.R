setwd("~/winData/multiomics_networks_simulation")
rm(list = ls(all = T))
#setwd("~/GitHub/multiomics_networks_simulation")

source("network_generation.R")
library(reshape2)


## TEMPLATE
#' when getting a list of dataframes with simulated expression profiles
#' Plot for each molecule the expression profile of the different individuals

plotexpprof = function(res){
  mols = setdiff(colnames(res[[1]]), c("time", "trial"))
  
  molsList = lapply(mols, function(m){
    molsummary = vector()
    for(ind in names(res)){
      resInd = res[[ind]]
      timexTrial = matrix(resInd[, m], ncol = max(resInd$trial), nrow = length(unique(resInd$time)))
      meanTrials = rowMeans(timexTrial)
      quantile1Trials = sapply(1:nrow(timexTrial), function(x){quantile(timexTrial[x, ], probs = 0.25)})
      quantile3Trials = sapply(1:nrow(timexTrial), function(x){quantile(timexTrial[x, ], probs = 0.975)})
      #molsummary = rbind(molsummary, data.frame("time" = resInd$time[1:nrow(timexTrial)], "ind" = rep(ind, nrow(timexTrial)), "mol" = rep(m, nrow(timexTrial)), "mean" = meanTrials, "quantile1" = quantile1Trials, "quantile3" = quantile3Trials))
      molsummary = rbind(molsummary, data.frame("time" = resInd$time[1:nrow(timexTrial)], "ind" = rep(ind, nrow(timexTrial)), "mean" = meanTrials, "quantile1" = quantile1Trials, "quantile3" = quantile3Trials))
    }
    return(molsummary)
  })
  names(molsList) = mols
  
  mylinetype = c("mean" = "solid", "quantile1" = "longdash", "quantile3" = "longdash")
  
  molsplot = lapply(mols, function(m){
    toplot = melt(molsList[[m]], id.vars = c("time", "ind"), measure.vars = c("mean", "quantile1", "quantile3"))
    myplot = ggplot() + geom_line(data = toplot, aes(x = time, y = value, color = ind, linetype = variable)) +
      geom_ribbon(data = molsList[[m]], aes(x = time, ymin=quantile1, ymax=quantile3, fill = ind), alpha=0.3) +
      scale_linetype_manual(values = mylinetype, guide = F) +
      ggtitle(m) + ylab("Molecule abundance")
    print(myplot)
    return(myplot)
  })
  
  return(molsplot)
}


getlasttimepoint = function(res){
  
  mols = setdiff(colnames(res[[1]]), c("time", "trial"))
  
  molsList.hist = lapply(mols, function(m){
    molsummary = vector()
    for(ind in names(res)){
      resInd = res[[ind]]
      timexTrial = matrix(resInd[, m], ncol = max(resInd$trial), nrow = length(unique(resInd$time)))
      molsummary = cbind(molsummary, timexTrial[nrow(timexTrial),])
    }
    colnames(molsummary) = names(res)
    return(molsummary)
  })
  names(molsList.hist) = mols
  
  return(molsList.hist)
  
}

# ------------------------------------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------- #
## 1: Create an system with only 1 protein-coding gene, no regulation
##    In the population each individual has only 1 QTL value !=1 to show the impact
##    of the different QTLs on the expression profiles
# ------------------------------------------------------------------------------------- #


mysystemargs = insilicosystemargs(G = 1, PC.p = 1)
insilicosystem = createInSilicoSystem(mysystemargs, empty = T)

myindivargs = insilicoindividualargs(ploidy = 1, ngenevariants = 1)
insilicopopulation = createPopulation(5, insilicosystem, myindivargs, sameInit = T)

## Manually change the value of the QTLs for each individual
qtlvals = c("qtlTCrate", "qtlTLrate", "qtlRDrate", "qtlPDrate")

for(i in 2:5){
  insilicopopulation$individualsList[[i]]$QTLeffects$GCN1[[qtlvals[i-1]]] = 0.5 ## set the corresponding QTL value to 0.5
}

res = simulateSystemStochasticParallel(insilicosystem, insilicopopulation, simtime = 5000, nepochs = 20, ntrialsPerInd = 50, simalgorithm = "SSA", returnStochModel = F)

molsplot = plotexpprof(res)

ggsave("/media/sf_data/RNAplot.png", plot = molsplot[[1]])
ggsave("/media/sf_data/Protplot.png", plot = molsplot[[2]])

save(insilicosystem, insilicopopulation, file = "/home/oangelin/Documents/plots05_07_2018.RData")

# ------------------------------------------------------------------------------------- #
## 2: Create an system with only 1 protein-coding gene, no regulation, ploidy of 2
##    2 alleles: 1 with default transcription rate, second with 0.5 transcription rate
# ------------------------------------------------------------------------------------- #

mysystemargs = insilicosystemargs(G = 1, PC.p = 1)
insilicosystem = createInSilicoSystem(mysystemargs, empty = T)

myindivargs = insilicoindividualargs(ploidy = 2, ngenevariants = 1)
insilicopopulation = createPopulation(3, insilicosystem, myindivargs, sameInit = T)

## Manually change the value of the QTLs for each individual

insilicopopulation$individualsList[[2]]$QTLeffects$GCN2[["qtlTCrate"]] = 0.5 ## set the corresponding QTL value to 0.5
insilicopopulation$individualsList[[3]]$QTLeffects$GCN1[["qtlTCrate"]] = 0.5 ## set the corresponding QTL value to 0.5
insilicopopulation$individualsList[[3]]$QTLeffects$GCN2[["qtlTCrate"]] = 0.5 ## set the corresponding QTL value to 0.5


restemp = simulateSystemStochasticParallel(insilicosystem, insilicopopulation, simtime = 1000, nepochs = 20, ntrialsPerInd = 50, simalgorithm = "SSA", returnStochModel = F)
#res2 = simulateSystemStochastic(insilicosystem, insilicopopulation, simtime = 1000, nepochs = 20, ntrialsPerInd = 50, simalgorithm = "SSA", returnStochModel = F)

res = lapply(restemp, mergeAllelesAbundance)

molsplot = plotexpprof(res)

ggsave("/media/sf_data/RNAplot2.png", plot = molsplot[[1]])
ggsave("/media/sf_data/Protplot2.png", plot = molsplot[[2]])

# ---------------
# plot histograms
# ---------------


molsList.hist = getlasttimepoint(res)

mols = setdiff(colnames(res[[1]]), c("time", "trial"))

molsplot.hist = lapply(mols, function(m){
  toplot = melt(molsList.hist[[m]])
  mybinwidth = min(sapply(levels(toplot$Var2), function(x){KernSmooth::dpih(toplot[toplot$Var2 == x, "value"])}))
  mymeans = data.frame("ind" = levels(toplot$Var2), "mean" = sapply(levels(toplot$Var2), function(x){mean(toplot[toplot$Var2 == x, "value"])}))
  myplot = ggplot() + geom_histogram(data = toplot, aes(x = value, fill = Var2), binwidth = mybinwidth, alpha = 0.8, position = "identity") + 
    geom_vline(data = mymeans, aes(xintercept = mean, color = ind), linetype="dashed", size=1) + 
    scale_color_discrete(guide = F) +
    ggtitle(m)
  print(myplot)
  return(myplot)
})

ggsave("/media/sf_data/RNAplot2_hist.png", plot = molsplot.hist[[1]])
ggsave("/media/sf_data/Protplot2_hist.png", plot = molsplot.hist[[2]])

save(insilicosystem, insilicopopulation, file = "/home/oangelin/Documents/plots2_05_07_2018.RData")


# ---------------------------------------------------------------------------------------- #
#' 3: One protein-coding gene, 5 variants: one orginal and 4 mutated version each affecting
#' a different QTL
#' Diploid individuals, 25 individuals (each possible combinations of the variants)
# ---------------------------------------------------------------------------------------- #

# source("network_generation.R")

mysystemargs = insilicosystemargs(G = 1, PC.p = 1)
insilicosystem = createInSilicoSystem(mysystemargs, empty = T)

myindivargs = insilicoindividualargs(ploidy = 2, ngenevariants = 5)

myvariants = list("1" = matrix(1.0, nrow = 9, ncol = 5, dimnames = list(c("qtlTCrate", "qtlRDrate", "qtlTCregbind", "qtlRDbindreg", "qtlactivity", "qtlTLrate", "qtlPDrate", "qtlTLregbind", "qtlPDregbind"), 1:5)))

## Manually change the value of the QTLs for each individual
qtlvals = c("qtlTCrate", "qtlTLrate", "qtlRDrate", "qtlPDrate")

for(i in 1:4){
  myvariants$`1`[qtlvals[i], i+1] = 0.5
}

## function to create an individual with the given allele combination
createmyindiv = function(allelecomb, variantsList, indargs){
  G = 1
  QTLeffects = vector("list", indargs$ploidy)
  names(QTLeffects) = indargs$gcnList
  individualvariants = as.data.frame(matrix(allelecomb, nrow = G, ncol = indargs$ploidy))
  names(individualvariants) = indargs$gcnList
  qtlnames = c("qtlTCrate", "qtlRDrate", "qtlTCregbind", "qtlRDbindreg", "qtlactivity", "qtlTLrate", "qtlPDrate", "qtlTLregbind", "qtlPDregbind")
  for(gcn in indargs$gcnList){
    QTLeffects[[gcn]] = vector("list", length(qtlnames))
    names(QTLeffects[[gcn]]) = qtlnames
    for(q in qtlnames){
      for(g in 1:G){
        QTLeffects[[gcn]][[q]][g] = variantsList[[g]][q, individualvariants[g, gcn]]
      }
    }
  }
  
  InitVar = vector("list", indargs$ploidy)
  names(InitVar) = indargs$gcnList
  for(gcn in indargs$gcnList){
    InitVar[[gcn]] = list("R" = rep(1.0, G),
                          "P" = rep(1.0, G))
  }
  value = list("QTLeffects" = QTLeffects, "haplotype" = individualvariants, "InitVar" = InitVar)
  attr(value, "class") = "insilicoindividual"
  
  return(value)
}

createmypop = function(myvariants, allallelecomb, indargs){
  nind = nrow(allallelecomb)
  genvariants = myvariants
  indnames = sapply(1:nind, function(x){paste0("Ind", x)})
  individualsList = vector("list", nind)
  names(individualsList) = indnames
  
  for(i in indnames){
    individualsList[[i]] = createmyindiv(allallelecomb[i,], genvariants, indargs)
  }
  
  value = list("GenesVariants" = genvariants, "individualsList" = individualsList, "indargs" = indargs)
  return(value)
}

allallelecomb = expand.grid("all1" = 1:5, "all2" = 1:5)
allallelecomb = t(apply(allallelecomb, 1, sort))
allallelecomb = allallelecomb[!duplicated(allallelecomb), ]
rownames(allallelecomb) = sapply(1:nrow(allallelecomb), function(x){paste0("Ind", x)})
insilicopopulation = createmypop(myvariants, allallelecomb, myindivargs)

restemp = simulateSystemStochasticParallel(insilicosystem, insilicopopulation, simtime = 5000, nepochs = 20, ntrialsPerInd = 50, simalgorithm = "SSA", returnStochModel = F)
res = lapply(restemp, mergeAllelesAbundance)

plotexpprof(res)

mols = setdiff(colnames(res[[1]]), c("time", "trial"))

variantseffects = c("original", "reduced\nTC rate", "reduced\nTL rate", "reduced\nRD rate", "reduced\nPD rate")

molsList = lapply(mols, function(m){
  molsummary = vector()
  for(ind in names(res)){
    resInd = res[[ind]]
    timexTrial = matrix(resInd[, m], ncol = max(resInd$trial), nrow = length(unique(resInd$time)))
    meanTrials = rowMeans(timexTrial)
    quantile1Trials = sapply(1:nrow(timexTrial), function(x){quantile(timexTrial[x, ], probs = 0.25)})
    quantile3Trials = sapply(1:nrow(timexTrial), function(x){quantile(timexTrial[x, ], probs = 0.975)})
    molsummary = rbind(molsummary, data.frame("time" = resInd$time[1:nrow(timexTrial)], "ind" = rep(ind, nrow(timexTrial)), "mean" = meanTrials, "quantile1" = quantile1Trials, "quantile3" = quantile3Trials, "GCN1" = rep(variantseffects[allallelecomb[ind, 1]], nrow(timexTrial)), "GCN2" = rep(variantseffects[allallelecomb[ind, 2]], nrow(timexTrial))))
    levels(molsummary$GCN1) = variantseffects = c("original", "reduced\nTC rate", "reduced\nTL rate", "reduced\nRD rate", "reduced\nPD rate")
    levels(molsummary$GCN2) = variantseffects = c("original", "reduced\nTC rate", "reduced\nTL rate", "reduced\nRD rate", "reduced\nPD rate")
  }
  return(molsummary)
})
names(molsList) = mols

mylinetype = c("mean" = "solid", "quantile1" = "longdash", "quantile3" = "longdash")

molsplot = lapply(mols, function(m){
  toplot = melt(molsList[[m]], id.vars = c("time", "ind", "GCN1", "GCN2"), measure.vars = c("mean", "quantile1", "quantile3"))
  myplot = ggplot() + geom_line(data = toplot, aes(x = time, y = value, color = ind, linetype = variable)) +
    geom_ribbon(data = molsList[[m]], aes(x = time, ymin=quantile1, ymax=quantile3, fill = ind), alpha=0.3) +
    scale_linetype_manual(values = mylinetype, guide = F) + facet_grid(GCN1 ~ GCN2) +
    ggtitle(m) + ylab("Molecule abundance")
  print(myplot)
  return(myplot)
})


# ---------------------------------------------------------------------------------------- #
#' 4: One protein-coding gene, 2 variants: one orginal and 1 mutated version 
#' with reduced transcription rate
#' Tetraploid individuals, 25 individuals (each possible combinations of the variants)
# ---------------------------------------------------------------------------------------- #


mysystemargs = insilicosystemargs(G = 1, PC.p = 1)
insilicosystem = createInSilicoSystem(mysystemargs, empty = T)

myindivargs = insilicoindividualargs(ploidy = 4, ngenevariants = 2)

myvariants = list("1" = matrix(1.0, nrow = 9, ncol = 2, dimnames = list(c("qtlTCrate", "qtlRDrate", "qtlTCregbind", "qtlRDbindreg", "qtlactivity", "qtlTLrate", "qtlPDrate", "qtlTLregbind", "qtlPDregbind"), 1:2)))

## Manually change the value of the QTLs for each individual
myvariants$`1`["qtlTCrate", 2] = 0.5

## function to create an individual with the given allele combination
createmyindiv = function(allelecomb, variantsList, indargs){
  G = 1
  QTLeffects = vector("list", indargs$ploidy)
  names(QTLeffects) = indargs$gcnList
  individualvariants = as.data.frame(matrix(allelecomb, nrow = G, ncol = indargs$ploidy))
  names(individualvariants) = indargs$gcnList
  qtlnames = c("qtlTCrate", "qtlRDrate", "qtlTCregbind", "qtlRDbindreg", "qtlactivity", "qtlTLrate", "qtlPDrate", "qtlTLregbind", "qtlPDregbind")
  for(gcn in indargs$gcnList){
    QTLeffects[[gcn]] = vector("list", length(qtlnames))
    names(QTLeffects[[gcn]]) = qtlnames
    for(q in qtlnames){
      for(g in 1:G){
        QTLeffects[[gcn]][[q]][g] = variantsList[[g]][q, individualvariants[g, gcn][[1]]]
      }
    }
  }
  
  InitVar = vector("list", indargs$ploidy)
  names(InitVar) = indargs$gcnList
  for(gcn in indargs$gcnList){
    InitVar[[gcn]] = list("R" = rep(1.0, G),
                          "P" = rep(1.0, G))
  }
  value = list("QTLeffects" = QTLeffects, "haplotype" = individualvariants, "InitVar" = InitVar)
  attr(value, "class") = "insilicoindividual"
  
  return(value)
}

createmypop = function(myvariants, allallelecomb, indargs){
  nind = nrow(allallelecomb)
  genvariants = myvariants
  indnames = sapply(1:nind, function(x){paste0("Ind", x)})
  individualsList = vector("list", nind)
  names(individualsList) = indnames
  
  for(i in indnames){
    individualsList[[i]] = createmyindiv(allallelecomb[i,], genvariants, indargs)
  }
  
  value = list("GenesVariants" = genvariants, "individualsList" = individualsList, "indargs" = indargs)
  return(value)
}

allallelecomb = as.data.frame(matrix(c(1,1,1,1,
                         1,1,1,2,
                         1,1,2,2,
                         1,2,2,2,
                         2,2,2,2), byrow = T, ncol = 4))
names(allallelecomb) = c("GCN1", "GCN2", "GCN3", "GCN4")
rownames(allallelecomb) = sapply(1:nrow(allallelecomb), function(x){paste0("Ind", x)})
insilicopopulation = createmypop(myvariants, allallelecomb, myindivargs)

restemp = simulateSystemStochasticParallel(insilicosystem, insilicopopulation, simtime = 5000, nepochs = 20, ntrialsPerInd = 100, simalgorithm = "SSA", returnStochModel = F)
res = lapply(restemp, mergeAllelesAbundance)

plotexpprof(res)

molsList.hist = getlasttimepoint(res)

mols = setdiff(colnames(res[[1]]), c("time", "trial"))

molsplot.hist = lapply(mols, function(m){
  toplot = melt(molsList.hist[[m]])
  mybinwidth = min(sapply(levels(toplot$Var2), function(x){KernSmooth::dpih(toplot[toplot$Var2 == x, "value"])}))
  mymeans = data.frame("ind" = levels(toplot$Var2), "mean" = sapply(levels(toplot$Var2), function(x){mean(toplot[toplot$Var2 == x, "value"])}))
  myplot = ggplot() + geom_histogram(data = toplot, aes(x = value, fill = Var2), binwidth = mybinwidth, alpha = 0.8, position = "identity") + 
    geom_vline(data = mymeans, aes(xintercept = mean, color = ind), linetype="dashed", size=1) + 
    scale_color_discrete(guide = F) +
  ggtitle(m)
  print(myplot)
  return(myplot)
})

