setwd("~/winData/multiomics_networks_simulation")
rm(list = ls(all = T))
#setwd("~/GitHub/multiomics_networks_simulation")

source("network_generation.R")
library(ggridges)

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
      quantile1Trials = sapply(1:nrow(timexTrial), function(x){quantile(timexTrial[x, ], probs = 0.025)})
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


getlasttimepoint = function(res){6

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



#### -------------------------------------------------------------------------------------------- ####
##                            Illustrating the regulatory interactions:
## 2 genes, one regulates the other, we only show how the regulation affects the expression profiles
#### -------------------------------------------------------------------------------------------- ####

#### 1 --> 2  - TC, positive regulation ####

mysystemargs = insilicosystemargs(G = 2, PC.p = 1)
insilicosystem = createInSilicoSystem(mysystemargs, empty = T)

## change the kinetics of gene 1 (regulatory gene)
insilicosystem$genes[1, "TCrate"] = 0.01
insilicosystem$genes[1, "RDrate"] = 0.001
insilicosystem$genes[1, "TLrate"] = 0.0005
insilicosystem$genes[1, "PDrate"] = 0.0001

## change the kinetics of gene 2 (target gene)
insilicosystem$genes[2, "TCrate"] = 0.05
insilicosystem$genes[2, "RDrate"] = 0.001
insilicosystem$genes[2, "TLrate"] = 0.001
insilicosystem$genes[2, "PDrate"] = 0.0002

insilicosystem = addEdg(insilicosystem, regulator = 1, target = 2, targetreaction = "TC", regsign = "1")
insilicosystem$mosystem$TCRN.edg[1, "TCbindingrate"] = 0.1
insilicosystem$mosystem$TCRN.edg[1, "TCunbindingrate"] = 5
insilicosystem$mosystem$TCRN.edg[1, "TCfoldchange"] = 10

myindivargs = insilicoindividualargs(ploidy = 1, ngenevariants = 1)
insilicopopulation = createPopulation(1, insilicosystem, myindivargs, sameInit = T)
res1pos = simulateSystemStochastic(insilicosystem, insilicopopulation, simtime = 3600, nepochs = 3600, ntrialsPerInd = 1000, simalgorithm = "SSA", returnStochModel = F)
#load("/media/sf_data/confirmation_results/conf_rep_res1.RData")

mols = setdiff(colnames(res1pos$resTable[[1]]), c("time", "trial"))
resInd = res1pos$resTable[["Ind1"]]

molsList = lapply(mols, function(m){
  molsummary = vector()
  timexTrial = matrix(resInd[, m], ncol = max(resInd$trial), nrow = length(unique(resInd$time)))
  meanTrials = rowMeans(timexTrial)
  quantile1Trials = sapply(1:nrow(timexTrial), function(x){quantile(timexTrial[x, ], probs = 0.025)})
  quantile3Trials = sapply(1:nrow(timexTrial), function(x){quantile(timexTrial[x, ], probs = 0.975)})
  molsummary = rbind(molsummary, data.frame("time" = resInd$time[1:nrow(timexTrial)], "moltype" = rep(substr(m, 1, 1), nrow(timexTrial)),  "molnumber" = rep(substr(m, 2, 2), nrow(timexTrial)), "mean" = meanTrials, "quantile1" = quantile1Trials, "quantile3" = quantile3Trials))
  return(molsummary)
})

names(molsList) = mols
molsall = do.call("rbind", molsList)
mylinetype = c("mean" = "solid", "quantile1" = "longdash", "quantile3" = "longdash")
molcol = c("1" = "red", "2" = "green3")

toplot = melt(molsall, id.vars = c("time", "moltype", "molnumber"), measure.vars = c("mean", "quantile1", "quantile3"))

plot1pos = ggplot() + geom_line(data = toplot, aes(x = time, y = value, color = molnumber, linetype = variable)) +
  geom_ribbon(data = molsall, aes(x = time, ymin=quantile1, ymax=quantile3, fill = molnumber), alpha=0.3) +
  scale_color_manual(values = molcol, name = "Gene ID", labels = c("1" = "Gene 1", "2" = "Gene 2")) +
  scale_fill_manual(values = molcol, name = "Gene ID", labels = c("1" = "Gene 1", "2" = "Gene 2")) +
  scale_linetype_manual(values = mylinetype, guide = F) +
  facet_wrap(~moltype, ncol = 1, scales = "free_y", labeller = as_labeller(c("R" = "Transcript abundance profiles", "P" = "Protein abundance profiles"))) +
  ggtitle("") + ylab("Absolute abundance (# molecules)") + xlab("Time (s)") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.x = element_text(size=12), strip.background = element_rect(colour = "white", fill = "white"))

print(plot1pos)

ggsave("/media/sf_data/confirmation_results/plot1_TCpos.pdf", plot = plot1pos, device = "pdf")


#### 1 --| 2 - TC, negative regulation ####

mysystemargs = insilicosystemargs(G = 2, PC.p = 1)
insilicosystem = createInSilicoSystem(mysystemargs, empty = T)


## change the kinetics of gene 1 (regulatory gene)
insilicosystem$genes[1, "TCrate"] = 0.01
insilicosystem$genes[1, "RDrate"] = 0.001
insilicosystem$genes[1, "TLrate"] = 0.0005
insilicosystem$genes[1, "PDrate"] = 0.0001

## change the kinetics of gene 2 (target gene)
insilicosystem$genes[2, "TCrate"] = 0.05
insilicosystem$genes[2, "RDrate"] = 0.001
insilicosystem$genes[2, "TLrate"] = 0.001
insilicosystem$genes[2, "PDrate"] = 0.0002

insilicosystem = addEdg(insilicosystem, regulator = 1, target = 2, targetreaction = "TC", regsign = "-1")
insilicosystem$mosystem$TCRN.edg[1, "TCbindingrate"] = 0.1
insilicosystem$mosystem$TCRN.edg[1, "TCunbindingrate"] = 5
insilicosystem$mosystem$TCRN.edg[1, "TCfoldchange"] = 0

myindivargs = insilicoindividualargs(ploidy = 1, ngenevariants = 1)
insilicopopulation = createPopulation(1, insilicosystem, myindivargs, sameInit = T)
res1neg = simulateSystemStochastic(insilicosystem, insilicopopulation, simtime = 3600, nepochs = 3600, ntrialsPerInd = 1000, simalgorithm = "SSA", returnStochModel = F)


mols = setdiff(colnames(res1neg$resTable[[1]]), c("time", "trial"))

resInd = res1neg$resTable[["Ind1"]]

molsList = lapply(mols, function(m){
  molsummary = vector()
  timexTrial = matrix(resInd[, m], ncol = max(resInd$trial), nrow = length(unique(resInd$time)))
  meanTrials = rowMeans(timexTrial)
  quantile1Trials = sapply(1:nrow(timexTrial), function(x){quantile(timexTrial[x, ], probs = 0.025)})
  quantile3Trials = sapply(1:nrow(timexTrial), function(x){quantile(timexTrial[x, ], probs = 0.975)})
  molsummary = rbind(molsummary, data.frame("time" = resInd$time[1:nrow(timexTrial)], "moltype" = rep(substr(m, 1, 1), nrow(timexTrial)),  "molnumber" = rep(substr(m, 2, 2), nrow(timexTrial)), "mean" = meanTrials, "quantile1" = quantile1Trials, "quantile3" = quantile3Trials))
  return(molsummary)
})

names(molsList) = mols
molsall = do.call("rbind", molsList)
mylinetype = c("mean" = "solid", "quantile1" = "longdash", "quantile3" = "longdash")
molcol = c("1" = "red", "2" = "green3")

toplot = melt(molsall, id.vars = c("time", "moltype", "molnumber"), measure.vars = c("mean", "quantile1", "quantile3"))

plot1neg = ggplot() + geom_line(data = toplot, aes(x = time, y = value, color = molnumber, linetype = variable)) +
  geom_ribbon(data = molsall, aes(x = time, ymin=quantile1, ymax=quantile3, fill = molnumber), alpha=0.3) +
  scale_color_manual(values = molcol, name = "Gene ID", labels = c("1" = "Gene 1", "2" = "Gene 2")) +
  scale_fill_manual(values = molcol, name = "Gene ID", labels = c("1" = "Gene 1", "2" = "Gene 2")) +
  scale_linetype_manual(values = mylinetype, guide = F) +
  facet_wrap(~moltype, ncol = 1, scales = "free_y", labeller = as_labeller(c("R" = "Transcript abundance profiles", "P" = "Protein abundance profiles"))) +
  ggtitle("") + ylab("Absolute abundance (# molecules)") + xlab("Time (s)") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.x = element_text(size=12), strip.background = element_rect(colour = "white", fill = "white"))

print(plot1neg)

ggsave("/media/sf_data/confirmation_results/plot1_TCneg.pdf", plot = plot1neg, device = "pdf")

save(res1pos, res1neg, file = "/media/sf_data/confirmation_results/conf_rep_res1.RData")

# plot1tot = ggarrange(plot1pos, plot1neg, nrow = 1, ncol = 2)
# ggsave("/media/sf_data/confirmation_results/plot1.pdf", plot = plot1tot, device = "pdf")


#### 1 --> 2  - TL, positive regulation ####

mysystemargs = insilicosystemargs(G = 2, PC.p = 1)
insilicosystem = createInSilicoSystem(mysystemargs, empty = T)

## change the kinetics of gene 1 (regulatory gene)
insilicosystem$genes[1, "TCrate"] = 0.01
insilicosystem$genes[1, "RDrate"] = 0.001
insilicosystem$genes[1, "TLrate"] = 0.0005
insilicosystem$genes[1, "PDrate"] = 0.0001
#insilicosystem$genes[1, "TargetReaction"] = "TL"

## change the kinetics of gene 2 (target gene)
insilicosystem$genes[2, "TCrate"] = 0.05
insilicosystem$genes[2, "RDrate"] = 0.001
insilicosystem$genes[2, "TLrate"] = 0.001
insilicosystem$genes[2, "PDrate"] = 0.0002

insilicosystem = addEdg(insilicosystem, regulator = 1, target = 2, targetreaction = "TL", regsign = "1")
insilicosystem$mosystem$TLRN.edg[1, "TLbindingrate"] = 0.1
insilicosystem$mosystem$TLRN.edg[1, "TLunbindingrate"] = 5
insilicosystem$mosystem$TLRN.edg[1, "TLfoldchange"] = 10

myindivargs = insilicoindividualargs(ploidy = 1, ngenevariants = 1)
insilicopopulation = createPopulation(1, insilicosystem, myindivargs, sameInit = T)
res1TLpos = simulateSystemStochastic(insilicosystem, insilicopopulation, simtime = 3600, nepochs = 3600, ntrialsPerInd = 1000, simalgorithm = "SSA", returnStochModel = F)
#load("/media/sf_data/confirmation_results/conf_rep_res1TL.RData")

mols = setdiff(colnames(res1TLpos$resTable[[1]]), c("time", "trial"))
resInd = res1TLpos$resTable[["Ind1"]]

molsList = lapply(mols, function(m){
  molsummary = vector()
  timexTrial = matrix(resInd[, m], ncol = max(resInd$trial), nrow = length(unique(resInd$time)))
  meanTrials = rowMeans(timexTrial)
  quantile1Trials = sapply(1:nrow(timexTrial), function(x){quantile(timexTrial[x, ], probs = 0.025)})
  quantile3Trials = sapply(1:nrow(timexTrial), function(x){quantile(timexTrial[x, ], probs = 0.975)})
  molsummary = rbind(molsummary, data.frame("time" = resInd$time[1:nrow(timexTrial)], "moltype" = rep(substr(m, 1, 1), nrow(timexTrial)),  "molnumber" = rep(substr(m, 2, 2), nrow(timexTrial)), "mean" = meanTrials, "quantile1" = quantile1Trials, "quantile3" = quantile3Trials))
  return(molsummary)
})

names(molsList) = mols
molsall = do.call("rbind", molsList)
mylinetype = c("mean" = "solid", "quantile1" = "longdash", "quantile3" = "longdash")
molcol = c("1" = "red", "2" = "green3")

toplot = melt(molsall, id.vars = c("time", "moltype", "molnumber"), measure.vars = c("mean", "quantile1", "quantile3"))

plot1TLpos = ggplot() + geom_line(data = toplot, aes(x = time, y = value, color = molnumber, linetype = variable)) +
  geom_ribbon(data = molsall, aes(x = time, ymin=quantile1, ymax=quantile3, fill = molnumber), alpha=0.3) +
  scale_color_manual(values = molcol, name = "Gene ID", labels = c("1" = "Gene 1", "2" = "Gene 2")) +
  scale_fill_manual(values = molcol, name = "Gene ID", labels = c("1" = "Gene 1", "2" = "Gene 2")) +
  scale_linetype_manual(values = mylinetype, guide = F) +
  facet_wrap(~moltype, ncol = 1, scales = "free_y", labeller = as_labeller(c("R" = "Transcript abundance profiles", "P" = "Protein abundance profiles"))) +
  ggtitle("") + ylab("Absolute abundance (# molecules)") + xlab("Time (s)") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.x = element_text(size=12), strip.background = element_rect(colour = "white", fill = "white"))

print(plot1TLpos)

ggsave("/media/sf_data/confirmation_results/plot1_TLpos.pdf", plot = plot1TLpos, device = "pdf")


#### 1 --> 2  - TL, negative regulation ####

mysystemargs = insilicosystemargs(G = 2, PC.p = 1)
insilicosystem = createInSilicoSystem(mysystemargs, empty = T)

## change the kinetics of gene 1 (regulatory gene)
insilicosystem$genes[1, "TCrate"] = 0.01
insilicosystem$genes[1, "RDrate"] = 0.001
insilicosystem$genes[1, "TLrate"] = 0.0005
insilicosystem$genes[1, "PDrate"] = 0.0001
#insilicosystem$genes[1, "TargetReaction"] = "TL"

## change the kinetics of gene 2 (target gene)
insilicosystem$genes[2, "TCrate"] = 0.05
insilicosystem$genes[2, "RDrate"] = 0.001
insilicosystem$genes[2, "TLrate"] = 0.001
insilicosystem$genes[2, "PDrate"] = 0.0002

insilicosystem = addEdg(insilicosystem, regulator = 1, target = 2, targetreaction = "TL", regsign = "1")
insilicosystem$mosystem$TLRN.edg[1, "TLbindingrate"] = 0.1
insilicosystem$mosystem$TLRN.edg[1, "TLunbindingrate"] = 5
insilicosystem$mosystem$TLRN.edg[1, "TLfoldchange"] = 0

myindivargs = insilicoindividualargs(ploidy = 1, ngenevariants = 1)
insilicopopulation = createPopulation(1, insilicosystem, myindivargs, sameInit = T)
res1TLneg = simulateSystemStochastic(insilicosystem, insilicopopulation, simtime = 3600, nepochs = 3600, ntrialsPerInd = 1000, simalgorithm = "SSA", returnStochModel = F)
#load("/media/sf_data/confirmation_results/conf_rep_res1.RData")

mols = setdiff(colnames(res1TLneg$resTable[[1]]), c("time", "trial"))
resInd = res1TLneg$resTable[["Ind1"]]

molsList = lapply(mols, function(m){
  molsummary = vector()
  timexTrial = matrix(resInd[, m], ncol = max(resInd$trial), nrow = length(unique(resInd$time)))
  meanTrials = rowMeans(timexTrial)
  quantile1Trials = sapply(1:nrow(timexTrial), function(x){quantile(timexTrial[x, ], probs = 0.025)})
  quantile3Trials = sapply(1:nrow(timexTrial), function(x){quantile(timexTrial[x, ], probs = 0.975)})
  molsummary = rbind(molsummary, data.frame("time" = resInd$time[1:nrow(timexTrial)], "moltype" = rep(substr(m, 1, 1), nrow(timexTrial)),  "molnumber" = rep(substr(m, 2, 2), nrow(timexTrial)), "mean" = meanTrials, "quantile1" = quantile1Trials, "quantile3" = quantile3Trials))
  return(molsummary)
})

names(molsList) = mols
molsall = do.call("rbind", molsList)
mylinetype = c("mean" = "solid", "quantile1" = "longdash", "quantile3" = "longdash")
molcol = c("1" = "red", "2" = "green3")

toplot = melt(molsall, id.vars = c("time", "moltype", "molnumber"), measure.vars = c("mean", "quantile1", "quantile3"))

plot1TLneg = ggplot() + geom_line(data = toplot, aes(x = time, y = value, color = molnumber, linetype = variable)) +
  geom_ribbon(data = molsall, aes(x = time, ymin=quantile1, ymax=quantile3, fill = molnumber), alpha=0.3) +
  scale_color_manual(values = molcol, name = "Gene ID", labels = c("1" = "Gene 1", "2" = "Gene 2")) +
  scale_fill_manual(values = molcol, name = "Gene ID", labels = c("1" = "Gene 1", "2" = "Gene 2")) +
  scale_linetype_manual(values = mylinetype, guide = F) +
  facet_wrap(~moltype, ncol = 1, scales = "free_y", labeller = as_labeller(c("R" = "Transcript abundance profiles", "P" = "Protein abundance profiles"))) +
  ggtitle("") + ylab("Absolute abundance (# molecules)") + xlab("Time (s)") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.x = element_text(size=12), strip.background = element_rect(colour = "white", fill = "white"))

print(plot1TLneg)

ggsave("/media/sf_data/confirmation_results/plot1_TLneg.pdf", plot = plot1TLneg, device = "pdf")

save(res1TLpos, res1TLneg, file = "/media/sf_data/confirmation_results/conf_rep_res1TL.RData")

#### 1 --> 2  -  RNA decay ####

mysystemargs = insilicosystemargs(G = 2, PC.p = 1)
insilicosystem = createInSilicoSystem(mysystemargs, empty = T)

## change the kinetics of gene 1 (regulatory gene)0
insilicosystem$genes[1, "TCrate"] = 0.001
insilicosystem$genes[1, "RDrate"] = 0.0001
insilicosystem$genes[1, "TLrate"] = 0.0005
insilicosystem$genes[1, "PDrate"] = 0.0001

## change the kinetics of gene 2 (target gene)
insilicosystem$genes[2, "TCrate"] = 0.005
insilicosystem$genes[2, "RDrate"] = 0.0001
insilicosystem$genes[2, "TLrate"] = 0.001
insilicosystem$genes[2, "PDrate"] = 0.0002

insilicosystem = addEdg(insilicosystem, regulator = 1, target = 2, targetreaction = "RD", regsign = "1")
insilicosystem$mosystem$RDRN.edg[1, "RDbindingrate"] = 0.00001

myindivargs = insilicoindividualargs(ploidy = 1, ngenevariants = 1)
insilicopopulation = createPopulation(1, insilicosystem, myindivargs, sameInit = T)
res1rna = simulateSystemStochastic(insilicosystem, insilicopopulation, simtime = 3600, nepochs = 3600, ntrialsPerInd = 1000, simalgorithm = "SSA", returnStochModel = F)
# load("/media/sf_data/confirmation_results/conf_rep_res1_decay.RData")

mols = setdiff(colnames(res1rna$resTable[[1]]), c("time", "trial"))
resInd = res1rna$resTable[["Ind1"]]

molsList = lapply(mols, function(m){
  molsummary = vector()
  timexTrial = matrix(resInd[, m], ncol = max(resInd$trial), nrow = length(unique(resInd$time)))
  meanTrials = rowMeans(timexTrial)
  quantile1Trials = sapply(1:nrow(timexTrial), function(x){quantile(timexTrial[x, ], probs = 0.025)})
  quantile3Trials = sapply(1:nrow(timexTrial), function(x){quantile(timexTrial[x, ], probs = 0.975)})
  molsummary = rbind(molsummary, data.frame("time" = resInd$time[1:nrow(timexTrial)], "moltype" = rep(substr(m, 1, 1), nrow(timexTrial)),  "molnumber" = rep(substr(m, 2, 2), nrow(timexTrial)), "mean" = meanTrials, "quantile1" = quantile1Trials, "quantile3" = quantile3Trials))
  return(molsummary)
})

names(molsList) = mols
molsall = do.call("rbind", molsList)
mylinetype = c("mean" = "solid", "quantile1" = "longdash", "quantile3" = "longdash")
molcol = c("1" = "red", "2" = "green3")

toplot = melt(molsall, id.vars = c("time", "moltype", "molnumber"), measure.vars = c("mean", "quantile1", "quantile3"))

plot1rna = ggplot() + geom_line(data = toplot, aes(x = time, y = value, color = molnumber, linetype = variable)) +
  geom_ribbon(data = molsall, aes(x = time, ymin=quantile1, ymax=quantile3, fill = molnumber), alpha=0.3) +
  scale_color_manual(values = molcol, name = "Gene ID", labels = c("1" = "Gene 1", "2" = "Gene 2")) +
  scale_fill_manual(values = molcol, name = "Gene ID", labels = c("1" = "Gene 1", "2" = "Gene 2")) +
  scale_linetype_manual(values = mylinetype, guide = F) +
  facet_wrap(~moltype, ncol = 1, scales = "free_y", labeller = as_labeller(c("R" = "Transcript abundance profiles", "P" = "Protein abundance profiles"))) +
  ggtitle("") + ylab("Absolute abundance (# molecules)") + xlab("Time (s)") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.x = element_text(size=12), strip.background = element_rect(colour = "white", fill = "white"))

print(plot1rna)

ggsave("/media/sf_data/confirmation_results/plot1_RD.pdf", plot = plot1rna, device = "pdf")


#### 1 --> 2   -  Protein decay ####

mysystemargs = insilicosystemargs(G = 2, PC.p = 1)
insilicosystem = createInSilicoSystem(mysystemargs, empty = T)

## change the kinetics of gene 1 (regulatory gene)0
insilicosystem$genes[1, "TCrate"] = 0.001
insilicosystem$genes[1, "RDrate"] = 0.0001
insilicosystem$genes[1, "TLrate"] = 0.0005
insilicosystem$genes[1, "PDrate"] = 0.0001

## change the kinetics of gene 2 (target gene)
insilicosystem$genes[2, "TCrate"] = 0.005
insilicosystem$genes[2, "RDrate"] = 0.0001
insilicosystem$genes[2, "TLrate"] = 0.001
insilicosystem$genes[2, "PDrate"] = 0.0002

insilicosystem = addEdg(insilicosystem, regulator = 1, target = 2, targetreaction = "PD", regsign = "1")
insilicosystem$mosystem$PDRN.edg[1, "PDbindingrate"] = 0.0001

myindivargs = insilicoindividualargs(ploidy = 1, ngenevariants = 1)
insilicopopulation = createPopulation(1, insilicosystem, myindivargs, sameInit = T)
res1prot = simulateSystemStochastic(insilicosystem, insilicopopulation, simtime = 3600, nepochs = 3600, ntrialsPerInd = 1000, simalgorithm = "SSA", returnStochModel = F)

mols = setdiff(colnames(res1prot$resTable[[1]]), c("time", "trial"))
#rnas = grep("^R", mols, value = TRUE)
#prots = grep("^P", mols, value = TRUE)

resInd = res1prot$resTable[["Ind1"]]

molsList = lapply(mols, function(m){
  molsummary = vector()
  timexTrial = matrix(resInd[, m], ncol = max(resInd$trial), nrow = length(unique(resInd$time)))
  meanTrials = rowMeans(timexTrial)
  quantile1Trials = sapply(1:nrow(timexTrial), function(x){quantile(timexTrial[x, ], probs = 0.025)})
  quantile3Trials = sapply(1:nrow(timexTrial), function(x){quantile(timexTrial[x, ], probs = 0.975)})
  molsummary = rbind(molsummary, data.frame("time" = resInd$time[1:nrow(timexTrial)], "moltype" = rep(substr(m, 1, 1), nrow(timexTrial)),  "molnumber" = rep(substr(m, 2, 2), nrow(timexTrial)), "mean" = meanTrials, "quantile1" = quantile1Trials, "quantile3" = quantile3Trials))
  return(molsummary)
})

names(molsList) = mols
molsall = do.call("rbind", molsList)
mylinetype = c("mean" = "solid", "quantile1" = "longdash", "quantile3" = "longdash")
molcol = c("1" = "red", "2" = "green3")

toplot = melt(molsall, id.vars = c("time", "moltype", "molnumber"), measure.vars = c("mean", "quantile1", "quantile3"))

plot1prot = ggplot() + geom_line(data = toplot, aes(x = time, y = value, color = molnumber, linetype = variable)) +
  geom_ribbon(data = molsall, aes(x = time, ymin=quantile1, ymax=quantile3, fill = molnumber), alpha=0.3) +
  scale_color_manual(values = molcol, name = "Gene ID", labels = c("1" = "Gene 1", "2" = "Gene 2")) +
  scale_fill_manual(values = molcol, name = "Gene ID", labels = c("1" = "Gene 1", "2" = "Gene 2")) +
  scale_linetype_manual(values = mylinetype, guide = F) +
  facet_wrap(~moltype, ncol = 1, scales = "free_y", labeller = as_labeller(c("R" = "Transcript abundance profiles", "P" = "Protein abundance profiles"))) +
  ggtitle("") + ylab("Absolute abundance (# molecules)") + xlab("Time (s)") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.x = element_text(size=12), strip.background = element_rect(colour = "white", fill = "white"))

print(plot1prot)

ggsave("/media/sf_data/confirmation_results/plot1_PD.pdf", plot = plot1prot, device = "pdf")

save(res1rna, res1prot, file = "/media/sf_data/confirmation_results/conf_rep_res1_decay.RData")

# plot1tot_decay = ggarrange(plot1rna, plot1prot, nrow = 1, ncol = 2)
# ggsave("/media/sf_data/confirmation_results/plot1_decay.pdf", plot = plot1tot_decay, device = "pdf")

#### -------------------------------------------------------------------------------------------- ####
##                               Illustrating the genetic variants:
## System 1 of simulation_tests: In the population each individual has only 1 QTL value !=1 to
##  show the impact of the different QTLs on the expression profiles
#### -------------------------------------------------------------------------------------------- ####
#----

mysystemargs = insilicosystemargs(G = 1, PC.p = 1)
insilicosystem = createInSilicoSystem(mysystemargs, empty = T)

## change the kinetics of gene 1
insilicosystem$genes[1, "TCrate"] = 0.1
insilicosystem$genes[1, "RDrate"] = 0.001
insilicosystem$genes[1, "TLrate"] = 0.01
insilicosystem$genes[1, "PDrate"] = 0.001

myindivargs = insilicoindividualargs(ploidy = 1, ngenevariants = 1)
insilicopopulation = createPopulation(5, insilicosystem, myindivargs, sameInit = T)

## Manually change the value of the QTLs for each individual
qtlvals = c("qtlTCrate", "qtlTLrate", "qtlRDrate", "qtlPDrate")

for(i in 2:5){
  insilicopopulation$individualsList[[i]]$QTLeffects$GCN1[[qtlvals[i-1]]] = 0.5 ## set the corresponding QTL value to 0.5
}

resmut = simulateSystemStochasticParallel(insilicosystem, insilicopopulation, simtime = 3600, nepochs = 2, ntrialsPerInd = 5000, simalgorithm = "SSA", returnStochModel = F)

mols = setdiff(colnames(resmut[[1]]), c("time", "trial"))
molsList.hist = getlasttimepoint(resmut)

molsList = lapply(mols, function(m){
  melted = melt(molsList.hist[[m]])
  return(data.frame("mol" = rep(substr(m, 1, 1), nrow(melted)), "Ind" = melted$Var2, "Value" = melted$value))
})
names(molsList) = mols
molsall = do.call("rbind", molsList)

plotmut = ggplot(molsall, aes(x = Value, y = Ind, fill = mol, colour = mol)) + geom_density_ridges(alpha = 0.8, rel_min_height = 0.01, scale = 1.5) +
  facet_wrap(~mol, scale = "free_x", labeller = as_labeller(c("R" = "Transcript abundance density\nat t = 1 hour (5000 simulations)", "P" = "Protein abundance density\nat t = 1 hour (5000 simulations)"))) +
  scale_colour_manual(values = c("R" = "darkgoldenrod1", "P" = "dodgerblue"), guide = F) +
  scale_fill_manual(values = c("R" = "darkgoldenrod1", "P" = "dodgerblue"), guide = F) +
  scale_y_discrete(limits = rev(levels(molsall$Ind)), labels = c("Ind1" = "Ind1\n(Original\nallele)", "Ind2" = "Ind2\n(Reduced\ntranscription rate)", "Ind3" = "Ind3\n(Reduced\ntranslation rate)",
                              "Ind4" = "Ind4\n(Reduced\nRNA decay rate)", "Ind5" = "Ind5\n(Reduced\nprotein decay rate)")) +
  xlab("Molecule abundance (# molecules)") + ylab("Simulated individuals") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.border = element_blank(), panel.spacing = unit(1, "lines"),
                                                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey45"),
                                                       strip.text.x = element_text(size=12), strip.background = element_rect(colour = "white", fill = "white"))
print(plotmut)

ggsave("/media/sf_data/confirmation_results/plot2mut.pdf", plot = plotmut, device = "pdf")

save(resmut, file = "/media/sf_data/confirmation_results/conf_rep_res2mut.RData")



#### -------------------------------------------------------------------------------------------- ####
##                               Illustrating the ploidy level:
## System 4 and 4' of simulation_tests: One protein-coding gene, 2 variants: one original and 1
## mutated version with reduced transcription rate. Tetraploid individuals, 5 individuals
## (each possible allele dosage)
#### -------------------------------------------------------------------------------------------- ####
#----

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

## Part 1: QTL effect of 0.5 ----

mysystemargs = insilicosystemargs(G = 1, PC.p = 1)
insilicosystem = createInSilicoSystem(mysystemargs, empty = T)

## change the kinetics of gene 1
insilicosystem$genes[1, "TCrate"] = 0.1
insilicosystem$genes[1, "RDrate"] = 0.001
insilicosystem$genes[1, "TLrate"] = 0.01
insilicosystem$genes[1, "PDrate"] = 0.001

myindivargs = insilicoindividualargs(ploidy = 4, ngenevariants = 2)

myvariants = list("1" = matrix(1.0, nrow = 9, ncol = 2, dimnames = list(c("qtlTCrate", "qtlRDrate", "qtlTCregbind", "qtlRDbindreg", "qtlactivity", "qtlTLrate", "qtlPDrate", "qtlTLregbind", "qtlPDregbind"), 1:2)))

## Manually change the value of the QTLs for each individual
myvariants$`1`["qtlTCrate", 2] = 0.5

allallelecomb = as.data.frame(matrix(c(1,1,1,1,
                                       1,1,1,2,
                                       1,1,2,2,
                                       1,2,2,2,
                                       2,2,2,2), byrow = T, ncol = 4))
names(allallelecomb) = c("GCN1", "GCN2", "GCN3", "GCN4")
rownames(allallelecomb) = sapply(1:nrow(allallelecomb), function(x){paste0("Ind", x)})
insilicopopulation = createmypop(myvariants, allallelecomb, myindivargs)

resploidy = simulateSystemStochasticParallel(insilicosystem, insilicopopulation, simtime = 3600, nepochs = 2, ntrialsPerInd = 5000, simalgorithm = "SSA", returnStochModel = F)
resploidymerged = lapply(resploidy, mergeAllelesAbundance)

molsList.hist = getlasttimepoint(resploidymerged)
mols = setdiff(colnames(resploidymerged[[1]]), c("time", "trial"))

temp = lapply(mols, function(m){
  melted = melt(molsList.hist[[m]])
  return(data.frame("mol" = rep(m, nrow(melted)), "ind" = melted$Var2, "value" = melted$value))
})

molsall = do.call(rbind, temp)
mycols = c("#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8","#253494")
names(mycols) = c("Ind5", "Ind4", "Ind3", "Ind2", "Ind1")
  # c("aaaa", "Aaaa", "AAaa", "AAAa", "AAAA")

temp = lapply(mols, function(m){
  mymeans = colMeans(molsList.hist[[m]])
  return(data.frame("mol" = rep(m, length(mymeans)), "ind" = names(mymeans), "mean" = mymeans))
})
meanstoplot = do.call(rbind, temp)

plotploidy = ggplot() + geom_histogram(data = molsall, aes(x = value, color = ind, fill = ind), alpha = 0.8, position = "identity") +
  geom_vline(data = meanstoplot, aes(xintercept = mean, color = ind), linetype="dashed", size=1) +
  scale_color_manual(name = "Simulated individuals", values = mycols, labels = c("Ind1" = "Ind1\n(AAAA)", "Ind2" = "Ind2\n(AAAa)", "Ind3" = "Ind3\n(AAaa)", "Ind4" = "Ind4\n(Aaaa)", "Ind5" = "Ind5\n(aaaa)"), guide = F) +
  scale_fill_manual(name = "Simulated individuals",values = mycols, labels = c("Ind1" = "Ind1\n(AAAA)", "Ind2" = "Ind2\n(AAAa)", "Ind3" = "Ind3\n(AAaa)", "Ind4" = "Ind4\n(Aaaa)", "Ind5" = "Ind5\n(aaaa)")) +
  facet_wrap(~mol, scale = "free", labeller = as_labeller(c("R1" = "Transcript abundance density\nat t = 1 hour (5000 simulations)", "P1" = "Protein abundance density\nat t = 1 hour (5000 simulations)"))) +
  theme_bw() + theme(legend.position = "top", legend.direction = "horizontal", plot.title = element_text(hjust = 0.5), panel.border = element_blank(), panel.spacing = unit(1, "lines"),
                                                                                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey45"),
                                                                                         strip.text.x = element_text(size=12), strip.background = element_rect(colour = "white", fill = "white")) +

  xlab("Molecule abundance (# molecules)")

print(plotploidy)

ggsave("/media/sf_data/confirmation_results/plot3ploidy.pdf", plot = plotploidy, device = "pdf")

save(resploidy, file = "/media/sf_data/confirmation_results/conf_rep_res3ploidy.RData")


## Part2: vary the value of QTL effect ----

mysystemargs = insilicosystemargs(G = 1, PC.p = 1)
insilicosystem = createInSilicoSystem(mysystemargs, empty = T)

## change the kinetics of gene 1
insilicosystem$genes[1, "TCrate"] = 0.1
insilicosystem$genes[1, "RDrate"] = 0.001
insilicosystem$genes[1, "TLrate"] = 0.01
insilicosystem$genes[1, "PDrate"] = 0.001

qtleffectrange = seq(from = 1, to = 0.1, by = -0.1)
myvariants = list("1" = matrix(1.0, nrow = 9, ncol = length(qtleffectrange) + 1, dimnames = list(c("qtlTCrate", "qtlRDrate", "qtlTCregbind", "qtlRDbindreg", "qtlactivity", "qtlTLrate", "qtlPDrate", "qtlTLregbind", "qtlPDregbind"), 1:(length(qtleffectrange) + 1))))

## Manually change the value of the QTLs for each individual

for(i in 1:length(qtleffectrange)){
  myvariants$`1`["qtlTCrate", i+1] = qtleffectrange[i]
}

myindivargs = insilicoindividualargs(ploidy = 4, ngenevariants = length(qtleffectrange) + 1)


allallelecomb = vector()

for(i in 2:(length(qtleffectrange)+1)){
  allallelecomb = rbind(allallelecomb, as.data.frame(matrix(c(1,1,1,1,
                                                              1,1,1,i,
                                                              1,1,i,i,
                                                              1,i,i,i,
                                                              i,i,i,i), byrow = T, ncol = 4)))
}

names(allallelecomb) = c("GCN1", "GCN2", "GCN3", "GCN4")
rownames(allallelecomb) = sapply(1:nrow(allallelecomb), function(x){paste0("Ind", x)})
insilicopopulation = createmypop(myvariants, allallelecomb, myindivargs)

resploidy2 = simulateSystemStochasticParallel(insilicosystem, insilicopopulation, simtime = 3600, nepochs = 2, ntrialsPerInd = 5000, simalgorithm = "SSA", returnStochModel = F)
resploidymerged2 = lapply(resploidy2, mergeAllelesAbundance)

molsList.hist = getlasttimepoint(resploidymerged2)

mols = setdiff(colnames(resploidymerged2[[1]]), c("time", "trial"))

qtleffectval = rep(qtleffectrange, each = 5)
dosage = rep(c("AAAA", "AAAa", "AAaa", "Aaaa", "aaaa"), length(qtleffectrange))
names(qtleffectval) = names(dosage) = colnames(molsList.hist[[1]])

#mycols = brewer.pal(5, "YlOrRd")
mycols = c("#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8","#253494")
names(mycols) = c("aaaa", "Aaaa", "AAaa", "AAAa", "AAAA")


pretoplot = melt(molsList.hist[["R1"]])
toplot = data.frame("Rep" = pretoplot$Var1, "Ind" = pretoplot$Var2, "Value" = pretoplot$value, "QTLeffect" = factor(qtleffectval[pretoplot$Var2], levels = sort(qtleffectrange), labels = sapply(sort(qtleffectrange), as.character)), "dosage" = dosage[pretoplot$Var2])

plotploidy2 = ggplot(toplot, aes(x = Value, y = QTLeffect, fill = dosage, color = dosage)) + geom_density_ridges(alpha = 0.8, rel_min_height = 0.01, scale = 1.5) +
  scale_color_manual(limits = rev(levels(toplot$dosage)), name = "Simulated individuals", values = mycols, labels = c("AAAA" = "Ind1 (AAAA)", "AAAa" = "Ind2 (AAAa)", "AAaa" = "Ind3 (AAaa)", "Aaaa" = "Ind4 (Aaaa)", "aaaa" = "Ind5 (aaaa)")) +
  scale_fill_manual(limits = rev(levels(toplot$dosage)), name = "Simulated individuals",values = mycols, labels = c("AAAA" = "Ind1 (AAAA)", "AAAa" = "Ind2 (AAAa)", "AAaa" = "Ind3 (AAaa)", "Aaaa" = "Ind4 (Aaaa)", "aaaa" = "Ind5 (aaaa)")) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.border = element_blank(), panel.spacing = unit(1, "lines"),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey45"),
                     strip.text.x = element_text(size=12), strip.background = element_rect(colour = "white", fill = "white"), legend.title = element_text( vjust = 0.5)) +
  xlab("RNA abundance (# molecules)") + ylab("QTL effect coefficient") + ggtitle("RNA abundance at t = 1 hour (5000 simulations)")
print(plotploidy2)

ggsave("/media/sf_data/confirmation_results/plot3ploidy2.pdf", plot = plotploidy2, device = "pdf")

save(resploidy2, file = "/media/sf_data/confirmation_results/conf_rep_res23ploidy2.RData")



##########################################################################################################################################################
#                                                     OTHER PLOTS
##########################################################################################################################################################

outdeg = 1:100
l = 5
g = 1/l
exp = (1/l)*exp(-outdeg/l)
exp = exp/sum(exp)
pow = outdeg^(-g)
pow = pow/sum(pow)
toplot = data.frame("distribution" = rep(c("exp", "pow"), each = length(outdeg)), "kout" = rep(outdeg, 2), "proba" = c(exp, pow))

plotoutdeg = ggplot(toplot, aes(x = kout, y = proba, colour = distribution)) + geom_line() +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.border = element_blank(), panel.spacing = unit(1, "lines"),
                                                                                                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey45"),
                                                                                                    strip.text.x = element_text(size=12), strip.background = element_rect(colour = "white", fill = "white"), legend.title = element_text( vjust = 0.5)) +
  xlab("Out-degree") + ylab("Probability") + scale_colour_manual(name = "Distribution", limits = c("pow", "exp"), labels = c("pow" = "Power-law", "exp" = "Exponential"), values = c("pow" = "red", "exp" = "green3"))
print(plotoutdeg)

ggsave("/media/sf_data/confirmation_results/plotoutdegs.pdf", plot = plotoutdeg, device = "pdf")

# ------

outdeg = 1:50
l = 1:5
toplot = data.frame("kout" = numeric(), "l" = numeric(), "proba" = numeric())

for(coef in l){
  exp = (1/coef)*exp(-outdeg/coef)
  exp = exp/sum(exp)
  toplot = rbind(toplot, data.frame("kout" = outdeg, "l" = rep(coef, length(outdeg)), "proba" = exp))
}
toplot$l = factor(toplot$l, levels = paste0(l))
ggplot(toplot, aes(x = kout, y = proba, colour = l)) + geom_line()


toplot2 = data.frame("kout" = numeric(), "l" = numeric(), "proba" = numeric())

for(coef in l){
  pow = outdeg^(-coef)
  pow = pow/sum(pow)
  toplot2 = rbind(toplot2, data.frame("kout" = outdeg, "l" = rep(coef, length(outdeg)), "proba" = pow))
}
toplot2$l = factor(toplot2$l, levels = paste0(l))
ggplot(toplot2, aes(x = kout, y = proba, colour = l)) + geom_line()


# _____________________________________________________
## Simulate stochastic expression profiles
res = list()
t = seq(1, 50, by = 0.1)
x1 = 50 + 100*t/(t + 3)
x1rand = x1 + rnorm(length(x2), mean = 0, sd = 2)
res[[1]] = data.frame("t" = t, "mol" = paste0(1), "level" = x1rand, stringsAsFactors = F)

x2 = 80 - 20*t/(t + 8)
x2rand = x2 + rnorm(length(x2), mean = 0, sd = 1)
res[[2]] = data.frame("t" = t, "mol" = paste0(2), "level" = x2rand, stringsAsFactors = F)

x3rand = 110 + rnorm(length(x2), mean = 0, sd = 3)
res[[3]] = data.frame("t" = t, "mol" = paste0(3), "level" = x3rand, stringsAsFactors = F)

toplot = do.call(rbind, res)

exprplot = ggplot(toplot, aes(x = t, y = level, colour = mol)) + geom_line(size = 1) + scale_colour_discrete(guide = F) +
  xlab("Time") + ylab("Abundance") +
  theme_bw() + theme(panel.border = element_blank(), panel.spacing = unit(1, "lines"),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey45"), axis.title = element_text(size = 18), axis.text	= element_text(size = 14))
print(exprplot)

ggsave("C:\\Users\\oangelin\\OneDrive - Massey University\\Documents\\Confirmation\\Presentation\\exprprofiles.png", exprplot, device = "png")


# _____________________________________________________
## Simulate deterministic vs stochastic expression profile

t = seq(1, 50, by = 0.1)
x1 = 50 + 100*t/(t + 3)
x1rand = x1 + rnorm(length(x1), mean = 0, sd = 2)

toplotdet = data.frame("t" = t, "level" = x1, stringsAsFactors = F)

exprplotdet = ggplot(toplotdet, aes(x = t, y = level)) + geom_line(size = 1, colour = "red") +
  xlab("Time") + ylab("Abundance") +
  theme_bw() + theme(panel.border = element_blank(), panel.spacing = unit(1, "lines"),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey45"), axis.title = element_text(size = 18), axis.text	= element_text(size = 14))
print(exprplotdet)

ggsave("C:\\Users\\oangelin\\OneDrive - Massey University\\Documents\\Confirmation\\Presentation\\detprofile.png", exprplotdet, device = "png")

toplotsto = data.frame("t" = t, "level" = x1rand, stringsAsFactors = F)

exprplotsto = ggplot(toplotsto, aes(x = t, y = level)) + geom_line(size = 1, colour = "red") +
  xlab("Time") + ylab("Abundance") +
  theme_bw() + theme(panel.border = element_blank(), panel.spacing = unit(1, "lines"),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey45"), axis.title = element_text(size = 18), axis.text	= element_text(size = 14))
print(exprplotsto)

ggsave("C:\\Users\\oangelin\\OneDrive - Massey University\\Documents\\Confirmation\\Presentation\\stoprofile.png", exprplotsto, device = "png")
