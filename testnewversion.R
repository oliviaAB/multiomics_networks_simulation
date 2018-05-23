library(tictoc)

setwd("~/winData/multiomics_networks_simulation")
#setwd("~/GitHub/multiomics_networks_simulation")



source("network_generation.R")


mysystemargs = insilicosystemargs(G = 50)
myinsilicosystem = createInSilicoSystem(mysystemargs)

myindivargs = insilicoindividualargs()
mypopulation = createPopulation(20, myinsilicosystem, myindivargs)





##############################################

source("network_generation.R")

mysystemargs = insilicosystemargs(G = 50)
myindividualargs = insilicoindividualargs(ploidy = 4)

mysystemgenes = createGenes(mysystemargs)

mymosystem = createMultiOmicNetwork(mysystemgenes, mysystemargs)

myinsilicosystem = createInSilicoSystem(mysystemargs)

sysargs = mysystemargs
indargs = myindividualargs
mosystem = mymosystem
nod = mysystemgenes

tic()
for(i in 1:100){
  mysystemargs = insilicosystemargs(G = 100)
  myinsilicosystem = createInSilicoSystem(mysystemargs)
  }
toc()

##############################################
library(RColorBrewer)
library(gridExtra)
library(tictoc)

setwd("~/winData/multiomics_networks_simulation")
#setwd("~/GitHub/multiomics_networks_simulation")
source("network_generation.R")

mysystemargs = insilicosystemargs(G = 500)
tic(); myinsilicosystem = createInSilicoSystem(mysystemargs); toc()

## plotmosystem

plotGlobalSystem(myinsilicosystem)

## ---------------------------- ##
## Plot each regulatory network ##
## ---------------------------- ##

mycolsCS = c("PC" = "#b13e25",  "NC" = "#602377", "Tot" = "#31161F")

mycolsGF = brewer.pal(6, "RdYlBu")
names(mycolsGF) = c("TC", "TL", "RD", "PD", "PTM", "MR")

reactionsnames = c("TC" = "transcription", "TL" = "translation", "RD" = "RNA decay", "PD" = "protein decay", "PTM" = "protein post-translational modification", "MR" = "metabolic reaction")


temp = c("TC", "TL", "RD", "PD", "PTM")
concnod = list("TC" = c("PC", "NC"), "TL" = c("PC"), "RD" = c("PC", "NC"), "PD" = c("PC"), "PTM" = c("PC"))


for(t in temp){
  edgtot = insilicosystem$mosystem$edg[insilicosystem$mosystem$edg$TargetReaction == t, ]
  edgcomp = insilicosystem$mosystem[[paste0(t, "RN.edg")]]
  genesid = insilicosystem$genes[, c("id", "coding", "TargetReaction")]
  
  ## Plot out-degree distribution
  
  regs = unique(edgtot$from)
  outdegreedf = data.frame("Regid" = regs, "RegCoding" = genesid[regs, "coding"], "Outdegree" = sapply(regs, function(x){sum(edgtot$from == x)}))
  if(nrow(outdegreedf) == 0)   outdegreedf = data.frame("Regid" = numeric(0), "RegCoding" = character(0), "Outdegree" = numeric(0))
  
  gODT = ggplot(outdegreedf, aes(x = Outdegree)) + 
    geom_histogram(binwidth = 10, center = 5, fill = mycolsCS["Tot"]) + 
    theme( plot.title = element_text(hjust = 0.5)) + 
    ggtitle(paste("Out-degree distribution of", reactionsnames[t], "regulators")) + xlab("Number of targets") + ylab("Frequency")
  
  gODPC = ggplot(outdegreedf[outdegreedf$RegCoding == "PC", ], aes(x = Outdegree)) + 
    geom_histogram(binwidth = 5, center = 2.5, fill = mycolsCS["PC"]) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    scale_fill_discrete(guide=FALSE) +
    ggtitle(paste("Out-degree distribution of protein-coding", reactionsnames[t], "regulators")) + xlab("Number of targets") + ylab("Frequency")
  
  gODNC = ggplot(outdegreedf[outdegreedf$RegCoding == "NC", ], aes(x = Outdegree)) + 
    geom_histogram(binwidth = 5, center = 2.5, fill = mycolsCS["NC"]) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    scale_fill_discrete(guide=FALSE) +
    ggtitle(paste("Out-degree distribution of noncoding", reactionsnames[t], "regulators")) + xlab("Number of targets") + ylab("Frequency")
  
  ggarrange(gODT, gODPC, gODNC, ncol = 3)
  
  
  ## Plot out-degree distribution
  tars = genesid[genesid$coding %in% concnod[[t]], "id"]
  indegreedf = data.frame("Tarid" = tars, "IndegreePC" = sapply(tars, function(x){sum(edgtot$to[edgtot$RegBy == "PC"] == x)}), "IndegreeNC" = sapply(tars, function(x){sum(edgtot$to[edgtot$RegBy == "NC"] == x)}), "IndegreeTot" = rep(0, length(tars)))
  indegreedf$IndegreeTot = indegreedf$IndegreePC + indegreedf$IndegreeNC
  
  ## to test whether there is any regulation of this type we refer to the data frame outdegreedf (because indegreedf will never have 0 rows unless the target id list is empty)
  if(nrow(outdegreedf) == 0)   indegreedf = data.frame("Tarid" = numeric(0), "IndegreePC" = numeric(0), "IndegreeNC" = numeric(0), "IndegreeTot" = numeric(0))
  
  
  gIDT = ggplot(indegreedf, aes(x = IndegreeTot)) + 
    geom_histogram(binwidth = 1, center = 0.5) + 
    theme( plot.title = element_text(hjust = 0.5)) + 
    ggtitle(paste("In-degree distribution of genes -\n", reactionsnames[t], "regulation")) + xlab("Number of regulators") + ylab("Frequency")
  
  gIDPC = ggplot(indegreedf, aes(x = IndegreePC)) + 
    geom_histogram(binwidth = 1, center = 0.5, fill = mycolsCS["PC"]) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    scale_fill_discrete(guide=FALSE) +
    ggtitle(paste("In-degree distribution of genes -\n only protein coding", reactionsnames[t], "regulators")) + xlab("Number of protein coding regulators") + ylab("Frequency")
  
  gIDNC = ggplot(indegreedf, aes(x = IndegreeNC)) + 
    geom_histogram(binwidth = 1, center = 0.5, fill = mycolsCS["NC"]) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    scale_fill_discrete(guide=FALSE) +
    ggtitle(paste("In-degree distribution of genes -\n only noncoding", reactionsnames[t], "regulators")) + xlab("Number of noncoding regulators") + ylab("Frequency")
  
  # ggarrange(gIDT, gIDPC, gIDNC, ncol = 3)
  
  ggarrange(gODT, gODPC, gODNC, gIDT, gIDPC, gIDNC, nrow = 2, ncol = 3)
  
  
  ## Plot kinetic parameters of the interactions
  kineticplotList = list()
  for(k in grep(paste0("^", t), names(edgcomp), value = T)){
    myplot = ggplot(edgcomp, aes_string(x = k)) +
      geom_histogram() + 
      ggtitle(sub(paste0("^",t), "", k)) + xlab(sub(paste0("^",t), "", k)) + ylab("Frequency")
      
    assign(paste0("plot", k), myplot)
    kineticplotList[[k]] = myplot
  }
  
  if(length(kineticplotList) > 0) ggarrange(plots = kineticplotList, ncol = length(kineticplotList))
  
}


table(test$RegSign)
table(test$RegBy)
table(test[,c("RegSign", "RegBy")])
test = data.frame("edgid" = 1:nrow(test), test)


testfondue = melt(test, id.vars = "edgid", measure.vars = c("RegSign", "RegBy"))
ggplot(testfondue, aes(x = variable, fill = value)) + geom_bar() + coord_flip()
