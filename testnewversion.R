library(tictoc)

setwd("~/winData/multiomics_networks_simulation")
#setwd("~/GitHub/multiomics_networks_simulation")
source("network_generation.R")

mysystemargs = insilicosystemargs(G = 50)
myindividualargs = insilicoindividualargs(ploidy = 4)

mysystemgenes = createGenes(mysystemargs)

mymosystem = createMultiOmicNetwork(mysystemgenes, mysystemargs)


sysargs = mysystemargs
indargs = myindividualargs
mosystem = mymosystem
nod = mysystemgenes

tic()
for(i in 1:1){
  mysystemargs = insilicosystemargs(G = 100)
  mysystemgenes = createGenes(mysystemargs)
  mymosystem = createMultiOmicNetwork(mysystemgenes, mysystemargs)
  }
toc()


##############################################
library(RColorBrewer)
library(gridExtra)

setwd("~/winData/multiomics_networks_simulation")
#setwd("~/GitHub/multiomics_networks_simulation")
source("network_generation.R")

mysystemargs = insilicosystemargs(G = 50)
myindividualargs = insilicoindividualargs(ploidy = 4)

mysystemgenes = createGenes(mysystemargs)

mymosystem = createMultiOmicNetwork(mysystemgenes, mysystemargs)

## plotmosystem

## Define the colour palettes
mycolsCS = c("darkorchid4", "darkgreen")
names(mycolsCS) = c("PC", "NC")

mycolsGF = brewer.pal(6, "RdYlBu")
names(mycolsGF) = c("TC", "TL", "RD", "PD", "PTM", "MR")

## Plot composition of system
gCS = ggplot(mysystemgenes, aes(x = "1", fill = coding)) +
  geom_bar() + 
  scale_fill_manual(values = mycolsCS, drop = F, name = "Coding status", labels = c("PC" = "Protein coding", "NC" = "Noncoding")) + 
  coord_flip() + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Coding status of genes in the system") + ylab("Number of genes") 


gGF = ggplot(mysystemgenes, aes(x = coding, fill = factor(TargetReaction, levels = rev(c("TC", "TL", "RD", "PD", "PTM", "MR"))))) +
  geom_bar() + 
  scale_fill_manual(values = mycolsGF, breaks = c("TC", "TL", "RD", "PD", "PTM", "MR"), drop = F, name = "Gene function ") + 
  scale_x_discrete(limits = c("PC", "NC"), labels = c("Protein-coding", "Noncoding")) +
  coord_flip() +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5), plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Coding status and function of genes in the system") + xlab("Coding status") + ylab("Number of genes")

grid.arrange(gCS, gGF, ncol = 1, heights = c(0.3, 0.7))


## Plot overview of regulators

genindex = 1:mysystemargs$G
indegtot = sapply(genindex, function(x){sum(mymosystem$edg$to == x)})
indegGF = sapply(genindex, function(x){sum(mymosystem$edg$to == x & mymosystem$edg$TargetReaction == "TC")})
genorder = genindex[order(indegGF, decreasing = T)]
indegtot = indegtot[order(indegGF, decreasing = T)]
genorder = genorder[order(indegtot, decreasing = T)]

gIDTGF = ggplot(mymosystem$edg, aes(x = factor(to, as.character(genorder)), fill = factor(TargetReaction, levels = rev(c("TC", "TL", "RD", "PD"))))) + 
  geom_bar() + 
  scale_fill_manual(values = mycolsGF, breaks = c("TC", "TL", "RD", "PD"), name = "Type of\nregulation") + 
  theme(axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Number and type of regulators for each gene") + xlab("Genes in the system") + ylab("Number regulators")

# grid.arrange(gCS, gGF, gIDTGF, layout_matrix = rbind(c(1, 1, 1, 3, 3, 3), c(2, 2, 2, 3, 3, 3), c(2, 2, 2, 3, 3, 3)))

genindex = 1:mysystemargs$G
indegtot = sapply(genindex, function(x){sum(mymosystem$edg$to == x)})
indegCS = sapply(genindex, function(x){sum(mymosystem$edg$to == x & mymosystem$edg$RegBy == "PC")})
genorder = genindex[order(indegCS, decreasing = T)]
indegtot = indegtot[order(indegCS, decreasing = T)]
genorder = genorder[order(indegtot, decreasing = T)]

gIDTCS = ggplot(mymosystem$edg, aes(x = factor(to, as.character(genorder)), fill = factor(RegBy, levels = c("NC", "PC")))) + 
  geom_bar() + 
  scale_fill_manual(values = mycolsCS, breaks = c("PC", "NC"), name = "Type of\nregulator") + 
  theme(axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Number and type of regulators for each gene") + xlab("Genes in the system") + ylab("Number regulators")

grid.arrange(gIDTCS, gIDTGF, ncol = 2)

## Plot kinetic parameters of the genes

concnod = list("TCrate" = c("PC", "NC"), "TLrate" = c("PC"), "RDrate" = c("PC", "NC"), "PDrate" = c("PC"))
graphtitles = list("TCrate" = "Transcription rate", "TLrate" = "Translation rate", "RDrate" = "RNA decay rate", "PDrate" = "Protein decay")

for(g in names(concnod)){
  myplot = ggplot(mysystemgenes[mysystemgenes$coding %in% concnod[[g]],], aes_string(x = g)) + geom_histogram() +
    ggplot(graphtitles[[g]]) + xlab(paste(graphtitles[[g]], "(1/s)", sep = " "))
    
  assign(paste0("plot", g), myplot)
}

grid.arrange(plotTCrate, plotTLrate, plotRDrate, plotPDrate, ncol = 2)

## Plot each regulatory network

temp = c("TC", "TL", "RD", "PD", "PTM")
concnod = list("TC" = c("PC", "NC"), "TL" = c("PC"), "RD" = c("PC", "NC"), "PD" = c("PC"), "PTM" = c("PC"))

for(t in temp){
  edgtot = mymosystem$edg[mymosystem$edg$TargetReaction == t, ]
  edgcomp = mymosystem[[paste0(t, "RN.edg")]]
  
  
}
