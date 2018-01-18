rm(list = ls(all = T))
cat("\n \n \n \n \n \n \n \n")

source("simulation.R")

G = 5
P = 3
M = 3
MR = 1

network = rand_network(G, P, M, MR)
indiv = rand_indiv(network)
