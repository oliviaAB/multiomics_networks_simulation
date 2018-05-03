##########################################################################################################################
###                                JULIA FUNCTIONS FOR NETWORK GENERATION                                              ###
##########################################################################################################################
#=if !haskey(Pkg.installed(), "StatsBase") 
	Pkg.add("StatsBase")
end=#
using StatsBase

# ------------------------------------------- #
## FUNCTIONS FOR SAMPLING FROM DISTRIBUTIONS ##
# ------------------------------------------- #

## Sample from an discrete exponential distribution
function sampleexpon(n, lambda, max)
  prob = (1/float(lambda))*exp.(-(1:max)/float(lambda))
  cumprob = cumsum(prob / sum(prob))
  rnb = rand(Int(n))
  res = Int64[]
  for nb in rnb
    push!(res, findfirst(x -> x >= nb, cumprob))
  end 
  return res
end


## Sample from an discrete power-law distribution
function samplepowerlaw(n, gamma, max)
  prob = (1:max).^(-float(gamma))
  cumprob = cumsum(prob / sum(prob))
  rnb = rand(Int(n))
  res = Int64[]
  for nb in rnb
    push!(res, findfirst(x -> x >= nb, cumprob))
  end 
  return res
end


# ------------------------------------------------------------------------------------------------ #
## FUNCTIONS FOR COMPUTING PROBA FOR PREFERENTIAL ATTACHMENT OR "INVERSE PREFERENTIAL ATTACHMENT" ##
# ------------------------------------------------------------------------------------------------ #

# Linear preferential attachment from Barabasi-Albert
function probaPA(nodes, edg)
  ki = getInDeg(nodes, edg) + 1 # in-degree of each node, add 1 so that nodes with 0 in-degree can still be picked
  res = ki/sum(ki)
  res = res/sum(res)
  return res
end

# "Inverse" preferential attachment from Lachgar
function probaInvPA(nodes, edg)
  ki = getInDeg(nodes, edg) # in-degree of each node 
  if sum(ki) == 0
    res = fill(1, length(nodes))
  else
    res = 1 - ki/sum(ki)
  end
  res = res / sum(res)
  return res
end



# ------------------------------------------------------- #
## FUNCTIONS FOR GETTING THE IN- AND OUT-DEGREE OF NODES ##
# ------------------------------------------------------- #
## Edges of a graph are stored in a mx2 array (where m is the number of edges)

function getInDeg(nodes, edg)
  to = sort(edg[:,2])
  res = [length(searchsorted(to, n)) for n in nodes]
  return res
end


function getOutDeg(nodes, edg)
  from = sort(edg[:,1])
  res = [length(searchsorted(from, n)) for n in nodes]
  return res
end

# Function for checking if an edge between two nodes exist
# from is an array of nodes, to is a single nodes
function isEdge(from, to, edg)
  edgS = edg[sortperm(edg[:,2]),:]
  toedg = edgS[searchsorted(edgS[:,2], to),:]
  toedg = sort(toedg[:,1])
  [length(searchsorted(toedg, f))>0 for f in from]
end



# ------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------- #



## MAIN FUNCTION - GENERATE A GRAPH WITH SPECIFIED IN- AND OUT- DEGREE DISTRIBUTION


function nwgeneration(reg, target, indeg, outdeg, outdegexp, autoregproba, twonodesloop, edg = Array{Int64}(0,2), targetweight = []) 
## Input:
##    - reg: list of regulator nodes
##    - target: list of target nodes
##    - indeg: string variable (either "exponential" or "powerlaw") specifying the type of preferential attachment used to construct the network
##    - outdeg: string variable (either "exponential" or "powerlaw") specifying the type of distribution from which the out-degree of regulators are sampled
##    - outdegexp: the exponent of the out-degree distribution
##    - autoregproba: probability that a regulatory molecule regulates itself
##    - twonodesloop: do we allow 2-nodes loops? can be true or false
##    - optional argument edg: a 2-D array of edges already existing, to take into account; if not given creates an empty array
##    - optional argument targetweight: an array of weight for the targets, defining the probability of each target to be selected (multiplied to the probability computed from in-degree of target)

## Output:
##    - edg: A 2D array of edges, 1st column: from, 2nd column: to

  # Ensure that reg and target are arrays
  if typeof(reg) == Int64
    reg = [reg]
  end
  if typeof(target) == Int64
    target = [target]
  end

  ## If a target weight vector is provided check that its length matches the number of target
  if length(targetweight)>0 & length(targetweight)!=length(target)
  	error("targetweight must match the length of target")
  end
  ## If no target weight is provided targetweight will simply be one for all targets
  if length(targetweight) == 0
  	targetweight = fill(1.0, length(target))
  end


  ## Get the function for sampling from the desired out- degree distribution
  if outdeg == "exponential"
    foutdeg = getfield(current_module(), Symbol("sampleexpon"))
  elseif outdeg == "powerlaw"
    foutdeg = getfield(current_module(), Symbol("samplepowerlaw"))  
  else
    error("Argument outdeg non-valid: must be exponential or powerlaw")
  end

  ## Get the function for computing the probability of each target node to receive a new incoming edge for the desired in-degree distribution
  ##  If in-degree is power law, use the model of preferential attachment from Barabasi-Albert
  ##  If in-degree is exponential, use the model of preferential attachment from Lachgar
  if indeg == "exponential"
    findeg = getfield(current_module(), Symbol("probaInvPA"))
  elseif indeg == "powerlaw"
    findeg = getfield(current_module(), Symbol("probaPA"))
  else
    error("Argument indeg non-valid: must be exponential or powerlaw")
  end

  ## Sample the number of target (out-degree) for each regulator
  out = foutdeg(length(reg), outdegexp, length(target))
  sort!(out, rev = true)

  ## Create the mx2 array of edges, 1st column = from, 2nd column = to
  # edg = Array{Int64}(0,2)

  ## For each regulator, sample from the target list its target according to its number of targets specified in the out variable
  for r in eachindex(reg)

    probTar = findeg(target, edg) # compute for each target the probability of being regulated by r
    # probTar = targetweight .* findeg(target, edg) # compute for each target the probability of being regulated by r, weighted by the weight of each target

    ## we exclude the regulator r from the list of potential targets (autoregulation is treated later)
    if reg[r] in target
      probTar[findfirst(y -> y == reg[r], target)] = 0
    end

    ## How to deal with loops, i.e. if one or more target(s) already control the regulator r
    ## if twonodesloop == false we don't authorise the regulator to target a node that controls it
    if !twonodesloop
      exEdg = isEdge(target, reg[r], edg)
      probTar[exEdg] = 0
    end

    ## Make sure that the out-degree of regulator r doesn't exceed the number of targets with non-null proba 
    out[r] = min(out[r], sum(probTar .> 0))
    if out[r] == 0
      println("OUT-DEGREE SET TO 0")
    end

    ## Sample targets of regulator
    sa = StatsBase.sample(target, Weights(probTar), out[r], replace = false)

    ## Add the created edges in edg
    edg = vcat(edg, [fill(reg[r], out[r]) sa])

    ## Add an autoregulatory edge with probability autoregproba
    if rand() <= autoregproba
      edg = vcat(edg, [reg[r] reg[r]])
    end

  end

  return edg

end




function nwgenerationFromInDegree(reg, target, indeg, outdeg, outdegexp, autoregproba, twonodesloop, edg = Array{Int64}(0,2), targetweight = []) 
## Input:
##    - reg: list of regulator nodes
##    - target: list of target nodes
##    - indeg: string variable (either "exponential" or "powerlaw") specifying the type of preferential attachment used to construct the network
##    - outdeg: string variable (either "exponential" or "powerlaw") specifying the type of distribution from which the out-degree of regulators are sampled
##    - outdegexp: the exponent of the out-degree distribution
##    - autoregproba: probability that a regulatory molecule regulates itself
##    - twonodesloop: do we allow 2-nodes loops? can be true or false
##    - optional argument edg: a 2-D array of edges already existing, to take into account; if not given creates an empty array
##    - optional argument targetweight: an array of weight for the targets, defining the probability of each target to be selected (multiplied to the probability computed from in-degree of target)

## Output:
##    - edg: A 2D array of edges, 1st column: from, 2nd column: to


  # Ensure that reg and target are arrays
  if typeof(reg) == Int64
    reg = [reg]
  end
  if typeof(target) == Int64
    target = [target]
  end

  ## If a target weight vector is provided check that its length matches the number of target
  if length(targetweight)>0 & length(targetweight)!=length(target)
  	error("targetweight must match the length of target")
  end
  ## If no target weight is provided targetweight will simply be one for all targets
  if length(targetweight) == 0
  	targetweight = fill(1.0, length(target))
  end


  ## Get the function for sampling from the desired out- degree distribution
  if outdeg == "exponential"
    foutdeg = getfield(current_module(), Symbol("sampleexpon"))
  elseif outdeg == "powerlaw"
    foutdeg = getfield(current_module(), Symbol("samplepowerlaw"))  
  else
    error("Argument outdeg non-valid: must be exponential or powerlaw")
  end

  ## Get the function for computing the probability of each target node to receive a new incoming edge for the desired in-degree distribution
  ##  If in-degree is power law, use the model of preferential attachment from Barabasi-Albert
  ##  If in-degree is exponential, use the model of preferential attachment from Lachgar
  if indeg == "exponential"
    findeg = getfield(current_module(), Symbol("probaInvPA"))
  elseif indeg == "powerlaw"
    findeg = getfield(current_module(), Symbol("probaPA"))
  else
    error("Argument indeg non-valid: must be exponential or powerlaw")
  end

  ## Sample the number of target (out-degree) for each regulator
  out = foutdeg(length(reg), outdegexp, length(target))
  sort!(out, rev = true)

  ## Create the mx2 array of edges, 1st column = from, 2nd column = to
  # edg = Array{Int64}(0,2)

  ## For each regulator, sample from the target list its target according to its number of targets specified in the out variable
  for r in eachindex(reg)

    probTar = findeg(target, edg) # compute for each target the probability of being regulated by r
    # probTar = targetweight .* findeg(target, edg) # compute for each target the probability of being regulated by r, weighted by the weight of each target

    ## we exclude the regulator r from the list of potential targets (autoregulation is treated later)
    if reg[r] in target
      probTar[findfirst(y -> y == reg[r], target)] = 0
    end

    ## How to deal with loops, i.e. if one or more target(s) already control the regulator r
    ## if twonodesloop == false we don't authorise the regulator to target a node that controls it
    if !twonodesloop
      exEdg = isEdge(target, reg[r], edg)
      probTar[exEdg] = 0
    end

    ## Make sure that the out-degree of regulator r doesn't exceed the number of targets with non-null proba 
    out[r] = min(out[r], sum(probTar .> 0))
    if out[r] == 0
      println("OUT-DEGREE SET TO 0")
    end

    ## Sample targets of regulator
    sa = StatsBase.sample(target, Weights(probTar), out[r], replace = false)

    ## Add the created edges in edg
    edg = vcat(edg, [fill(reg[r], out[r]) sa])

    ## Add an autoregulatory edge with probability autoregproba
    if rand() <= autoregproba
      edg = vcat(edg, [reg[r] reg[r]])
    end

  end

  return edg

end





function nwgenerationSR(regList, target, indeg, outdegList, outdegexpList, autoregprobaList, twonodesloop) 
## Input:
##    - reg: list of lists of regulator nodes (e.g. the first list represents the TFs, and the second list the lncRNAs). Each regulator list can have a different out-degree distribution
##            but the edges is added considering the "global" in-degree of the targets, i.e. without distinction between the different types of regulators
##    - target: list of target nodes
##    - indeg: string variable (either "exponential" or "powerlaw") specifying the type of preferential attachment used to construct the network
##    - outdeg: list of string variables (either "exponential" or "powerlaw") specifying the type of distribution from which the out-degree of regulators are sampled for each regulator list
##    - outdegexp: list of the exponents of the out-degree distribution for each of the regulator list
##    - autoregprobaList: list of probability that a regulatory molecule regulates itself for each regulator list
##    - twonodesloop: do we allow 2-nodes loops? can be true or false


## Output:
##    - edg: A 2D array of edges, 1st column: from, 2nd column: to

  ## Check that the length of the the regulator list matches the length of the out-degree distribution list and out degree distribution exponent list
  if length(regList)!=length(outdegList) | length(regList)!=length(outdegexpList) | length(regList)!=length(autoregprobaList)
    error("Make sure the length of regList, outdegList, outdegexpList and autoregprobaList is the same \n")
  end

  ## For each list of regulators, sample an out-degree from specified distribution for each regulator
  out = []
  reg = [] # gives the ID of regulators of the different lists
  autoregproba = [] # gives the probability of doing autoregulation for each regulator 

  for l in 1:length(regList)

    ## Get the function for sampling from the desired out- degree distribution
    if outdegList[l] == "exponential"
      foutdeg = getfield(current_module(), Symbol("sampleexpon"))
    elseif outdegList[l] == "powerlaw"
      foutdeg = getfield(current_module(), Symbol("samplepowerlaw"))  
    else
      error("Argument outdeg non-valid: must be exponential or powerlaw")
    end

    ## Sample the number of target (out-degree) for each regulator in the regulator list l
    append!(out, foutdeg(length(regList[l]), outdegexpList[l], length(target)))
    append!(reg, regList[l])  ## keep in memory the id of each regulator associated to the out-degree just sampled
    append!(autoregproba, fill(autoregprobaList[l], length(regList[l]))) ## keep in memory for each regulator its proba of creating an autoregulatory edge
  end

  reg = reg[sortperm(out, rev = true)] ## sort regulator id according to their out-degree
  autoregproba = autoregproba[sortperm(out, rev = true)] ## sort accordingly the probability of each regulator to perform autoregulation
  sort!(out, rev = true) ## sort the out-degree of the different regulators, without considering the list they are from

  ## Get the function for computing the probability of each target node to receive a new incoming edge for the desired in-degree distribution
  ##  If in-degree is power law, use the model of preferential attachment from Barabasi-Albert
  ##  If in-degree is exponential, use the model of preferential attachment from Lachgar
  if indeg == "exponential"
    findeg = getfield(current_module(), Symbol("probaInvPA"))
  elseif indeg == "powerlaw"
    findeg = getfield(current_module(), Symbol("probaPA"))
  else
    error("Argument indeg non-valid: must be exponential or powerlaw")
  end



  ## Create the mx2 array of edges (m = number of edges), 1st column = from, 2nd column = to
  edg = Array{Int64}(0,2)

  ## For each regulator, sample from the target list its target according to its number of targets specified in the out variable
  for r in eachindex(reg)

    probTar = findeg(target, edg) # compute for each target the probability of being regulated by r


    ## we exclude the regulator r from the list of potential targets (autoregulation is treated later)
    if reg[r] in target
      probTar[findfirst(y -> y == reg[r], target)] = 0
    end

    ## How to deal with loops, i.e. if one or more target(s) already control the regulator r
    ## if twonodesloop == false we don't authorise the regulator to target a node that controls it
    if !twonodesloop
      exEdg = isEdge(target, reg[r], edg)
      probTar[exEdg] = 0
    end

    ## Make sure that the out-degree of regulator r doesn't exceed the number of targets with non-null proba 
    out[r] = min(out[r], sum(probTar .> 0))
    if out[r] == 0
      println("OUT-DEGREE SET TO 0")
    end

    ## Sample targets of regulator
    sa = StatsBase.sample(target, Weights(probTar), out[r], replace = false)

    ## Add the created edges in edg
    edg = vcat(edg, [fill(reg[r], out[r]) sa])

    ## Add an autoregulatory edge with probability autoregproba
    if rand() <= autoregproba[r]
      edg = vcat(edg, [reg[r] reg[r]])
    end

  end

  return edg

end


function whatisit(x)
  println(x)
  println(typeof(x))
end







### TEMPORARY - TRIALS

#=
reg = collect(1:3)
target = collect(1:10)
indeg = "exponential"
outdeg = "powerlaw"
outdegexp = 0.8
autoregproba = 0.1
twonodesloop = false



regList = [collect(1:5), collect(6:12)]
target = collect(1:50)
indeg = "exponential"
outdegList = ["powerlaw", "powerlaw"]
outdegexpList = [2.2, 1]
autoregprobaList = [0.2, 0.9]
twonodesloop = false
=#