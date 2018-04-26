##########################################################################################################################
###                                JULIA FUNCTIONS FOR NETWORK GENERATION                                              ###
##########################################################################################################################

using StatsBase

# ------------------------------------------- #
## FUNCTIONS FOR SAMPLING FROM DISTRIBUTIONS ##
# ------------------------------------------- #

## Sample from an discrete exponential distribution
function sampleexpon(n, lambda, max)
  prob = (1/lambda)*exp.(-(1:max)/lambda)
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
  prob = (1:max).^(-gamma)
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
  return [sum(map(i -> i == to, edg[:,2])) for to in nodes]
end


function getOutDeg(nodes, edg)
  return [sum(map(i -> i == from, edg[:,1])) for from in nodes]
end


# Function for checking if an edge between two nodes exist
# from is an array of nodes, to is a single nodes
function isEdge(from, to, edg)
  [any(all(edg .== [f to], 2)) for f in from]
end



# ------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------- #



## MAIN FUNCTION - GENERATE A GRAPH WITH SPECIFIED IN- AND OUT- DEGREE DISTRIBUTION

function nwgeneration(reg, target, indeg, outdeg, outdegexp) 

  # Ensure that reg and target are arrays
  if typeof(reg) == String
    reg = [reg]
  end
  if typeof(target) == String
    target = [target]
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
  #out = sort!(out, rev = true)

  ## Create the mx2 array of edges (m = number of edges), 1st column = from, 2nd column = to
  edg = fill("", (sum(out),2))

  ## For each regulator, sample from the target list its target according to its number of targets specified in the out variable
  for r in eachindex(reg)

    println(r)
    probTar = findeg(target, edg) # compute for each target the probability of being regulated by r

    ## How to deal with self-regulatory edges: if regulator r is in 
    ##   Here we try to limit the number of self-regulatory edges
    if reg[r] in target
      probTar[findfirst(y -> y == reg[r], target)] = probTar[findfirst(y -> y == reg[r], target)] / 2
    end

    ## How to deal with loops, i.e. if one or more target(s) already control the regulator r
    ## Here we try to reduce the number of loops
    probTar[isEdge(target, reg[r], edg)] = probTar[isEdge(target, reg[r], edg)]/2

    ## Make sure that the out-degree of regulator r doesn't exceed the number of targets with non-null proba 
    out[r] = min(out[r], sum(probTar .> 0))
    if out[r] == 0
      println("OUT-DEGREE SET TO 0")
    end

    ## Sample targets of regulator
    sa = StatsBase.sample(target, Weights(probTar), out[r], replace = false)

    ## Add the created edges in edg
    edg[(cumsum(out)-out+1)[r]:cumsum(out)[r], 1] = reg[r] ## edges come from the regulator ...
    edg[(cumsum(out)-out+1)[r]:cumsum(out)[r], 2] = sa ## ... and go to each target

  end

  return edg

end


function whatisit(x)
  println(typeof(x))
end


### TEMPORARY - TRIALS

#=target = ["G$(i)" for i=1:10]
reg = target[1:4]
indeg = "exponential"
outdeg = "powerlaw"
outdegexp = 0.8=#