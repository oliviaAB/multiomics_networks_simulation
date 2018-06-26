##########################################################################################################################
###                                JULIA FUNCTIONS FOR NETWORK GENERATION                                              ###
##########################################################################################################################

## If not already installed, install required packages
#if !haskey(Pkg.installed(), "StatsBase") 
#	Pkg.add("StatsBase")
#end

#= if !haskey(Pkg.installed(), "Combinatorics") 
  Pkg.add("Combinatorics")
end
=#

using StatsBase

using Combinatorics

using BioSimulator


# using Combinatorics

# --------------------------------------------------- #
## FUNCTIONS FOR GENERATING QTL VALUES FOR EACH GENE ##
# --------------------------------------------------- #



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


function nwgeneration(reg, target, indeg, outdeg, outdegexp, autoregproba, twonodesloop, edg = Array{Int64}(0,2)) #, targetweight = []) 
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

#=
  ## If a target weight vector is provided check that its length matches the number of target
  if length(targetweight)>0 & length(targetweight)!=length(target)
  	error("targetweight must match the length of target")
  end
  ## If no target weight is provided targetweight will simply be one for all targets
  if length(targetweight) == 0
  	targetweight = fill(1.0, length(target))
  end
=#

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

    ## we exclude the regulator r from the list of potential targets (autoregulation is addressed later)
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



function combreg(edgfrom, edgto, edgsign, p, complexsize, reac)

  edg = [edgfrom edgto edgsign]
  edg = edg[sortperm(edg[:,2]),:] ## sort edges according to the target ids
  complexes = Dict()
  complexesid = 1 ## to give a unique ID to each created complex
  rowstoremove = []
  edgtoadd = Array{Any}(0,3)
  complexsize = Int(complexsize)

  target = unique(edgto)

  for tar in target
    temp = searchsorted(edg[:, 2], tar)
    regs = edg[temp, :] ## look for the edges coming to tar
    
    # regspos = regs[regs[:,3] .== "1", :] ## look for the edges coming to tar representing positive regulation
    edgpos = temp[regs[:,3] .== "1"] ## stock the row number of the positive edges
    # regsneg = regs[regs[:,3] .== "-1", :] ## look for the edges coming to tar representing negative regulation
    edgneg = temp[regs[:,3] .== "-1"] ## stock the row number of the negative edges

    ntrypos = div(length(edgpos), complexsize) ## max number of complexes you can simultaneously create from the positive regulators of tar
    ntryneg = div(length(edgneg), complexsize) ## max number of complexes you can simultaneously create from the negative regulators of tar

    for i in 1:ntrypos
      if rand()<= p ## with probability p form a complex
        compo = sample(edgpos, complexsize, replace=false) ## sample the components of the complex (in fact sample rows in the edg matrix)
        edgpos = setdiff(edgpos, compo) ## remove the selected rows (regulators) from the list of possible components (future sampling)
        append!(rowstoremove, compo) ## the edges correpsponding to the selected components of the complex will be removed
        compid = string("C", reac, complexesid) ## create the new complex ID
        complexesid +=1
        complexes[compid] = edg[compo, 1] ## in the dictionary of complexes add the composition (ie array of components) of the new complex
        edgtoadd = vcat(edgtoadd, [compid tar "1"]) ## create a new regulatory edge from the complex to the target, with positive regulation
      end
    end

    ## repeat for the negative regulation
    for i in 1:ntryneg
      if rand()<= p ## with probability p form a complex
        compo = sample(edgneg, complexsize, replace=false) ## sample the components of the complex (in fact sample rows in the edg matrix)
        edgneg = setdiff(edgneg, compo) ## remove the selected rows (regulators) from the list of possible components (future sampling)
        append!(rowstoremove, compo) ## the edges correpsponding to the selected components of the complex will be removed
        compid = string("C", reac, complexesid) ## create the new complex ID
        complexesid +=1
        complexes[compid] = edg[compo, 1] ## in the dictionary of complexes add the composition (ie array of components) of the new complex
        edgtoadd = vcat(edgtoadd, [compid tar "-1"]) ## create a new regulatory edge from the complex to the target, with positive regulation
      end
    end

  end


  ## remove the individual edges from regulators forming a complex
  edg = edg[setdiff(collect(1:size(edg)[1]), rowstoremove), :]
  ## add the new edges corresponding to regulation by complexes
  edg = vcat(edg, edgtoadd)

  edg[:,1] = [string(i) for i in edg[:,1]] ## Transform the integer ID of regulators into String ID 

  return Dict("newedg" => edg, "Complexes" => complexes)

end







# ------------------------------------------------------------------------------------------------ #
##              FUNCTIONS FOR CREATING THE REACTION LIST AND PROPENSITY LIST                      ##
# ------------------------------------------------------------------------------------------------ #

## Input:
##    - promList: an 1D array, each elemnt being a 2D array. Each elements correspond to the possible active states of a given promoter binding site (1st column) and the associated fold change (2nd column)
function combinactivepromstates(promList)
  proms = [t[:, 1] for t in promList]
  fcs = [t[:, 2] for t in promList]
  elsize = [length(t) for t in proms]
  combproms = Array{String}(prod(elsize), length(elsize))
  combfcs = Array{Float64}(prod(elsize), length(elsize))

  for i in eachindex(elsize)
    combproms[:, i] = repeat(proms[i], inner = prod(elsize[(i+1):end]), outer = prod(elsize[1:(i-1)]))
    combfcs[:, i] = repeat(fcs[i], inner = prod(elsize[(i+1):end]), outer = prod(elsize[1:(i-1)]))
  end

  return Dict("proms" => combproms, "fcs" => combfcs)

end


## Input:
##    - promList: an 1D array, each elemnt being an array. Each elements correspond to all possible states of the binding site
function combinallpromstates(proms)
  elsize = [length(t) for t in proms]
  combproms = Array{String}(prod(elsize), length(elsize))

  for i in eachindex(elsize)
    combproms[:, i] = repeat(proms[i], inner = prod(elsize[(i+1):end]), outer = prod(elsize[1:(i-1)]))
  end

  return combproms

end



## Gives all the possible combinations of length elem of the vector vect
function allposscomb(vect, elem)

  nb = length(vect)
  comb = Array{eltype(vect)}(nb^elem, elem)

  for i in 1:elem
    comb[:,i] = repeat(vect, inner = [nb^(elem-i)], outer = [nb^(i-1)])
  end

  return comb
end

## Generate all possible combinations of promoter states given the number of promoter and their regulatory effect (activtor or repressor)
## Input:
##    - edgSign: array with for each regulator the sign of its regulation ("1" for activation and "-1" for repression)
function allactivepromstates(edgSign)

  act = find(edgSign .== "1")
  R = length(act)
  Rtot = length(edgSign)
  comb = zeros(Int64, 2^R, Rtot)

  temp = [0, 1]
  for i in 1:R
    comb[:,act[i]] = repeat(temp, inner = [2^(R-i)], outer = [2^(i-1)])
  end

  return comb
end


function reactBioSim(reacList, prodList)
  return join(reacList, " + ") * " --> " * join(prodList, " + ")
end





function createbindingpromreactions(edg, exprstep, promPrefix, nod, functform, complexes, complexeskinetics, complexsize, complexvariants, gcnList)

  spec = []
  inicond = []
  react = []
  reactnames = []
  propens = []

  regsingl = find( x -> !in('C', x), edg["from"]) ## identify single-molecule regulators (i.e. not regulatory complexes)
  regcompl = find( x -> in('C', x), edg["from"]) ## identify regulatory complexes

  promActiveStates = Dict(string(i)*j => [] for i in nod["id"], j in gcnList) ## dictionary, 1 element for each allele version of each gene. Value is an array with m elements, 1 for each regulator. 
                                                                              ## The j-th element of this array is a matrix: rows: all possible active states of the binding site of regulator j (column 1: name of the binding site state, column 2 fold-change associated with this state)
  promAllStates = Dict(string(i)*j => [] for i in nod["id"], j in gcnList) ## dictionary, 1 element for each allele version of each gene. Value is an array with m elements, 1 for each regulator. 
                                                                              ## The j-th element of this array is an array, listing all the possible states of the binding site of regulator j (even if not active)

  ## Generate the binding of single regulators on binding sites of targets
  for r in regsingl, gcn in gcnList
    tarid = edg["to"][r]
    tar = string(tarid) * gcn
    reg = edg["from"][r]

    ## Create the binding site for regulator
    prom = promPrefix * tar * "reg" * reg
    
    ## tempactivestates is a matrix where each row correspond to an active state of the promoter: 
    ##    1st column = name of species corresponding to the promoter state, 2nd column = fold change associated with the promoter state
    tempactivestates = [prom*"F" 1]
    boundactive = edg["RegSign"][r] == "1" ## are the bound states of the promoter actives i.e. able to transcribe? No if the regulator is a repressor
    tempallstates = [prom*"F"]

    ## Add the free form of the promoter to the list of species (the bound form is specific to the bound molecule)
    push!(spec, prom*"F")
    ## Add the initial abundance of the free form of the promoter to the list of initial conditions. 
    ## The initial abundance of free promoter set to 1 if it is a binding site on the DNA (exprstep = "TC")
    ## If we are working on translation (exprstep = "TL"), each free binding site has an initial abundance of TCrate*qtlTCrate/(RDrate*qtlRDrate) (RNA initial abundance)
    if exprstep == "TC"
      #push!(inicond, :(1))
      push!(inicond, "1")
    else
      #push!(inicond, :(($(nod["TCrate"][tarid]/nod["RDrate"][tarid])*QTLeffects[$(gcn)]["qtlTCrate"][$(tarid)]/QTLeffects[$(gcn)]["qtlRDrate"][$(tarid)])*InitVar[$(gcn)]["R"][$(tarid)]))
        push!(inicond, """($(nod["TCrate"][tarid]/nod["RDrate"][tarid])*QTLeffects["$(gcn)"]["qtlTCrate"][$(tarid)]/QTLeffects["$(gcn)"]["qtlRDrate"][$(tarid)])*InitVar["$(gcn)"]["R"][$(tarid)]""")
    end    

    for gcnreg in gcnList
      prombound = prom*gcnreg*"B"
      ## add the binding reaction to the list
      push!(spec, prombound) ## add the bound promoter to the list of species
      #push!(inicond, :(0)) ## add its initial abundance to the list of initial conditions. Initial abundance of bound promoter set to 0
      push!(inicond, "0") ## add its initial abundance to the list of initial conditions. Initial abundance of bound promoter set to 0
      if boundactive
        tempactivestates = vcat(tempactivestates, [prombound edg[exprstep*"foldchange"][r]])
      end
      push!(tempallstates, prombound)
      push!(react, reactBioSim([prom*"F", functform[reg]*gcnreg], [prombound])) ## promF + reg -> promregB
      push!(reactnames, "binding"*prom*gcnreg)
      #push!(propens, :($(edg[exprstep*"bindingrate"][r])*QTLeffects[$(gcn)][$("qtl"*exprstep*"regbind")][$(tarid)]*QTLeffects[$(gcnreg)]["qtlactivity"][$(parse(Int64, reg))]))
      push!(propens, """$(edg[exprstep*"bindingrate"][r])*QTLeffects["$(gcn)"]["$("qtl"*exprstep*"regbind")"][$(tarid)]*QTLeffects["$(gcnreg)"]["qtlactivity"][$(parse(Int64, reg))]""")
      ## add the unbinding reaction to the list 
      push!(react, reactBioSim([prombound], [prom*"F", functform[reg]*gcnreg])) ## promregB -> promF + reg
      push!(reactnames, "unbinding"*prom*gcnreg)
      #push!(propens, :($(edg[exprstep*"unbindingrate"][r])))
      push!(propens, """$(edg[exprstep*"unbindingrate"][r])""")
    end

    push!(promActiveStates[tar], tempactivestates)
    push!(promAllStates[tar], tempallstates)

  end

  ## Generate the binding of TL regulatory complexes on RNA-binding sites of targets
  for r in regcompl, gcn in gcnList
    tarid = edg["to"][r]
    tar = string(tarid) * gcn
    compl = edg["from"][r]

    ## Create the binding site for regulator
    prom = promPrefix * tar* "reg" * compl

    ## tempactivestates is a matrix where each lign correspond to an active state of the promoter: 
    ##    1st column = name of species corresponding to the promoter state, 2nd column = fold change associated with the promoter state
    tempactivestates = [prom*"F" 1]
    boundactive = edg["RegSign"][r] == "1" ## are the bound states of the promoter actives i.e. able to transcribe? No if the regulator is a repressor
    tempallstates = [prom*"F"]

    ## Add the free and bound forms of the promoter to the list of species
    push!(spec, prom*"F")
    ## Add the initial abundance of the free form of the promoter to the list of initial conditions. 
    ## The initial abundance of free promoter set to 1 if it is a binding site on the DNA (exprstep = "TC")
    ## If we are working on translation (exprstep = "TL"), each free binding site has an initial abundance of TCrate*qtlTCrate/(RDrate*qtlRDrate) (RNA initial abundance)
    if exprstep == "TC"
      #push!(inicond, :(1))
      push!(inicond, "1")
    else
      #push!(inicond, :(($(nod["TCrate"][tarid]/nod["RDrate"][tarid])*QTLeffects[$(gcn)]["qtlTCrate"][$(tarid)]/QTLeffects[$(gcn)]["qtlRDrate"][$(tarid)])*InitVar[$(gcn)]["R"][$(tarid)]))
      push!(inicond, """($(nod["TCrate"][tarid]/nod["RDrate"][tarid])*QTLeffects["$(gcn)"]["qtlTCrate"][$(tarid)]/QTLeffects["$(gcn)"]["qtlRDrate"][$(tarid)])*InitVar["$(gcn)"]["R"][$(tarid)]""")
    end 

    for t in 1:size(complexvariants)[1]
        complvar = compl*"comp"*join([string(complexes[compl][i])*complexvariants[t, i] for i in 1:complexsize])
        ## Create the binding of the complex on the binding site
        prombound = prom*join(complexvariants[t, :])*"B"
        push!(spec, prombound)
        #push!(inicond, :(0)) ## add its initial abundance to the list of initial conditions. Initial abundance of bound promoter set to 0
        push!(inicond, "0") ## add its initial abundance to the list of initial conditions. Initial abundance of bound promoter set to 0
        if boundactive
            tempactivestates = vcat(tempactivestates, [prombound edg[exprstep*"foldchange"][r]])
        end
        push!(tempallstates, prombound)
        push!(react, reactBioSim([prom*"F", complvar], [prombound])) ## promF + complvar -> promB
        push!(reactnames, "binding"*prom*complvar)
        prop = """$(edg[exprstep*"bindingrate"][r])*QTLeffects["$(gcn)"]["$("qtl"*exprstep*"regbind")"][$(tarid)]"""*join(["""*QTLeffects["$(complexvariants[t, i])"]["qtlactivity"][$(complexes[compl][i])]""" for i in 1:complexsize])
        #push!(propens, parse(prop))
        push!(propens, prop)
        ## add the unbinding of the complex from the binding site
        push!(react, reactBioSim([prombound], [prom*"F", complvar])) ## promB -> promF + complvar
        push!(reactnames, "unbinding"*prom*complvar)
        #push!(propens, :($(edg[exprstep*"unbindingrate"][r])))
        push!(propens, """$(edg[exprstep*"unbindingrate"][r])""")
    end

    push!(promActiveStates[tar], tempactivestates)
    push!(promAllStates[tar], tempallstates)

  end

  return Dict("species" => spec, "initialconditions" => inicond, "reactions" => react, "reactionsnames" => reactnames, "propensities" => propens, "promActiveStates" => promActiveStates, "promAllStates" => promAllStates)

end



function generateReactionList(nod, edgTCRN, edgTLRN, edgRDRN, edgPDRN, edgPTMRN, complexes, complexeskinetics, complexsize, gcnList)
  
    ## Output of the function
    species = []
    initialconditions = []
    reactions = []
    reactionsnames = []
    propensities = []


    ## Specify the active form of each regulator, ie the species that performs the regulation
    functform = Dict(zip(map(string, nod["id"]), nod["ActiveForm"]))

    ## Create all possible allele combinations of the components of a complex
    complexvariants = allposscomb(gcnList, complexsize)

    ## Generate the formation and dissociation reactions of the different regulatory complexes
    for compl in keys(complexes), t in 1:size(complexvariants)[1]
        complvar = compl*"comp"*join([string(complexes[compl][i])*complexvariants[t, i] for i in 1:complexsize])
        push!(species, complvar) ## add the complex form to the list of species
        #push!(initialconditions, :(0)) ## add the initial abundance of the complex form to the list of initial conditions. For complexes, initial abundance set to 0
        push!(initialconditions, "0") ## add the initial abundance of the complex form to the list of initial conditions. For complexes, initial abundance set to 0
        ## Create the reaction of complex formation
        push!(reactions, reactBioSim([functform[string(complexes[compl][i])]*complexvariants[t, i] for i in 1:complexsize], [complvar])) ## sum of complex components -> compl
        push!(reactionsnames, "formation"*complvar) 
        #push!(propensities, :($(complexeskinetics[compl]["formationrate"])))
        push!(propensities, """$(complexeskinetics[compl]["formationrate"])""")
        ## Create the reaction of complex dissociation
        push!(reactions, reactBioSim([complvar], [functform[string(complexes[compl][i])]*complexvariants[t, i] for i in 1:complexsize])) ## compl -> sum of complex components
        push!(reactionsnames, "dissociation"*complvar) 
        #push!(propensities, :($(complexeskinetics[compl]["dissociationrate"])))
        push!(propensities, """$(complexeskinetics[compl]["dissociationrate"])""")
    end



    TLbinding = createbindingpromreactions(edgTLRN, "TL", "RBS", nod, functform, complexes, complexeskinetics, complexsize, complexvariants, gcnList)
    promactiveTL = TLbinding["promActiveStates"]
    promallTL = TLbinding["promAllStates"]

    TCbinding = createbindingpromreactions(edgTCRN, "TC", "Pr", nod, functform, complexes, complexeskinetics, complexsize, complexvariants, gcnList)
    promactiveTC = TCbinding["promActiveStates"]
    promallTC = TCbinding["promAllStates"]

    species  = vcat(species, TCbinding["species"], TLbinding["species"])
    initialconditions  = vcat(initialconditions, TCbinding["initialconditions"], TLbinding["initialconditions"])
    reactions  = vcat(reactions, TCbinding["reactions"], TLbinding["reactions"])
    reactionsnames  = vcat(reactionsnames, TCbinding["reactionsnames"], TLbinding["reactionsnames"])
    propensities  = vcat(propensities, TCbinding["propensities"], TLbinding["propensities"])


    ## CREATE TRANSCRIPTION AND RNA DECAY REACTIONS

    for g in nod["id"], gcn in gcnList
        gname = string(g) * gcn
        TCreg = promactiveTC[gname]
        TLreg = promallTL[gname]

        ## What is the RNA form of the gene? 
        ##    R[gene] if the gene is not controlled at the translation level
        ##    otherwise transcription produces the free (unbound) form of each binding site present on the RNA (1 site per TL regulator)
        if length(TLreg) == 0
            RNAform = ["R"*gname]
            push!(species, "R"*gname)
            #push!(initialconditions, :(($(nod["TCrate"][g]/nod["RDrate"][g])*QTLeffects[$(gcn)]["qtlTCrate"][$(g)]/QTLeffects[$(gcn)]["qtlRDrate"][$(g)])*InitVar[$(gcn)]["R"][$(g)]))
            push!(initialconditions, """($(nod["TCrate"][g]/nod["RDrate"][g])*QTLeffects["$(gcn)"]["qtlTCrate"][$(g)]/QTLeffects["$(gcn)"]["qtlRDrate"][$(g)])*InitVar["$(gcn)"]["R"][$(g)]""")

            transcriptforms = ["R"*gname] ## for RNA decay
        else
            RNAform = [t[1] for t in TLreg] ## The TLreg dictionary is constructed such that each free binding site name is at the position 1
            transcriptforms = combinallpromstates(TLreg) ## for RNA decay
        end

        if length(TCreg) == 0
            ## Create transcription of the gene 
            push!(reactions, reactBioSim([0], RNAform))
            push!(reactionsnames, "transcription"*gname)
            #push!(propensities, :($(nod["TCrate"][g])*QTLeffects[$(gcn)]["qtlTCrate"][$(g)]))
            push!(propensities, """$(nod["TCrate"][g])*QTLeffects["$(gcn)"]["qtlTCrate"][$(g)]""")
        else
            promcomb = combinactivepromstates(TCreg) ## gives all possible active combinations of the different binding site states
            for t in 1:size(promcomb["proms"])[1]
                push!(reactions, reactBioSim(promcomb["proms"][t,:], vcat(promcomb["proms"][t,:], RNAform)))
                push!(reactionsnames, "transcription"*join(promcomb["proms"][t,:]))
                #push!(propensities, :($(nod["TCrate"][g])*QTLeffects[$(gcn)]["qtlTCrate"][$(g)]*$(prod(promcomb["fcs"][t,:]))))
                push!(propensities, """$(nod["TCrate"][g])*QTLeffects["$(gcn)"]["qtlTCrate"][$(g)]*$(prod(promcomb["fcs"][t,:]))""")
            end
        end


        ## Generate RNA decay reactions
        edgreg = find(edgRDRN["to"] .== g)
        edgsingl = edgreg[map(x -> !in('C',x), edgRDRN["from"][edgreg])]
        edgcompl = edgreg[map(x -> in('C',x), edgRDRN["from"][edgreg])]

        for t in 1:size(transcriptforms)[1]
            ## Basal decay rate
            push!(reactions, reactBioSim(transcriptforms[t,:], [0]))
            push!(reactionsnames, "RNAdecay"*join(transcriptforms[t,:]))
            #push!(propensities, :($(nod["RDrate"][g])*QTLeffects[$(gcn)]["qtlRDrate"][$(g)]))
            push!(propensities, """$(nod["RDrate"][g])*QTLeffects["$(gcn)"]["qtlRDrate"][$(g)]""")

            for r in edgsingl, gcnreg in gcnList
                reg = edgRDRN["from"][r]
                push!(reactions, reactBioSim(vcat(transcriptforms[t,:], functform[reg]*gcnreg), [functform[reg]*gcnreg]))
                push!(reactionsnames, "RNAdecay"*join(transcriptforms[t,:])*"reg"*reg*gcnreg)
                #push!(propensities, :($(edgRDRN["RDbindingrate"][r])*QTLeffects[$(gcnreg)]["qtlactivity"][$(parse(Int64, reg))]))
                push!(propensities, """$(edgRDRN["RDbindingrate"][r])*QTLeffects["$(gcnreg)"]["qtlactivity"][$(parse(Int64, reg))]""")
            end

            for r in edgcompl, j in 1:size(complexvariants)[1]
                compl = edgRDRN["from"][r]
                complvar = compl*"comp"*join([string(complexes[compl][i])*complexvariants[j, i] for i in 1:complexsize])
                push!(reactions, reactBioSim(vcat(transcriptforms[t,:], complvar), [complvar]))
                push!(reactionsnames, "RNAdecay"*join(transcriptforms[t,:])*"reg"*complvar)
                propens = """$(edgRDRN["RDbindingrate"][r])"""*join(["*"*"""QTLeffects["$(complexvariants[j, i])"]["qtlactivity"][$(complexes[compl][i])]""" for i in 1:complexsize])
                #push!(propensities, parse(propens))
                push!(propensities, propens)
            end
        end
    end


    ## CREATE TRANSLATION AND PROTEIN DECAY REACTIONS
    for g in nod["id"][nod["coding"] .== "PC"], gcn in gcnList
        gname = string(g) * gcn
        TLreg = promactiveTL[gname]

        push!(species, "P"*gname)
        ## Initial abundance of a protein = TCrate/RDrate * TLrate/PDrate (given the QTL effects as well)
        #push!(initialconditions, :(($(nod["TCrate"][g]*nod["TLrate"][g]/(nod["RDrate"][g]*nod["PDrate"][g]))*QTLeffects[$(gcn)]["qtlTCrate"][$(g)]*QTLeffects[$(gcn)]["qtlTLrate"][$(g)]/(QTLeffects[$(gcn)]["qtlRDrate"][$(g)]*QTLeffects[$(gcn)]["qtlPDrate"][$(g)]))*InitVar[$(gcn)]["P"][$(g)]))
        push!(initialconditions, """($(nod["TCrate"][g]*nod["TLrate"][g]/(nod["RDrate"][g]*nod["PDrate"][g]))*QTLeffects["$(gcn)"]["qtlTCrate"][$(g)]*QTLeffects["$(gcn)"]["qtlTLrate"][$(g)]/(QTLeffects["$(gcn)"]["qtlRDrate"][$(g)]*QTLeffects["$(gcn)"]["qtlPDrate"][$(g)]))*InitVar["$(gcn)"]["P"][$(g)]""")


        if length(TLreg) == 0
            ## Create transltion of the gene
            push!(reactions, reactBioSim(["R"*gname], ["R"*gname, "P"*gname]))
            push!(reactionsnames, "translation"*gname)
            #push!(propensities, :($(nod["TLrate"][g])*QTLeffects[$(gcn)]["qtlTLrate"][$(g)]))
            push!(propensities, """$(nod["TLrate"][g])*QTLeffects["$(gcn)"]["qtlTLrate"][$(g)]""")

        else
            promcomb = combinactivepromstates(TLreg) ## gives all possible active combinations of the different binding site states
            for t in 1:size(promcomb["proms"])[1]
                push!(reactions, reactBioSim(promcomb["proms"][t,:], vcat(promcomb["proms"][t,:], "P"*gname)))
                push!(reactionsnames, "translation"*join(promcomb["proms"][t,:]))
                #push!(propensities, :($(nod["TLrate"][g])*QTLeffects[$(gcn)]["qtlTLrate"][$(g)]*$(prod(promcomb["fcs"][t,:]))))
                push!(propensities, """$(nod["TLrate"][g])*QTLeffects["$(gcn)"]["qtlTLrate"][$(g)]*$(prod(promcomb["fcs"][t,:]))""")
            end
        end

        ## Create basal protein decay rate
        push!(reactions, reactBioSim(["P"*gname], [0]))
        push!(reactionsnames, "proteindecay"*"P"*gname)
        #push!(propensities, :($(nod["PDrate"][g])*QTLeffects[$(gcn)]["qtlPDrate"][$(g)]))
        push!(propensities, """$(nod["PDrate"][g])*QTLeffects["$(gcn)"]["qtlPDrate"][$(g)]""")
        if nod["PTMform"][g] =="1"

            ## Add to the list of species the PTM form of the protein
            push!(species, "Pm"*string(g)*gcn)
            #push!(initialconditions, :(0)) ## initial abundance of modified proteins set to 0
            push!(initialconditions, "0") ## initial abundance of modified proteins set to 0

            push!(reactions, reactBioSim(["Pm"*gname], [0]))
            push!(reactionsnames, "proteindecay"*"Pm"*gname)
            #push!(propensities, :($(nod["PDrate"][g])*QTLeffects[$(gcn)]["qtlPDrate"][$(g)]))
            push!(propensities, """$(nod["PDrate"][g])*QTLeffects["$(gcn)"]["qtlPDrate"][$(g)]""")
        end
    end

      ## CREATE REGULATED PROTEIN DECAY

    regsingl = find( x -> !in('C', x), edgPDRN["from"]) ## identify single-molecule regulators (i.e. not regulatory complexes)
    regcompl = find( x -> in('C', x), edgPDRN["from"]) ## identify regulatory complexes

    for r in regsingl, gcn in gcnList
        tarid = edgPDRN["to"][r]
        tar = string(tarid) * gcn
        reg = edgPDRN["from"][r]

        pform = ["P"*tar]
        if nod["PTMform"][tarid] == "1"
          push!(pform, "Pm"*tar)
        end

        for p in pform, gcnreg in gcnList
            push!(reactions, reactBioSim([p, functform[reg]*gcnreg], [functform[reg]*gcnreg]))
            push!(reactionsnames, "proteindecay"*p*"reg"*reg*gcnreg)
            #push!(propensities, :($(edgPDRN["PDbindingrate"][r])*QTLeffects[$(gcnreg)]["qtlactivity"][$(parse(Int64, reg))]))
            push!(propensities, """$(edgPDRN["PDbindingrate"][r])*QTLeffects["$(gcnreg)"]["qtlactivity"][$(reg)]""")
        end
    end

    for r in regcompl, gcn in gcnList
        tarid = edgPDRN["to"][r]
        tar = string(tarid) * gcn
        compl = edgPDRN["from"][r]

        pform = ["P"*tar]
        if nod["PTMform"][tarid] == "1"
            push!(pform, "Pm"*tar)
        end

        for p in pform, t in 1:size(complexvariants)[1]
            complvariant = compl*"comp"*join([string(complexes[compl][i])*complexvariants[t, i] for i in 1:complexsize])
            push!(reactions, reactBioSim([p, complvariant], [complvariant]))
            push!(reactionsnames, "proteindecay"*p*"reg"*complvariant)
            propens = """$(edgPDRN["PDbindingrate"][r])"""*join(["*"*"""QTLeffects["$(complexvariants[t, i])"]["qtlactivity"][$(complexes[compl][i])]""" for i in 1:complexsize])
            #push!(propensities, parse(propens))
            push!(propensities, propens)
        end
    end


  ## CREATE PROTEIN POST-TRANSLATIONAL MODIFICATION


    regsingl = find( x -> !in('C', x), edgPTMRN["from"]) ## identify single-molecule regulators (i.e. not regulatory complexes)
    regcompl = find( x -> in('C', x), edgPTMRN["from"]) ## identify regulatory complexes

    for r in regsingl, gcn in gcnList
        tarid = edgPTMRN["to"][r]
        tar = string(tarid) * gcn
        reg = edgPTMRN["from"][r]
        ispos = edgPTMRN["RegSign"][r] == "1" ## is the regulator transforming the orginal protein into its modified form (RegSign = "1") or the opposite (RegSign = "-1")

        for gcnreg in gcnList
            if ispos
                push!(reactions, reactBioSim([ "P"*tar, functform[reg]*gcnreg], ["Pm"*tar, functform[reg]*gcnreg]))
                push!(reactionsnames, "PTM"*tar*"reg"*reg*gcnreg)
                #push!(propensities, :($(edgPTMRN["PTMbindingrate"][r])*QTLeffects[$(gcnreg)]["qtlactivity"][$(parse(Int64, reg))]))
                push!(propensities, """$(edgPTMRN["PTMbindingrate"][r])*QTLeffects["$(gcnreg)"]["qtlactivity"][$(reg)]""")
            else
                push!(reactions, reactBioSim([ "Pm"*tar, functform[reg]*gcnreg], ["P"*tar, functform[reg]*gcnreg]))
                push!(reactionsnames, "de-PTM"*tar*"reg"*reg*gcnreg)
                #push!(propensities, :($(edgPTMRN["PTMbindingrate"][r])*QTLeffects[$(gcnreg)]["qtlactivity"][$(parse(Int64, reg))])) 
                push!(propensities, """$(edgPTMRN["PTMbindingrate"][r])*QTLeffects["$(gcnreg)"]["qtlactivity"][$(reg)]""") 
            end
        end
    end

    for r in regcompl, gcn in gcnList
        tarid = edgPTMRN["to"][r]
        tar = string(tarid) * gcn
        compl = edgPTMRN["from"][r]
        ispos = edgPTMRN["RegSign"][r] == "1" ## is the regulator transforming the orginal protein into its modified form (RegSign = "1") or the opposite (RegSign = "-1")


        for t in 1:size(complexvariants)[1]
            complvariant = compl*"comp"*join([string(complexes[compl][i])*complexvariants[t, i] for i in 1:complexsize])
          
            if ispos
                push!(reactions, reactBioSim([ "P"*tar, complvariant], ["Pm"*tar, complvariant]))
                push!(reactionsnames, "PTM"*tar*"reg"*complvariant)
                propens = """$(edgPTMRN["PTMbindingrate"][r])"""*join(["*"*"""QTLeffects["$(complexvariants[t, i])"]["qtlactivity"][$(complexes[compl][i])]""" for i in 1:complexsize])
                #push!(propensities, parse(propens))
                push!(propensities, propens)
            else
                push!(reactions, reactBioSim([ "Pm"*tar, complvariant], ["P"*tar, complvariant]))
                push!(reactionsnames, "de-PTM"*tar*"reg"*complvariant)
                propens = """$(edgPTMRN["PTMbindingrate"][r])"""*join(["*"*"""QTLeffects["$(complexvariants[t, i])"]["qtlactivity"][$(complexes[compl][i])]""" for i in 1:complexsize])
                #push!(propensities, parse(propens))
                push!(propensities,propens)
            end
        end
    end

    return Dict("species" => species, "initialconditions" => initialconditions, "reactions" => reactions, "reactionsnames" => reactionsnames, "propensities" => propensities, "TCproms" => promallTC, "TLproms" => promallTL)
end

function getDictfromKey(mydict, mykey)
    return mydict[mykey]
end




# ------------------------------------------------------------------------------------------------ #
##              FUNCTIONS FOR CREATING THE REACTION LIST AND PROPENSITY LIST                      ##
# ------------------------------------------------------------------------------------------------ #



function res2df(output :: BioSimulator.SimData)
    t, data = get_data(output)
    id2ind  = output.id2ind

    d1, d2, d3 = size(data)

    ids = collect(keys(id2ind))

    df = DataFrame()
    df[:time] = repeat(collect(t), outer=[d3])
    for i = ids
        df[i] = vec(convert(Array, view(data, id2ind[i], :, :)))
    end
    df[:trial] = repeat(collect(1:d3), inner=[d2])
    return df
end


#=
function smallmodel()

    model = BioSimulator.Network("coucou")

    model <= BioSimulator.Species("X", 25)
    model <= BioSimulator.Reaction("transcription", 0.02, "0 --> X")
    model <= BioSimulator.Reaction("decay", 0.00069, "X --> 0")

    result = BioSimulator.simulate(model, algorithm = SSA, time = 1000.0, epochs = 1000)
    resdf = res2df(result)

    return resdf
end
=#


function allequal(x)
    all(y->y==x[1], x)
end

function stochasticsimulation(stochmodel, QTLeffects, InitVar, nod, simtime; modelname = "MySimulation", ntrials = 1, nepochs = -1, simalgorithm = "SSA")

    try
        ## If no nepochs provided, record abundance at each time units of the simulation
        if nepochs == -1
            nepochs = simtime
        end

        ## Convert the String name of the simulator to use into a BioSimulator object of type inheriting from Algorithm
        if in(simalgorithm, ["SSA", "FRM", "NRM", "ODM", "SAL"])
            simalgorithm = eval(parse("BioSimulator."*simalgorithm))
        else
            error("Specified algorithm is not implemented in module BioSimulator")
        end

        model = BioSimulator.Network(modelname)


        ## Add the species in the model, with their initial abundance
        for i in 1:length(stochmodel["species"])
            i0 = replace(stochmodel["initialconditions"][i], "QTLeffects", "$QTLeffects")
            i0 = replace(i0, "InitVar", "$InitVar")
            i0 = eval(parse(i0))
            #println(stochmodel["species"][i]* "\t"*string(i0))
            model <= BioSimulator.Species(stochmodel["species"][i], round(Int, i0))

            if !isa(i0, Number)
                println(stochmodel["initialconditions"][i])
                error("Pb i0")
            end
            if typeof(stochmodel["species"][i]) != String
                println(stochmodel["species"][i])
                error("Pb species")
            end
        end


        ## Add the reactions in the model, with their name and rate
        for i in eachindex(stochmodel["reactions"])
            prop = replace(stochmodel["propensities"][i], "QTLeffects", "$QTLeffects")
            #println(prop)
            prop = eval(parse(prop))
            model <= BioSimulator.Reaction(stochmodel["reactionsnames"][i], prop, stochmodel["reactions"][i])

            if !isa(prop, Number)
                println(stochmodel["propensities"][i])
                error("Pb prop")
            end
            if typeof(stochmodel["reactionsnames"][i]) != String
                println(stochmodel["reactionsnames"][i])
                error("Pb reactionsnames")
            end
            if typeof(stochmodel["reactions"][i]) != String
                println(stochmodel["reactions"][i])
                error("Pb reactions")
            end

        end

        #println("JULIA> Running simulation ...")
        #tic()
        #JSON.print(getfield(Main, :RJuliaSocket), toR("please wait"))    
        result = simulate(model, algorithm = simalgorithm, time = convert(Float64, simtime), epochs = round(Int64, nepochs), trials = convert(Int64, ntrials))
        #toc()
        #println("JULIA> Done.")

        resultdf = res2df(result)

        abundancedf = resultdf[:, [:time, :trial]]

        for g in collect(keys(stochmodel["TCproms"]))

            gid = parse(Int64, replace(g, r"GCN.+$","")) ## gives the gene id
            
            ## Check that for each binding site on the promoter of each gene, at each time the sum of the abundance of all species corresponding to all possible binding site states equals 1 (bc only 1 binding site per gene)
            for i in eachindex(stochmodel["TCproms"][g])
                prabundance = resultdf[:, map(Symbol, stochmodel["TCproms"][g][i])]
                sumprom = [sum(convert(Array, row)) for row in eachrow(prabundance)]
                if any(sumprom .!=1)
                    error("The sum of promoter states for "*stochmodel["TCproms"][g][i][1]*"not equal to 1.")
                end
            end

            if length(stochmodel["TLproms"][g]) > 0
                
                ## Check for each RNA that the different binding sites on the RNA are in equal abundance at each time of the simulation
                rbsabundance = [sum(convert(Array, resultdf[t, map(Symbol, x)]))for t in 1:size(resultdf)[1], x in stochmodel["TLproms"][g]]
                
                if !all(mapslices(allequal, rbsabundance, 2))
                    error("The abundance of the different binding sites on the RNA "*g*"are not equal.")
                end

                ## Add to abundancedf a column corresponding to the abundance of the RNA associated with g
                abundancedf[Symbol("R"*g)] = rbsabundance[:,1]
            else
                abundancedf[Symbol("R"*g)] = resultdf[:, Symbol("R"*g)]
            end


            ## MAYBE to change if we don't make the disctinction between original and modified protein
            if nod["coding"][gid] == "PC"
                abundancedf[Symbol("P"*g)] = resultdf[:, Symbol("P"*g)]

                if nod["PTMform"][gid] == "1"
                    abundancedf[Symbol("Pm"*g)] = resultdf[:, Symbol("Pm"*g)]
                end

            end
        end

        return abundancedf


    catch err
        isa(err, InterruptException) || rethrow(err)
        return nothing
    end

end

#=
function runstochsim(stochmodel, QTLeffects, InitVar, nod, simtime; modelname = "MySimulation", ntrials = 1, nepochs = -1, simalgorithm = "SSA")
    try
        res = stochasticsimulation(stochmodel, QTLeffects, InitVar, nod, simtime, modelname = modelname, ntrials = ntrials, nepochs = nepochs, simalgorithm = simalgorithm)
        return res
    catch err
        isa(err, InterruptException) || rethrow(err)
        return nothing
    end
end
=#


function whatisit(x)
  println(x)
  println(typeof(x))
end



# ------------------------------------------------------------------------------------------------ #
##              SOLUTION TO RUN A SIMULATION AND STOP IT IF RUNNING TIME IS                       ##
##                         LARGER THAN A DEFINED TIME LIMIT                                       ##
##                    CURRENTLY ONLY WORKS FROM THE JULIA TERMINAL                                ##
# ------------------------------------------------------------------------------------------------ #

##   !!!WARNING!!! DO NOT UNCOMMENT IN THE julia_functions.jl FILE OR IT WILL BUG

#=
p = addprocs(1)[1]
@everywhere include("winData/multiomics_networks_simulation/julia_functions.jl")


@everywhere function stochsimtimelimit(p, stochmodel, QTLeffects, InitVar, nod, simtime; modelname = "MySimulation", ntrials = 1, nepochs = -1, simalgorithm = "SSA", time_limit = Inf)

    output = Channel(1) ## create a Channel, that will store the simulation results
    @async put!(output, remotecall_fetch(stochasticsimulation, p, stochmodel, QTLeffects, InitVar, nod, 1)) ## start the simulation on the process 2 of the Julia evaluator

    start=time() ## start the timer

    while !isready(output) & (time() - start < time_limit) ## While the channel is not ready (i.e. computation not done) and the timer is < to the defined value (default no time limit)
      sleep(0.1)
    end

    if !isready(output) ## As soon as the Channel is ready or the time is out, check if the computation is done
      interrupt(2) ## if not interrupt the second process i.e. the computation
      data = "timeout" ## return a timeout message
      #close(output)
    else
      data = take!(output) ## if the computation is done take the result from the Channel
    end

    return data
end


## Loading an example
using JLD
gentil = load("/home/oangelin/Documents/goodmodel.jld")
QTLeffects = gentil["QTLeffects"]
InitVar = gentil["InitVar"]
stochmodel = gentil["stochmodel"]
nod = load("/home/oangelin/Documents/goodmodelnod.jld")["nod"]

## Testing on the example

@time stochsimtimelimit(p, stochmodel, QTLeffects, InitVar, nod, 1, time_limit = 10)
@time stochsimtimelimit(p, stochmodel, QTLeffects, InitVar, nod, 1, time_limit = 70)
@time stochsimtimelimit(p, stochmodel, QTLeffects, InitVar, nod, 1)
=#


