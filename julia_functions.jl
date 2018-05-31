##########################################################################################################################
###                                JULIA FUNCTIONS FOR NETWORK GENERATION                                              ###
##########################################################################################################################

## If not already installed, instsall required packages
if !haskey(Pkg.installed(), "StatsBase") 
	Pkg.add("StatsBase")
end

#= if !haskey(Pkg.installed(), "Combinatorics") 
  Pkg.add("Combinatorics")
end
=#

using StatsBase
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
        edgpos = setdiff(edgpos, compo) ## remove the selected rows (regulators) for possible components (future sampling)
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
        edgneg = setdiff(edgneg, compo) ## remove the selected rows (regulators) for possible components (future sampling)
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

  edg[:,1] = [string(i) for i in edg[:,1]] ## Transform the integer ID into String ID 

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
      push!(inicond, :(1))
    else
      push!(inicond, :(($(nod["TCrate"][tarid]/nod["RDrate"][tarid])*QTLeffects[$(gcn)]["qtlTCrate"][$(tarid)]/QTLeffects[$(gcn)]["qtlRDrate"][$(tarid)])*InitVar[$(gcn)]["R"][$(tarid)]))
    end    

    for gcnreg in gcnList
      prombound = prom*gcnreg*"B"
      ## add the binding reaction to the list
      push!(spec, prombound) ## add the bound promoter to the list of species
      push!(inicond, :(0)) ## add its initial abundance to the list of initial conditions. Initial abundance of bound promoter set to 0
      if boundactive
        tempactivestates = vcat(tempactivestates, [prombound edg[exprstep*"foldchange"][r]])
      end
      push!(tempallstates, prombound)
      push!(react, reactBioSim([prom*"F", functform[reg]*gcnreg], [prombound])) ## promF + reg -> promregB
      push!(reactnames, "binding"*prom*gcnreg)
      push!(propens, :($(edg[exprstep*"bindingrate"][r])*QTLeffects[$(gcn)][$("qtl"*exprstep*"regbind")][$(tarid)]*QTLeffects[$(gcnreg)]["qtlactivity"][$(parse(Int64, reg))]))
      ## add the unbinding reaction to the list 
      push!(react, reactBioSim([prombound], [prom*"F", functform[reg]*gcnreg])) ## promregB -> promF + reg
      push!(reactnames, "unbinding"*prom*gcnreg)
      push!(propens, :($(edg[exprstep*"unbindingrate"][r])))
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
      push!(inicond, :(1))
    else
      push!(inicond, :(($(nod["TCrate"][tarid]/nod["RDrate"][tarid])*QTLeffects[$(gcn)]["qtlTCrate"][$(tarid)]/QTLeffects[$(gcn)]["qtlRDrate"][$(tarid)])*InitVar[$(gcn)]["R"][$(tarid)]))
    end 

    for t in 1:size(complexvariants)[1]
      complvar = compl*"comp"*join([string(complexes[compl][i])*complexvariants[t, i] for i in 1:complexsize])
      ## Create the binding of the complex on the binding site
      prombound = prom*join(complexvariants[t, :])*"B"
      push!(spec, prombound)
      push!(inicond, :(0)) ## add its initial abundance to the list of initial conditions. Initial abundance of bound promoter set to 0
      if boundactive
        tempactivestates = vcat(tempactivestates, [prombound edg[exprstep*"foldchange"][r]])
      end
      push!(tempallstates, prombound)
      push!(react, reactBioSim([prom*"F", complvar], [prombound])) ## promF + complvar -> promB
      push!(reactnames, "binding"*prom*complvar)
      prop = """$(edg[exprstep*"bindingrate"][r])*QTLeffects["$(gcn)"]["$("qtl"*exprstep*"regbind")"][$(tarid)]"""*join(["""*QTLeffects["$(complexvariants[t, i])"]["qtlactivity"][$(complexes[compl][i])]""" for i in 1:complexsize])
      push!(propens, parse(prop))
      ## add the unbinding of the complex from the binding site
      push!(react, reactBioSim([prombound], [prom*"F", complvar])) ## promB -> promF + complvar
      push!(reactnames, "unbinding"*prom*complvar)
      push!(propens, :($(edg[exprstep*"unbindingrate"][r])))
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
        push!(initialconditions, :(0)) ## add the initial abundance of the complex form to the list of initial conditions. For complexes, initial abundance set to 0
        ## Create the reaction of complex formation
        push!(reactions, reactBioSim([functform[string(complexes[compl][i])]*complexvariants[t, i] for i in 1:complexsize], [complvar])) ## sum of complex components -> compl
        push!(reactionsnames, "formation"*complvar) 
        push!(propensities, :($(complexeskinetics[compl]["formationrate"])))
        ## Create the reaction of complex dissociation
        push!(reactions, reactBioSim([complvar], [functform[string(complexes[compl][i])]*complexvariants[t, i] for i in 1:complexsize])) ## compl -> sum of complex components
        push!(reactionsnames, "dissociation"*complvar) 
        push!(propensities, :($(complexeskinetics[compl]["dissociationrate"])))
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
            push!(initialconditions, :(($(nod["TCrate"][g]/nod["RDrate"][g])*QTLeffects[$(gcn)]["qtlTCrate"][$(g)]/QTLeffects[$(gcn)]["qtlRDrate"][$(g)])*InitVar[$(gcn)]["R"][$(g)]))
            transcriptforms = ["R"*gname] ## for RNA decay
        else
            RNAform = [t[1] for t in TLreg] ## The TLreg dictionary is constructed such that each free binding site name is at the position 1
            transcriptforms = combinallpromstates(TLreg) ## for RNA decay
        end

        if length(TCreg) == 0
            ## Create transcription of the gene 
            push!(reactions, reactBioSim([0], RNAform))
            push!(reactionsnames, "transcription"*gname)
            push!(propensities, :($(nod["TCrate"][g])*QTLeffects[$(gcn)]["qtlTCrate"][$(g)]))
        else
            promcomb = combinactivepromstates(TCreg) ## gives all possible active combinations of the different binding site states
            for t in 1:size(promcomb["proms"])[1]
                push!(reactions, reactBioSim(promcomb["proms"][t,:], vcat(promcomb["proms"][t,:], RNAform)))
                push!(reactionsnames, "transcription"*join(promcomb["proms"][t,:]))
                push!(propensities, :($(nod["TCrate"][g])*QTLeffects[$(gcn)]["qtlTCrate"][$(g)]*$(prod(promcomb["fcs"][t,:]))))
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
            push!(propensities, :($(nod["RDrate"][g])*QTLeffects[$(gcn)]["qtlRDrate"][$(g)]))

            for r in edgsingl, gcnreg in gcnList
                reg = edgRDRN["from"][r]
                push!(reactions, reactBioSim(vcat(transcriptforms[t,:], functform[reg]*gcnreg), [functform[reg]*gcnreg]))
                push!(reactionsnames, "RNAdecay"*join(transcriptforms[t,:])*"reg"*reg*gcnreg)
                push!(propensities, :($(edgRDRN["RDbindingrate"][r])*QTLeffects[$(gcnreg)]["qtlactivity"][$(parse(Int64, reg))]))
            end

            for r in edgcompl, j in 1:size(complexvariants)[1]
                compl = edgRDRN["from"][r]
                complvar = compl*"comp"*join([string(complexes[compl][i])*complexvariants[j, i] for i in 1:complexsize])
                push!(reactions, reactBioSim(vcat(transcriptforms[t,:], complvar), [complvar]))
                push!(reactionsnames, "RNAdecay"*join(transcriptforms[t,:])*"reg"*complvar)
                propens = """$(edgRDRN["RDbindingrate"][r])"""*join(["*"*"""QTLeffects["$(complexvariants[j, i])"]["qtlactivity"][$(complexes[compl][i])]""" for i in 1:complexsize])
                push!(propensities, parse(propens))
            end
        end
    end


    ## CREATE TRANSLATION AND PROTEIN DECAY REACTIONS
    for g in nod["id"][nod["coding"] .== "PC"], gcn in gcnList
        gname = string(g) * gcn
        TLreg = promactiveTL[gname]

        push!(species, "P"*gname)
        ## Initial abundance of a protein = TCrate/RDrate * TLrate/PDrate (given the QTL effects as well)
        push!(initialconditions, :(($(nod["TCrate"][g]*nod["TLrate"][g]/(nod["RDrate"][g]*nod["PDrate"][g]))*QTLeffects[$(gcn)]["qtlTCrate"][$(g)]*QTLeffects[$(gcn)]["qtlTLrate"][$(g)]/(QTLeffects[$(gcn)]["qtlRDrate"][$(g)]*QTLeffects[$(gcn)]["qtlPDrate"][$(g)]))*InitVar[$(gcn)]["P"][$(g)]))


        if length(TLreg) == 0
            ## Create transltion of the gene
            push!(reactions, reactBioSim(["R"*gname], ["R"*gname, "P"*gname]))
            push!(reactionsnames, "translation"*gname)
            push!(propensities, :($(nod["TLrate"][g])*QTLeffects[$(gcn)]["qtlTLrate"][$(g)]))

        else
            promcomb = combinactivepromstates(TLreg) ## gives all possible active combinations of the different binding site states
            for t in 1:size(promcomb["proms"])[1]
                push!(reactions, reactBioSim(promcomb["proms"][t,:], vcat(promcomb["proms"][t,:], "P"*gname)))
                push!(reactionsnames, "translation"*join(promcomb["proms"][t,:]))
                push!(propensities, :($(nod["TLrate"][g])*QTLeffects[$(gcn)]["qtlTLrate"][$(g)]*$(prod(promcomb["fcs"][t,:]))))
            end
        end

        ## Create basal protein decay rate
        push!(reactions, reactBioSim(["P"*gname], [0]))
        push!(reactionsnames, "proteindecay"*"P"*gname)
        push!(propensities, :($(nod["PDrate"][g])*QTLeffects[$(gcn)]["qtlPDrate"][$(g)]))
        if nod["PTMform"][g] =="1"

            ## Add to the list of species the PTM form of the protein
            push!(species, "Pm"*string(g)*gcn)
            push!(initialconditions, :(0)) ## initial abundance of modified proteins set to 0

            push!(reactions, reactBioSim(["Pm"*gname], [0]))
            push!(reactionsnames, "proteindecay"*"Pm"*gname)
            push!(propensities, :($(nod["PDrate"][g])*QTLeffects[$(gcn)]["qtlPDrate"][$(g)]))
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
            push!(propensities, :($(edgPDRN["PDbindingrate"][r])*QTLeffects[$(gcnreg)]["qtlactivity"][$(parse(Int64, reg))]))
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
            push!(propensities, parse(propens))
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
                push!(propensities, :($(edgPTM["PTMbindingrate"][r])*QTLeffects[$(gcnreg)]["qtlactivity"][$(parse(Int64, reg))]))
            else
                push!(reactions, reactBioSim([ "Pm"*tar, functform[reg]*gcnreg], ["P"*tar, functform[reg]*gcnreg]))
                push!(reactionsnames, "de-PTM"*tar*"reg"*reg*gcnreg)
                push!(propensities, :($(edgPTM["PTMbindingrate"][r])*QTLeffects[$(gcnreg)]["qtlactivity"][$(parse(Int64, reg))])) 
            end
        end
    end

    for r in regcompl, gcn in gcnList
        tarid = edgPTMRN["to"][r]
        tar = string(tarid) * gcn
        compl = edgPTMRN["from"][r]
        ispos = edgPTMRN["RegSign"][r] == "1" ## is the regulator transforming the orginal protein into its modified form (RegSign = "1") or the opposite (RegSign = "-1")


        for p in pform, t in 1:size(complexvariants)[1]
            complvariant = compl*"comp"*join([string(complexes[compl][i])*complexvariants[t, i] for i in 1:complexsize])
          
            if ispos
                push!(reactions, reactBioSim([ "P"*tar, complvariant], ["Pm"*tar, complvariant]))
                push!(reactionsnames, "PTM"*tar*"reg"*reg*complvariant)
                propens = """$(edgPTM["PTMbindingrate"][r])"""*join(["*"*"""QTLeffects["$(complexvariants[t, i])"]["qtlactivity"][$(complexes[compl][i])]""" for i in 1:complexsize])
                push!(propensities, parse(propens))
            else
                push!(reactions, reactBioSim([ "Pm"*tar, complvariant], ["P"*tar, complvariant]))
                push!(reactionsnames, "de-PTM"*tar*"reg"*reg*complvariant)
                propens = """$(edgPTM["PTMbindingrate"][r])"""*join(["*"*"""QTLeffects["$(complexvariants[t, i])"]["qtlactivity"][$(complexes[compl][i])]""" for i in 1:complexsize])
                push!(propensities, parse(propens))
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


using BioSimulator




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


function smallmodel()


    model = BioSimulator.Network("coucou")

    model <= BioSimulator.Species("X", 25)
    model <= BioSimulator.Reaction("transcription", 0.02, "0 --> X")
    model <= BioSimulator.Reaction("decay", 0.00069, "X --> 0")

    result = BioSimulator.simulate(model, algorithm = SSA, time = 1000.0, epochs = 1000)
    resdf = res2df(result)

    return resdf
end

function allequal(x)
    all(y->y==x[1], x)
end

function stochasticsimulation(stochmodel, QTLeffectslocal, InitVarlocal, nod, simtime; modelname = "MySimulation", ntrials = 1, nepochs = -1, simalgorithm = "SSA")

    ## Only way to use eval inside the function
    global QTLeffects = QTLeffectslocal
    global InitVar = InitVarlocal

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
        i0 = round(Int, eval(stochmodel["initialconditions"][i]))
        #println(stochmodel["species"][i]* "\t"*string(i0))
        model <= BioSimulator.Species(stochmodel["species"][i], i0)
    end


    ## Add the reactions in the model, with their name and rate
    for i in eachindex(stochmodel["reactions"])
        model <= BioSimulator.Reaction(stochmodel["reactionsnames"][i], eval(stochmodel["propensities"][i]), stochmodel["reactions"][i])
    end

    println("Running simulation ...")
    tic();result = simulate(model, algorithm = simalgorithm, time = convert(Float64, simtime), epochs = round(Int64, nepochs), trials = convert(Int64, ntrials));toc()

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
end

#=

## system of 20 nodes
workspace()
include("julia_functions.jl")
nod = Dict{String,Any}(Pair{String,Any}("TargetReaction", String["TL", "TL", "TL", "TC", "TC", "RD", "TC", "TL", "TC", "TC", "TL", "TC", "TC", "TC", "TL", "RD", "RD", "MR", "TL", "TL", "RD", "RD", "RD", "TC", "TL", "TL", "TC", "TC", "TL", "TC", "TL", "TL", "TC", "TL", "TC", "TL", "TC", "TC", "RD", "TL", "TL", "TC", "PTM", "RD", "TC", "TL", "TL", "PD", "RD", "RD"]),Pair{String,Any}("TCrate", [0.0383058, 0.0354193, 0.069086, 0.0617387, 0.0369238, 0.0624596, 0.0773182, 0.0833146, 0.0502687, 0.0313804, 0.0527055, 0.0432903, 0.0856088, 0.0585721, 0.0452048, 0.0356704, 0.0208695, 0.0670213, 0.0835078, 0.0779509, 0.0776592, 0.0541112, 0.0537629, 0.050601, 0.0935788, 0.0785906, 0.0299127, 0.0495242, 0.0799262, 0.0830136, 0.0184689, 0.0875529, 0.024178, 0.0702895, 0.0728294, 0.0295397, 0.0367362, 0.0986029, 0.011211, 0.0592016, 0.0767439, 0.0478291, 0.04266, 0.0607304, 0.0426131, 0.0921907, 0.0305072, 0.0778754, 0.0362273, 0.0529586]),Pair{String,Any}("RDrate", [0.000413907, 0.000396511, 0.000516262, 0.000287356, 0.000577367, 0.000364166, 0.000644745, 0.000420521, 0.000540833, 0.00100604, 0.00361011, 0.000694444, 0.000315856, 0.000350877, 0.00041425, 0.003663, 0.000322789, 0.000588582, 0.01, 0.00040016, 0.000402576, 0.000288517, 0.00215517, 0.000422654, 0.000745156, 0.0108696, 0.000290951, 0.00157729, 0.00030003, 0.000286287, 0.000786782, 0.000357782, 0.000281611, 0.00038432, 0.000336587, 0.000567215, 0.000361141, 0.000417188, 0.000470367, 0.000331675, 0.000787402, 0.00102354, 0.000865801, 0.000551268, 0.000450248, 0.000684463, 0.00077821, 0.00137174, 0.0025974, 0.000648508]),Pair{String,Any}("id", [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50]),Pair{String,Any}("ActiveForm", String["Pm1", "R2", "R3", "Pm4", "Pm5", "P6", "Pm7", "P8", "R9", "P10", "Pm11", "P12", "Pm13", "Pm14", "Pm15", "P16", "Pm17", "P18", "P19", "R20", "R21", "R22", "Pm23", "P24", "R25", "P26", "R27", "Pm28", "Pm29", "R30", "Pm31", "P32", "Pm33", "R34", "R35", "R36", "Pm37", "R38", "R39", "Pm40", "R41", "Pm42", "Pm43", "Pm44", "P45", "Pm46", "Pm47", "P48", "P49", "R50"]),Pair{String,Any}("coding", String["PC", "NC", "NC", "PC", "PC", "PC", "PC", "PC", "NC", "PC", "PC", "PC", "PC", "PC", "PC", "PC", "PC", "PC", "PC", "NC", "NC", "NC", "PC", "PC", "NC", "PC", "NC", "PC", "PC", "NC", "PC", "PC", "PC", "NC", "NC", "NC", "PC", "NC", "NC", "PC", "NC", "PC", "PC", "PC", "PC", "PC", "PC", "PC", "PC", "NC"]),Pair{String,Any}("TLrate", [4.26806, 0.0, 0.0, 0.620594, 2.80055, 4.59256, 1.04708, 2.53738, 0.0, 1.42765, 0.910486, 0.591064, 4.53948, 2.94196, 1.39797, 3.28256, 0.998617, 4.0215, 2.9922, 0.0, 0.0, 0.0, 1.75604, 4.17647, 0.0, 4.12805, 0.0, 2.73211, 3.54725, 0.0, 4.38978, 3.06126, 4.33202, 0.0, 0.0, 0.0, 3.14571, 0.0, 0.0, 4.79517, 0.0, 0.658, 2.52966, 2.13485, 1.68262, 3.07637, 3.53223, 3.55141, 3.73643, 0.0]),Pair{String,Any}("nameid", String["G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "G_11", "G_12", "G_13", "G_14", "G_15", "G_16", "G_17", "G_18", "G_19", "G_20", "G_21", "G_22", "G_23", "G_24", "G_25", "G_26", "G_27", "G_28", "G_29", "G_30", "G_31", "G_32", "G_33", "G_34", "G_35", "G_36", "G_37", "G_38", "G_39", "G_40", "G_41", "G_42", "G_43", "G_44", "G_45", "G_46", "G_47", "G_48", "G_49", "G_50"]),Pair{String,Any}("PTMform", String["1", "0", "0", "1", "1", "0", "1", "0", "0", "0", "1", "0", "1", "1", "1", "0", "1", "0", "0", "0", "0", "0", "1", "0", "0", "0", "0", "1", "1", "0", "1", "0", "1", "0", "0", "0", "1", "0", "0", "1", "0", "1", "1", "1", "0", "1", "1", "0", "0", "0"]),Pair{String,Any}("PDrate", [8.76424e-5, 0.0, 0.0, 0.000133529, 0.000127551, 7.11136e-5, 9.59049e-5, 8.21085e-5, 0.0, 0.000140568, 0.000113109, 9.72479e-5, 8.62069e-5, 0.000105809, 9.37207e-5, 8.05737e-5, 0.00017304, 8.85426e-5, 7.62195e-5, 0.0, 0.0, 0.0, 0.000134517, 0.000146327, 0.0, 8.09323e-5, 0.0, 0.000121197, 0.00014453, 0.0, 7.36485e-5, 0.000121374, 8.54555e-5, 0.0, 0.0, 0.0, 7.60572e-5, 0.0, 0.0, 0.000140115, 0.0, 8.66551e-5, 0.000111532, 0.000106849, 7.32172e-5, 0.000160179, 0.000104232, 0.000107504, 7.64117e-5, 0.0]))

edgTCRN = Dict{String,Any}(Pair{String,Any}("TCbindingrate", [0.00321405, 0.00473435, 0.00298111, 0.00997359, 0.00291898, 0.00353278, 0.00401311, 0.00918377, 0.00882659, 0.00991747, 0.00986185, 0.00773114, 0.00469209, 0.00833767, 0.00769984, 0.00308694, 0.0092009, 0.00546315, 0.00678325, 0.00528333, 0.00498977, 0.00628615, 0.00391206, 0.00753273, 0.00601903, 0.0036557, 0.0084864, 0.00780943, 0.00222953, 0.00341767, 0.00862045, 0.00686601, 0.00132944, 0.00947477, 0.0052604, 0.00279102, 0.00649986, 0.00375578, 0.00460729, 0.00232223, 0.00508685, 0.00581595, 0.00475835]),Pair{String,Any}("TCunbindingrate", [0.00953332, 0.00153581, 0.00438116, 0.00426973, 0.00776465, 0.00440608, 0.00537871, 0.00877117, 0.00303946, 0.00433087, 0.00435096, 0.00107793, 0.00760973, 0.00725071, 0.00517657, 0.00694052, 0.00655797, 0.00328065, 0.00940147, 0.0037251, 0.00747392, 0.00848613, 0.00760885, 0.0058852, 0.00403472, 0.00363345, 0.00683965, 0.00949643, 0.00515598, 0.00686772, 0.00853534, 0.00479064, 0.0057276, 0.00617271, 0.0040537, 0.0011888, 0.00258087, 0.00234804, 0.00195996, 0.00905108, 0.00372552, 0.00246426, 0.00836741]),Pair{String,Any}("TargetReaction", String["TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC"]),Pair{String,Any}("RegSign", String["1", "1", "1", "1", "1", "1", "1", "-1", "-1", "-1", "1", "1", "1", "-1", "1", "1", "-1", "-1", "-1", "-1", "1", "1", "1", "-1", "1", "1", "1", "1", "-1", "1", "-1", "-1", "-1", "-1", "1", "-1", "-1", "-1", "1", "1", "-1", "-1", "1"]),Pair{String,Any}("to", [46, 33, 4, 31, 26, 37, 44, 14, 1, 23, 45, 28, 32, 18, 17, 8, 6, 24, 16, 26, 47, 8, 19, 24, 43, 46, 14, 7, 7, 9, 22, 23, 27, 31, 33, 36, 40, 43, 44, 48, 49, 24, 40]),Pair{String,Any}("from", String["9", "9", "9", "9", "9", "9", "9", "9", "9", "27", "27", "27", "27", "27", "27", "27", "27", "30", "30", "30", "30", "30", "30", "35", "35", "35", "38", "5", "42", "10", "33", "5", "14", "12", "4", "4", "37", "13", "4", "28", "7", "CTC1", "CTC2"]),Pair{String,Any}("TCfoldchange", [20.0, 2.0, 29.0, 10.0, 24.0, 20.0, 13.0, 0.0, 0.0, 0.0, 5.0, 13.0, 12.0, 0.0, 21.0, 9.0, 0.0, 0.0, 0.0, 0.0, 23.0, 8.0, 2.0, 0.0, 3.0, 11.0, 6.0, 24.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 14.0, 0.0, 0.0, 0.0, 23.0, 19.0, 0.0, 0.0, 4.0]))
edgTLRN = Dict{String,Any}(Pair{String,Any}("TLbindingrate", [0.0075634, 0.00899251, 0.00914503, 0.00456027, 0.00866226, 0.00804967, 0.00706122, 0.00706322, 0.00907255, 0.00944675, 0.00577241, 0.00354473, 0.00517906, 0.00873183, 0.00656527, 0.00351437, 0.00115263, 0.00213262, 0.00395243, 0.00406976, 0.0048827, 0.00555369, 0.00323992, 0.00846262, 0.00549503, 0.00416131, 0.00727439, 0.00981221, 0.00233314, 0.00859894, 0.00359648, 0.00142693, 0.00249592, 0.00712456, 0.00172458, 0.00782158, 0.00392808]),Pair{String,Any}("TargetReaction", String["TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL", "TL"]),Pair{String,Any}("TLunbindingrate", [0.00487906, 0.00538211, 0.00634111, 0.00801036, 0.00602788, 0.00736443, 0.00988095, 0.0019663, 0.00963572, 0.00709988, 0.00963074, 0.00309585, 0.00739528, 0.00336303, 0.00119043, 0.008372, 0.00958686, 0.00171852, 0.00974889, 0.00771954, 0.00842176, 0.00285371, 0.00383953, 0.00978057, 0.00262683, 0.00275653, 0.00125189, 0.00166714, 0.00178441, 0.00221861, 0.00592483, 0.00340203, 0.00730519, 0.00985448, 0.00845491, 0.0023016, 0.00690244]),Pair{String,Any}("TLfoldchange", [0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 29.0, 0.0, 28.0, 0.0, 20.0, 16.0, 11.0, 25.0, 0.0, 25.0, 0.0, 0.0, 0.0, 19.0]),Pair{String,Any}("RegSign", String["-1", "-1", "1", "-1", "-1", "-1", "-1", "-1", "-1", "-1", "-1", "-1", "-1", "-1", "-1", "-1", "-1", "-1", "-1", "-1", "-1", "1", "-1", "1", "-1", "1", "-1", "1", "1", "1", "1", "-1", "1", "-1", "-1", "-1", "1"]),Pair{String,Any}("to", [49, 45, 26, 23, 33, 19, 6, 16, 23, 49, 33, 18, 16, 32, 17, 29, 31, 17, 12, 28, 40, 7, 7, 10, 11, 12, 15, 15, 17, 28, 29, 29, 37, 40, 44, 47, 10]),Pair{String,Any}("from", String["2", "2", "2", "2", "2", "2", "3", "3", "3", "3", "3", "3", "20", "20", "20", "25", "25", "25", "34", "36", "41", "1", "19", "32", "26", "31", "11", "46", "40", "15", "1", "29", "29", "8", "1", "47", "CTL1"]))
edgRDRN = Dict{String,Any}(Pair{String,Any}("TargetReaction", String["RD", "RD", "RD", "RD", "RD", "RD", "RD", "RD", "RD", "RD", "RD", "RD", "RD", "RD", "RD", "RD", "RD", "RD", "RD", "RD", "RD", "RD", "RD", "RD", "RD", "RD", "RD"]),Pair{String,Any}("RDbindingrate", [0.00697804, 0.00774964, 0.00449371, 0.00691126, 0.00959166, 0.00921492, 0.00198896, 0.00931222, 0.00596598, 0.00868057, 0.00212075, 0.00382286, 0.00546214, 0.00956411, 0.00476764, 0.00695973, 0.00762754, 0.00366035, 0.00382369, 0.00257641, 0.0045423, 0.00148735, 0.00206899, 0.00265115, 0.00646582, 0.00915574, 0.00710227]),Pair{String,Any}("RDunbindingrate", [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00255529, 0.0, 0.0, 0.0, 0.0, 0.00489535, 0.00824349, 0.0, 0.0, 0.00915262, 0.00921839, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00705485, 0.0]),Pair{String,Any}("RegSign", String["1", "1", "1", "1", "1", "1", "1", "1", "-1", "1", "1", "1", "1", "-1", "-1", "1", "1", "-1", "-1", "1", "1", "1", "1", "1", "1", "-1", "1"]),Pair{String,Any}("to", [4, 28, 22, 10, 46, 11, 29, 5, 7, 16, 17, 18, 19, 21, 25, 26, 27, 29, 30, 30, 33, 39, 46, 49, 50, 24, 42]),Pair{String,Any}("from", String["21", "21", "21", "22", "22", "39", "50", "16", "6", "16", "6", "16", "6", "17", "6", "6", "6", "44", "6", "49", "17", "23", "17", "49", "23", "CRD1", "CRD2"]))
edgPDRN = Dict{String,Any}(Pair{String,Any}("TargetReaction", String["PD", "PD", "PD"]),Pair{String,Any}("RegBy", String["PCreg", "PCreg", "PCreg"]),Pair{String,Any}("RegSign", String["-1", "1", "1"]),Pair{String,Any}("to", [31, 14, 48]),Pair{String,Any}("from", String["48", "48", "48"]),Pair{String,Any}("PDunbindingrate", [0.00570818, 0.0, 0.0]),Pair{String,Any}("PDbindingrate", [0.00289001, 0.00727593, 0.00850916]))
edgPTMRN = Dict("from" => [], "to" => [], "TargetReaction" => [], "RegSign" => [], "PTMbindingrate" => [])

complexes = Dict("CRD2"=>[6, 23],"CRD1"=>[6, 16],"CTC1"=>[24, 5],"CTC2"=>[24, 45],"CTL1"=>[8, 47])
complexeskinetics = Dict("CRD2"=>Dict("formationrate"=>0.00604369,"dissociationrate"=>0.00756089),"CRD1"=>Dict("formationrate"=>0.00493894,"dissociationrate"=>0.00610697),"CTC1"=>Dict("formationrate"=>0.00410344,"dissociationrate"=>0.00283723),"CTC2"=>Dict("formationrate"=>0.00346069,"dissociationrate"=>0.00521991),"CTL1"=>Dict("formationrate"=>0.00449745,"dissociationrate"=>0.00485848))
complexsize = 2

xploidy = 4
gcnList = ["GCN"*string(i) for i in 1:xploidy]

tic(); test = generateReactionList(nod, edgTCRN, edgTLRN, edgRDRN, edgPDRN, edgPTMRN, complexes, complexeskinetics, complexsize, gcnList); toc()
=#



#=
workspace()
using StatsBase
N= 20
nod = Dict("id" => collect(1:N), "coding" => sample(["PC" "NC"], N), "PTMform" => sample(["0" "1"], N), "ActiveForm" => fill("P", N) .*map(string, collect(1:N)),
 "TCrate" => rand(N), "TLrate" => rand(N), "RDrate" => rand(N), "PDrate" => rand(N))
functform = Dict(zip(map(string, nod["id"]), nod["ActiveForm"]))

E = 60
tempCr = sample(collect(1:N), E)
edgTCRN = Dict("from" => map(string, sample(collect(1:N), E)), "to" => sort(tempCr), "TargetReaction" => fill("TC", E), "RegSign" => sample(["1", "-1"], E), "TCbindingrate" => rand(E), "TCunbindingrate" => rand(E), "TCfoldchange" => sample(1:10, E))

E = 40
tempCr = sample(collect(1:N), E)
edgTLRN = Dict("from" => map(string, sample(collect(1:N), E)), "to" => sort(tempCr), "TargetReaction" => fill("TL", E), "RegSign" => sample(["1", "-1"], E), "TLbindingrate" => rand(E), "TLunbindingrate" => rand(E), "TLfoldchange" => sample(1:10, E))

E = 20
tempCr = sample(collect(1:N), E)
edgPTMRN = Dict("from" => map(string, sample(collect(1:N), E)), "to" => sort(tempCr), "TargetReaction" => fill("TL", E), "RegSign" => sample(["1", "-1"], E), "PTMbindingrate" => rand(E))


#----
edgTCRN = Dict("from" => [], "to" => [], "TargetReaction" => [], "RegSign" => [], "TCbindingrate" => [], "TCunbindingrate" => [], "TCfoldchange" => [])
edgTLRN = Dict("from" => [], "to" => [], "TargetReaction" => [], "RegSign" => [], "TLbindingrate" => [], "TLunbindingrate" => [], "TLfoldchange" => [])
#----

complexes = Dict()

xploidy = 4
gcnList = ["GCN"*string(i) for i in 1:xploidy]



# --------------------------------------------

using BioSimulator
include("winData/multiomics_networks_simulation/julia_functions.jl")

stochmodel = Dict("reactionsnames"=>Any["transcription1GCN1", "RNAdecayR1GCN1", "transcription1GCN2", "RNAdecayR1GCN2", "transcription1GCN3", "RNAdecayR1GCN3", "transcription1GCN4", "RNAdecayR1GCN4", "transcription2GCN1", "RNAdecayR2GCN1", "transcription2GCN2", "RNAdecayR2GCN2", "transcription2GCN3", "RNAdecayR2GCN3", "transcription2GCN4", "RNAdecayR2GCN4", "transcription3GCN1", "RNAdecayR3GCN1", "transcription3GCN2", "RNAdecayR3GCN2", "transcription3GCN3", "RNAdecayR3GCN3", "transcription3GCN4", "RNAdecayR3GCN4", "transcription4GCN1", "RNAdecayR4GCN1", "transcription4GCN2", "RNAdecayR4GCN2", "transcription4GCN3", "RNAdecayR4GCN3", "transcription4GCN4", "RNAdecayR4GCN4", "transcription5GCN1", "RNAdecayR5GCN1", "transcription5GCN2", "RNAdecayR5GCN2", "transcription5GCN3", "RNAdecayR5GCN3", "transcription5GCN4", "RNAdecayR5GCN4", "transcription6GCN1", "RNAdecayR6GCN1", "transcription6GCN2", "RNAdecayR6GCN2", "transcription6GCN3", "RNAdecayR6GCN3", "transcription6GCN4", "RNAdecayR6GCN4", "transcription7GCN1", "RNAdecayR7GCN1", "transcription7GCN2", "RNAdecayR7GCN2", "transcription7GCN3", "RNAdecayR7GCN3", "transcription7GCN4", "RNAdecayR7GCN4", "transcription8GCN1", "RNAdecayR8GCN1", "transcription8GCN2", "RNAdecayR8GCN2", "transcription8GCN3", "RNAdecayR8GCN3", "transcription8GCN4", "RNAdecayR8GCN4", "transcription9GCN1", "RNAdecayR9GCN1", "transcription9GCN2", "RNAdecayR9GCN2", "transcription9GCN3", "RNAdecayR9GCN3", "transcription9GCN4", "RNAdecayR9GCN4", "transcription10GCN1", "RNAdecayR10GCN1", "transcription10GCN2", "RNAdecayR10GCN2", "transcription10GCN3", "RNAdecayR10GCN3", "transcription10GCN4", "RNAdecayR10GCN4", "translation1GCN1", "proteindecayP1GCN1", "translation1GCN2", "proteindecayP1GCN2", "translation1GCN3", "proteindecayP1GCN3", "translation1GCN4", "proteindecayP1GCN4", "translation2GCN1", "proteindecayP2GCN1", "translation2GCN2", "proteindecayP2GCN2", "translation2GCN3", "proteindecayP2GCN3", "translation2GCN4", "proteindecayP2GCN4", "translation3GCN1", "proteindecayP3GCN1", "translation3GCN2", "proteindecayP3GCN2", "translation3GCN3", "proteindecayP3GCN3", "translation3GCN4", "proteindecayP3GCN4", "translation4GCN1", "proteindecayP4GCN1", "translation4GCN2", "proteindecayP4GCN2", "translation4GCN3", "proteindecayP4GCN3", "translation4GCN4", "proteindecayP4GCN4", "translation5GCN1", "proteindecayP5GCN1", "translation5GCN2", "proteindecayP5GCN2", "translation5GCN3", "proteindecayP5GCN3", "translation5GCN4", "proteindecayP5GCN4", "translation7GCN1", "proteindecayP7GCN1", "translation7GCN2", "proteindecayP7GCN2", "translation7GCN3", "proteindecayP7GCN3", "translation7GCN4", "proteindecayP7GCN4"],"propensities"=>Any[:(0.020032 * ((QTLeffects["GCN1"])["qtlTCrate"])[1]), :(0.000694444 * ((QTLeffects["GCN1"])["qtlRDrate"])[1]), :(0.020032 * ((QTLeffects["GCN2"])["qtlTCrate"])[1]), :(0.000694444 * ((QTLeffects["GCN2"])["qtlRDrate"])[1]), :(0.020032 * ((QTLeffects["GCN3"])["qtlTCrate"])[1]), :(0.000694444 * ((QTLeffects["GCN3"])["qtlRDrate"])[1]), :(0.020032 * ((QTLeffects["GCN4"])["qtlTCrate"])[1]), :(0.000694444 * ((QTLeffects["GCN4"])["qtlRDrate"])[1]), :(0.0435213 * ((QTLeffects["GCN1"])["qtlTCrate"])[2]), :(0.000534759 * ((QTLeffects["GCN1"])["qtlRDrate"])[2]), :(0.0435213 * ((QTLeffects["GCN2"])["qtlTCrate"])[2]), :(0.000534759 * ((QTLeffects["GCN2"])["qtlRDrate"])[2]), :(0.0435213 * ((QTLeffects["GCN3"])["qtlTCrate"])[2]), :(0.000534759 * ((QTLeffects["GCN3"])["qtlRDrate"])[2]), :(0.0435213 * ((QTLeffects["GCN4"])["qtlTCrate"])[2]), :(0.000534759 * ((QTLeffects["GCN4"])["qtlRDrate"])[2]), :(0.0766692 * ((QTLeffects["GCN1"])["qtlTCrate"])[3]), :(0.00310559 * ((QTLeffects["GCN1"])["qtlRDrate"])[3]), :(0.0766692 * ((QTLeffects["GCN2"])["qtlTCrate"])[3]), :(0.00310559 * ((QTLeffects["GCN2"])["qtlRDrate"])[3]), :(0.0766692 * ((QTLeffects["GCN3"])["qtlTCrate"])[3]), :(0.00310559 * ((QTLeffects["GCN3"])["qtlRDrate"])[3]), :(0.0766692 * ((QTLeffects["GCN4"])["qtlTCrate"])[3]), :(0.00310559 * ((QTLeffects["GCN4"])["qtlRDrate"])[3]), :(0.0897525 * ((QTLeffects["GCN1"])["qtlTCrate"])[4]), :(0.000286451 * ((QTLeffects["GCN1"])["qtlRDrate"])[4]), :(0.0897525 * ((QTLeffects["GCN2"])["qtlTCrate"])[4]), :(0.000286451 * ((QTLeffects["GCN2"])["qtlRDrate"])[4]), :(0.0897525 * ((QTLeffects["GCN3"])["qtlTCrate"])[4]), :(0.000286451 * ((QTLeffects["GCN3"])["qtlRDrate"])[4]), :(0.0897525 * ((QTLeffects["GCN4"])["qtlTCrate"])[4]), :(0.000286451 * ((QTLeffects["GCN4"])["qtlRDrate"])[4]), :(0.071785 * ((QTLeffects["GCN1"])["qtlTCrate"])[5]), :(0.00247525 * ((QTLeffects["GCN1"])["qtlRDrate"])[5]), :(0.071785 * ((QTLeffects["GCN2"])["qtlTCrate"])[5]), :(0.00247525 * ((QTLeffects["GCN2"])["qtlRDrate"])[5]), :(0.071785 * ((QTLeffects["GCN3"])["qtlTCrate"])[5]), :(0.00247525 * ((QTLeffects["GCN3"])["qtlRDrate"])[5]), :(0.071785 * ((QTLeffects["GCN4"])["qtlTCrate"])[5]), :(0.00247525 * ((QTLeffects["GCN4"])["qtlRDrate"])[5]), :(0.0514678 * ((QTLeffects["GCN1"])["qtlTCrate"])[6]), :(0.000312891 * ((QTLeffects["GCN1"])["qtlRDrate"])[6]), :(0.0514678 * ((QTLeffects["GCN2"])["qtlTCrate"])[6]), :(0.000312891 * ((QTLeffects["GCN2"])["qtlRDrate"])[6]), :(0.0514678 * ((QTLeffects["GCN3"])["qtlTCrate"])[6]), :(0.000312891 * ((QTLeffects["GCN3"])["qtlRDrate"])[6]), :(0.0514678 * ((QTLeffects["GCN4"])["qtlTCrate"])[6]), :(0.000312891 * ((QTLeffects["GCN4"])["qtlRDrate"])[6]), :(0.0220153 * ((QTLeffects["GCN1"])["qtlTCrate"])[7]), :(0.000307692 * ((QTLeffects["GCN1"])["qtlRDrate"])[7]), :(0.0220153 * ((QTLeffects["GCN2"])["qtlTCrate"])[7]), :(0.000307692 * ((QTLeffects["GCN2"])["qtlRDrate"])[7]), :(0.0220153 * ((QTLeffects["GCN3"])["qtlTCrate"])[7]), :(0.000307692 * ((QTLeffects["GCN3"])["qtlRDrate"])[7]), :(0.0220153 * ((QTLeffects["GCN4"])["qtlTCrate"])[7]), :(0.000307692 * ((QTLeffects["GCN4"])["qtlRDrate"])[7]), :(0.0721572 * ((QTLeffects["GCN1"])["qtlTCrate"])[8]), :(0.000376506 * ((QTLeffects["GCN1"])["qtlRDrate"])[8]), :(0.0721572 * ((QTLeffects["GCN2"])["qtlTCrate"])[8]), :(0.000376506 * ((QTLeffects["GCN2"])["qtlRDrate"])[8]), :(0.0721572 * ((QTLeffects["GCN3"])["qtlTCrate"])[8]), :(0.000376506 * ((QTLeffects["GCN3"])["qtlRDrate"])[8]), :(0.0721572 * ((QTLeffects["GCN4"])["qtlTCrate"])[8]), :(0.000376506 * ((QTLeffects["GCN4"])["qtlRDrate"])[8]), :(0.0658958 * ((QTLeffects["GCN1"])["qtlTCrate"])[9]), :(0.000824402 * ((QTLeffects["GCN1"])["qtlRDrate"])[9]), :(0.0658958 * ((QTLeffects["GCN2"])["qtlTCrate"])[9]), :(0.000824402 * ((QTLeffects["GCN2"])["qtlRDrate"])[9]), :(0.0658958 * ((QTLeffects["GCN3"])["qtlTCrate"])[9]), :(0.000824402 * ((QTLeffects["GCN3"])["qtlRDrate"])[9]), :(0.0658958 * ((QTLeffects["GCN4"])["qtlTCrate"])[9]), :(0.000824402 * ((QTLeffects["GCN4"])["qtlRDrate"])[9]), :(0.0930883 * ((QTLeffects["GCN1"])["qtlTCrate"])[10]), :(0.000316957 * ((QTLeffects["GCN1"])["qtlRDrate"])[10]), :(0.0930883 * ((QTLeffects["GCN2"])["qtlTCrate"])[10]), :(0.000316957 * ((QTLeffects["GCN2"])["qtlRDrate"])[10]), :(0.0930883 * ((QTLeffects["GCN3"])["qtlTCrate"])[10]), :(0.000316957 * ((QTLeffects["GCN3"])["qtlRDrate"])[10]), :(0.0930883 * ((QTLeffects["GCN4"])["qtlTCrate"])[10]), :(0.000316957 * ((QTLeffects["GCN4"])["qtlRDrate"])[10]), :(3.34576 * ((QTLeffects["GCN1"])["qtlTLrate"])[1]), :(7.62369e-5 * ((QTLeffects["GCN1"])["qtlPDrate"])[1]), :(3.34576 * ((QTLeffects["GCN2"])["qtlTLrate"])[1]), :(7.62369e-5 * ((QTLeffects["GCN2"])["qtlPDrate"])[1]), :(3.34576 * ((QTLeffects["GCN3"])["qtlTLrate"])[1]), :(7.62369e-5 * ((QTLeffects["GCN3"])["qtlPDrate"])[1]), :(3.34576 * ((QTLeffects["GCN4"])["qtlTLrate"])[1]), :(7.62369e-5 * ((QTLeffects["GCN4"])["qtlPDrate"])[1]), :(3.00529 * ((QTLeffects["GCN1"])["qtlTLrate"])[2]), :(0.00012364 * ((QTLeffects["GCN1"])["qtlPDrate"])[2]), :(3.00529 * ((QTLeffects["GCN2"])["qtlTLrate"])[2]), :(0.00012364 * ((QTLeffects["GCN2"])["qtlPDrate"])[2]), :(3.00529 * ((QTLeffects["GCN3"])["qtlTLrate"])[2]), :(0.00012364 * ((QTLeffects["GCN3"])["qtlPDrate"])[2]), :(3.00529 * ((QTLeffects["GCN4"])["qtlTLrate"])[2]), :(0.00012364 * ((QTLeffects["GCN4"])["qtlPDrate"])[2]), :(1.15585 * ((QTLeffects["GCN1"])["qtlTLrate"])[3]), :(0.00010269 * ((QTLeffects["GCN1"])["qtlPDrate"])[3]), :(1.15585 * ((QTLeffects["GCN2"])["qtlTLrate"])[3]), :(0.00010269 * ((QTLeffects["GCN2"])["qtlPDrate"])[3]), :(1.15585 * ((QTLeffects["GCN3"])["qtlTLrate"])[3]), :(0.00010269 * ((QTLeffects["GCN3"])["qtlPDrate"])[3]), :(1.15585 * ((QTLeffects["GCN4"])["qtlTLrate"])[3]), :(0.00010269 * ((QTLeffects["GCN4"])["qtlPDrate"])[3]), :(3.5602 * ((QTLeffects["GCN1"])["qtlTLrate"])[4]), :(0.000100725 * ((QTLeffects["GCN1"])["qtlPDrate"])[4]), :(3.5602 * ((QTLeffects["GCN2"])["qtlTLrate"])[4]), :(0.000100725 * ((QTLeffects["GCN2"])["qtlPDrate"])[4]), :(3.5602 * ((QTLeffects["GCN3"])["qtlTLrate"])[4]), :(0.000100725 * ((QTLeffects["GCN3"])["qtlPDrate"])[4]), :(3.5602 * ((QTLeffects["GCN4"])["qtlTLrate"])[4]), :(0.000100725 * ((QTLeffects["GCN4"])["qtlPDrate"])[4]), :(0.567232 * ((QTLeffects["GCN1"])["qtlTLrate"])[5]), :(9.17936e-5 * ((QTLeffects["GCN1"])["qtlPDrate"])[5]), :(0.567232 * ((QTLeffects["GCN2"])["qtlTLrate"])[5]), :(9.17936e-5 * ((QTLeffects["GCN2"])["qtlPDrate"])[5]), :(0.567232 * ((QTLeffects["GCN3"])["qtlTLrate"])[5]), :(9.17936e-5 * ((QTLeffects["GCN3"])["qtlPDrate"])[5]), :(0.567232 * ((QTLeffects["GCN4"])["qtlTLrate"])[5]), :(9.17936e-5 * ((QTLeffects["GCN4"])["qtlPDrate"])[5]), :(4.44623 * ((QTLeffects["GCN1"])["qtlTLrate"])[7]), :(0.00010015 * ((QTLeffects["GCN1"])["qtlPDrate"])[7]), :(4.44623 * ((QTLeffects["GCN2"])["qtlTLrate"])[7]), :(0.00010015 * ((QTLeffects["GCN2"])["qtlPDrate"])[7]), :(4.44623 * ((QTLeffects["GCN3"])["qtlTLrate"])[7]), :(0.00010015 * ((QTLeffects["GCN3"])["qtlPDrate"])[7]), :(4.44623 * ((QTLeffects["GCN4"])["qtlTLrate"])[7]), :(0.00010015 * ((QTLeffects["GCN4"])["qtlPDrate"])[7])],"reactions"=>Any["0 --> R1GCN1", "R1GCN1 --> 0", "0 --> R1GCN2", "R1GCN2 --> 0", "0 --> R1GCN3", "R1GCN3 --> 0", "0 --> R1GCN4", "R1GCN4 --> 0", "0 --> R2GCN1", "R2GCN1 --> 0", "0 --> R2GCN2", "R2GCN2 --> 0", "0 --> R2GCN3", "R2GCN3 --> 0", "0 --> R2GCN4", "R2GCN4 --> 0", "0 --> R3GCN1", "R3GCN1 --> 0", "0 --> R3GCN2", "R3GCN2 --> 0", "0 --> R3GCN3", "R3GCN3 --> 0", "0 --> R3GCN4", "R3GCN4 --> 0", "0 --> R4GCN1", "R4GCN1 --> 0", "0 --> R4GCN2", "R4GCN2 --> 0", "0 --> R4GCN3", "R4GCN3 --> 0", "0 --> R4GCN4", "R4GCN4 --> 0", "0 --> R5GCN1", "R5GCN1 --> 0", "0 --> R5GCN2", "R5GCN2 --> 0", "0 --> R5GCN3", "R5GCN3 --> 0", "0 --> R5GCN4", "R5GCN4 --> 0", "0 --> R6GCN1", "R6GCN1 --> 0", "0 --> R6GCN2", "R6GCN2 --> 0", "0 --> R6GCN3", "R6GCN3 --> 0", "0 --> R6GCN4", "R6GCN4 --> 0", "0 --> R7GCN1", "R7GCN1 --> 0", "0 --> R7GCN2", "R7GCN2 --> 0", "0 --> R7GCN3", "R7GCN3 --> 0", "0 --> R7GCN4", "R7GCN4 --> 0", "0 --> R8GCN1", "R8GCN1 --> 0", "0 --> R8GCN2", "R8GCN2 --> 0", "0 --> R8GCN3", "R8GCN3 --> 0", "0 --> R8GCN4", "R8GCN4 --> 0", "0 --> R9GCN1", "R9GCN1 --> 0", "0 --> R9GCN2", "R9GCN2 --> 0", "0 --> R9GCN3", "R9GCN3 --> 0", "0 --> R9GCN4", "R9GCN4 --> 0", "0 --> R10GCN1", "R10GCN1 --> 0", "0 --> R10GCN2", "R10GCN2 --> 0", "0 --> R10GCN3", "R10GCN3 --> 0", "0 --> R10GCN4", "R10GCN4 --> 0", "R1GCN1 --> P1GCN1", "P1GCN1 --> 0", "R1GCN2 --> P1GCN2", "P1GCN2 --> 0", "R1GCN3 --> P1GCN3", "P1GCN3 --> 0", "R1GCN4 --> P1GCN4", "P1GCN4 --> 0", "R2GCN1 --> P2GCN1", "P2GCN1 --> 0", "R2GCN2 --> P2GCN2", "P2GCN2 --> 0", "R2GCN3 --> P2GCN3", "P2GCN3 --> 0", "R2GCN4 --> P2GCN4", "P2GCN4 --> 0", "R3GCN1 --> P3GCN1", "P3GCN1 --> 0", "R3GCN2 --> P3GCN2", "P3GCN2 --> 0", "R3GCN3 --> P3GCN3", "P3GCN3 --> 0", "R3GCN4 --> P3GCN4", "P3GCN4 --> 0", "R4GCN1 --> P4GCN1", "P4GCN1 --> 0", "R4GCN2 --> P4GCN2", "P4GCN2 --> 0", "R4GCN3 --> P4GCN3", "P4GCN3 --> 0", "R4GCN4 --> P4GCN4", "P4GCN4 --> 0", "R5GCN1 --> P5GCN1", "P5GCN1 --> 0", "R5GCN2 --> P5GCN2", "P5GCN2 --> 0", "R5GCN3 --> P5GCN3", "P5GCN3 --> 0", "R5GCN4 --> P5GCN4", "P5GCN4 --> 0", "R7GCN1 --> P7GCN1", "P7GCN1 --> 0", "R7GCN2 --> P7GCN2", "P7GCN2 --> 0", "R7GCN3 --> P7GCN3", "P7GCN3 --> 0", "R7GCN4 --> P7GCN4", "P7GCN4 --> 0"],"initialconditions"=>Any[:(((28.8461 * ((QTLeffects["GCN1"])["qtlTCrate"])[1]) / ((QTLeffects["GCN1"])["qtlRDrate"])[1]) * ((InitVar["GCN1"])["R"])[1]), :(((28.8461 * ((QTLeffects["GCN2"])["qtlTCrate"])[1]) / ((QTLeffects["GCN2"])["qtlRDrate"])[1]) * ((InitVar["GCN2"])["R"])[1]), :(((28.8461 * ((QTLeffects["GCN3"])["qtlTCrate"])[1]) / ((QTLeffects["GCN3"])["qtlRDrate"])[1]) * ((InitVar["GCN3"])["R"])[1]), :(((28.8461 * ((QTLeffects["GCN4"])["qtlTCrate"])[1]) / ((QTLeffects["GCN4"])["qtlRDrate"])[1]) * ((InitVar["GCN4"])["R"])[1]), :(((81.3848 * ((QTLeffects["GCN1"])["qtlTCrate"])[2]) / ((QTLeffects["GCN1"])["qtlRDrate"])[2]) * ((InitVar["GCN1"])["R"])[2]), :(((81.3848 * ((QTLeffects["GCN2"])["qtlTCrate"])[2]) / ((QTLeffects["GCN2"])["qtlRDrate"])[2]) * ((InitVar["GCN2"])["R"])[2]), :(((81.3848 * ((QTLeffects["GCN3"])["qtlTCrate"])[2]) / ((QTLeffects["GCN3"])["qtlRDrate"])[2]) * ((InitVar["GCN3"])["R"])[2]), :(((81.3848 * ((QTLeffects["GCN4"])["qtlTCrate"])[2]) / ((QTLeffects["GCN4"])["qtlRDrate"])[2]) * ((InitVar["GCN4"])["R"])[2]), :(((24.6875 * ((QTLeffects["GCN1"])["qtlTCrate"])[3]) / ((QTLeffects["GCN1"])["qtlRDrate"])[3]) * ((InitVar["GCN1"])["R"])[3]), :(((24.6875 * ((QTLeffects["GCN2"])["qtlTCrate"])[3]) / ((QTLeffects["GCN2"])["qtlRDrate"])[3]) * ((InitVar["GCN2"])["R"])[3]), :(((24.6875 * ((QTLeffects["GCN3"])["qtlTCrate"])[3]) / ((QTLeffects["GCN3"])["qtlRDrate"])[3]) * ((InitVar["GCN3"])["R"])[3]), :(((24.6875 * ((QTLeffects["GCN4"])["qtlTCrate"])[3]) / ((QTLeffects["GCN4"])["qtlRDrate"])[3]) * ((InitVar["GCN4"])["R"])[3]), :(((313.326 * ((QTLeffects["GCN1"])["qtlTCrate"])[4]) / ((QTLeffects["GCN1"])["qtlRDrate"])[4]) * ((InitVar["GCN1"])["R"])[4]), :(((313.326 * ((QTLeffects["GCN2"])["qtlTCrate"])[4]) / ((QTLeffects["GCN2"])["qtlRDrate"])[4]) * ((InitVar["GCN2"])["R"])[4]), :(((313.326 * ((QTLeffects["GCN3"])["qtlTCrate"])[4]) / ((QTLeffects["GCN3"])["qtlRDrate"])[4]) * ((InitVar["GCN3"])["R"])[4]), :(((313.326 * ((QTLeffects["GCN4"])["qtlTCrate"])[4]) / ((QTLeffects["GCN4"])["qtlRDrate"])[4]) * ((InitVar["GCN4"])["R"])[4]), :(((29.0011 * ((QTLeffects["GCN1"])["qtlTCrate"])[5]) / ((QTLeffects["GCN1"])["qtlRDrate"])[5]) * ((InitVar["GCN1"])["R"])[5]), :(((29.0011 * ((QTLeffects["GCN2"])["qtlTCrate"])[5]) / ((QTLeffects["GCN2"])["qtlRDrate"])[5]) * ((InitVar["GCN2"])["R"])[5]), :(((29.0011 * ((QTLeffects["GCN3"])["qtlTCrate"])[5]) / ((QTLeffects["GCN3"])["qtlRDrate"])[5]) * ((InitVar["GCN3"])["R"])[5]), :(((29.0011 * ((QTLeffects["GCN4"])["qtlTCrate"])[5]) / ((QTLeffects["GCN4"])["qtlRDrate"])[5]) * ((InitVar["GCN4"])["R"])[5]), :(((164.491 * ((QTLeffects["GCN1"])["qtlTCrate"])[6]) / ((QTLeffects["GCN1"])["qtlRDrate"])[6]) * ((InitVar["GCN1"])["R"])[6]), :(((164.491 * ((QTLeffects["GCN2"])["qtlTCrate"])[6]) / ((QTLeffects["GCN2"])["qtlRDrate"])[6]) * ((InitVar["GCN2"])["R"])[6]), :(((164.491 * ((QTLeffects["GCN3"])["qtlTCrate"])[6]) / ((QTLeffects["GCN3"])["qtlRDrate"])[6]) * ((InitVar["GCN3"])["R"])[6]), :(((164.491 * ((QTLeffects["GCN4"])["qtlTCrate"])[6]) / ((QTLeffects["GCN4"])["qtlRDrate"])[6]) * ((InitVar["GCN4"])["R"])[6]), :(((71.5497 * ((QTLeffects["GCN1"])["qtlTCrate"])[7]) / ((QTLeffects["GCN1"])["qtlRDrate"])[7]) * ((InitVar["GCN1"])["R"])[7]), :(((71.5497 * ((QTLeffects["GCN2"])["qtlTCrate"])[7]) / ((QTLeffects["GCN2"])["qtlRDrate"])[7]) * ((InitVar["GCN2"])["R"])[7]), :(((71.5497 * ((QTLeffects["GCN3"])["qtlTCrate"])[7]) / ((QTLeffects["GCN3"])["qtlRDrate"])[7]) * ((InitVar["GCN3"])["R"])[7]), :(((71.5497 * ((QTLeffects["GCN4"])["qtlTCrate"])[7]) / ((QTLeffects["GCN4"])["qtlRDrate"])[7]) * ((InitVar["GCN4"])["R"])[7]), :(((191.649 * ((QTLeffects["GCN1"])["qtlTCrate"])[8]) / ((QTLeffects["GCN1"])["qtlRDrate"])[8]) * ((InitVar["GCN1"])["R"])[8]), :(((191.649 * ((QTLeffects["GCN2"])["qtlTCrate"])[8]) / ((QTLeffects["GCN2"])["qtlRDrate"])[8]) * ((InitVar["GCN2"])["R"])[8]), :(((191.649 * ((QTLeffects["GCN3"])["qtlTCrate"])[8]) / ((QTLeffects["GCN3"])["qtlRDrate"])[8]) * ((InitVar["GCN3"])["R"])[8]), :(((191.649 * ((QTLeffects["GCN4"])["qtlTCrate"])[8]) / ((QTLeffects["GCN4"])["qtlRDrate"])[8]) * ((InitVar["GCN4"])["R"])[8]), :(((79.9315 * ((QTLeffects["GCN1"])["qtlTCrate"])[9]) / ((QTLeffects["GCN1"])["qtlRDrate"])[9]) * ((InitVar["GCN1"])["R"])[9]), :(((79.9315 * ((QTLeffects["GCN2"])["qtlTCrate"])[9]) / ((QTLeffects["GCN2"])["qtlRDrate"])[9]) * ((InitVar["GCN2"])["R"])[9]), :(((79.9315 * ((QTLeffects["GCN3"])["qtlTCrate"])[9]) / ((QTLeffects["GCN3"])["qtlRDrate"])[9]) * ((InitVar["GCN3"])["R"])[9]), :(((79.9315 * ((QTLeffects["GCN4"])["qtlTCrate"])[9]) / ((QTLeffects["GCN4"])["qtlRDrate"])[9]) * ((InitVar["GCN4"])["R"])[9]), :(((293.694 * ((QTLeffects["GCN1"])["qtlTCrate"])[10]) / ((QTLeffects["GCN1"])["qtlRDrate"])[10]) * ((InitVar["GCN1"])["R"])[10]), :(((293.694 * ((QTLeffects["GCN2"])["qtlTCrate"])[10]) / ((QTLeffects["GCN2"])["qtlRDrate"])[10]) * ((InitVar["GCN2"])["R"])[10]), :(((293.694 * ((QTLeffects["GCN3"])["qtlTCrate"])[10]) / ((QTLeffects["GCN3"])["qtlRDrate"])[10]) * ((InitVar["GCN3"])["R"])[10]), :(((293.694 * ((QTLeffects["GCN4"])["qtlTCrate"])[10]) / ((QTLeffects["GCN4"])["qtlRDrate"])[10]) * ((InitVar["GCN4"])["R"])[10]), :(((1.26595e6 * ((QTLeffects["GCN1"])["qtlTCrate"])[1] * ((QTLeffects["GCN1"])["qtlTLrate"])[1]) / (((QTLeffects["GCN1"])["qtlRDrate"])[1] * ((QTLeffects["GCN1"])["qtlPDrate"])[1])) * ((InitVar["GCN1"])["P"])[1]), :(((1.26595e6 * ((QTLeffects["GCN2"])["qtlTCrate"])[1] * ((QTLeffects["GCN2"])["qtlTLrate"])[1]) / (((QTLeffects["GCN2"])["qtlRDrate"])[1] * ((QTLeffects["GCN2"])["qtlPDrate"])[1])) * ((InitVar["GCN2"])["P"])[1]), :(((1.26595e6 * ((QTLeffects["GCN3"])["qtlTCrate"])[1] * ((QTLeffects["GCN3"])["qtlTLrate"])[1]) / (((QTLeffects["GCN3"])["qtlRDrate"])[1] * ((QTLeffects["GCN3"])["qtlPDrate"])[1])) * ((InitVar["GCN3"])["P"])[1]), :(((1.26595e6 * ((QTLeffects["GCN4"])["qtlTCrate"])[1] * ((QTLeffects["GCN4"])["qtlTLrate"])[1]) / (((QTLeffects["GCN4"])["qtlRDrate"])[1] * ((QTLeffects["GCN4"])["qtlPDrate"])[1])) * ((InitVar["GCN4"])["P"])[1]), :(((1.97821e6 * ((QTLeffects["GCN1"])["qtlTCrate"])[2] * ((QTLeffects["GCN1"])["qtlTLrate"])[2]) / (((QTLeffects["GCN1"])["qtlRDrate"])[2] * ((QTLeffects["GCN1"])["qtlPDrate"])[2])) * ((InitVar["GCN1"])["P"])[2]), :(((1.97821e6 * ((QTLeffects["GCN2"])["qtlTCrate"])[2] * ((QTLeffects["GCN2"])["qtlTLrate"])[2]) / (((QTLeffects["GCN2"])["qtlRDrate"])[2] * ((QTLeffects["GCN2"])["qtlPDrate"])[2])) * ((InitVar["GCN2"])["P"])[2]), :(((1.97821e6 * ((QTLeffects["GCN3"])["qtlTCrate"])[2] * ((QTLeffects["GCN3"])["qtlTLrate"])[2]) / (((QTLeffects["GCN3"])["qtlRDrate"])[2] * ((QTLeffects["GCN3"])["qtlPDrate"])[2])) * ((InitVar["GCN3"])["P"])[2]), :(((1.97821e6 * ((QTLeffects["GCN4"])["qtlTCrate"])[2] * ((QTLeffects["GCN4"])["qtlTLrate"])[2]) / (((QTLeffects["GCN4"])["qtlRDrate"])[2] * ((QTLeffects["GCN4"])["qtlPDrate"])[2])) * ((InitVar["GCN4"])["P"])[2]), :(((2.77873e5 * ((QTLeffects["GCN1"])["qtlTCrate"])[3] * ((QTLeffects["GCN1"])["qtlTLrate"])[3]) / (((QTLeffects["GCN1"])["qtlRDrate"])[3] * ((QTLeffects["GCN1"])["qtlPDrate"])[3])) * ((InitVar["GCN1"])["P"])[3]), :(((2.77873e5 * ((QTLeffects["GCN2"])["qtlTCrate"])[3] * ((QTLeffects["GCN2"])["qtlTLrate"])[3]) / (((QTLeffects["GCN2"])["qtlRDrate"])[3] * ((QTLeffects["GCN2"])["qtlPDrate"])[3])) * ((InitVar["GCN2"])["P"])[3]), :(((2.77873e5 * ((QTLeffects["GCN3"])["qtlTCrate"])[3] * ((QTLeffects["GCN3"])["qtlTLrate"])[3]) / (((QTLeffects["GCN3"])["qtlRDrate"])[3] * ((QTLeffects["GCN3"])["qtlPDrate"])[3])) * ((InitVar["GCN3"])["P"])[3]), :(((2.77873e5 * ((QTLeffects["GCN4"])["qtlTCrate"])[3] * ((QTLeffects["GCN4"])["qtlTLrate"])[3]) / (((QTLeffects["GCN4"])["qtlRDrate"])[3] * ((QTLeffects["GCN4"])["qtlPDrate"])[3])) * ((InitVar["GCN4"])["P"])[3]), :(((1.10747e7 * ((QTLeffects["GCN1"])["qtlTCrate"])[4] * ((QTLeffects["GCN1"])["qtlTLrate"])[4]) / (((QTLeffects["GCN1"])["qtlRDrate"])[4] * ((QTLeffects["GCN1"])["qtlPDrate"])[4])) * ((InitVar["GCN1"])["P"])[4]), :(((1.10747e7 * ((QTLeffects["GCN2"])["qtlTCrate"])[4] * ((QTLeffects["GCN2"])["qtlTLrate"])[4]) / (((QTLeffects["GCN2"])["qtlRDrate"])[4] * ((QTLeffects["GCN2"])["qtlPDrate"])[4])) * ((InitVar["GCN2"])["P"])[4]), :(((1.10747e7 * ((QTLeffects["GCN3"])["qtlTCrate"])[4] * ((QTLeffects["GCN3"])["qtlTLrate"])[4]) / (((QTLeffects["GCN3"])["qtlRDrate"])[4] * ((QTLeffects["GCN3"])["qtlPDrate"])[4])) * ((InitVar["GCN3"])["P"])[4]), :(((1.10747e7 * ((QTLeffects["GCN4"])["qtlTCrate"])[4] * ((QTLeffects["GCN4"])["qtlTLrate"])[4]) / (((QTLeffects["GCN4"])["qtlRDrate"])[4] * ((QTLeffects["GCN4"])["qtlPDrate"])[4])) * ((InitVar["GCN4"])["P"])[4]), :(((1.7921e5 * ((QTLeffects["GCN1"])["qtlTCrate"])[5] * ((QTLeffects["GCN1"])["qtlTLrate"])[5]) / (((QTLeffects["GCN1"])["qtlRDrate"])[5] * ((QTLeffects["GCN1"])["qtlPDrate"])[5])) * ((InitVar["GCN1"])["P"])[5]), :(((1.7921e5 * ((QTLeffects["GCN2"])["qtlTCrate"])[5] * ((QTLeffects["GCN2"])["qtlTLrate"])[5]) / (((QTLeffects["GCN2"])["qtlRDrate"])[5] * ((QTLeffects["GCN2"])["qtlPDrate"])[5])) * ((InitVar["GCN2"])["P"])[5]), :(((1.7921e5 * ((QTLeffects["GCN3"])["qtlTCrate"])[5] * ((QTLeffects["GCN3"])["qtlTLrate"])[5]) / (((QTLeffects["GCN3"])["qtlRDrate"])[5] * ((QTLeffects["GCN3"])["qtlPDrate"])[5])) * ((InitVar["GCN3"])["P"])[5]), :(((1.7921e5 * ((QTLeffects["GCN4"])["qtlTCrate"])[5] * ((QTLeffects["GCN4"])["qtlTLrate"])[5]) / (((QTLeffects["GCN4"])["qtlRDrate"])[5] * ((QTLeffects["GCN4"])["qtlPDrate"])[5])) * ((InitVar["GCN4"])["P"])[5]), :(((3.17649e6 * ((QTLeffects["GCN1"])["qtlTCrate"])[7] * ((QTLeffects["GCN1"])["qtlTLrate"])[7]) / (((QTLeffects["GCN1"])["qtlRDrate"])[7] * ((QTLeffects["GCN1"])["qtlPDrate"])[7])) * ((InitVar["GCN1"])["P"])[7]), :(((3.17649e6 * ((QTLeffects["GCN2"])["qtlTCrate"])[7] * ((QTLeffects["GCN2"])["qtlTLrate"])[7]) / (((QTLeffects["GCN2"])["qtlRDrate"])[7] * ((QTLeffects["GCN2"])["qtlPDrate"])[7])) * ((InitVar["GCN2"])["P"])[7]), :(((3.17649e6 * ((QTLeffects["GCN3"])["qtlTCrate"])[7] * ((QTLeffects["GCN3"])["qtlTLrate"])[7]) / (((QTLeffects["GCN3"])["qtlRDrate"])[7] * ((QTLeffects["GCN3"])["qtlPDrate"])[7])) * ((InitVar["GCN3"])["P"])[7]), :(((3.17649e6 * ((QTLeffects["GCN4"])["qtlTCrate"])[7] * ((QTLeffects["GCN4"])["qtlTLrate"])[7]) / (((QTLeffects["GCN4"])["qtlRDrate"])[7] * ((QTLeffects["GCN4"])["qtlPDrate"])[7])) * ((InitVar["GCN4"])["P"])[7])],"species"=>Any["R1GCN1", "R1GCN2", "R1GCN3", "R1GCN4", "R2GCN1", "R2GCN2", "R2GCN3", "R2GCN4", "R3GCN1", "R3GCN2", "R3GCN3", "R3GCN4", "R4GCN1", "R4GCN2", "R4GCN3", "R4GCN4", "R5GCN1", "R5GCN2", "R5GCN3", "R5GCN4", "R6GCN1", "R6GCN2", "R6GCN3", "R6GCN4", "R7GCN1", "R7GCN2", "R7GCN3", "R7GCN4", "R8GCN1", "R8GCN2", "R8GCN3", "R8GCN4", "R9GCN1", "R9GCN2", "R9GCN3", "R9GCN4", "R10GCN1", "R10GCN2", "R10GCN3", "R10GCN4", "P1GCN1", "P1GCN2", "P1GCN3", "P1GCN4", "P2GCN1", "P2GCN2", "P2GCN3", "P2GCN4", "P3GCN1", "P3GCN2", "P3GCN3", "P3GCN4", "P4GCN1", "P4GCN2", "P4GCN3", "P4GCN4", "P5GCN1", "P5GCN2", "P5GCN3", "P5GCN4", "P7GCN1", "P7GCN2", "P7GCN3", "P7GCN4"])

QTLeffects = Dict("GCN3"=>Dict("qtlPDregbind"=>[1.0, 0.903064, 1.0, 1.01958, 0.908109, 0.0, 0.985851, 0.0, 0.0, 0.0],"qtlTLrate"=>[1.0, 0.98387, 1.12259, 0.947592, 1.02447, 0.0, 0.897265, 0.0, 0.0, 0.0],"qtlTCrate"=>[1.0, 1.17192, 1.0, 0.973951, 1.07116, 1.11191, 1.0003, 1.04185, 0.96825, 1.0],"qtlRDbindreg"=>[1.0, 0.969047, 0.883935, 1.16487, 1.18382, 0.926073, 1.13709, 1.16358, 1.20015, 1.0],"qtlTCregbind"=>[1.0, 0.939463, 1.0, 0.996446, 1.11144, 1.0, 0.92738, 1.0, 1.11756, 1.0],"qtlactivity"=>[1.0, 1.07163, 1.0, 0.948701, 1.00975, 1.0, 1.00126, 0.875307, 1.03423, 1.0],"qtlRDrate"=>[1.0, 1.0, 1.0, 1.08387, 1.02613, 1.0, 1.20753, 1.06082, 1.0, 1.0],"qtlTLregbind"=>[1.0, 0.969293, 1.0, 0.9764, 0.932668, 0.0, 0.918416, 0.0, 0.0, 0.0],"qtlPDrate"=>[1.0, 1.14071, 1.0, 0.93931, 1.02331, 0.0, 1.05153, 0.0, 0.0, 0.0]),"GCN2"=>Dict("qtlPDregbind"=>[1.0, 0.903064, 1.0, 0.904002, 1.0, 0.0, 0.947351, 0.0, 0.0, 0.0],"qtlTLrate"=>[0.888289, 0.98387, 1.0, 0.633331, 0.927585, 0.0, 1.03055, 0.0, 0.0, 0.0],"qtlTCrate"=>[1.0, 1.17192, 1.0, 1.12524, 1.0, 1.0, 0.993596, 0.99752, 1.0, 1.0],"qtlRDbindreg"=>[1.0, 0.969047, 1.0, 1.04387, 1.0, 0.888253, 1.0, 0.843412, 1.0, 1.0],"qtlTCregbind"=>[1.0, 0.939463, 1.0, 1.15848, 1.0, 0.919143, 1.0, 0.813784, 0.975534, 1.0],"qtlactivity"=>[1.0, 1.07163, 1.0, 0.864578, 1.0, 1.13318, 0.918793, 1.08566, 1.0, 1.0],"qtlRDrate"=>[1.0, 1.0, 1.0, 1.01624, 1.0, 1.0, 1.0, 0.926551, 0.941337, 1.0],"qtlTLregbind"=>[1.0, 0.969293, 1.0, 1.05993, 1.0, 0.0, 0.992744, 0.0, 0.0, 0.0],"qtlPDrate"=>[1.0, 1.14071, 1.0, 1.00784, 1.0, 0.0, 0.990836, 0.0, 0.0, 0.0]),"GCN4"=>Dict("qtlPDregbind"=>[1.0, 0.903064, 1.0, 1.0, 1.0, 0.0, 0.985851, 0.0, 0.0, 0.0],"qtlTLrate"=>[0.888289, 0.98387, 1.04376, 1.0, 1.0, 0.0, 0.897265, 0.0, 0.0, 0.0],"qtlTCrate"=>[1.0, 1.17192, 1.0, 1.0, 1.0, 1.0, 1.0003, 1.0, 0.934712, 1.09404],"qtlRDbindreg"=>[1.0, 0.969047, 1.08212, 1.0, 1.0, 1.0, 1.13709, 1.0, 0.953415, 1.08835],"qtlTCregbind"=>[1.0, 0.939463, 1.0, 1.0, 1.0, 1.0, 0.92738, 1.0, 1.0, 0.874129],"qtlactivity"=>[1.0, 1.07163, 0.985766, 1.0, 1.0, 1.0, 1.00126, 1.0, 0.87951, 0.774213],"qtlRDrate"=>[1.0, 1.0, 0.923331, 1.0, 1.0, 1.0, 1.20753, 1.0, 1.0, 1.0],"qtlTLregbind"=>[1.0, 0.969293, 0.947754, 0.790886, 1.0, 0.0, 0.918416, 0.0, 0.0, 0.0],"qtlPDrate"=>[1.0, 1.14071, 1.0, 1.0, 1.0, 0.0, 1.05153, 0.0, 0.0, 0.0]),"GCN1"=>Dict("qtlPDregbind"=>[1.15029, 1.16139, 1.0, 1.0, 0.965733, 0.0, 1.0, 0.0, 0.0, 0.0],"qtlTLrate"=>[1.0, 0.961051, 1.12259, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0],"qtlTCrate"=>[1.0, 1.24792, 1.0, 0.894408, 1.0, 1.0, 1.0, 0.99752, 0.96825, 1.0],"qtlRDbindreg"=>[1.0, 0.83549, 0.883935, 1.07553, 1.0, 0.870481, 1.0, 0.843412, 1.20015, 1.0],"qtlTCregbind"=>[1.0, 1.08633, 1.0, 1.06498, 1.03092, 1.0, 1.0, 0.813784, 1.11756, 1.0],"qtlactivity"=>[1.0, 0.791183, 1.0, 1.0, 0.948967, 1.0, 1.0, 1.08566, 1.03423, 1.0],"qtlRDrate"=>[1.0, 0.818898, 1.0, 1.16617, 1.09283, 1.0, 1.0, 0.926551, 1.0, 1.0],"qtlTLregbind"=>[1.0, 0.9524, 1.0, 0.970783, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0],"qtlPDrate"=>[1.11645, 0.940561, 1.0, 0.949767, 0.837376, 0.0, 1.0, 0.0, 0.0, 0.0]))
InitVar = Dict("GCN3"=>Dict("P"=>[0.922647, 0.827737, 0.962663, 1.07635, 0.848211, 1.01387, 1.15418, 0.997701, 0.989415, 1.06221],"R"=>[1.01076, 0.903517, 0.931234, 0.948392, 1.08051, 1.07374, 0.843925, 0.956626, 1.1761, 1.09528]),"GCN2"=>Dict("P"=>[0.853842, 0.93087, 1.12092, 0.964489, 1.08834, 1.05661, 0.833048, 1.06602, 1.16125, 1.10232],"R"=>[0.97554, 0.932472, 1.0511, 0.898645, 1.05261, 1.05966, 0.7716, 0.90272, 0.807083, 1.06358]),"GCN4"=>Dict("P"=>[0.990594, 1.09941, 1.0421, 0.930396, 0.81033, 0.90876, 0.86762, 1.13135, 0.93453, 1.04141],"R"=>[1.00492, 1.1564, 1.14431, 0.93374, 0.970863, 1.13996, 1.1174, 0.978497, 1.10544, 1.00269]),"GCN1"=>Dict("P"=>[1.07192, 1.15348, 1.06973, 0.926324, 1.04767, 0.955735, 1.08041, 0.796943, 0.988385, 0.919596],"R"=>[0.855099, 0.818183, 0.857508, 0.910632, 1.00853, 0.99178, 0.934635, 0.879292, 0.843702, 0.908478]))

model = Network("coucou")

for i in 1:length(stochmodel["species"])
  i0 = round(Int, eval(stochmodel["initialconditions"][i]))
  println(i0)
  model <= BioSimulator.Species(stochmodel["species"][i], i0)
  println(model.species_list)
end

for i in eachindex(stochmodel["reactions"])
  model <= BioSimulator.Reaction(stochmodel["reactionsnames"][i], eval(stochmodel["propensities"][i]), stochmodel["reactions"][i])
end

result = simulate(model, algorithm=SSA, time = 4.0, trials = 1)

=#





function whatisit(x)
  println(x)
  println(typeof(x))
end


function superfunction(ar, ind)
  println(size(ar))
  println(size(ar[ind]))
  println(typeof(ar[ind]))
  return ar[ind]
end




