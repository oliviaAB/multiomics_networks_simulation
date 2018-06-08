
using BioSimulator
include("winData/multiomics_networks_simulation/julia_functions.jl")


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


model = BioSimulator.Network("coucou")

model <= BioSimulator.Species("X", 25)
model <= BioSimulator.Reaction("transcription", 0.02, "0 --> X")
model <= BioSimulator.Reaction("decay", 0.00069, "X --> 0")

result = BioSimulator.simulate(model, algorithm = SSA, time = 1000.0, epochs = 1000)
resdf = res2df(result)

plot(resdf[:, :time], resdf[:, :X])

# -----------------

stochmodel = Dict("reactionsnames"=>Any["transcription1GCN1", "RNAdecayR1GCN1", "transcription1GCN2", "RNAdecayR1GCN2", "transcription1GCN3", "RNAdecayR1GCN3", "transcription1GCN4", "RNAdecayR1GCN4", "transcription2GCN1", "RNAdecayR2GCN1", "transcription2GCN2", "RNAdecayR2GCN2", "transcription2GCN3", "RNAdecayR2GCN3", "transcription2GCN4", "RNAdecayR2GCN4", "transcription3GCN1", "RNAdecayR3GCN1", "transcription3GCN2", "RNAdecayR3GCN2", "transcription3GCN3", "RNAdecayR3GCN3", "transcription3GCN4", "RNAdecayR3GCN4", "transcription4GCN1", "RNAdecayR4GCN1", "transcription4GCN2", "RNAdecayR4GCN2", "transcription4GCN3", "RNAdecayR4GCN3", "transcription4GCN4", "RNAdecayR4GCN4", "transcription5GCN1", "RNAdecayR5GCN1", "transcription5GCN2", "RNAdecayR5GCN2", "transcription5GCN3", "RNAdecayR5GCN3", "transcription5GCN4", "RNAdecayR5GCN4", "transcription6GCN1", "RNAdecayR6GCN1", "transcription6GCN2", "RNAdecayR6GCN2", "transcription6GCN3", "RNAdecayR6GCN3", "transcription6GCN4", "RNAdecayR6GCN4", "transcription7GCN1", "RNAdecayR7GCN1", "transcription7GCN2", "RNAdecayR7GCN2", "transcription7GCN3", "RNAdecayR7GCN3", "transcription7GCN4", "RNAdecayR7GCN4", "transcription8GCN1", "RNAdecayR8GCN1", "transcription8GCN2", "RNAdecayR8GCN2", "transcription8GCN3", "RNAdecayR8GCN3", "transcription8GCN4", "RNAdecayR8GCN4", "transcription9GCN1", "RNAdecayR9GCN1", "transcription9GCN2", "RNAdecayR9GCN2", "transcription9GCN3", "RNAdecayR9GCN3", "transcription9GCN4", "RNAdecayR9GCN4", "transcription10GCN1", "RNAdecayR10GCN1", "transcription10GCN2", "RNAdecayR10GCN2", "transcription10GCN3", "RNAdecayR10GCN3", "transcription10GCN4", "RNAdecayR10GCN4", "translation1GCN1", "proteindecayP1GCN1", "translation1GCN2", "proteindecayP1GCN2", "translation1GCN3", "proteindecayP1GCN3", "translation1GCN4", "proteindecayP1GCN4", "translation2GCN1", "proteindecayP2GCN1", "translation2GCN2", "proteindecayP2GCN2", "translation2GCN3", "proteindecayP2GCN3", "translation2GCN4", "proteindecayP2GCN4", "translation3GCN1", "proteindecayP3GCN1", "translation3GCN2", "proteindecayP3GCN2", "translation3GCN3", "proteindecayP3GCN3", "translation3GCN4", "proteindecayP3GCN4", "translation4GCN1", "proteindecayP4GCN1", "translation4GCN2", "proteindecayP4GCN2", "translation4GCN3", "proteindecayP4GCN3", "translation4GCN4", "proteindecayP4GCN4", "translation5GCN1", "proteindecayP5GCN1", "translation5GCN2", "proteindecayP5GCN2", "translation5GCN3", "proteindecayP5GCN3", "translation5GCN4", "proteindecayP5GCN4", "translation7GCN1", "proteindecayP7GCN1", "translation7GCN2", "proteindecayP7GCN2", "translation7GCN3", "proteindecayP7GCN3", "translation7GCN4", "proteindecayP7GCN4"],
    "propensities"=>Any[:(0.020032 * ((QTLeffects["GCN1"])["qtlTCrate"])[1]), :(0.000694444 * ((QTLeffects["GCN1"])["qtlRDrate"])[1]), :(0.020032 * ((QTLeffects["GCN2"])["qtlTCrate"])[1]), :(0.000694444 * ((QTLeffects["GCN2"])["qtlRDrate"])[1]), :(0.020032 * ((QTLeffects["GCN3"])["qtlTCrate"])[1]), :(0.000694444 * ((QTLeffects["GCN3"])["qtlRDrate"])[1]), :(0.020032 * ((QTLeffects["GCN4"])["qtlTCrate"])[1]), :(0.000694444 * ((QTLeffects["GCN4"])["qtlRDrate"])[1]), :(0.0435213 * ((QTLeffects["GCN1"])["qtlTCrate"])[2]), :(0.000534759 * ((QTLeffects["GCN1"])["qtlRDrate"])[2]), :(0.0435213 * ((QTLeffects["GCN2"])["qtlTCrate"])[2]), :(0.000534759 * ((QTLeffects["GCN2"])["qtlRDrate"])[2]), :(0.0435213 * ((QTLeffects["GCN3"])["qtlTCrate"])[2]), :(0.000534759 * ((QTLeffects["GCN3"])["qtlRDrate"])[2]), :(0.0435213 * ((QTLeffects["GCN4"])["qtlTCrate"])[2]), :(0.000534759 * ((QTLeffects["GCN4"])["qtlRDrate"])[2]), :(0.0766692 * ((QTLeffects["GCN1"])["qtlTCrate"])[3]), :(0.00310559 * ((QTLeffects["GCN1"])["qtlRDrate"])[3]), :(0.0766692 * ((QTLeffects["GCN2"])["qtlTCrate"])[3]), :(0.00310559 * ((QTLeffects["GCN2"])["qtlRDrate"])[3]), :(0.0766692 * ((QTLeffects["GCN3"])["qtlTCrate"])[3]), :(0.00310559 * ((QTLeffects["GCN3"])["qtlRDrate"])[3]), :(0.0766692 * ((QTLeffects["GCN4"])["qtlTCrate"])[3]), :(0.00310559 * ((QTLeffects["GCN4"])["qtlRDrate"])[3]), :(0.0897525 * ((QTLeffects["GCN1"])["qtlTCrate"])[4]), :(0.000286451 * ((QTLeffects["GCN1"])["qtlRDrate"])[4]), :(0.0897525 * ((QTLeffects["GCN2"])["qtlTCrate"])[4]), :(0.000286451 * ((QTLeffects["GCN2"])["qtlRDrate"])[4]), :(0.0897525 * ((QTLeffects["GCN3"])["qtlTCrate"])[4]), :(0.000286451 * ((QTLeffects["GCN3"])["qtlRDrate"])[4]), :(0.0897525 * ((QTLeffects["GCN4"])["qtlTCrate"])[4]), :(0.000286451 * ((QTLeffects["GCN4"])["qtlRDrate"])[4]), :(0.071785 * ((QTLeffects["GCN1"])["qtlTCrate"])[5]), :(0.00247525 * ((QTLeffects["GCN1"])["qtlRDrate"])[5]), :(0.071785 * ((QTLeffects["GCN2"])["qtlTCrate"])[5]), :(0.00247525 * ((QTLeffects["GCN2"])["qtlRDrate"])[5]), :(0.071785 * ((QTLeffects["GCN3"])["qtlTCrate"])[5]), :(0.00247525 * ((QTLeffects["GCN3"])["qtlRDrate"])[5]), :(0.071785 * ((QTLeffects["GCN4"])["qtlTCrate"])[5]), :(0.00247525 * ((QTLeffects["GCN4"])["qtlRDrate"])[5]), :(0.0514678 * ((QTLeffects["GCN1"])["qtlTCrate"])[6]), :(0.000312891 * ((QTLeffects["GCN1"])["qtlRDrate"])[6]), :(0.0514678 * ((QTLeffects["GCN2"])["qtlTCrate"])[6]), :(0.000312891 * ((QTLeffects["GCN2"])["qtlRDrate"])[6]), :(0.0514678 * ((QTLeffects["GCN3"])["qtlTCrate"])[6]), :(0.000312891 * ((QTLeffects["GCN3"])["qtlRDrate"])[6]), :(0.0514678 * ((QTLeffects["GCN4"])["qtlTCrate"])[6]), :(0.000312891 * ((QTLeffects["GCN4"])["qtlRDrate"])[6]), :(0.0220153 * ((QTLeffects["GCN1"])["qtlTCrate"])[7]), :(0.000307692 * ((QTLeffects["GCN1"])["qtlRDrate"])[7]), :(0.0220153 * ((QTLeffects["GCN2"])["qtlTCrate"])[7]), :(0.000307692 * ((QTLeffects["GCN2"])["qtlRDrate"])[7]), :(0.0220153 * ((QTLeffects["GCN3"])["qtlTCrate"])[7]), :(0.000307692 * ((QTLeffects["GCN3"])["qtlRDrate"])[7]), :(0.0220153 * ((QTLeffects["GCN4"])["qtlTCrate"])[7]), :(0.000307692 * ((QTLeffects["GCN4"])["qtlRDrate"])[7]), :(0.0721572 * ((QTLeffects["GCN1"])["qtlTCrate"])[8]), :(0.000376506 * ((QTLeffects["GCN1"])["qtlRDrate"])[8]), :(0.0721572 * ((QTLeffects["GCN2"])["qtlTCrate"])[8]), :(0.000376506 * ((QTLeffects["GCN2"])["qtlRDrate"])[8]), :(0.0721572 * ((QTLeffects["GCN3"])["qtlTCrate"])[8]), :(0.000376506 * ((QTLeffects["GCN3"])["qtlRDrate"])[8]), :(0.0721572 * ((QTLeffects["GCN4"])["qtlTCrate"])[8]), :(0.000376506 * ((QTLeffects["GCN4"])["qtlRDrate"])[8]), :(0.0658958 * ((QTLeffects["GCN1"])["qtlTCrate"])[9]), :(0.000824402 * ((QTLeffects["GCN1"])["qtlRDrate"])[9]), :(0.0658958 * ((QTLeffects["GCN2"])["qtlTCrate"])[9]), :(0.000824402 * ((QTLeffects["GCN2"])["qtlRDrate"])[9]), :(0.0658958 * ((QTLeffects["GCN3"])["qtlTCrate"])[9]), :(0.000824402 * ((QTLeffects["GCN3"])["qtlRDrate"])[9]), :(0.0658958 * ((QTLeffects["GCN4"])["qtlTCrate"])[9]), :(0.000824402 * ((QTLeffects["GCN4"])["qtlRDrate"])[9]), :(0.0930883 * ((QTLeffects["GCN1"])["qtlTCrate"])[10]), :(0.000316957 * ((QTLeffects["GCN1"])["qtlRDrate"])[10]), :(0.0930883 * ((QTLeffects["GCN2"])["qtlTCrate"])[10]), :(0.000316957 * ((QTLeffects["GCN2"])["qtlRDrate"])[10]), :(0.0930883 * ((QTLeffects["GCN3"])["qtlTCrate"])[10]), :(0.000316957 * ((QTLeffects["GCN3"])["qtlRDrate"])[10]), :(0.0930883 * ((QTLeffects["GCN4"])["qtlTCrate"])[10]), :(0.000316957 * ((QTLeffects["GCN4"])["qtlRDrate"])[10]), :(3.34576 * ((QTLeffects["GCN1"])["qtlTLrate"])[1]), :(7.62369e-5 * ((QTLeffects["GCN1"])["qtlPDrate"])[1]), :(3.34576 * ((QTLeffects["GCN2"])["qtlTLrate"])[1]), :(7.62369e-5 * ((QTLeffects["GCN2"])["qtlPDrate"])[1]), :(3.34576 * ((QTLeffects["GCN3"])["qtlTLrate"])[1]), :(7.62369e-5 * ((QTLeffects["GCN3"])["qtlPDrate"])[1]), :(3.34576 * ((QTLeffects["GCN4"])["qtlTLrate"])[1]), :(7.62369e-5 * ((QTLeffects["GCN4"])["qtlPDrate"])[1]), :(3.00529 * ((QTLeffects["GCN1"])["qtlTLrate"])[2]), :(0.00012364 * ((QTLeffects["GCN1"])["qtlPDrate"])[2]), :(3.00529 * ((QTLeffects["GCN2"])["qtlTLrate"])[2]), :(0.00012364 * ((QTLeffects["GCN2"])["qtlPDrate"])[2]), :(3.00529 * ((QTLeffects["GCN3"])["qtlTLrate"])[2]), :(0.00012364 * ((QTLeffects["GCN3"])["qtlPDrate"])[2]), :(3.00529 * ((QTLeffects["GCN4"])["qtlTLrate"])[2]), :(0.00012364 * ((QTLeffects["GCN4"])["qtlPDrate"])[2]), :(1.15585 * ((QTLeffects["GCN1"])["qtlTLrate"])[3]), :(0.00010269 * ((QTLeffects["GCN1"])["qtlPDrate"])[3]), :(1.15585 * ((QTLeffects["GCN2"])["qtlTLrate"])[3]), :(0.00010269 * ((QTLeffects["GCN2"])["qtlPDrate"])[3]), :(1.15585 * ((QTLeffects["GCN3"])["qtlTLrate"])[3]), :(0.00010269 * ((QTLeffects["GCN3"])["qtlPDrate"])[3]), :(1.15585 * ((QTLeffects["GCN4"])["qtlTLrate"])[3]), :(0.00010269 * ((QTLeffects["GCN4"])["qtlPDrate"])[3]), :(3.5602 * ((QTLeffects["GCN1"])["qtlTLrate"])[4]), :(0.000100725 * ((QTLeffects["GCN1"])["qtlPDrate"])[4]), :(3.5602 * ((QTLeffects["GCN2"])["qtlTLrate"])[4]), :(0.000100725 * ((QTLeffects["GCN2"])["qtlPDrate"])[4]), :(3.5602 * ((QTLeffects["GCN3"])["qtlTLrate"])[4]), :(0.000100725 * ((QTLeffects["GCN3"])["qtlPDrate"])[4]), :(3.5602 * ((QTLeffects["GCN4"])["qtlTLrate"])[4]), :(0.000100725 * ((QTLeffects["GCN4"])["qtlPDrate"])[4]), :(0.567232 * ((QTLeffects["GCN1"])["qtlTLrate"])[5]), :(9.17936e-5 * ((QTLeffects["GCN1"])["qtlPDrate"])[5]), :(0.567232 * ((QTLeffects["GCN2"])["qtlTLrate"])[5]), :(9.17936e-5 * ((QTLeffects["GCN2"])["qtlPDrate"])[5]), :(0.567232 * ((QTLeffects["GCN3"])["qtlTLrate"])[5]), :(9.17936e-5 * ((QTLeffects["GCN3"])["qtlPDrate"])[5]), :(0.567232 * ((QTLeffects["GCN4"])["qtlTLrate"])[5]), :(9.17936e-5 * ((QTLeffects["GCN4"])["qtlPDrate"])[5]), :(4.44623 * ((QTLeffects["GCN1"])["qtlTLrate"])[7]), :(0.00010015 * ((QTLeffects["GCN1"])["qtlPDrate"])[7]), :(4.44623 * ((QTLeffects["GCN2"])["qtlTLrate"])[7]), :(0.00010015 * ((QTLeffects["GCN2"])["qtlPDrate"])[7]), :(4.44623 * ((QTLeffects["GCN3"])["qtlTLrate"])[7]), :(0.00010015 * ((QTLeffects["GCN3"])["qtlPDrate"])[7]), :(4.44623 * ((QTLeffects["GCN4"])["qtlTLrate"])[7]), :(0.00010015 * ((QTLeffects["GCN4"])["qtlPDrate"])[7])],
    "reactions"=>Any["0 --> R1GCN1", "R1GCN1 --> 0", "0 --> R1GCN2", "R1GCN2 --> 0", "0 --> R1GCN3", "R1GCN3 --> 0", "0 --> R1GCN4", "R1GCN4 --> 0", "0 --> R2GCN1", "R2GCN1 --> 0", "0 --> R2GCN2", "R2GCN2 --> 0", "0 --> R2GCN3", "R2GCN3 --> 0", "0 --> R2GCN4", "R2GCN4 --> 0", "0 --> R3GCN1", "R3GCN1 --> 0", "0 --> R3GCN2", "R3GCN2 --> 0", "0 --> R3GCN3", "R3GCN3 --> 0", "0 --> R3GCN4", "R3GCN4 --> 0", "0 --> R4GCN1", "R4GCN1 --> 0", "0 --> R4GCN2", "R4GCN2 --> 0", "0 --> R4GCN3", "R4GCN3 --> 0", "0 --> R4GCN4", "R4GCN4 --> 0", "0 --> R5GCN1", "R5GCN1 --> 0", "0 --> R5GCN2", "R5GCN2 --> 0", "0 --> R5GCN3", "R5GCN3 --> 0", "0 --> R5GCN4", "R5GCN4 --> 0", "0 --> R6GCN1", "R6GCN1 --> 0", "0 --> R6GCN2", "R6GCN2 --> 0", "0 --> R6GCN3", "R6GCN3 --> 0", "0 --> R6GCN4", "R6GCN4 --> 0", "0 --> R7GCN1", "R7GCN1 --> 0", "0 --> R7GCN2", "R7GCN2 --> 0", "0 --> R7GCN3", "R7GCN3 --> 0", "0 --> R7GCN4", "R7GCN4 --> 0", "0 --> R8GCN1", "R8GCN1 --> 0", "0 --> R8GCN2", "R8GCN2 --> 0", "0 --> R8GCN3", "R8GCN3 --> 0", "0 --> R8GCN4", "R8GCN4 --> 0", "0 --> R9GCN1", "R9GCN1 --> 0", "0 --> R9GCN2", "R9GCN2 --> 0", "0 --> R9GCN3", "R9GCN3 --> 0", "0 --> R9GCN4", "R9GCN4 --> 0", "0 --> R10GCN1", "R10GCN1 --> 0", "0 --> R10GCN2", "R10GCN2 --> 0", "0 --> R10GCN3", "R10GCN3 --> 0", "0 --> R10GCN4", "R10GCN4 --> 0", "R1GCN1 --> R1GCN1 + P1GCN1", "P1GCN1 --> 0", "R1GCN2 --> R1GCN2 + P1GCN2", "P1GCN2 --> 0", "R1GCN3 --> R1GCN3 + P1GCN3", "P1GCN3 --> 0", "R1GCN4 --> R1GCN4 + P1GCN4", "P1GCN4 --> 0", "R2GCN1 --> R2GCN1 + P2GCN1", "P2GCN1 --> 0", "R2GCN2 --> R2GCN2 + P2GCN2", "P2GCN2 --> 0", "R2GCN3 --> R2GCN3 + P2GCN3", "P2GCN3 --> 0", "R2GCN4 --> R2GCN4 + P2GCN4", "P2GCN4 --> 0", "R3GCN1 --> R3GCN1 + P3GCN1", "P3GCN1 --> 0", "R3GCN2 --> R3GCN2 + P3GCN2", "P3GCN2 --> 0", "R3GCN3 --> R3GCN3 + P3GCN3", "P3GCN3 --> 0", "R3GCN4 --> R3GCN4 + P3GCN4", "P3GCN4 --> 0", "R4GCN1 --> R4GCN1 + P4GCN1", "P4GCN1 --> 0", "R4GCN2 --> R4GCN2 + P4GCN2", "P4GCN2 --> 0", "R4GCN3 --> R4GCN3 + P4GCN3", "P4GCN3 --> 0", "R4GCN4 --> R4GCN4 + P4GCN4", "P4GCN4 --> 0", "R5GCN1 --> R5GCN1 + P5GCN1", "P5GCN1 --> 0", "R5GCN2 --> R5GCN2 + P5GCN2", "P5GCN2 --> 0", "R5GCN3 --> R5GCN3 + P5GCN3", "P5GCN3 --> 0", "R5GCN4 --> R5GCN4 + P5GCN4", "P5GCN4 --> 0", "R7GCN1 --> R7GCN1 + P7GCN1", "P7GCN1 --> 0", "R7GCN2 --> R7GCN2 + P7GCN2", "P7GCN2 --> 0", "R7GCN3 --> R7GCN3 + P7GCN3", "P7GCN3 --> 0", "R7GCN4 --> R7GCN4 + P7GCN4", "P7GCN4 --> 0"],
    "initialconditions"=>Any[:(((28.8461 * ((QTLeffects["GCN1"])["qtlTCrate"])[1]) / ((QTLeffects["GCN1"])["qtlRDrate"])[1]) * ((Initvar["GCN1"])["R"])[1]), :(((28.8461 * ((QTLeffects["GCN2"])["qtlTCrate"])[1]) / ((QTLeffects["GCN2"])["qtlRDrate"])[1]) * ((Initvar["GCN2"])["R"])[1]), :(((28.8461 * ((QTLeffects["GCN3"])["qtlTCrate"])[1]) / ((QTLeffects["GCN3"])["qtlRDrate"])[1]) * ((Initvar["GCN3"])["R"])[1]), :(((28.8461 * ((QTLeffects["GCN4"])["qtlTCrate"])[1]) / ((QTLeffects["GCN4"])["qtlRDrate"])[1]) * ((Initvar["GCN4"])["R"])[1]), :(((81.3848 * ((QTLeffects["GCN1"])["qtlTCrate"])[2]) / ((QTLeffects["GCN1"])["qtlRDrate"])[2]) * ((Initvar["GCN1"])["R"])[2]), :(((81.3848 * ((QTLeffects["GCN2"])["qtlTCrate"])[2]) / ((QTLeffects["GCN2"])["qtlRDrate"])[2]) * ((Initvar["GCN2"])["R"])[2]), :(((81.3848 * ((QTLeffects["GCN3"])["qtlTCrate"])[2]) / ((QTLeffects["GCN3"])["qtlRDrate"])[2]) * ((Initvar["GCN3"])["R"])[2]), :(((81.3848 * ((QTLeffects["GCN4"])["qtlTCrate"])[2]) / ((QTLeffects["GCN4"])["qtlRDrate"])[2]) * ((Initvar["GCN4"])["R"])[2]), :(((24.6875 * ((QTLeffects["GCN1"])["qtlTCrate"])[3]) / ((QTLeffects["GCN1"])["qtlRDrate"])[3]) * ((Initvar["GCN1"])["R"])[3]), :(((24.6875 * ((QTLeffects["GCN2"])["qtlTCrate"])[3]) / ((QTLeffects["GCN2"])["qtlRDrate"])[3]) * ((Initvar["GCN2"])["R"])[3]), :(((24.6875 * ((QTLeffects["GCN3"])["qtlTCrate"])[3]) / ((QTLeffects["GCN3"])["qtlRDrate"])[3]) * ((Initvar["GCN3"])["R"])[3]), :(((24.6875 * ((QTLeffects["GCN4"])["qtlTCrate"])[3]) / ((QTLeffects["GCN4"])["qtlRDrate"])[3]) * ((Initvar["GCN4"])["R"])[3]), :(((313.326 * ((QTLeffects["GCN1"])["qtlTCrate"])[4]) / ((QTLeffects["GCN1"])["qtlRDrate"])[4]) * ((Initvar["GCN1"])["R"])[4]), :(((313.326 * ((QTLeffects["GCN2"])["qtlTCrate"])[4]) / ((QTLeffects["GCN2"])["qtlRDrate"])[4]) * ((Initvar["GCN2"])["R"])[4]), :(((313.326 * ((QTLeffects["GCN3"])["qtlTCrate"])[4]) / ((QTLeffects["GCN3"])["qtlRDrate"])[4]) * ((Initvar["GCN3"])["R"])[4]), :(((313.326 * ((QTLeffects["GCN4"])["qtlTCrate"])[4]) / ((QTLeffects["GCN4"])["qtlRDrate"])[4]) * ((Initvar["GCN4"])["R"])[4]), :(((29.0011 * ((QTLeffects["GCN1"])["qtlTCrate"])[5]) / ((QTLeffects["GCN1"])["qtlRDrate"])[5]) * ((Initvar["GCN1"])["R"])[5]), :(((29.0011 * ((QTLeffects["GCN2"])["qtlTCrate"])[5]) / ((QTLeffects["GCN2"])["qtlRDrate"])[5]) * ((Initvar["GCN2"])["R"])[5]), :(((29.0011 * ((QTLeffects["GCN3"])["qtlTCrate"])[5]) / ((QTLeffects["GCN3"])["qtlRDrate"])[5]) * ((Initvar["GCN3"])["R"])[5]), :(((29.0011 * ((QTLeffects["GCN4"])["qtlTCrate"])[5]) / ((QTLeffects["GCN4"])["qtlRDrate"])[5]) * ((Initvar["GCN4"])["R"])[5]), :(((164.491 * ((QTLeffects["GCN1"])["qtlTCrate"])[6]) / ((QTLeffects["GCN1"])["qtlRDrate"])[6]) * ((Initvar["GCN1"])["R"])[6]), :(((164.491 * ((QTLeffects["GCN2"])["qtlTCrate"])[6]) / ((QTLeffects["GCN2"])["qtlRDrate"])[6]) * ((Initvar["GCN2"])["R"])[6]), :(((164.491 * ((QTLeffects["GCN3"])["qtlTCrate"])[6]) / ((QTLeffects["GCN3"])["qtlRDrate"])[6]) * ((Initvar["GCN3"])["R"])[6]), :(((164.491 * ((QTLeffects["GCN4"])["qtlTCrate"])[6]) / ((QTLeffects["GCN4"])["qtlRDrate"])[6]) * ((Initvar["GCN4"])["R"])[6]), :(((71.5497 * ((QTLeffects["GCN1"])["qtlTCrate"])[7]) / ((QTLeffects["GCN1"])["qtlRDrate"])[7]) * ((Initvar["GCN1"])["R"])[7]), :(((71.5497 * ((QTLeffects["GCN2"])["qtlTCrate"])[7]) / ((QTLeffects["GCN2"])["qtlRDrate"])[7]) * ((Initvar["GCN2"])["R"])[7]), :(((71.5497 * ((QTLeffects["GCN3"])["qtlTCrate"])[7]) / ((QTLeffects["GCN3"])["qtlRDrate"])[7]) * ((Initvar["GCN3"])["R"])[7]), :(((71.5497 * ((QTLeffects["GCN4"])["qtlTCrate"])[7]) / ((QTLeffects["GCN4"])["qtlRDrate"])[7]) * ((Initvar["GCN4"])["R"])[7]), :(((191.649 * ((QTLeffects["GCN1"])["qtlTCrate"])[8]) / ((QTLeffects["GCN1"])["qtlRDrate"])[8]) * ((Initvar["GCN1"])["R"])[8]), :(((191.649 * ((QTLeffects["GCN2"])["qtlTCrate"])[8]) / ((QTLeffects["GCN2"])["qtlRDrate"])[8]) * ((Initvar["GCN2"])["R"])[8]), :(((191.649 * ((QTLeffects["GCN3"])["qtlTCrate"])[8]) / ((QTLeffects["GCN3"])["qtlRDrate"])[8]) * ((Initvar["GCN3"])["R"])[8]), :(((191.649 * ((QTLeffects["GCN4"])["qtlTCrate"])[8]) / ((QTLeffects["GCN4"])["qtlRDrate"])[8]) * ((Initvar["GCN4"])["R"])[8]), :(((79.9315 * ((QTLeffects["GCN1"])["qtlTCrate"])[9]) / ((QTLeffects["GCN1"])["qtlRDrate"])[9]) * ((Initvar["GCN1"])["R"])[9]), :(((79.9315 * ((QTLeffects["GCN2"])["qtlTCrate"])[9]) / ((QTLeffects["GCN2"])["qtlRDrate"])[9]) * ((Initvar["GCN2"])["R"])[9]), :(((79.9315 * ((QTLeffects["GCN3"])["qtlTCrate"])[9]) / ((QTLeffects["GCN3"])["qtlRDrate"])[9]) * ((Initvar["GCN3"])["R"])[9]), :(((79.9315 * ((QTLeffects["GCN4"])["qtlTCrate"])[9]) / ((QTLeffects["GCN4"])["qtlRDrate"])[9]) * ((Initvar["GCN4"])["R"])[9]), :(((293.694 * ((QTLeffects["GCN1"])["qtlTCrate"])[10]) / ((QTLeffects["GCN1"])["qtlRDrate"])[10]) * ((Initvar["GCN1"])["R"])[10]), :(((293.694 * ((QTLeffects["GCN2"])["qtlTCrate"])[10]) / ((QTLeffects["GCN2"])["qtlRDrate"])[10]) * ((Initvar["GCN2"])["R"])[10]), :(((293.694 * ((QTLeffects["GCN3"])["qtlTCrate"])[10]) / ((QTLeffects["GCN3"])["qtlRDrate"])[10]) * ((Initvar["GCN3"])["R"])[10]), :(((293.694 * ((QTLeffects["GCN4"])["qtlTCrate"])[10]) / ((QTLeffects["GCN4"])["qtlRDrate"])[10]) * ((Initvar["GCN4"])["R"])[10]), :(((1.26595e6 * ((QTLeffects["GCN1"])["qtlTCrate"])[1] * ((QTLeffects["GCN1"])["qtlTLrate"])[1]) / (((QTLeffects["GCN1"])["qtlRDrate"])[1] * ((QTLeffects["GCN1"])["qtlPDrate"])[1])) * ((Initvar["GCN1"])["P"])[1]), :(((1.26595e6 * ((QTLeffects["GCN2"])["qtlTCrate"])[1] * ((QTLeffects["GCN2"])["qtlTLrate"])[1]) / (((QTLeffects["GCN2"])["qtlRDrate"])[1] * ((QTLeffects["GCN2"])["qtlPDrate"])[1])) * ((Initvar["GCN2"])["P"])[1]), :(((1.26595e6 * ((QTLeffects["GCN3"])["qtlTCrate"])[1] * ((QTLeffects["GCN3"])["qtlTLrate"])[1]) / (((QTLeffects["GCN3"])["qtlRDrate"])[1] * ((QTLeffects["GCN3"])["qtlPDrate"])[1])) * ((Initvar["GCN3"])["P"])[1]), :(((1.26595e6 * ((QTLeffects["GCN4"])["qtlTCrate"])[1] * ((QTLeffects["GCN4"])["qtlTLrate"])[1]) / (((QTLeffects["GCN4"])["qtlRDrate"])[1] * ((QTLeffects["GCN4"])["qtlPDrate"])[1])) * ((Initvar["GCN4"])["P"])[1]), :(((1.97821e6 * ((QTLeffects["GCN1"])["qtlTCrate"])[2] * ((QTLeffects["GCN1"])["qtlTLrate"])[2]) / (((QTLeffects["GCN1"])["qtlRDrate"])[2] * ((QTLeffects["GCN1"])["qtlPDrate"])[2])) * ((Initvar["GCN1"])["P"])[2]), :(((1.97821e6 * ((QTLeffects["GCN2"])["qtlTCrate"])[2] * ((QTLeffects["GCN2"])["qtlTLrate"])[2]) / (((QTLeffects["GCN2"])["qtlRDrate"])[2] * ((QTLeffects["GCN2"])["qtlPDrate"])[2])) * ((Initvar["GCN2"])["P"])[2]), :(((1.97821e6 * ((QTLeffects["GCN3"])["qtlTCrate"])[2] * ((QTLeffects["GCN3"])["qtlTLrate"])[2]) / (((QTLeffects["GCN3"])["qtlRDrate"])[2] * ((QTLeffects["GCN3"])["qtlPDrate"])[2])) * ((Initvar["GCN3"])["P"])[2]), :(((1.97821e6 * ((QTLeffects["GCN4"])["qtlTCrate"])[2] * ((QTLeffects["GCN4"])["qtlTLrate"])[2]) / (((QTLeffects["GCN4"])["qtlRDrate"])[2] * ((QTLeffects["GCN4"])["qtlPDrate"])[2])) * ((Initvar["GCN4"])["P"])[2]), :(((2.77873e5 * ((QTLeffects["GCN1"])["qtlTCrate"])[3] * ((QTLeffects["GCN1"])["qtlTLrate"])[3]) / (((QTLeffects["GCN1"])["qtlRDrate"])[3] * ((QTLeffects["GCN1"])["qtlPDrate"])[3])) * ((Initvar["GCN1"])["P"])[3]), :(((2.77873e5 * ((QTLeffects["GCN2"])["qtlTCrate"])[3] * ((QTLeffects["GCN2"])["qtlTLrate"])[3]) / (((QTLeffects["GCN2"])["qtlRDrate"])[3] * ((QTLeffects["GCN2"])["qtlPDrate"])[3])) * ((Initvar["GCN2"])["P"])[3]), :(((2.77873e5 * ((QTLeffects["GCN3"])["qtlTCrate"])[3] * ((QTLeffects["GCN3"])["qtlTLrate"])[3]) / (((QTLeffects["GCN3"])["qtlRDrate"])[3] * ((QTLeffects["GCN3"])["qtlPDrate"])[3])) * ((Initvar["GCN3"])["P"])[3]), :(((2.77873e5 * ((QTLeffects["GCN4"])["qtlTCrate"])[3] * ((QTLeffects["GCN4"])["qtlTLrate"])[3]) / (((QTLeffects["GCN4"])["qtlRDrate"])[3] * ((QTLeffects["GCN4"])["qtlPDrate"])[3])) * ((Initvar["GCN4"])["P"])[3]), :(((1.10747e7 * ((QTLeffects["GCN1"])["qtlTCrate"])[4] * ((QTLeffects["GCN1"])["qtlTLrate"])[4]) / (((QTLeffects["GCN1"])["qtlRDrate"])[4] * ((QTLeffects["GCN1"])["qtlPDrate"])[4])) * ((Initvar["GCN1"])["P"])[4]), :(((1.10747e7 * ((QTLeffects["GCN2"])["qtlTCrate"])[4] * ((QTLeffects["GCN2"])["qtlTLrate"])[4]) / (((QTLeffects["GCN2"])["qtlRDrate"])[4] * ((QTLeffects["GCN2"])["qtlPDrate"])[4])) * ((Initvar["GCN2"])["P"])[4]), :(((1.10747e7 * ((QTLeffects["GCN3"])["qtlTCrate"])[4] * ((QTLeffects["GCN3"])["qtlTLrate"])[4]) / (((QTLeffects["GCN3"])["qtlRDrate"])[4] * ((QTLeffects["GCN3"])["qtlPDrate"])[4])) * ((Initvar["GCN3"])["P"])[4]), :(((1.10747e7 * ((QTLeffects["GCN4"])["qtlTCrate"])[4] * ((QTLeffects["GCN4"])["qtlTLrate"])[4]) / (((QTLeffects["GCN4"])["qtlRDrate"])[4] * ((QTLeffects["GCN4"])["qtlPDrate"])[4])) * ((Initvar["GCN4"])["P"])[4]), :(((1.7921e5 * ((QTLeffects["GCN1"])["qtlTCrate"])[5] * ((QTLeffects["GCN1"])["qtlTLrate"])[5]) / (((QTLeffects["GCN1"])["qtlRDrate"])[5] * ((QTLeffects["GCN1"])["qtlPDrate"])[5])) * ((Initvar["GCN1"])["P"])[5]), :(((1.7921e5 * ((QTLeffects["GCN2"])["qtlTCrate"])[5] * ((QTLeffects["GCN2"])["qtlTLrate"])[5]) / (((QTLeffects["GCN2"])["qtlRDrate"])[5] * ((QTLeffects["GCN2"])["qtlPDrate"])[5])) * ((Initvar["GCN2"])["P"])[5]), :(((1.7921e5 * ((QTLeffects["GCN3"])["qtlTCrate"])[5] * ((QTLeffects["GCN3"])["qtlTLrate"])[5]) / (((QTLeffects["GCN3"])["qtlRDrate"])[5] * ((QTLeffects["GCN3"])["qtlPDrate"])[5])) * ((Initvar["GCN3"])["P"])[5]), :(((1.7921e5 * ((QTLeffects["GCN4"])["qtlTCrate"])[5] * ((QTLeffects["GCN4"])["qtlTLrate"])[5]) / (((QTLeffects["GCN4"])["qtlRDrate"])[5] * ((QTLeffects["GCN4"])["qtlPDrate"])[5])) * ((Initvar["GCN4"])["P"])[5]), :(((3.17649e6 * ((QTLeffects["GCN1"])["qtlTCrate"])[7] * ((QTLeffects["GCN1"])["qtlTLrate"])[7]) / (((QTLeffects["GCN1"])["qtlRDrate"])[7] * ((QTLeffects["GCN1"])["qtlPDrate"])[7])) * ((Initvar["GCN1"])["P"])[7]), :(((3.17649e6 * ((QTLeffects["GCN2"])["qtlTCrate"])[7] * ((QTLeffects["GCN2"])["qtlTLrate"])[7]) / (((QTLeffects["GCN2"])["qtlRDrate"])[7] * ((QTLeffects["GCN2"])["qtlPDrate"])[7])) * ((Initvar["GCN2"])["P"])[7]), :(((3.17649e6 * ((QTLeffects["GCN3"])["qtlTCrate"])[7] * ((QTLeffects["GCN3"])["qtlTLrate"])[7]) / (((QTLeffects["GCN3"])["qtlRDrate"])[7] * ((QTLeffects["GCN3"])["qtlPDrate"])[7])) * ((Initvar["GCN3"])["P"])[7]), :(((3.17649e6 * ((QTLeffects["GCN4"])["qtlTCrate"])[7] * ((QTLeffects["GCN4"])["qtlTLrate"])[7]) / (((QTLeffects["GCN4"])["qtlRDrate"])[7] * ((QTLeffects["GCN4"])["qtlPDrate"])[7])) * ((Initvar["GCN4"])["P"])[7])],"species"=>Any["R1GCN1", "R1GCN2", "R1GCN3", "R1GCN4", "R2GCN1", "R2GCN2", "R2GCN3", "R2GCN4", "R3GCN1", "R3GCN2", "R3GCN3", "R3GCN4", "R4GCN1", "R4GCN2", "R4GCN3", "R4GCN4", "R5GCN1", "R5GCN2", "R5GCN3", "R5GCN4", "R6GCN1", "R6GCN2", "R6GCN3", "R6GCN4", "R7GCN1", "R7GCN2", "R7GCN3", "R7GCN4", "R8GCN1", "R8GCN2", "R8GCN3", "R8GCN4", "R9GCN1", "R9GCN2", "R9GCN3", "R9GCN4", "R10GCN1", "R10GCN2", "R10GCN3", "R10GCN4", "P1GCN1", "P1GCN2", "P1GCN3", "P1GCN4", "P2GCN1", "P2GCN2", "P2GCN3", "P2GCN4", "P3GCN1", "P3GCN2", "P3GCN3", "P3GCN4", "P4GCN1", "P4GCN2", "P4GCN3", "P4GCN4", "P5GCN1", "P5GCN2", "P5GCN3", "P5GCN4", "P7GCN1", "P7GCN2", "P7GCN3", "P7GCN4"])

QTLeffects = Dict("GCN3"=>Dict("qtlPDregbind"=>[1.0, 0.903064, 1.0, 1.01958, 0.908109, 0.0, 0.985851, 0.0, 0.0, 0.0],"qtlTLrate"=>[1.0, 0.98387, 1.12259, 0.947592, 1.02447, 0.0, 0.897265, 0.0, 0.0, 0.0],"qtlTCrate"=>[1.0, 1.17192, 1.0, 0.973951, 1.07116, 1.11191, 1.0003, 1.04185, 0.96825, 1.0],"qtlRDbindreg"=>[1.0, 0.969047, 0.883935, 1.16487, 1.18382, 0.926073, 1.13709, 1.16358, 1.20015, 1.0],"qtlTCregbind"=>[1.0, 0.939463, 1.0, 0.996446, 1.11144, 1.0, 0.92738, 1.0, 1.11756, 1.0],"qtlactivity"=>[1.0, 1.07163, 1.0, 0.948701, 1.00975, 1.0, 1.00126, 0.875307, 1.03423, 1.0],"qtlRDrate"=>[1.0, 1.0, 1.0, 1.08387, 1.02613, 1.0, 1.20753, 1.06082, 1.0, 1.0],"qtlTLregbind"=>[1.0, 0.969293, 1.0, 0.9764, 0.932668, 0.0, 0.918416, 0.0, 0.0, 0.0],"qtlPDrate"=>[1.0, 1.14071, 1.0, 0.93931, 1.02331, 0.0, 1.05153, 0.0, 0.0, 0.0]),"GCN2"=>Dict("qtlPDregbind"=>[1.0, 0.903064, 1.0, 0.904002, 1.0, 0.0, 0.947351, 0.0, 0.0, 0.0],"qtlTLrate"=>[0.888289, 0.98387, 1.0, 0.633331, 0.927585, 0.0, 1.03055, 0.0, 0.0, 0.0],"qtlTCrate"=>[1.0, 1.17192, 1.0, 1.12524, 1.0, 1.0, 0.993596, 0.99752, 1.0, 1.0],"qtlRDbindreg"=>[1.0, 0.969047, 1.0, 1.04387, 1.0, 0.888253, 1.0, 0.843412, 1.0, 1.0],"qtlTCregbind"=>[1.0, 0.939463, 1.0, 1.15848, 1.0, 0.919143, 1.0, 0.813784, 0.975534, 1.0],"qtlactivity"=>[1.0, 1.07163, 1.0, 0.864578, 1.0, 1.13318, 0.918793, 1.08566, 1.0, 1.0],"qtlRDrate"=>[1.0, 1.0, 1.0, 1.01624, 1.0, 1.0, 1.0, 0.926551, 0.941337, 1.0],"qtlTLregbind"=>[1.0, 0.969293, 1.0, 1.05993, 1.0, 0.0, 0.992744, 0.0, 0.0, 0.0],"qtlPDrate"=>[1.0, 1.14071, 1.0, 1.00784, 1.0, 0.0, 0.990836, 0.0, 0.0, 0.0]),"GCN4"=>Dict("qtlPDregbind"=>[1.0, 0.903064, 1.0, 1.0, 1.0, 0.0, 0.985851, 0.0, 0.0, 0.0],"qtlTLrate"=>[0.888289, 0.98387, 1.04376, 1.0, 1.0, 0.0, 0.897265, 0.0, 0.0, 0.0],"qtlTCrate"=>[1.0, 1.17192, 1.0, 1.0, 1.0, 1.0, 1.0003, 1.0, 0.934712, 1.09404],"qtlRDbindreg"=>[1.0, 0.969047, 1.08212, 1.0, 1.0, 1.0, 1.13709, 1.0, 0.953415, 1.08835],"qtlTCregbind"=>[1.0, 0.939463, 1.0, 1.0, 1.0, 1.0, 0.92738, 1.0, 1.0, 0.874129],"qtlactivity"=>[1.0, 1.07163, 0.985766, 1.0, 1.0, 1.0, 1.00126, 1.0, 0.87951, 0.774213],"qtlRDrate"=>[1.0, 1.0, 0.923331, 1.0, 1.0, 1.0, 1.20753, 1.0, 1.0, 1.0],"qtlTLregbind"=>[1.0, 0.969293, 0.947754, 0.790886, 1.0, 0.0, 0.918416, 0.0, 0.0, 0.0],"qtlPDrate"=>[1.0, 1.14071, 1.0, 1.0, 1.0, 0.0, 1.05153, 0.0, 0.0, 0.0]),"GCN1"=>Dict("qtlPDregbind"=>[1.15029, 1.16139, 1.0, 1.0, 0.965733, 0.0, 1.0, 0.0, 0.0, 0.0],"qtlTLrate"=>[1.0, 0.961051, 1.12259, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0],"qtlTCrate"=>[1.0, 1.24792, 1.0, 0.894408, 1.0, 1.0, 1.0, 0.99752, 0.96825, 1.0],"qtlRDbindreg"=>[1.0, 0.83549, 0.883935, 1.07553, 1.0, 0.870481, 1.0, 0.843412, 1.20015, 1.0],"qtlTCregbind"=>[1.0, 1.08633, 1.0, 1.06498, 1.03092, 1.0, 1.0, 0.813784, 1.11756, 1.0],"qtlactivity"=>[1.0, 0.791183, 1.0, 1.0, 0.948967, 1.0, 1.0, 1.08566, 1.03423, 1.0],"qtlRDrate"=>[1.0, 0.818898, 1.0, 1.16617, 1.09283, 1.0, 1.0, 0.926551, 1.0, 1.0],"qtlTLregbind"=>[1.0, 0.9524, 1.0, 0.970783, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0],"qtlPDrate"=>[1.11645, 0.940561, 1.0, 0.949767, 0.837376, 0.0, 1.0, 0.0, 0.0, 0.0]))
Initvar = Dict("GCN3"=>Dict("P"=>[0.922647, 0.827737, 0.962663, 1.07635, 0.848211, 1.01387, 1.15418, 0.997701, 0.989415, 1.06221],"R"=>[1.01076, 0.903517, 0.931234, 0.948392, 1.08051, 1.07374, 0.843925, 0.956626, 1.1761, 1.09528]),"GCN2"=>Dict("P"=>[0.853842, 0.93087, 1.12092, 0.964489, 1.08834, 1.05661, 0.833048, 1.06602, 1.16125, 1.10232],"R"=>[0.97554, 0.932472, 1.0511, 0.898645, 1.05261, 1.05966, 0.7716, 0.90272, 0.807083, 1.06358]),"GCN4"=>Dict("P"=>[0.990594, 1.09941, 1.0421, 0.930396, 0.81033, 0.90876, 0.86762, 1.13135, 0.93453, 1.04141],"R"=>[1.00492, 1.1564, 1.14431, 0.93374, 0.970863, 1.13996, 1.1174, 0.978497, 1.10544, 1.00269]),"GCN1"=>Dict("P"=>[1.07192, 1.15348, 1.06973, 0.926324, 1.04767, 0.955735, 1.08041, 0.796943, 0.988385, 0.919596],"R"=>[0.855099, 0.818183, 0.857508, 0.910632, 1.00853, 0.99178, 0.934635, 0.879292, 0.843702, 0.908478]))

model = Network("coucou")

for i in 1:length(stochmodel["species"])
    i0 = round(Int, eval(stochmodel["initialconditions"][i]))
    #println(stochmodel["species"][i]* "\t"*string(i0))
    model <= BioSimulator.Species(stochmodel["species"][i], i0)
end

for i in eachindex(stochmodel["reactions"])
    #println(join([stochmodel["reactionsnames"][i], eval(stochmodel["propensities"][i]), stochmodel["reactions"][i]], "\t", "\n"))
    model <= BioSimulator.Reaction(stochmodel["reactionsnames"][i], eval(stochmodel["propensities"][i]), stochmodel["reactions"][i])
end

result = simulate(model, algorithm=SSA, time = 1000.0, epochs = 1000, trials = 1)


test = res2df(result)
plot(test[:, :time], test[:, :R1GCN1])


genes = collect(1:10)
gcnList = ["GCN1", "GCN2", "GCN3", "GCN4"]

for g in genes
##    toplotRNA = [test[t, Symbol("R"*string(g).*gcn)] for t in 1:size(test)[1], gcn in gcnList]
##    toplotProt = [test[t, Symbol("P"*string(g).*gcn)] for t in 1:size(test)[1], gcn in gcnList]
##    gui(plot(test[:, :time], hcat(toplotRNA, toplotProt), show = true, reuse = false))

    toplotRNA = [test[t, Symbol("R"*string(g).*gcn)] for t in 1:size(test)[1], gcn in gcnList]
    if in(g, [1,2,3,4,5,7])
        toplotProt = [test[t, Symbol("P"*string(g).*gcn)] for t in 1:size(test)[1], gcn in gcnList]
        gui(Plots.plot(test[:, :time], hcat(toplotRNA, toplotProt)[:, [1,5,2,6,3,7,4,8]], layout = @layout([a b]), show = true, reuse = false))
    else
        gui(Plots.plot(test[:, :time], toplotRNA, show = true, reuse = false))
    end
end



#############################################################################



QTLeffects = Dict("GCN3"=>Dict("qtlPDregbind"=>[0.0, 1.0, 0.996807, 0.868407, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],"qtlTLrate"=>[0.0, 1.06523, 1.04637, 0.996554, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],"qtlTCrate"=>[1.00504, 1.0, 1.0, 1.0519, 1.0, 0.939048, 1.0, 1.0, 1.0, 1.0],"qtlRDbindreg"=>[1.01377, 1.0, 1.15503, 1.0, 1.16826, 1.01473, 1.0, 1.0, 1.0, 1.0],"qtlTCregbind"=>[0.980469, 1.0, 0.957135, 1.14658, 1.0, 0.987476, 1.0, 1.0, 1.0, 1.0],"qtlactivity"=>[0.914794, 1.0, 0.949882, 0.923203, 1.0, 1.01069, 1.0, 1.0, 1.0, 1.0],"qtlRDrate"=>[0.921899, 1.0, 0.988992, 0.899535, 1.0, 0.912394, 1.0, 1.0, 1.0, 1.0],"qtlTLregbind"=>[0.0, 1.0, 0.933088, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],"qtlPDrate"=>[0.0, 1.0, 1.00773, 1.0762, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]),"GCN2"=>Dict("qtlPDregbind"=>[0.0, 1.0, 0.996807, 0.868407, 0.0, 0.0, 1.0, 0.97277, 1.0, 1.0],"qtlTLrate"=>[0.0, 1.06523, 1.04637, 0.996554, 0.0, 0.0, 1.0, 1.12512, 1.0, 1.0],"qtlTCrate"=>[1.00504, 1.0, 1.0, 1.0519, 1.09237, 1.08641, 1.0, 0.901656, 1.04507, 1.0],"qtlRDbindreg"=>[1.01377, 1.0, 1.15503, 1.0, 0.959294, 0.818953, 1.0, 1.0482, 1.0, 1.0],"qtlTCregbind"=>[0.980469, 1.0, 0.957135, 1.14658, 1.22001, 1.02222, 1.0, 1.15073, 1.0, 1.0],"qtlactivity"=>[0.914794, 1.0, 0.949882, 0.923203, 1.07234, 1.01743, 1.0, 0.974383, 1.0, 1.0],"qtlRDrate"=>[0.921899, 1.0, 0.988992, 0.899535, 1.02111, 1.04598, 1.0, 1.06141, 1.01365, 1.0],"qtlTLregbind"=>[0.0, 1.0, 0.933088, 1.0, 0.0, 0.0, 1.0, 1.07071, 1.0, 1.0],"qtlPDrate"=>[0.0, 1.0, 1.00773, 1.0762, 0.0, 0.0, 1.0, 1.03252, 1.0, 1.0]),"GCN4"=>Dict("qtlPDregbind"=>[0.0, 1.0, 1.13171, 0.868407, 0.0, 0.0, 1.05037, 1.0, 1.0, 1.0],"qtlTLrate"=>[0.0, 1.0, 0.969644, 0.996554, 0.0, 0.0, 1.07574, 1.0, 0.956107, 1.0],"qtlTCrate"=>[1.0, 1.0, 1.0, 1.0519, 1.09237, 0.939048, 0.922822, 1.0, 1.0, 1.0],"qtlRDbindreg"=>[1.0, 1.0, 1.0, 1.0, 0.959294, 1.01473, 0.927262, 1.0, 0.853407, 1.0],"qtlTCregbind"=>[1.0, 1.0, 1.0, 1.14658, 1.22001, 0.987476, 0.931289, 1.0, 1.0, 1.0],"qtlactivity"=>[1.0, 1.0, 1.09505, 0.923203, 1.07234, 1.01069, 0.902183, 1.0, 1.03853, 1.0],"qtlRDrate"=>[1.0, 1.0, 1.0, 0.899535, 1.02111, 0.912394, 0.836054, 1.0, 1.0, 1.0],"qtlTLregbind"=>[0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],"qtlPDrate"=>[0.0, 1.0, 1.03975, 1.0762, 0.0, 0.0, 0.92047, 1.0, 1.0, 1.0]),"GCN1"=>Dict("qtlPDregbind"=>[0.0, 1.0, 1.13171, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],"qtlTLrate"=>[0.0, 1.05908, 0.969644, 1.06483, 0.0, 0.0, 1.0, 0.796604, 1.0, 1.0],"qtlTCrate"=>[1.08969, 1.0, 1.0, 1.0, 1.14363, 1.04624, 1.0, 1.01885, 1.04507, 1.0],"qtlRDbindreg"=>[1.08142, 0.991458, 1.0, 1.0, 0.959343, 1.0, 1.0, 0.980113, 1.0, 1.0],"qtlTCregbind"=>[1.0, 1.0, 1.0, 0.778718, 1.0, 1.01828, 1.0, 1.15897, 1.0, 1.0],"qtlactivity"=>[0.94031, 1.0, 1.09505, 0.946715, 1.0, 1.0, 1.0, 0.945421, 1.0, 1.0],"qtlRDrate"=>[1.0, 1.09692, 1.0, 1.0, 1.0, 1.04338, 1.0, 0.896263, 1.01365, 1.0],"qtlTLregbind"=>[0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.05649, 1.03906, 1.0, 1.0],"qtlPDrate"=>[0.0, 1.0, 1.03975, 1.0, 0.0, 0.0, 1.0, 1.01075, 1.0, 1.0]))

InitVar = Dict("GCN3"=>Dict("P"=>[0.80262, 0.96886, 0.928253, 0.922065, 0.879601, 0.990862, 1.05296, 0.888656, 1.17134, 1.09654],"R"=>[0.904391, 1.08354, 1.1342, 1.06773, 0.962775, 1.08633, 1.09988, 1.04768, 1.09627, 0.980728]),"GCN2"=>Dict("P"=>[1.03636, 0.852727, 0.907479, 1.13289, 0.689987, 1.1046, 1.072, 1.10135, 1.00242, 0.848639],"R"=>[0.889878, 1.03319, 1.03879, 0.960242, 0.953933, 0.95715, 0.939481, 0.883476, 0.851746, 1.15431]),"GCN4"=>Dict("P"=>[1.07095, 1.09463, 1.1466, 0.954728, 0.901106, 1.17743, 1.14382, 0.993785, 1.14216, 0.974802],"R"=>[1.0076, 0.996241, 1.21269, 0.931376, 1.0446, 0.95608, 0.964789, 1.12328, 1.05084, 1.08357]),"GCN1"=>Dict("P"=>[1.09452, 1.07919, 0.866409, 1.09506, 1.03763, 1.09511, 0.907524, 0.962947, 0.867508, 1.17224],"R"=>[1.16395, 0.96986, 0.988619, 0.992444, 0.940217, 0.985759, 0.777827, 0.929596, 0.832051, 0.883847]))

nod = Dict{String,Any}(Pair{String,Any}("TargetReaction", String["TL", "TL", "TC", "TC", "TL", "RD", "TC", "TL", "TL", "TL"]),Pair{String,Any}("TCrate", [0.0876407, 0.0590708, 0.0874417, 0.017104, 0.0708739, 0.0741194, 0.0804994, 0.0102489, 0.0217458, 0.0622714]),Pair{String,Any}("RDrate", [0.000417537, 0.000398565, 0.000361141, 0.00048852, 0.000350877, 0.000358551, 0.00458716, 0.00129199, 0.000702741, 0.000365364]),Pair{String,Any}("id", [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),Pair{String,Any}("ActiveForm", String["R1", "Pm2", "Pm3", "P4", "R5", "R6", "Pm7", "Pm8", "P9", "Pm10"]),Pair{String,Any}("coding", String["NC", "PC", "PC", "PC", "NC", "NC", "PC", "PC", "PC", "PC"]),Pair{String,Any}("TLrate", [0.0, 4.70877, 2.47364, 3.30288, 0.0, 0.0, 4.7373, 2.23753, 1.76261, 3.7403]),Pair{String,Any}("PTMform", String["0", "1", "1", "0", "0", "0", "1", "1", "0", "1"]),Pair{String,Any}("PDrate", [0.0, 8.8739e-5, 0.000151953, 7.41565e-5, 0.0, 0.0, 8.75197e-5, 9.34842e-5, 8.85897e-5, 0.000114916]))

simtime = 100
nepochs = -1
ntrials = 1
simalgorithm = "SSA"
modelname = "coucou"



genes = collect(1:10)
gcnList = ["GCN1", "GCN2", "GCN3", "GCN4"]

for g in genes
##    toplotRNA = [abundancedf[t, Symbol("R"*string(g).*gcn)] for t in 1:size(abundancedf)[1], gcn in gcnList]
##    toplotProt = [abundancedf[t, Symbol("P"*string(g).*gcn)] for t in 1:size(abundancedf)[1], gcn in gcnList]
##    gui(plot(abundancedf[:, :time], hcat(toplotRNA, toplotProt), show = true, reuse = false))

    toplotRNA = [abundancedf[t, Symbol("R"*string(g).*gcn)] for t in 1:size(abundancedf)[1], gcn in gcnList]
    if nod["coding"][g] == "PC"
        toplotProt = [abundancedf[t, Symbol("P"*string(g).*gcn)] for t in 1:size(abundancedf)[1], gcn in gcnList]
        gui(Plots.plot(abundancedf[:, :time], hcat(toplotRNA, toplotProt)[:, [1,5,2,6,3,7,4,8]], layout = @layout([a b]), show = true, reuse = false, title = "Gene "*string(g)))

    else
        gui(Plots.plot(abundancedf[:, :time], toplotRNA, show = true, reuse = false, title = "Gene "*string(g)))
    end

end


#############################################################################

using JLD
using BioSimulator
mechant = load("/home/oangelin/Documents/errmodel.jld")
QTLeffects = mechant["QTLeffects"]
InitVar = mechant["InitVar"]
stochmodel = mechant["stochmodel"]


gentil = load("/home/oangelin/Documents/goodmodel.jld")
QTLeffects = gentil["QTLeffects"]
InitVar = gentil["InitVar"]
stochmodel = gentil["stochmodel"]


model = BioSimulator.Network("coucouBad")


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
    reacname = replace(stochmodel["reactionsnames"][i], "de-", "de")
    #model <= BioSimulator.Reaction(stochmodel["reactionsnames"][i], prop, stochmodel["reactions"][i])
    model <= BioSimulator.Reaction(reacname, prop, stochmodel["reactions"][i])

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


#result = simulate(model, algorithm = SSA, time = 10.0, epochs = 1, trials = 1)



time = 10.0
epochs = 2
trials = 1
algorithm = SSA
D = Array
T = SSA

# extract model information
c = n_species(model)
d = n_reactions(model)

species = species_list(model)
reactions = reaction_list(model)

# create simulation data structures
x0, rxn, id, id2ind = BioSimulator.make_datastructs(species, reactions, c, d)
xt     = copy(x0)
alg    = algorithm(time)

output = BioSimulator.SimData(
    id2ind,
    linspace(0.0, time, epochs + 1),
    #SharedArray{eltype(x0)}(c, epochs + 1, trials),
    D{eltype(x0)}(c, epochs + 1, trials),
    alg.stats
)

BioSimulator.init!(alg, xt, rxn)

#simulate_wrapper!(output, xt, x0, alg, rxn)

data = output.data

i = 1
len = div(size(data, 3), Threads.nthreads())
domain = ((i-1)*len+1):i*len

#simulate_chunk!(output, xt, x0, alg, rxn, domain)
# here Xt = xt, X0 = x0, algorithm = alg, reactions = rxn, trial_set = domain

Xt = xt
X0 = x0


a = BioSimulator.propensities(rxn)

# for trial in trial_set
trial = 1

copy!(Xt, X0)

BioSimulator.update_all_propensities!(a, rxn, Xt)
BioSimulator.reset!(alg, a)

#simulate!(output, Xt, algorithm, reactions, trial)

epoch = 1



#props = zeros(Float64, length(a), 10000)
count = 0
tic(); while !done(alg)
#tic(); for count in 1:size(props, 2)
    t     = BioSimulator.get_time(alg)
    epoch = BioSimulator.update!(output, Xt, t, epoch, trial)

    BioSimulator.step!(alg, Xt, rxn)

#    props[:, count] = BioSimulator.propensities(rxn)

    count = count+1
    if count == 100000
        count = 0
        println(t)
        println(findmax(BioSimulator.propensities(rxn).cache))
    end

end; toc()


#BioSimulator.step!(alg, Xt, rxn)
