module MPGenbb

using Random
using Distributions
using Documenter
using DataFrames

include("MomentumFuncs.jl")
export  getpMag, 
        getpx, 
        getpy, 
        getpz

include("StringFuncs.jl")
export getParticleString, getHeaderString, getEventString

include("SamplingFuncs.jl")
export sampleEnergy, samplePhi, sampleTheta, sample_energies


end # module MPgenbb