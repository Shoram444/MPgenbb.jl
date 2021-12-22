module MPgenbb

using Random
using Distributions
using Documenter

include("MomentumFuncs.jl")
export getpMag, getpx, getpy, getpz

include("StringFuncs.jl")
export getParticleString, getHeaderString, getEventString

include("SamplingFuncs.jl")
export sampleEnergy, samplePhi, sampleTheta


end # module MPgenbb