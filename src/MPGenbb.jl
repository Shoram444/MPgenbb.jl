module MPGenbb

using Random
using Distributions
using Documenter
using DataFrames

include("MomentumFuncs.jl")
export  get_pMag, 
        get_px, 
        get_py, 
        get_pz

include("StringFuncs.jl")
export get_particle_string, get_header_string, get_event_string

include("SamplingFuncs.jl")
export sample_phi, sample_theta, sample_energies


end # module MPgenbb