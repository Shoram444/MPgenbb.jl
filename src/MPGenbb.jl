module MPGenbb

using Random
using Distributions
using Documenter
using DataFrames

include("MomentumFuncs.jl")
export  get_pMag,            # returns the magnitude of the momentum vector calculated from kinetic energy T
        get_px,              # x-component of momentum vector
        get_py,              # y-component of momentum vector
        get_pz,              # z-component of momentum vector
        get_first_vector,    # returns the momentum vector of the first electron (independent, direction generated uniformaly)
        get_second_vector,   # returns the momentum vector of the second electron (dependent on first, direction generated based on θdif distribuition)
        energy_to_momentum   # returns the calculated kinetic energy from momentum

include("StringFuncs.jl")
export get_particle_string,  
       get_header_string, 
       get_event_string

       
include("SamplingFuncs.jl")
export sample_phi,              # returns ϕ azimuthal angle uniformly
       sample_theta,            # returns θ polar angle uniformly
       sample_energies,         # returns !!tuple!! of electron kintetic energies T1, T2
       sample_theta_dif,        # returns θdif as the angle between the 2 electrons (based on 1 - cos(θ) angular distribution)
       solve_quadratic,         # returns x1, x2 of a quadratic equation 
       sample_energies_uniform, # returns tuple of electron kinetic energies T1, T2 given by uniform distribution
       sample_theta_dif_uniform


end # module MPgenbb