"""
    sampleEnergy(Q_VALUE)

Description of sampleEnergy
------------------------------
Returns (E₁, E₂) tuple of electron energies in MeV.

Currently used pdf is Beta(2,4) - see Distributions.Beta

Sum of electron kinetic energies is equal Q-value. => (E₁ is random, E₂ = Q - E₁)

"""
function sampleEnergy(Q_VALUE)             # Currently sampling random energy from Beta distribution. Should Change when spectra available
    T = rand(Q_VALUE*Distributions.Beta(2,4),1)[1]

    return (T, Q_VALUE-T)
end


"""
    samplePhi()

Description of samplePhi
------------------------------
Returns (ϕ₁, ϕ₂) tuple in radians Uniformly distributed. See Distributions.Uniform

"""
function samplePhi()
    return 2*π*rand(Distributions.Uniform(),2)
end

"""
    sampleTheta()

Description of sampleTheta
------------------------------
Returns (θ₁, θ₂) tuple in radians Uniformly distributed. See Distributions.Uniform

Sampling equation is:
``` acos(1 - 2*U) ; where U is uniform random number; U ∈ (0,1). ```

"""
function sampleTheta()
    return acos.(1 .-2*rand(Distributions.Uniform(),2))
end
