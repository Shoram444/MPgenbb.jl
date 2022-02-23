"""
    getParticleString(PART_TYPE, T, MASS, ϕ, θ)

Description of getParticleString
--------------------------------
Returns one-line string for particle in form "TYPE pₓ p𝔂 p𝐳 0"

example:

"electron px = 1.0 py = 1.0 pz = 1.0 dt = 0"

"3 1.0 1.0 1.0 0"
"""
function getParticleString(PART_TYPE, T, MASS, ϕ, θ)
    pMag = getpMag(T, MASS)
    px   = getpx(pMag, ϕ, θ)
    py   = getpy(pMag, ϕ, θ)
    pz   = getpz(pMag,  θ)

    return string(PART_TYPE)*" "*string(px)*" "*string(py)*" "*string(pz)*" 0"
end

"""
    getHeaderString(ID)

Description of getHeaderString
------------------------------
Returns one-line string for header in form "ID 0 2" (2 stands for number of particles)

example:

"ID = 0 t = 0 particles = 2"

"0 0 2"
"""
function getHeaderString(ID)
    return string(ID)*" 0 2"
end


"""
    getEventString(ID, PART_TYPE, T, MASS, ϕ, θ)

Description of getEventString
------------------------------
Returns one event, separated at the end by #.
"""
function getEventString(ID, PART_TYPE, T, MASS, ϕ, θ)
    sHeader= getHeaderString(ID)
    sPart1 = getParticleString(PART_TYPE, T[1], MASS, ϕ[1], θ[1])
    sPart2 = getParticleString(PART_TYPE, T[2], MASS, ϕ[2], θ[2])

    return sHeader*"\n"*sPart1*"\n"*sPart2*"\n#\n"
end
