"""
    get_pMag(T::Real, MASS::Real)

Description of get_pMag
--------------------------
Computes magnitude of momentum as:
``` |p| = √((MASS + T)² - MASS²) ```.

Input arguments are:
+ T    : kinetic energy in [MeV]
+ MASS : particle mass in [MeV/c²]
"""
function get_pMag(T::Real, MASS::Real)
	return sqrt( (MASS+T)^2 - MASS^2)
end

"""
    get_px(pMag, ϕ, θ)

Description of get_px
--------------------------
Computes x-component of the momentum as:
``` pₓ = |p|*cos(ϕ)*sin(θ) ```

Input arguments are:
+ pMag : momentum magnitude in [MeV]
+ ϕ    : azimuth angle in [rad]
+ θ    : polar angle in [rad]
"""
function get_px(pMag, ϕ, θ)
    return pMag*cos(ϕ)*sin(θ)
end

"""
    get_py(pMag, ϕ, θ)

Description of get_py
--------------------------
Computes y-component of the momentum as:
``` p𝔂 = |p|*sin(ϕ)*sin(θ) ```

Input arguments are:
+ pMag : momentum magnitude in [MeV]
+ ϕ    : azimuth angle in [rad]
+ θ    : polar angle in [rad]
"""
function get_py(pMag, ϕ, θ)
    return pMag*sin(ϕ)*sin(θ)
end

"""
    ```get_pz(pMag, θ)```

Description of get_pz
--------------------------
Computes z-component of the momentum as:
`` p𝐳 = |p|*cos(θ) ``

Input arguments are:
+ pMag : momentum magnitude in [MeV]
+ θ    : polar angle in [rad]
"""
function get_pz(pMag, θ)
    return pMag*cos(θ)
end



"""
```get_first_vector(T::Real, MASS::Real)```

Description of ```get_first_vector```
-------------------------------
Returns the momentum vector of the first electron. Momentum magnitude is calculated from the electron kinetic energy (T). 
The direction is generated uniformly in a sphere.

The momentum components are (in order):
    + pMag - Magnitude of the momentum vector (given by the electron Energy)
    + px   - x-component
    + py   - y-component
    + pz   - z-component
    + θ    - polar angle  (used in calculation of the second momentum vector)
    + ϕ    - azimuthal angle (used in calculation of the second momentum vector)
"""
function get_first_vector(T::Real, MASS::Real)
    pMag = get_pMag(T, MASS)
    ϕ = sample_phi()
    θ = sample_theta()
    px = get_px(pMag, ϕ, θ)
    py = get_py(pMag, ϕ, θ)
    pz = get_pz(pMag, θ)
    return [pMag, px, py, pz, θ, ϕ]
end


"""
```get_pPrime(T::Real, MASS::Real)```

Description of get_pPrime
-------------------------------
Returns the momentum (prime) vector of the second electron. The vector is generated as a random point from a circle projected onto the x-y plane. 
The circle radius is given by ``r = \\sin(θdif)|p1|``, where ``\\theta_{dif}`` is the generated angle between the two electrons and ``|p1|`` 
is the magnitude of the momentum vector of the first electron. The process assumes p1 to be directed along the z-axis (hence projection onto x-y plane).

The vector components are (in order):
    + ``\\sin(\\theta_{dif})\\cos(\\phi_2)pMag`` - x-component
    + ``\\sin(\\theta_{dif})\\sin(\\phi_2)pMag`` - y-component
    + ``\\cos(\\theta_{dif})pMag``               - z-component


```get_pPrime(T::Real, MASS::Real, θdif)```

Generates momentum (prime) vector of the second electron, where θdif is predefined. 

"""
function get_pPrime(T::Real, MASS::Real)
    pMag = get_pMag(T, MASS)
    θdif = sample_theta_dif()
    ϕ2   = sample_phi()
    
    return [ sin(θdif)*cos(ϕ2)*pMag ; sin(θdif)*sin(ϕ2)*pMag ; cos(θdif)*pMag ] # returns column vector of p2'
end

function get_pPrime(T::Real, MASS::Real, θdif)
    pMag = get_pMag(T, MASS)
    ϕ2   = sample_phi()
    
    return [ sin(θdif)*cos(ϕ2)*pMag ; sin(θdif)*sin(ϕ2)*pMag ; cos(θdif)*pMag ] # returns column vector of p2'
end



"""
```get_second_vector(T::Real, MASS::Real, p1::DataFrameRow) ```   p1 given as dataframe row 
```get_second_vector(T::Real, MASS::Real, p1::Vector{<:Real})``` p1 given as vector

Description of get_second_vector
--------------------------------
Returns the momentum vector of the second electron. Momentum magnitude is calculated from the electron kinetic energy (T). 
The direction is generated from sampling θdif angle (see ```sample_theta_dif()```). In order to calculate the direction: first a p2Prime vector is generated
(see ```get_pPrime(T::Real, MASS::Real)```) by sampling uniformly from a circle projected onto x-y plane, 
where the circle is determined by the p1 vector directed along z-axis. Afterwards, p2Prime is rotated along 1. y-axis, 2. z-axis to obtain p2 in the direction 
of a cone around p1, where the cone is given by θdif. The equation for p2 is:

`` p2 = Rz*Ry*p2Prime ``

The momentum components are (in order):
    + pMag - Magnitude of the momentum vector (given by the electron Energy)
    + px   - x-component
    + py   - y-component
    + pz   - z-component
    + θ    - polar angle  (used in calculation of the second momentum vector)
    + ϕ    - azimuthal angle (used in calculation of the second momentum vector)
"""
function get_second_vector(T::Real, MASS::Real, p1::DataFrameRow)
    p2Prime = get_pPrime(T, MASS) # vector generated by randomly chosing point on circle projected onto xy plane

    p2 = [ 
        cos(p1.phi)*( p2Prime[1]*cos(p1.theta ) + p2Prime[3]*sin(p1.theta) ) - p2Prime[2]*sin(p1.phi)
        sin(p1.phi)*( p2Prime[1]*cos(p1.theta ) + p2Prime[3]*sin(p1.theta) ) + p2Prime[2]*cos(p1.phi)
        p2Prime[3]*cos(p1.theta ) - p2Prime[1]*sin(p1.theta )
    ]
    
    p2Mag = get_pMag(T, MASS)   
    θ2 = acos(p2[3]/p2Mag)
    ϕ2 = atan(p2[2],p2[1])

    if (p2[1] == NaN)
        get_second_vector(T, MASS, p1)
    end
    
    return [ p2Mag, p2[1], p2[2], p2[3], θ2, ϕ2 ]
end

function get_second_vector(T::Real, MASS::Real, p1::Vector{<:Real})
    p2Prime = get_pPrime(T, MASS) # vector generated by randomly chosing point on circle projected onto xy plane

    p2 = [ 
        cos(p1[6] )*( p2Prime[1]*cos(p1[5] ) + p2Prime[3]*sin(p1[5]) ) - p2Prime[2]*sin(p1[6])
        sin(p1[6] )*( p2Prime[1]*cos(p1[5] ) + p2Prime[3]*sin(p1[5]) ) + p2Prime[2]*cos(p1[6])
        p2Prime[3]cos(p1[5] ) - p2Prime[1]*sin(p1[5] )
    ]
    
    p2Mag = get_pMag(T, MASS)   
    θ2 = acos(p2[3]/p2Mag)
    ϕ2 = atan(p2[2],p2[1])
    
    if (p2[1] == NaN)
        get_second_vector(T, MASS, p1)
    end

    return [ p2Mag, p2[1], p2[2], p2[3], θ2, ϕ2 ]
end

function get_second_vector(T::Real, MASS::Real, p1::DataFrameRow, θdif::Real)
    p2Prime = get_pPrime(T, MASS, θdif) # vector generated by randomly chosing point on circle projected onto xy plane
    
    p2 = [ 
        cospi(p1.phi /π)*( p2Prime[1]*cospi(p1.theta /π) + p2Prime[3]*sinpi(p1.theta/π) ) - p2Prime[2]*sinpi(p1.phi/π)
        sinpi(p1.phi /π)*( p2Prime[1]*cospi(p1.theta /π) + p2Prime[3]*sinpi(p1.theta/π) ) + p2Prime[2]*cospi(p1.phi/π)
        p2Prime[3]cospi(p1.theta /π) - p2Prime[1]*sinpi(p1.theta /π)
    ]
    
    p2Mag = get_pMag(T, MASS)
    θ2 = acos(p2[3]/p2Mag)
    ϕ2 = atan(p2[2],p2[1])
    
    if (p2[1] == NaN)
        get_second_vector(T, MASS, p1, θdif)
    end

    return [ p2Mag, p2[1], p2[2], p2[3], θ2, ϕ2 ]
end

function get_second_vector(T::Real, MASS::Real, p1::Vector{<:Real},  θdif::Real)
    p2Prime = get_pPrime(T, MASS, θdif) # vector generated by randomly chosing point on circle projected onto xy plane

    p2 = [ 
        cospi(p1[6] /π)*( p2Prime[1]*cospi(p1[5] /π) + p2Prime[3]*sinpi(p1[5]/π) ) - p2Prime[2]*sinpi(p1[6]/π)
        sinpi(p1[6] /π)*( p2Prime[1]*cospi(p1[5] /π) + p2Prime[3]*sinpi(p1[5]/π) ) + p2Prime[2]*cospi(p1[6]/π)
        p2Prime[3]cospi(p1[5] /π) - p2Prime[1]*sinpi(p1[5] /π)
    ]
    
    p2Mag = get_pMag(T, MASS)   
    θ2 = acos(p2[3]/p2Mag)
    ϕ2 = atan(p2[2],p2[1])  

        
    if (p2[1] == NaN)
        get_second_vector(T, MASS, p1, θdif)
    end
    
    return [ p2Mag, p2[1], p2[2], p2[3], θ2, ϕ2 ]
end