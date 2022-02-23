"""
    getpMag(T::Real, MASS::Real)

Description of getpMag
--------------------------
Computes magnitude of momentum as:
``` |p| = √((MASS + T)² - MASS²) ```.

Input arguments are:
+ T    : kinetic energy in [MeV]
+ MASS : particle mass in [MeV/c²]
"""
function getpMag(T::Real, MASS::Real)
	return sqrt( (MASS+T)^2 - MASS^2)
end

"""
    getpx(pMag, ϕ, θ)

Description of getpMag
--------------------------
Computes x-component of the momentum as:
``` pₓ = |p|*cos(ϕ)*sin(θ) ```

Input arguments are:
+ pMag : momentum magnitude in [MeV]
+ ϕ    : azimuth angle in [rad]
+ θ    : zenith angle in [rad]
"""
function getpx(pMag, ϕ, θ)
    return pMag*cos(ϕ)*sin(θ)
end

"""
    getpy(pMag, ϕ, θ)

Description of getpMag
--------------------------
Computes y-component of the momentum as:
``` p𝔂 = |p|*sin(ϕ)*sin(θ) ```

Input arguments are:
+ pMag : momentum magnitude in [MeV]
+ ϕ    : azimuth angle in [rad]
+ θ    : zenith angle in [rad]
"""
function getpy(pMag, ϕ, θ)
    return pMag*sin(ϕ)*sin(θ)
end

"""
    getpz(pMag, ϕ, θ)

Description of getpMag
--------------------------
Computes z-component of the momentum as:
``` p𝐳 = |p|*cos(θ) ```

Input arguments are:
+ pMag : momentum magnitude in [MeV]
+ θ    : zenith angle in [rad]
"""
function getpz(pMag, θ)
    return pMag*cos(θ)
end
