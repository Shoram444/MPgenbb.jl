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
+ θ    : zenith angle in [rad]
"""
function get_px(pMag, ϕ, θ)
    return pMag*cos(θ)*sin(ϕ)
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
+ θ    : zenith angle in [rad]
"""
function get_py(pMag, ϕ, θ)
    return pMag*sin(ϕ)*sin(θ)
end

"""
    get_pz(pMag, ϕ, θ)

Description of get_pz
--------------------------
Computes z-component of the momentum as:
``` p𝐳 = |p|*cos(θ) ```

Input arguments are:
+ pMag : momentum magnitude in [MeV]
+ θ    : zenith angle in [rad]
"""
function get_pz(pMag, θ)
    return pMag*cos(ϕ)
end
