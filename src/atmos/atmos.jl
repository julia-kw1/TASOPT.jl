module atmosphere
using Roots
export atmos, find_altitude_from_density

"""
    atmos(h, ΔT)
    
Atmospheric functions ` T(h)`, `ρ(h)` etc
valid to `h`=20km, `p(h)` valid to `h`=70km.

Also calculates viscosity using Sutherland's law. Non-standard sea-level temperatures are allowed
with an ISA + ΔT like model.

Units:
- [h]   = km ASL
- [T]   = Kelvin
- [p]   = Pa
- [ρ]   = kg/m^3
- [a]   = m/s
- [μ]   = kg/m-s 
"""
function atmos(h::Float64, ΔT::Float64 = 0.0)
    # convert h to m
    h_m = h * 1.0e3
    
    # constants
    cp  = 1004.0   # J/kg-K
    ɣ = 1.4 
    Tsuth = 110.4 # K, Sutherland's temperature

    # sea-level values
    pSL = 1.01325e5 # Pa
    μSL  = 1.7894e-5  # kg/m-s
    g0 = 9.80665 # m/s^2
    R = 287.05287 # J/kg-K

    # layer values for ISA
    h1 = 11.0e3 # m
    h2 = 20.0e3 
    h3 = 32.0e3
    h4 = 47.0e3

    # standard temperature at layer
    T0 = 288.15 # K
    T1 = 216.65
    T2 = 216.65
    T3 = 228.65

    # lapse rates
    L0 = -6.5e-3 # K/m
    L1 = 0.0
    L2 = 1.0e-3
    L3 = 2.8e-3

    # base pressures
    p1 = pSL * (T1/T0)^(-g0/(L0*R))
    p2 = p1 * exp(-g0*(h2-h1)/(R*T1))
    p3 = p2 * (T3/T2)^(-g0/(L2*R))
    T4 = T3 + L3*(h4-h3)

    Tstd = 0.0
    p = 0.0

    if 0.0 <= h_m < h1
        Tstd = T0 + L0 * h_m
        p = pSL * (Tstd/T0)^(-g0/(L0*R))
    elseif h1 <= h_m < h2
        Tstd = T1
        p = p1 * exp(-g0*(h_m-h1)/(R*T1))
    elseif h2 <= h_m < h3
        Tstd = T2 + L2 * (h_m - h2)
        p = p2 * (Tstd/T2)^(-g0/(L2*R))
    elseif h3 <= h_m <= h4
        Tstd = T3 + L3 * (h_m - h3)
        p = p3 * (Tstd/T3)^(-g0/(L3*R))
    else
        error("Altitude out of range for pressure calculation (h <= 47km)")
    end

    # apply temperature offset
    T = Tstd + ΔT

    # calculate density, speed of sound, and viscosity
    ρ = p / (R * T)
    a = sqrt(ɣ * R * T)
    μ = μSL * (T/T0)^(3/2) * (T0 + Tsuth) / (T + Tsuth)

 return T,p,ρ,a,μ

end # atmos

"""
    find_altitude_from_density(ρ::Float64, ΔT::Float64 = 0.0) 
    
Uses a non-linear solver to find the altitude corresponding to a given air density.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ρ::Float64`: air density (kg/m^3)
    - `ΔT::Float64`: temperature difference from standard atmosphere (K)
    
    **Outputs:**
    - `h::Float64`: altitude (km)
"""
function find_altitude_from_density(ρ::Float64, ΔT::Float64 = 0.0) 
    res(x) = atmos(x, ΔT)[3] - ρ #Residual for density
    h = find_zero(res, 0.0)
    return h
end

end