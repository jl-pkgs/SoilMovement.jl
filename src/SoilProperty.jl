const DENSITY_ICE = 917 # kg m-3 (Cutnell & Johnson. 1995. "Physics," 3rd ed.)
const DENSITY_WATER = 1000 # kg m-3
const DRAINAGE_DECAY_FACTOR = 2.5 # m-1 (CLM 4.5)
const MATRIC_POTENTIAL_ORGANIC = -10.3 # mm (CLM 4.0)
const PERCOLATION_THRESHOLD = 0.5 # From percolation theory and CLM 4.0
const SOCC_MAX = 130e3 # g m-3 (Lawrence and Slater 2008)
const SOIL_FREEZING = 273.15 # Temp. below which soil is frozen (deg K)


"""
The Clapp & Hornberger exponent, "B"
"""
function b(profile::SoilProfile)
  b_min = 2.91 .+ 0.159 .* (profile.frac_clay .* 100)
  return ((1 .- profile.frac_organic) .* b_min) .+ (profile.frac_organic .* 2.7)
end


"""
Fraction of soil layer allows percolation in organic material
"""
function frac_percolating(profile::SoilProfile)
  return ifelse.(profile.frac_organic .< PERCOLATION_THRESHOLD, 0, (1 - PERCOLATION_THRESHOLD) .^ -0.139 .* (profile.frac_organic .- PERCOLATION_THRESHOLD) .^ 0.139 .* profile.frac_organic)
end


"""
Bulk saturated hydraulic conductivity of the soil (mm s-1)
"""
function ksat(profile::SoilProfile)
  f_uncon = 1 .- frac_percolating(profile)
  return (f_uncon .* ksat_uncon(profile)) .+ ((1 .- f_uncon) .* profile.params[:ksat_om])
end

"""
Saturated hydraulic conductivity for mineral soil (mm s-1)
"""
function ksat_min(profile::SoilProfile)
  return 0.0070556 .* 10 .^ (-0.884 .+ 0.0153 .* profile.frac_sand .* 100)
end

"""
Hydraulic conductivity of the saturated, unconnected fraction (mm s-1)
"""
function ksat_uncon(profile::SoilProfile)
  f_uncon = 1 .- frac_percolating(profile)
  return f_uncon .* (1 ./ ((1 .- profile.frac_organic) ./ ksat_min(profile) .+ (profile.frac_organic .- frac_percolating(profile)) ./ profile.params[:ksat_om]))
end

"""
Saturated matric potential, in millimeters (mm)
"""
function psi_sat(profile::SoilProfile)
  return ((1 .- profile.frac_organic) .* psi_sat_min(profile)) .+ (profile.frac_organic .* MATRIC_POTENTIAL_ORGANIC)
end

"""
Saturated matric potential of mineral soil, in millimeters (mm)
"""
function psi_sat_min(profile::SoilProfile)
  return -10 .* 10 .^ (1.88 .- 0.0131 .* profile.frac_sand .* 100)
end

"""
Soil root fraction in each layer
"""
function root_fraction(profile::SoilProfile)
  beta = [NaN, 0.959, 0.962, 0.976, 0.966, 0.943, 0.964, 0.961, 0.961][profile.pft]
  depth_cm = -profile.depths_m .* 100
  return beta .^ vcat([0.0], depth_cm[1:end-1]) .- beta .^ depth_cm
end

"""
Saturation water content (or saturation porosity)
"""
function theta_sat(profile::SoilProfile)
  return fill(profile.porosity, size(profile.depths_m))
end

"""
Convert a matric potential to a corresponding VWC
"""
function potential_to_vwc(profile::SoilProfile, psi)
  theta_crit = (theta_sat(profile) .* (psi ./ psi_sat(profile))) .^ -(1 ./ b(profile))
  return clamp.(theta_crit, 0, 1)
end

"""
Critical point (VWC) or field capacity of the soil, conventionally
defined as -0.033 MPa. Returns the equivalent volumetric water content
(m3 m-3).
"""
function field_capacity(profile::SoilProfile)
  psi_crit = -0.033e6 / 9.8
  return potential_to_vwc(profile, psi_crit)
end

"""
Field capacity of the soil, from Balland et al. (2008).
"""
function field_capacity_balland(profile::SoilProfile)
  return theta_sat(profile) .* (0.565 .+ ((0.991 .- 0.565) .* sqrt.(profile.frac_clay))) .* exp.(-((0.103 .* profile.frac_sand) .- (0.785 .* profile.frac_organic)) ./ theta_sat(profile))
end

"""
Permanent wilting point of the soil, conventionally defined as
-1.5 MPa (Tolk et al. 2003). Returns the equivalent volumetric water
content (m3 m-3).
"""
function wilting_point(profile::SoilProfile)
  psi_wilt = -1.5e6 / 9.8
  return potential_to_vwc(profile, psi_wilt)
end

"""
Permanent wilting point of the soil, from Balland et al. (2008).
"""
function wilting_point_balland(profile::SoilProfile)
  fc = field_capacity_balland(profile)
  return fc .* (0.17 .+ ((0.832 .- 0.17) .* sqrt.(profile.frac_clay))) .* exp.(-1.4 .* profile.frac_organic ./ fc)
end

"""
The ice fraction of the combined liquid and ice water volumes, after
the empirical formulation by Decker and Zeng (2006, Geophys. Res.
Lett.).
"""
function f_ice(profile::SoilProfile, vwc, temp_k, alpha=2, beta=4)
  wetness = vwc ./ theta_sat(profile)
  f_ice = (1 .- exp.(alpha .* wetness .^ beta .* (temp_k .- SOIL_FREEZING))) ./ exp.(1 .- wetness)
  f_ice = clamp.(f_ice, 0, 1)
  return f_ice
end

"""
The ice fraction of the combined liquid and ice water volumes, after
the empirical formulation by from the European Centre for Medium Range
Weather Forecasting (ECMWF), as described by Decker and Zeng (2006,
Geophys. Res. Lett.).
"""
function f_ice2(profile::SoilProfile, vwc, temp_k, field_capacity=nothing)
  if field_capacity === nothing
    field_capacity = field_capacity_balland(profile)
  end
  vwc_ice = (field_capacity ./ 2) .* (1 .- abs.(sin.((π .* (temp_k .- SOIL_FREEZING .- 2)) ./ 4)))
  vwc_ice = clamp.(vwc_ice, 0, field_capacity)
  vwc_ice = ifelse.(temp_k .> SOIL_FREEZING + 1, 0, ifelse.(temp_k .< SOIL_FREEZING - 3, field_capacity, vwc_ice))
  return vwc_ice ./ vwc
end

"""
Ice impedance of the soil layers.
"""
function f_impedance(profile::SoilProfile, vwc, f_ice)
  return 10 .^ (-6 .* ((vwc .* f_ice) ./ theta_sat(profile)))
end

"""
Hydraulic conductivity (mm s-1) of each soil layer, as a function of
the soil water and ice volumes.
"""
function h2o_conductivity(profile::SoilProfile, vwc, f_ice)
  vwc_liq = vwc .* (1 .- f_ice)
  impedance_i = f_impedance(profile, vwc, f_ice)
  impedance_n = impedance_i[end]
  b_exp = 2 .* b(profile) .+ 3
  k = impedance_i[1:end-1] .* ksat(profile)[1:end-1] .* ((0.5 .* (vwc_liq[1:end-1] .+ vwc_liq[2:end])) ./ theta_sat(profile)[2:end]) .^ b_exp[1:end-1]
  kn = impedance_n .* ksat(profile)[end] .* (vwc_liq[end] ./ theta_sat(profile)[end]) .^ b_exp[end]
  return vcat(k, kn)
end

# end
# module SoilProfile
"""
The soil matric potential (mm), defined at the "node depth," or at the
midpoint of the soil layer.

```math
\\psi_i = \\psi_{sat,i}\\left(
  \\frac{\theta_i}{\theta_{sat,i}}
\\right)^{-B_i}
\\quad\\mbox{where}\\quad \\psi_i \\ge -1\\times 10^8;\\quad
0.01 \\le \\frac{\\theta_i}{\theta_{sat,i}} \\le 1
```

Where `psi_sat` is the saturated soil matric potential, `B` is the
Clapp & Hornberger exponent; see `SoilProfile.h2o_conductivity()`.

# Parameters
- `vwc::Array{Float64, 2}`: (Z x 1) array of soil volumetric water content (VWC)
- `f_ice::Array{Float64, 2}`: (Z x 1) array of the ice fraction

# Returns
- `Array{Float64, 2}`: (Z x 1) array of soil matric potential
"""
function matric_potential(profile::SoilProfile, vwc, f_ice)
  vwc_liq = vwc .* (1 .- f_ice)
  quo = vwc_liq ./ profile._theta_sat
  quo = clamp.(quo, 0.01, 1)
  psi0 = profile._psi_sat .* quo .^ -profile._b
  return clamp.(psi0, -1e8, Inf)
end


@doc raw"""
The maximum infiltration capacity of the (surface) soil layer.

```math
q_{max} = (1 - f_{sat}) \Theta_{ice} k_{sat}
```

Where `f_sat` is the fraction of the land surface that is saturated,
`Theta_ice` is the impedance due to ice, and `k_sat` is the saturated
hydraulic conductivity.

# Parameters
- `vwc::Array{Float64, 2}`: (Z x 1) array of soil volumetric water content (VWC)
- `temp_k::Array{Float64, 2}`: (Z x 1) array of soil temperatures in degrees K
- `f_saturated::Array{Float64, 2}`: (Z x 1) array of the fraction of the land surface that is saturated
- `f_ice::Union{Array{Float64, 2}, Nothing}`: (Optional) (Z x 1) array of the ice fraction; will be calculated based on VWC and temperature if None

# Returns
- `Array{Float64, 1}`: The maximum infiltration capacity (kg m-2 s-1); array of shape (N,)
"""
function max_infiltration(profile::SoilProfile, vwc, temp_k, f_saturated, f_ice=nothing)
  if f_ice === nothing
    f_ice = profile.f_ice(vwc, temp_k)
  end
  impedance_i = profile.f_impedance(vwc, f_ice)
  return (1 .- f_saturated) .* impedance_i[1] .* profile._ksat[1]
end

