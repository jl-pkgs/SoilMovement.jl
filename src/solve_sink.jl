"""
Calculates hydraulic sink term for each soil layer. Currently, this is limited
to transpiration from each soil layer. In the future, this may possibly include
evaporation from the soil surface. The soil water soil water stress constraint
is estimated as described in:

Verhoef, A., & Egea, G. (2014). Modeling plant transpiration under limited soil
  water: Comparison of different plant and soil hydraulic parameterizations and
  preliminary implications for their use in land surface models. *Agricultural
  and Forest Meteorology*, 191, 22–32.

NOTE: This function does not distinguish between liquid θ and the overall (ice
and water) θ. This is for simplicity, and because it is assumed that potential
transpiration is already limited by temperature. Increasing transpiration (e.g.,
by setting q to a value less than 1) will increase the magnitude of dry-down but
not necessarily decrease peak soil moisture in near-surface layers during
infiltration events.

# Parameters
- `θ::Array{Float64, 2}`: Current soil moisture (θ) in each soil layer
- `transpiration::Float64`: Total potential (unconstrained) transpiration rate
  (kg m-2 s-1), a daily scalar value
- `q::Float64`: Curvature exponent for the soil water stress factor (Default: 1)
- `use_balland::Bool`: True to use the self-consistent formulae for field
  capacity and wilting point from Balland et al. (2008); False to define those
  based on soil matric potentials of -0.033 MPa and -1.5 MPa, respectively
  (Default: True)
- `clip_to_saturation::Bool`: True to force field capacity to be no larger than
  the saturation porosity (takes the minimum of field capacity and saturation
  porosity) (Default: True)

# Returns
- `Tuple{Array{Float64, 2}, Nothing, Nothing}`: A 3-tuple of (transpiration,
  soil_evaporation, canopy_evaporation) arrays; transpiration is a (Z x 1) array
  and the other two elements are currently `None`.
"""
function solve_sink(profile::SoilProfile, θ, transpiration, q=1, use_balland=true, clip_to_saturation=true)

  if !use_balland
    fc = profile.field_capacity
    wp = profile.wilting_point
  else
    fc = profile.field_capacity_balland
    wp = profile.wilting_point_balland
  end

  clip_to_saturation && (clamp!(fc, 0, profile.θ_sat))
  # fc = min.(fc, profile.θ_sat)

  stress = (θ .- wp) ./ (fc .- wp)
  stress = clamp.(stress, 0, 1)

  return (stress .^ q) .* profile._root_fraction .* transpiration # trans_i
end
