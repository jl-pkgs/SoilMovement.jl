"""
Calculates hydraulic sink term for each soil layer. Currently, this is limited
to transpiration from each soil layer. In the future, this may possibly include
evaporation from the soil surface. The soil water soil water stress constraint
is estimated as described in:

Verhoef, A., & Egea, G. (2014). Modeling plant transpiration under limited soil
  water: Comparison of different plant and soil hydraulic parameterizations and
  preliminary implications for their use in land surface models. *Agricultural
  and Forest Meteorology*, 191, 22–32.

NOTE: This function does not distinguish between liquid VWC and the overall (ice
and water) VWC. This is for simplicity, and because it is assumed that potential
transpiration is already limited by temperature. Increasing transpiration (e.g.,
by setting q to a value less than 1) will increase the magnitude of dry-down but
not necessarily decrease peak soil moisture in near-surface layers during
infiltration events.

# Parameters
- `vwc::Array{Float64, 2}`: Current soil moisture (VWC) in each soil layer
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
function solve_sink(profile::SoilProfile, vwc, transpiration, q=1, use_balland=true, clip_to_saturation=true)
  fc = profile.field_capacity_balland
  wp = profile.wilting_point_balland

  if !use_balland
    fc = profile.field_capacity
    wp = profile.wilting_point
  end

  if clip_to_saturation
    fc = min.(fc, profile._theta_sat)
  end
  stress = (vwc .- wp) ./ (fc .- wp)
  stress = clamp.(stress, 0, 1)
  trans_i = (stress .^ q) .* profile._root_fraction .* transpiration
  return (trans_i, nothing, nothing)
end


@doc raw"""
Solves for volumetric water content (VWC) at a single time step for each depth
using a tridiagonal system of equations for the water balance. The free drainage
("flux") bottom boundary condition is always enforced because, otherwise, soil
layers will saturate quickly; the aquifer is uncoupled and assumed to lie below
the soil layer. Other considerations:

1. Below the soil profile there is no ice (ice-filled fraction is zero); this is
    out of necessity in calculating derivatives but is also consistent with a
    geothermal heat flux that maintains above-freezing conditions below the
    profile.
2. Sub-surface runoff is computed separately; it is one of the values returned
    by this function and should be subtracted from the soil VWC profile after
    updating with the change in VWC estimated by this function.
3. Similarly, if there is a perched, saturated layer above a frozen layer, the
    returned value of "runoff" includes lateral drainage from the perched
    layer(s), which should also be subtracted from the soil VWC profile.

At a minimum, runoff includes drainage according to the free-drainage or "flux"
bottom boundary condition of CLM 5.0:

```math
q_{drain} = k_i + \left[
\frac{\partial\, k}{\partial\, \theta_{liq}} \times \Delta \theta_{liq}
\right]_i
```

Where `k` is the hydraulic conductivity. If saturated conditions exist
within the soil column (saturated from the bottom-up), then additional,
lateral sub-surface runoff is calculated as described in CLM 4.0
Technical Note, Section 7.5:

```math
q_{drain} = \Theta_{ice}\, 10\,\mathrm{sin}(\beta )\,
    \mathrm{exp}(-f_{drain} z_{\nabla})
    \quad\mbox{where}\quad f_{drain} = 2.5\,\mathrm{m}^{-1}
```

Where `beta` is the slope, `z_nabla` is the depth to the top of the
saturated zone, and `Theta_ice` is the impedance due to ice. The
specific yield is calculated:

```math
S_y = \theta_{sat}\left(1 - \left(
1 + \frac{z_{\nabla}}{\Psi_{sat}}
\right)^{-1/B}\right)
```

# Parameters

- `influx::Union{Float64, Array{Float64, 1}}`: Scalar or N-dimensional array of
  water influx at the surface layer, in units of (kg m-2 s-1) or (mm s-1), as 1
  mm of water over an area of 1 m-2 weighs 1 kg.
- `vwc::Array{Float64, 2}`: (Z x 1) array of total soil volumetric water content
  (VWC) for the current time step, including both liquid and ice water content
- `temp_k::Array{Float64, 2}`: (Z x 1) array of soil temperatures in degrees K
  for the current time step
- `dt::Int`: Size of time step (secs)
- `transpiration::Union{Array{Float64, 2}, Nothing}`: (Optional) Transpiration
  in each soil layer (kg m-2 s-1), a (Z x 1) array
- `saturated_below::Bool`: True to invoke a virtual soil layer below the soil
  profile that is fully saturated (Default: False)

# Returns
- `Tuple{Array{Float64, 2}, Tuple{Array{Float64, 2}, Array{Float64, 2}},
  Array{Float64, 2}}`: Returns a 3-tuple of (`solution`, `flows`, `runoff`)

  + `solution`: the change in VWC in each layer; 
  + `flows`: a tuple of (`q_in`, `q_out`) where `q_in` is the flow into each layer from the layer
  above, and `q_out` is the flow out of each layer (all flows in units of mm
  s-1); 
  + `runoff`: the change in VWC due to lateral sub-surface runoff.
"""
function solve_vwc(profile::SoilProfile, influx, vwc, temp_k, dt, transpiration=nothing, saturated_below=false)
  function _dk_dliq(vwc_liq, temp_k, mean_impedance, bottom_vwc_liq)
    mean_vwc_liq = vcat(
      0.5 * (vwc_liq[1:end-1] .+ vwc_liq[2:end]),
      0.5 * (vwc_liq[end] .+ bottom_vwc_liq)
    )
    result = (2 * self._b .+ 3) .* mean_impedance .* self._ksat .* (mean_vwc_liq ./ self._theta_sat) .^ (2 * self._b .+ 2) .* (0.5 ./ self._theta_sat)
    return result
  end

  function _drain_runoff(vwc, solution, mean_impedance, dk_dliq)
    drain_runoff = copy(self._zeros)
    sp_yield = Inf
    if self._z_bedrock < self._z[end]
      drain_runoff[end] = k[end] + (dk_dliq[end] * solution[end])
    end
    sat_mask = falses(size(vwc))
    for i in 1:length(vwc)
      j = length(vwc) - i + 1
      if vwc[j] >= 0.9 * self._theta_sat[j]
        sat_mask[j] = true
      else
        break
      end
    end
    if any(sat_mask)
      table_depth = self._depths_m[1:end-1][sat_mask[2:end]]
      if length(table_depth) > 0
        drain_runoff[end] = mean(mean_impedance[sat_mask]) * 10 * sin(self._slope) * exp(-self.DRAINAGE_DECAY_FACTOR * maximum(table_depth))
        sp_yield = mean(self._theta_sat[sat_mask]) * (1 - (1 + maximum(table_depth) / mean(self._psi_sat[sat_mask]))^(-1 / mean(self._b[sat_mask])))
      end
    end
    return (drain_runoff, sp_yield)
  end

  function _drain_perched(vwc, temp_k, f_ice)
    drain_perched = copy(self._zeros)
    sp_yield = Inf
    perched = (vwc[1:end-1] .>= 0.9 * self._theta_sat[1:end-1]) .& (f_ice[2:end] .> 0)
    if any(perched)
      i_perch = findfirst(perched)
      i_frost = findlast(perched) + 1
      impedance = self.f_impedance(vwc, f_ice)
      ksat_perch = 10e-5 * sin(self._slope) * sum((impedance.*self._ksat.*self._thickness_mm)[i_perch:i_frost]) / sum(self._thickness_mm[i_perch:i_frost])
      drain_perched[i_frost-1] = ksat_perch * (self._depths_m[i_frost-1] - self._depths_m[i_perch-1]) * 1e3
      perch_mask = vcat(perched, false)
      sp_yield = mean(self._theta_sat[perch_mask]) * (1 - (1 + (self._depths_m[i_frost-1] * 1e3) / mean(self._psi_sat[perch_mask]))^(-1 / mean(self._b[perch_mask])))
    end
    return (drain_perched, sp_yield)
  end

  if all(temp_k .> 276)
    f_ice = zeros(size(vwc))
  else
    f_ice = profile.f_ice(vwc, temp_k)
  end
  vwc_liq = vwc .* (1 .- f_ice)
  vwc_ice = vwc .* f_ice
  mean_f_ice = 0.5 .* vcat(vwc_ice[2:end], zeros(size(vwc_ice)[1:1]))
  mean_impedance = profile.f_impedance(vwc, mean_f_ice)
  psi = profile.matric_potential(vwc, f_ice)
  k = profile.h2o_conductivity(vwc, f_ice)
  nans = fill(NaN, size(k)[1:1])

  if saturated_below
    dk_dliq = _dk_dliq(vwc_liq, temp_k, mean_impedance, mean([vwc_liq[end], profile._theta_sat[end]]))
  else
    dk_dliq = _dk_dliq(vwc_liq, temp_k, mean_impedance, vwc_liq[end])
  end
  dk0_dliq = vcat(NaN, dk_dliq[1:end-1])

  dpsi_dliq = clamp.(-profile._b .* (psi ./ vwc_liq), 0.01 * profile._theta_sat, profile._theta_sat)

  n_diff = profile._z_node .- vcat(NaN, profile._z_node[1:end-1])
  dqout0_dliq0 = -((vcat(nans, k[1:end-1]) ./ n_diff) .* vcat(nans, dpsi_dliq[1:end-1])) .- dk0_dliq .* ((vcat(nans, psi[1:end-1]) .- psi) .+ n_diff) ./ n_diff

  dqout0_dliq = ((vcat(nans, k[1:end-1]) ./ n_diff) .* dpsi_dliq) .- dk0_dliq .* ((vcat(nans, psi[1:end-1]) .- psi) .+ n_diff) ./ n_diff

  dqout_dliq = vcat(dqout0_dliq0[2:end], NaN)
  dqout_dliq[end] = dk_dliq[end] ./ profile._theta_sat[end]

  dqout_dliq1 = vcat(dqout0_dliq[2:end], NaN)

  contrast = (psi[1:end-1] .- psi[2:end]) .+ (profile._z_node[2:end] .- profile._z_node[1:end-1]) ./ (profile._z_node[2:end] .- profile._z_node[1:end-1])
  q_out = vcat(-k[1:end-1] .* contrast, -k[end])
  q_in = vcat(influx, -k[1:end-1] .* contrast)

  nn = length(self._depths_m)
  lhs = zeros(nn, nn)
  rhs = fill(NaN, nn)
  for z in 1:nn
    if z == 1
      b = dqout_dliq[z] - (self._thickness_mm[z] / dt)
      c = dqout_dliq1[z]
      lhs[z, 1:2] = [b, c]
    elseif z < self.nlayer
      a = -dqout0_dliq0[z]
      b = dqout_dliq[z] - dqout0_dliq[z] - (self._thickness_mm[z] / dt)
      c = dqout_dliq1[z]
      lhs[z, (z-1):(z+1)] = [a, b, c]
    else
      a = -dqout0_dliq0[z]
      b = dqout_dliq[z] - dqout0_dliq[z] - (self._thickness_mm[z] / dt)
      lhs[z, end-1:end] = [a, b]
    end
    if isnothing(transpiration)
      r = q_in[z] - q_out[z]
    else
      r = q_in[z] - q_out[z] - transpiration[z]
    end
    rhs[z] = r
  end

  banded = vcat(
    hcat(diagm(1 => lhs[2:end, :]), zeros(nn)),
    diagm(0 => lhs),
    hcat(zeros(nn), diagm(-1 => lhs[:, 2:end]))
  )
  solution = tridiag_solver(lhs, rhs, banded=banded)[:, 1]

  q_runoff, sp_yield = _drain_runoff(vwc, solution, mean_impedance, dk_dliq)
  runoff = abs.(q_runoff) .* (dt ./ thickness_mm)
  if isfinite(sp_yield)
    runoff = ifelse.(abs.(runoff) .> sp_yield, sp_yield, runoff)
  end

  # Compute lateral sub-surface drainage from perched zone; convert
  #   to (change in) VWC and compare to specific yield
  q_perched, sp_yield = _drain_perched(vwc, temp_k, f_ice)
  drainage = abs.(q_perched) .* (dt ./ thickness_mm)
  if isfinite(sp_yield)
    drainage = ifelse.(abs.(drainage) .> sp_yield, sp_yield, drainage)
  end

  solution, (q_runoff, q_perched), runoff .+ drainage
end


export solve_vwc
