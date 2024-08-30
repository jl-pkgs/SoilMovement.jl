function _drain_runoff(soil, θ, dθ, mean_impedance, dk_dliq)
  (; θ_sat, ψ_sat, b) = soil

  drain_runoff = copy(soil.zeros)
  sp_yield = Inf

  if soil.z_bedrock < soil.z[end]
    drain_runoff[end] = k[end] + (dk_dliq[end] * dθ[end])
  end

  sat_mask = falses(size(θ))
  n = length(θ)
  for i in 1:n
    j = n - i + 1
    if θ[j] >= 0.9 * θ_sat[j]
      sat_mask[j] = true
    else
      break
    end
  end

  if any(sat_mask)
    table_depth = θ_sat[1:end-1][sat_mask[2:end]]

    if length(table_depth) > 0
      drain_runoff[end] = mean(mean_impedance[sat_mask]) * 10 * sin(soil.slope) * exp(-DRAINAGE_DECAY_FACTOR * maximum(table_depth))
      sp_yield = mean(θ_sat[sat_mask]) * (1 - (1 + maximum(table_depth) / mean(ψ_sat[sat_mask]))^(-1 / mean(b[sat_mask])))
    end
  end
  return (drain_runoff, sp_yield)
end


function _drain_perched(soil, θ, temp_k, f_ice)
  (; Ksat, θ_sat, b) = soil
  drain_perched = copy(soil.zeros)
  sp_yield = Inf
  Δz = soil.thickness_mm

  perched = (θ[1:end-1] .>= 0.9 * θ_sat[1:end-1]) .& (f_ice[2:end] .> 0) # 饱和，且含冰

  if any(perched)
    i_perch = findfirst(perched)
    i_frost = findlast(perched) + 1
    inds = i_perch:i_frost

    Θ_ice = f_impedance(soil, θ, f_ice)
    ksat_perch = 10e-5 * sin(soil.slope) * sum((Θ_ice.*Ksat.*Δz)[inds]) / sum(Δz[inds]) # 2.7.107

    drain_perched[i_frost-1] = ksat_perch * (θ_sat[i_frost-1] - θ_sat[i_perch-1]) * 1e3 # error here

    perch_mask = vcat(perched, false)
    _θ_sat = mean(θ_sat[perch_mask])
    _ψ_sat = mean(ψ_sat[perch_mask])
    _b = mean(soil.b[perch_mask])
    sp_yield = _θ_sat * (1 - (1 + (θ_sat[i_frost-1] * 1e3) / _ψ_sat)^(-1 / _b))
  end

  return drain_perched, sp_yield
end


@doc raw"""
Solves for volumetric water content (θ) at a single time step for each depth
using a tridiagonal system of equations for the water balance. The free drainage
("flux") bottom boundary condition is always enforced because, otherwise, soil
layers will saturate quickly; the aquifer is uncoupled and assumed to lie below
the soil layer. Other considerations:

1. Below the soil profile there is no ice (ice-filled fraction is zero); this is
    out of necessity in calculating derivatives but is also consistent with a
    geothermal heat flux that maintains above-freezing conditions below the
    profile.
2. Sub-surface runoff is computed separately; it is one of the values returned
    by this function and should be subtracted from the soil θ profile after
    updating with the change in θ estimated by this function.
3. Similarly, if there is a perched, saturated layer above a frozen layer, the
    returned value of "runoff" includes lateral drainage from the perched
    layer(s), which should also be subtracted from the soil θ profile.

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
1 + \frac{z_{\nabla}}{\ψ_{sat}}
\right)^{-1/B}\right)
```

# Parameters

- `influx::Union{Float64, Array{Float64, 1}}`: Scalar or N-dimensional array of
  water influx at the surface layer, in units of (kg m-2 s-1) or (mm s-1), as 1
  mm of water over an area of 1 m-2 weighs 1 kg.
- `θ::Array{Float64, 2}`: (Z x 1) array of total soil volumetric water content
  (θ) for the current time step, including both liquid and ice water content
- `temp_k::Array{Float64, 2}`: (Z x 1) array of soil temperatures in degrees K
  for the current time step
- `dt::Int`: Size of time step (secs)
- `transpiration::Union{Array{Float64, 2}, Nothing}`: (Optional) Transpiration
  in each soil layer (kg m-2 s-1), a (Z x 1) array
- `saturated_below::Bool`: True to invoke a virtual soil layer below the soil
  profile that is fully saturated (Default: False)

# Returns
- `Tuple{Array{Float64, 2}, Tuple{Array{Float64, 2}, Array{Float64, 2}},
  Array{Float64, 2}}`: Returns a 3-tuple of (`dθ`, `flows`, `runoff`)

  + `dθ`: the change in θ in each layer; 
  + `flows`: a tuple of (`q_in`, `q_out`) where `q_in` is the flow into each layer from the layer
  above, and `q_out` is the flow out of each layer (all flows in units of mm
  s-1); 
  + `runoff`: the change in θ due to lateral sub-surface runoff.
"""
function solve_θ(soil::SoilProfile, influx, θ, temp_k, dt, transpiration=nothing, saturated_below=false)

  function _dk_dliq(soil, θ_liq, temp_k, mean_impedance, θ_liq_bot)
    (; Ksat, θ_sat, b) = soil
    mean_θ_liq = vcat(
      0.5 * (θ_liq[1:end-1] .+ θ_liq[2:end]),
      0.5 * (θ_liq[end] .+ θ_liq_bot)
    )
    (2 * soil.b .+ 3) .* mean_impedance .* Ksat .* (mean_θ_liq ./ θ_sat) .^ (2 * soil.b .+ 2) .* (0.5 ./ θ_sat)
  end

  if all(temp_k .> 276)
    f_ice = zeros(size(θ))
  else
    f_ice = soil.f_ice(θ, temp_k)
  end

  θ_liq = θ .* (1 .- f_ice)
  θ_ice = θ .* f_ice
  mean_f_ice = 0.5 .* vcat(θ_ice[2:end], zeros(size(θ_ice)[1:1]))
  mean_impedance = soil.f_impedance(θ, mean_f_ice)
  ψ = soil.matric_potential(θ, f_ice)
  k = soil.h2o_conductivity(θ, f_ice)
  nans = fill(NaN, size(k)[1:1])

  if saturated_below
    dk_dliq = _dk_dliq(soil, θ_liq, temp_k, mean_impedance, mean([θ_liq[end], soil.θ_sat[end]]))
  else
    dk_dliq = _dk_dliq(soil, θ_liq, temp_k, mean_impedance, θ_liq[end])
  end
  dk0_dliq = vcat(NaN, dk_dliq[1:end-1])

  dψ_dliq = clamp.(-soil._b .* (ψ ./ θ_liq), 0.01 * soil.θ_sat, soil.θ_sat)

  n_diff = soil._z_node .- vcat(NaN, soil._z_node[1:end-1])
  Δqout0_dliq0 = -((vcat(nans, k[1:end-1]) ./ n_diff) .* vcat(nans, dψ_dliq[1:end-1])) .- dk0_dliq .* ((vcat(nans, ψ[1:end-1]) .- ψ) .+ n_diff) ./ n_diff

  Δqout0_dliq = ((vcat(nans, k[1:end-1]) ./ n_diff) .* dψ_dliq) .- dk0_dliq .* ((vcat(nans, ψ[1:end-1]) .- ψ) .+ n_diff) ./ n_diff

  Δqout_dliq = vcat(Δqout0_dliq0[2:end], NaN)
  Δqout_dliq[end] = dk_dliq[end] ./ soil.θ_sat[end]

  Δqout_dliq1 = vcat(Δqout0_dliq[2:end], NaN)

  contrast = (ψ[1:end-1] .- ψ[2:end]) .+ (soil._z_node[2:end] .- soil._z_node[1:end-1]) ./ (soil._z_node[2:end] .- soil._z_node[1:end-1])
  q_out = vcat(-k[1:end-1] .* contrast, -k[end])
  q_in = vcat(influx, -k[1:end-1] .* contrast)

  nn = length(self.θ_sat)
  lhs = zeros(nn, nn)
  rhs = fill(NaN, nn)

  for z in 1:nn
    if z == 1
      b = Δqout_dliq[z] - (soil.thickness_mm[z] / dt)
      c = Δqout_dliq1[z]
      lhs[z, 1:2] = [b, c]
    elseif z < self.nlayer
      a = -Δqout0_dliq0[z]
      b = Δqout_dliq[z] - Δqout0_dliq[z] - (soil.thickness_mm[z] / dt)
      c = Δqout_dliq1[z]
      lhs[z, (z-1):(z+1)] = [a, b, c]
    else
      a = -Δqout0_dliq0[z]
      b = Δqout_dliq[z] - Δqout0_dliq[z] - (soil.thickness_mm[z] / dt)
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
  dθ = tridiag_solver(lhs, rhs, banded=banded)[:, 1]

  q_runoff, sp_yield = _drain_runoff(soil, θ, dθ, mean_impedance, dk_dliq)
  runoff = abs.(q_runoff) .* (dt ./ thickness_mm)
  if isfinite(sp_yield)
    runoff = ifelse.(abs.(runoff) .> sp_yield, sp_yield, runoff)
  end

  # Compute lateral sub-surface drainage from perched zone; convert to (change
  #   in) θ and compare to specific yield
  q_perched, sp_yield = _drain_perched(soil, θ, temp_k, f_ice)
  drainage = abs.(q_perched) .* (dt ./ thickness_mm)

  if isfinite(sp_yield)
    drainage = ifelse.(abs.(drainage) .> sp_yield, sp_yield, drainage)
  end

  dθ, (q_runoff, q_perched), runoff .+ drainage
end
