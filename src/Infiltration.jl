# module Infiltration
# end

"""
A soil water infiltration model, based on the SoilProfile class and
facilitating a maximum infiltration rate, transpiration loss, sub-surface
drainage, and with adaptive time stepping. Outstanding issues:

1. Frozen layers may exceed the saturation porosity because liqud water
    content can't be moved to or from those layers during rebalancing.

# Parameters
- `soil_model::SoilProfile`: Soil profile model
- `dt_min::Int`: Minimum number of seconds a sub-daily time step can take
- `f_ice_cutoff::Float64`: Ice fraction cutoff on the interval [0, 1], but the value should
    be >= 0.95. If the ice fraction exceeds this value, a rebalancing
    of soil moisture will not be performed. This can be necessary to
    avoid running into impossible balancing scenarios.
- `debug::Bool`: True to perform some (potentially expensive) validation checks at
    runtime (Default: False)
"""
@with_kw mutable struct InfiltrationModel
  soil_model::Any
  soil::Any
  dt_min::Int = 10
  f_ice_cutoff::Float64 = 0.96
  debug::Bool = false
end

"""
Runs the soil water infiltration model forward in time for a certain
number of days.

# Parameters
- `vwc::Array{Float64, 2}`: (Z x 1) array of the initial soil volumetric water content (VWC) profile
- `temp_profile::Array{Float64, 2}`: (Z x 1) array of soil temperatures in degrees K for the current time step
- `transpiration::Union{Array{Float64, 1}, Nothing}`: Sequence of the total daily potential (unconstrained) transpiration rate (kg m-2 s-1)
- `influx::Array{Float64, 1}`: Sequence of the daily water infiltration rate at the surface layer, in units of (kg m-2 s-1) or (mm s-1)
- `f_saturated::Array{Float64, 1}`: Sequence of the daily fraction of the land surface that is saturated
- `dt::Int`: Size of time step (secs)
- `n_days::Union{Int, Nothing}`: Number of days to run; defaults to the size of `influx`
- `ltol::Float64`: Lower bound for error tolerance in adaptive time stepping
- `utol::Float64`: Upper bound for error tolerance in adaptive time stepping
- `climatology::Bool`: True to run in climatology mode; i.e., if the input driver data are a 365-day climatology, the day index should be recycled
- `adaptive::Bool`: True to use adaptive time stepping: dynamic adjustment of the sub-daily time step based on the error in water balance (Default: True)

# Returns
- `Tuple{Array{Float64, 2}, Array{Float64, 2}, Union{Array{Float64, 2}, Nothing}}`: 3-tuple of `(vwc, err, psi)` where `vwc` is the soil moisture time series, a (Z x T) array; `err` is the estimated truncation error, a (Z x T) array; `psi` is the estimated soil matric potential, a (Z x T) array, where T is time and Z is the number of layers.
"""
function run(model::InfiltrationModel, vwc, temp_profile, transpiration, influx, f_saturated, dt, n_days::Union{Int,Nothing}=nothing, ltol::Float64=1e-2, utol::Float64=1e-1, climatology::Bool=false, adaptive::Bool=true)
  if n_days === nothing
    n_days = size(influx, 1)
  end

  est_vwc = fill(NaN, model.soil.nlayer, n_days)
  est_error = fill(NaN, model.soil.nlayer, n_days)
  est_psi = nothing

  if model.debug
    est_psi = fill(NaN, model.soil.nlayer, n_days)
    if !climatology
      @assert length(transpiration) == n_days
    end
  end

  iterations = 1:n_days
  if climatology
    iterations = iterations .% 365
  end

  for (i, d) in enumerate(iterations)
    successful = false
    while !successful && (model.dt_min <= dt < SECS_PER_DAY ÷ 2)
      _trans = transpiration === nothing ? nothing : transpiration[d]
      args = (vwc, temp_profile[:, d], _trans, influx[d], f_saturated[d], dt)
      if model.debug
        vwc, de = step_daily(model, args...)
      else
        try
          vwc, de = step_daily(model, args...)
        catch
          println("ERROR: Ending prematurely due to error in InfiltrationModel.step_daily()")
          return (est_vwc, est_error, est_psi)
        end
      end
      
      err = maximum(mean(abs.(hcat(de...)), dims=1))
      if !adaptive || (ltol < err <= utol)
        successful = true
        continue
      end
      if err <= ltol
        successful = true
      end
      d_dt = err <= ltol ? 2 : 0.5
      if model.dt_min <= (dt * d_dt) < SECS_PER_DAY ÷ 2
        dt = Int(dt * d_dt)
      end
    end
    est_vwc[:, i] = vwc[:]
    est_error[:, i] = mean(hcat(de...), dims=1)[:]
    if model.debug
      f_ice = model.soil.f_ice(vwc, temp_profile[:, d])
      psi = model.soil.matric_potential(vwc, f_ice)
      est_psi[:, i] = psi[:]
    end
  end
  return (est_vwc, est_error, est_psi)
end

"""
Executes a single daily time step of the soil water infiltration model.

# Parameters
- `vwc::Array{Float64, 2}`: (Z x 1) array of the initial soil volumetric water content (VWC) profile
- `temp_profile::Array{Float64, 2}`: (Z x 1) array of soil temperatures in degrees K for the current time step
- `transpiration::Float64`: Total potential (unconstrained) transpiration rate (kg m-2 s-1), a daily scalar value
- `influx::Float64`: Scalar or N-dimensional array of water infiltration at the surface layer, in units of (kg m-2 s-1) or (mm s-1)
- `f_saturated::Float64`: The fraction of the land surface that is saturated
- `dt::Int`: Size of time step (secs)

# Returns
- `Tuple{Array{Float64, 2}, Array{Float64, 1}}`: 2-tuple of `(vwc, error)` where `vwc` is the updated soil moisture profile and `error` is the estimated truncation error.
"""
function step_daily(model::InfiltrationModel, vwc, temp_profile, transpiration, influx, f_saturated, dt)

  function rebalance(vwc, temp_k, thickness_mm)
    if all(temp_k .> 276)
      f_ice = zeros(size(vwc))
    else
      f_ice = model.soil.f_ice(vwc, temp_k)
    end
    wliq = (vwc .- (vwc .* f_ice)) .* -thickness_mm
    wliq_max = (model.soil._theta_sat .- (vwc .* f_ice)) .* -thickness_mm
    i = 0
    while !all((0.01 .<= wliq .<= wliq_max) .| (wliq_max .< 0.01))
      @assert i < 1000
      excess = ifelse.(wliq .> wliq_max, wliq .- wliq_max, 0)
      deficit = ifelse.(wliq .< 0.01, 0.01 .- wliq, 0)
      excess = ifelse.(wliq_max .< 0.01, 0, excess)
      deficit = ifelse.(wliq_max .< 0.01, 0, deficit)
      wliq .-= excess .+ vcat(excess[2:end], 0)
      wliq .+= deficit .- vcat(0, deficit[1:end-1])
      i += 1
    end
    vwc = (wliq ./ -thickness_mm) .+ (vwc .* f_ice)
    return vwc
  end

  de = []
  iterations = 1:(SECS_PER_DAY÷dt)
  thickness_mm = model.soil._thickness_mm
  if model.debug
    @assert !haskey(transpiration, :length)
  end
  for t in iterations
    @assert all(0 .<= vwc .<= 1)
    actual_trans = zeros(size(vwc))
    if transpiration !== nothing
      actual_trans, _, _ = model.soil.solve_sink(vwc, transpiration)
    end
    max_influx = model.soil.max_infiltration(vwc, temp_profile, f_saturated)
    x, flows, runoff = model.soil.solve_vwc(min(influx, max_influx[1]), vwc, temp_profile, dt, actual_trans)
    q_in, q_out = flows
    vwc .+= x
    vwc .+= runoff
    vwc = rebalance(vwc, temp_profile, thickness_mm)
    if t > 1
      err = (dt / 2) * (((x .* thickness_mm) ./ dt) .- (q_in0 .- q_out0 .- actual_trans))
      push!(de, err)
    end
    q_in0, q_out0 = q_in, q_out
  end
  return (vwc, de)
end
