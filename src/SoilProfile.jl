# module SoilProfile
"""
Represents a soil profile. In this first version, the total porosity and the
sand and clay fractions are scalar fields that represent the entire soil column.
Total porosity stands in for the saturation porosity, which is a function of
organic matter and sand content in CLM 4.0. These scalar fields are also fixed
throughout the simulation; i.e., changes in soil organic carbon as part of some
coupled soil decomposition model do not propagate to changes in the organic
fraction.

Everything except the solution to the tridiagonal equation is vectorized so, for
now, input arrays must be (Z x N) for N = 1 only. As part of this limitation,
there is a check that the bedrock depth is a scalar.

NOTES:

1. Instead of tracking ice and liquid water content separately, VWC generally
    refers to the total (ice plus liquid) volumetric water content. The liquid
    water content is determined based on the ice fraction, when needed.
  
2. Hydraulic condutivity (k) is always measured at the bottom interface of a
    soil layer.

# Parameters
- `pft::Int`: Plant functional type (PFT) code
- `soc::Array{Float64, 2}`: Areal soil organic carbon (SOC) content (g C m-2) in
  each layer; should be a (Z x 1) array
- `sand::Union{Float64, Array{Float64, 2}}`: Sand content of soil, as proportion
  on [0,1], a (Z x 1) array
- `clay::Union{Float64, Array{Float64, 2}}`: Clay content of soil, as proportion
  on [0,1], a (Z x 1) array
- `porosity::Union{Float64, Array{Float64, 2}}`: Total porosity of the soil;
  should be a scalar or a (Z x 1) array on [0,1]
- `depths_m::Array{Float64, 2}`: Depths of each soil layer's bottom interface,
  in meters; should be a (Z x 1) array with negative values below the surface
"""
@with_kw mutable struct SoilProfile{FT<:Real}
  pft::Int = 1
  n_layers::Int = 10

  "Depths of each soil layer's bottom interface, in meters; with negative values below the surface"
  depths_m::Vector{FT} = zeros(n_layers)

  "Areal soil organic carbon (SOC) content (g C m-2) in each layer"
  soc::Vector{FT} = zeros(n_layers)
  sand::Vector{FT} = zeros(n_layers)
  clay::Vector{FT} = zeros(n_layers)
  porosity::Vector{FT} = zeros(n_layers)

  "Depth to bedrock (m); not currently used"
  bedrock::Float64 = NaN
  "Topographic slope (degrees)"
  slope::Float64 = NaN

  frac_clay::Vector{FT} = clay
  frac_sand::Vector{FT} = sand
  frac_organic::Vector{FT} = zeros(n_layers)
  
  "in mm"
  thickness_mm::Vector{FT} = zeros(n_layers)
  
  z::Vector{FT} = depth2z(depths_m) .* 1e3
  z_bedrock::Float64 = -abs(bedrock) .* 1e3
  z_node::Vector{FT} = zeros(n_layers)
  # zeros::Vector{FT} = zeros(n_layers)
  params::Dict{Symbol,Float64} = Dict(:ksat_om => 1e-1, :alpha => 3)
end


function SoilProfile(
  soc::Vector{FT},
  sand::Vector{FT},
  clay::Vector{FT},
  porosity::Vector{FT}, slope::FT, depths_m::Vector{FT}) where {FT}

  @assert all(depths_m .< 0) "Depths should be defined as negative downward from the soil surface"

  frac_organic = (soc ./ abs.(depths_m)) ./ SOCC_MAX
  @assert maximum(frac_organic) < 1 "Organic fraction > 1.0; check units of soil organic carbon"

  thickness_mm = FT(depths_m .- vcat(0.0, depths_m[1:end-1])) .* 1e3

  z_node = (depths_m .* 1e3) .- thickness_mm ./ 2 # not negative
  n_layers = length(depths_m)

  return SoilProfile{FT}(;pft, soc, sand, clay, porosity, bedrock, slope, depths_m, frac_clay, frac_organic, frac_sand, thickness_mm, z, z_bedrock, z_node, n_layers, params)
end

depth2z(depth::Vector{FT}) where {FT<:Real} = FT(depth .- vcat(FT(0.0), depth[1:end-1]))


export SoilProfile
