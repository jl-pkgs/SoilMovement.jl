module SoilSim

using Parameters
using LinearAlgebra
using Statistics

const SECS_PER_DAY = 86400

to_float(x::AbstractArray) = Float32.(x)


include("SoilProfile.jl")
include("SoilProperty.jl")
include("solve_vwc.jl")
include("solve_sink.jl")

include("Infiltration.jl")


export solve_θ, solve_sink, SoilProfile, SoilProperty, Infiltration


end # module SimSoil
