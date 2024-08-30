module SimSoil

using Parameters
using LinearAlgebra
using Statistics

const SECS_PER_DAY = 86400

to_float(x::AbstractArray) = Float32.(x)


include("SoilProfile.jl")
include("SoilProperty.jl")
include("solve_vwc.jl")
include("Infiltration.jl")


end # module SimSoil
