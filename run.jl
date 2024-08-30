using Pkg
Pkg.activate(".")
using SoilSim

depths_m = [-0.05, -0.1, -0.2, -0.4, -0.75, -1.5]
nlayer = length(depths_m)
s = SoilProfile{Float32}(; depths_m, nlayer)
display(s)
