using Test
using Random
using LinearAlgebra
using Statistics

# include("simsoil/core.jl")
const DEPTHS = [-0.05, -0.1, -0.2, -0.4, -0.75, -1.5]
const SOC_RATIOS = [1, 1.77, 2.86, 4.63, 6.46, 9.56] .* 3138

function setup()
    x = 0:364
    Random.seed!(9)
    influx = 1e-5 .+ (rand(length(x)) .* 5e-6) .+ 1e-5 .* sin.(x .* (π / 270))
    influx_low = 2e-6 .+ (rand(length(x)) .* 5e-7) .+ 1e-6 .* sin.(x .* (π / 270))
    cold = 270 .+ (rand(length(x)) .* 3) .+ 5 .* sin.(x .* (π / length(x)))
    temp_profile_cold = vcat([(-1, 0, 1, 1, 2, 2)[i] .+ cold for i in 1:length(DEPTHS)]...)
    warm = 290 .+ rand(length(x)) .+ 5 .* sin.(x .* (π / length(x)))
    temp_profile_warm = vcat([(-1, 0, 1, 1, 2, 2)[i] .+ warm for i in 1:length(DEPTHS)]...)
    potential_transpiration = 1e-4 .+ (4e-5 .* rand(365)) .+ 1e-4 .* sin.(x .* (π / 290))
    return influx, influx_low, temp_profile_cold, temp_profile_warm, potential_transpiration
end

influx, influx_low, temp_profile_cold, temp_profile_warm, potential_transpiration = setup()

defaults = Dict(
    :soc => SOC_RATIOS,
    :depths_m => DEPTHS,
    :bedrock => -3
)
SoilProfile(1, sand = 0.6, clay = 0.1, porosity = 0.3, slope = 0.01, defaults...)

soils = [
    SoilProfile(1, sand = 0.6, clay = 0.1, porosity = 0.3, slope = 0.01, defaults...),
    SoilProfile(1, sand = 0.6, clay = 0.1, porosity = 0.5, slope = 0.01, defaults...),
    SoilProfile(1, sand = 0.1, clay = 0.6, porosity = 0.3, slope = 0.01, defaults...),
    SoilProfile(1, sand = 0.1, clay = 0.6, porosity = 0.5, slope = 0.01, defaults...)
]

models = [InfiltrationModel(s) for s in soils]

models_drained = [
    InfiltrationModel(SoilProfile(1, sand = 0.1, clay = 0.6, porosity = 0.3, slope = 2, defaults...)),
    InfiltrationModel(SoilProfile(1, sand = 0.1, clay = 0.6, porosity = 0.3, slope = 2, defaults...))
]

models_transpiration = [
    InfiltrationModel(SoilProfile(4, sand = 0.6, clay = 0.1, porosity = 0.35, slope = 1, defaults...)),
    InfiltrationModel(SoilProfile(5, sand = 0.6, clay = 0.1, porosity = 0.35, slope = 1, defaults...))
]

@testset "SoilWaterInfiltrationTestSuite" begin
    @testset "test_soils_with_warm_temps" begin
        init_vwc = ones(length(DEPTHS)) .* 0.15
        f_wet = ones(length(influx)) .* 0.2
        for (i, model) in enumerate(models)
            results, _, _ = run(model, init_vwc, temp_profile_warm, nothing, influx, f_wet, 7200, adaptive = false)
            @test size(results) == (length(DEPTHS), length(f_wet))
            @test round(mean(results), digits=3) == (0.222, 0.315, 0.256, 0.323)[i]
            @test round(minimum(results), digits=3) == (0.150, 0.150, 0.150, 0.150)[i]
            @test round(maximum(results), digits=3) == (0.300, 0.500, 0.300, 0.422)[i]
            @test round(var(results), digits=4) == (0.0023, 0.0115, 0.0019, 0.0109)[i]
        end
    end

    @testset "test_soils_with_cold_temps" begin
        init_vwc = ones(length(DEPTHS)) .* 0.15
        f_wet = ones(length(influx)) .* 0.2
        for (i, model) in enumerate(models)
            if i > 0
                break
            end
            results, _, _ = run(model, init_vwc, temp_profile_cold, nothing, influx, f_wet, 1800, adaptive = true, ltol = 1e-3)
            @test size(results) == (length(DEPTHS), length(f_wet))
            @test round(mean(results), digits=3) == 0.217
            @test round(minimum(results), digits=3) == 0.148
            @test round(maximum(results), digits=3) == 0.336
            @test round(var(results), digits=4) == 0.0024
        end
    end

    @testset "test_soils_with_transpiration_with_warm_temps" begin
        init_vwc = ones(length(DEPTHS)) .* 0.15
        f_wet = ones(length(influx)) .* 0.1
        for (i, model) in enumerate(models_transpiration)
            results, _, _ = run(model, init_vwc, temp_profile_warm, potential_transpiration, influx, f_wet, 7200, adaptive = false)
            @test size(results) == (length(DEPTHS), length(f_wet))
            @test round(mean(results), digits=3) == (0.104, 0.107)[i]
            @test round(minimum(results), digits=3) == (0.078, 0.078)[i]
            @test round(maximum(results), digits=3) == (0.150, 0.150)[i]
            @test round(var(results), digits=5) == (0.00065, 0.0007)[i]
            if i == 0
                @test all(round.(results[:, end], digits=3) .== [0.119, 0.097, 0.079, 0.078, 0.086, 0.133])
            elseif i == 1
                @test all(round.(results[:, end], digits=3) .== [0.112, 0.087, 0.078, 0.081, 0.118, 0.143])
            end
        end
    end
end

@testset "InfiltrationModelTestSuite" begin
    depths = [-0.05, -0.1, -0.2, -0.4, -0.75, -1.5]
    soc = [1, 1.77, 2.86, 4.63, 6.46, 9.56] .* 3138
    sm_profile = [0.1374, 0.1368, 0.1357, 0.1336, 0.1297, 0.1270]
    temp_profile = [276, 277, 280, 284, 287, 290]

    @testset "test_infiltration_model_step_daily" begin
        model = InfiltrationModel(SoilProfile(1, soc, sand = 0.6, clay = 0.1, porosity = 0.4, bedrock = -1.8, slope = 0, depths_m = depths))
        final_vwc, de = step_daily(model, sm_profile, temp_profile, 0, 1e-6, 0.1, 3600)
        @test round(sum(final_vwc), digits=4) == 0.7989
        @test round(final_vwc[1], digits=4) == 0.1346
        @test length(de) == 23
        @test round(maximum(de), digits=5) == 0.00838
        @test round(mean(de), digits=5) == -0.00061

        final_vwc, de = step_daily(model, sm_profile, temp_profile, 0, 1e-6, 0.1, 450)
        @test round(sum(final_vwc), digits=4) == 0.7989
        @test round(final_vwc[1], digits=4) == 0.1346
        @test round(maximum(de), digits=5) == 0.00108
        @test round(mean(de), digits=5) == -8e-5
    end
end

@testset "SoilProfileTestSuite" begin
    depths = [-0.05, -0.1, -0.2, -0.4, -0.75, -1.5]
    soc = [1, 1.77, 2.86, 4.63, 6.46, 9.56] .* 3138
    sm_profile = [0.1374, 0.1368, 0.1357, 0.1336, 0.1297, 0.1270]
    temp_profile = [276, 277, 280, 284, 287, 290]

    @testset "test_constants" begin
        model = SoilProfile(1, soc, sand = 0.6
