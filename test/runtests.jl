using Test
using CSVFiles
using DataFrames
using MimiFAIR

@testset "FAIR" begin

#------------------------------------------------------------------------------
#   1. Carry out test to check the model runs.
#------------------------------------------------------------------------------

    @testset "FAIR-model" begin

    m = MimiFAIR.get_model()
    run(m)

    end # FAIR model run test.

#------------------------------------------------------------------------------
#   2. Carry out tests to make sure Mimi-FAIR matches the Python version.
#------------------------------------------------------------------------------
#=
    @testset "FAIR-Python" begin

    Precision = 1.0e-5

    m = MimiFAIR.get_model(rcp_scenario = "RCP85", start_year = 1765, end_year = 2500)
    run(m)

    # Load Python version of FAIR output.
    python_co2_conc      = DataFrame(load(joinpath(@__DIR__, "..", "data", "validation_data", "rcp85_python_fair_concentrations.csv"))).co2
    python_temperature   = DataFrame(load(joinpath(@__DIR__, "..", "data", "validation_data", "rcp85_python_fair_temperature.csv"))).temperature

    # Run tests for global temperature anomaly and atmospheric CO₂ concentrations.
    @test m[:co2_cycle, :C]            ≈ python_co2_conc atol = Precision
    @test m[:temperature, :T]          ≈ python_temperature atol = Precision

    end # FAIR Python comparison test.
=#
end # All FAIR tests.
