using Mimi
using DataFrames
using Test

include("../src/fair.jl")
using .Fair

@testset "FAIR" begin

#------------------------------------------------------------------------------
#   1. Run tests on the whole model
#------------------------------------------------------------------------------

@testset "FAIR-model" begin

global m = getfair()
run(m)

end #FAIR-model testset

#------------------------------------------------------------------------------
#   2. Run tests to make sure integration version (Mimi v0.5.0)
#   values match Mimi 0.4.0 values
#------------------------------------------------------------------------------

@testset "FAIR-integration" begin

Precision = 1.0e-11
nullvalue = -999.999

global m = getfair()
run(m)

for c in map(name, Mimi.compdefs(m)), v in Mimi.variable_names(m, c)
    
    #load data for comparison
    filepath = joinpath(@__DIR__, "../data/validation_data_v040/$c-$v.csv")        
    results = m[c, v]

    if typeof(results) <: Number
        validation_results = DataFrames.readtable(filepath)[1,1]
        
    else
        validation_results = convert(Array, DataFrames.readtable(filepath))

        #remove NaNs
        results[ismissing.(results)] .= nullvalue
        results[isnan.(results)] .= nullvalue
        validation_results[isnan.(validation_results)] .= nullvalue

        #match dimensions
        if size(validation_results,1) == 1
            validation_results = validation_results'
        end
        
    end
    @test results ≈ validation_results atol = Precision
    
end #for loop

end #FAIR-integration testset

end #FAIR testset
