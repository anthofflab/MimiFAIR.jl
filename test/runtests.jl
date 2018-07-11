using Mimi
using DataTables
using Base.Test

include(joinpath(@__DIR__, "../src/fair.jl"))
using fair

@testset "FAIR" begin

#------------------------------------------------------------------------------
#   1. Run tests on the whole model
#------------------------------------------------------------------------------

@testset "FAIR-model" begin

m = fair.FAIR
run(m)

end #FAIR-model testset

#------------------------------------------------------------------------------
#   2. Run tests to make sure integration version (Mimi v0.5.0)
#   values match Mimi 0.4.0 values
#------------------------------------------------------------------------------

@testset "FAIR-integration" begin

Mimi.reset_compdefs()

Precision = 1.0e-11
nullvalue = -999.999

m = fair.FAIR
run(m)

for c in map(name, Mimi.compdefs(m)), v in Mimi.variable_names(m, c)
    
    #load data for comparison
    filepath = "../data/validation_data_v040/$c-$v.csv"        
    results = m[c, v]

    if typeof(results) <: Number
        validation_results = DataFrames.readtable(filepath)[1,1]
        
    else
        validation_results = convert(Array, DataFrames.readtable(filepath))

        #match dimensions
        if size(validation_results,1) == 1
            validation_results = validation_results'
        end

        #remove NaNs
        results[isnan.(results)] = nullvalue
        validation_results[isnan.(validation_results)] = nullvalue
        
    end
    @test results ≈ validation_results atol = Precision
    
end #for loop

end #FAIR-integration testset

end #FAIR testset
