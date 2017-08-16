using Mimi
using DataFrames

include(joinpath(dirname(@__FILE__), "carboncycle.jl"))
include(joinpath(dirname(@__FILE__), "radiativeforcing.jl"))
include(joinpath(dirname(@__FILE__), "temperature.jl"))

function constructfair(;nsteps=736, scenario="rcp8.5", start_year = 1765)

    m = Model()
    setindex(m, :time, nsteps)

    # ---------------------------------------------
    # Add components to model
    # ---------------------------------------------
    addcomponent(m, carboncycle)
    addcomponent(m, radiativeforcing)
    addcomponent(m, temperature)

    # ---------------------------------------------
    # Read in data
    # ---------------------------------------------
    emissions_data  = readtable(joinpath(dirname(@__FILE__),"../data/rcp_scenarios/", scenario*"_emissions.csv"), allowcomments=true)
    forcing_data    = readtable(joinpath(dirname(@__FILE__),"../data/rcp_scenarios/", scenario*"_forcings.csv"), allowcomments=true)

    # Find index for start year and subset data
    start_index     = find(emissions_data[:Year] .== start_year)[1]
    emissions_data  = emissions_data[start_index:(start_index + nsteps-1), :]
    forcing_data    = forcing_data[start_index:(start_index + nsteps-1), :]

    # Create CO2 emissions variable (convert to ppm/year with 1ppm = 2.12 GtC) and non-CO2 radiative forcing variable
    E   = (emissions_data[:FossilCO2] + emissions_data[:OtherCO2])
    Fext= forcing_data[:SOLAR_RF] + forcing_data[:VOLCANIC_ANNUAL_RF] + forcing_data[:TOTAL_ANTHRO_RF] - forcing_data[:CO2_RF]

    # ---------------------------------------------
    # Set component parameters
    # ---------------------------------------------
    setparameter(m, :carboncycle, :C0, 278.0)
    setparameter(m, :carboncycle, :r0, 35.0)
    setparameter(m, :carboncycle, :rC, 0.02)
    setparameter(m, :carboncycle, :rT, 4.5)
    setparameter(m, :carboncycle, :a, [0.2173, 0.2240, 0.2824, 0.2763])
    setparameter(m, :carboncycle, :τ, [10.0^6, 394.4, 36.54, 4.304])
    setparameter(m, :carboncycle, :E, E)

    setparameter(m, :radiativeforcing, :C0, 278.0)
    setparameter(m, :radiativeforcing, :F2x, 3.74)
    setparameter(m, :radiativeforcing, :Fext, Fext)

    setparameter(m, :temperature, :d, [239.0, 4.1])
    setparameter(m, :temperature, :F2x, 3.74)

    # -----------------------------------------------
    # Create necessary connections between components
    # -----------------------------------------------
    connectparameter(m, :radiativeforcing, :C, :carboncycle, :C)
    connectparameter(m, :temperature, :F, :radiativeforcing, :F)
    connectparameter(m, :carboncycle, :T, :temperature, :T)

    # Return constructed model
    return m
end
