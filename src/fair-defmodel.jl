module fair

using Mimi
using DataFrames

include(joinpath(dirname(@__FILE__), "carboncycle.jl"))
include(joinpath(dirname(@__FILE__), "radiativeforcing.jl"))
include(joinpath(dirname(@__FILE__), "temperature.jl"))

export fair
const global scenario = "rcp8.5"
const global nsteps = 736
const global start_year = 1765

const global emissions_datafile = joinpath(dirname(@__FILE__),"../data/rcp_scenarios/", scenario*"_emissions.csv")
const global forcing_datafile = joinpath(dirname(@__FILE__),"../data/rcp_scenarios/", scenario*"_forcings.csv")

@defmodel fair begin

    index[time] = start_year:1:start_year + nsteps - 1

    # ---------------------------------------------
    # Read in data
    # ---------------------------------------------
    emissions_data  = readtable(emissions_datafile, allowcomments=true)
    forcing_data    = readtable(forcing_datafile, allowcomments=true)

    # Find index for start year and subset data
    start_index     = find(emissions_data[:Year] .== start_year)[1]
    emissions_data  = emissions_data[start_index:(start_index + nsteps-1), :]
    forcing_data    = forcing_data[start_index:(start_index + nsteps-1), :]

    # Create CO2 emissions variable and non-CO2 radiative forcing variable
    E   = (emissions_data[:FossilCO2] + emissions_data[:OtherCO2])
    Fext= forcing_data[:SOLAR_RF] + forcing_data[:VOLCANIC_ANNUAL_RF] + forcing_data[:TOTAL_ANTHRO_RF] - forcing_data[:CO2_RF]

    component(carboncycle)
    component(radiativeforcing)
    component(temperature)

    #  CARBON CYCLE COMPONENT
    carboncycle.C0  = 278.0
    carboncycle.r0  = 32.4
    carboncycle.rC  = 0.019
    carboncycle.rT  = 4.165
    carboncycle.a   = [0.2173, 0.2240, 0.2824, 0.2763]
    carboncycle.τ   = [10.0^6, 394.4, 36.54, 4.304]
    carboncycle.E   = E

    # RADIATIVE FORCING COMPONENT
    radiativeforcing.C0     = 278.0
    radiativeforcing.F2x    = 3.74
    radiativeforcing.Fext   = Fext

    # TEMPERATURE COMPONENT
    temperature.d   = [239.0, 4.1]
    temperature.q   = [0.33, 0.41]
    temperature.F2x = 3.74

    # CONNECTIONS
    carboncycle.C => radiativeforcing.C
    radiativeforcing.F => temperature.F
    # Note that dependence is on prior timestep ("[t-1]")
    temperature.T[t-1] => carboncycle.T

end

end #module
