module MimiFAIR

using Mimi
using DataFrames
using CSVFiles
using NLsolve

include("components/co2_cycle.jl")
include("components/radiativeforcing.jl")
include("components/temperature.jl")

export getfair

function getfair()

    global scenario = "rcp8.5"
    global nsteps = 736
    global start_year = 1765

    global emissions_datafile = joinpath(dirname(@__FILE__),"../data/rcp_scenarios/", scenario*"_emissions.csv")
    global forcing_datafile = joinpath(dirname(@__FILE__),"../data/rcp_scenarios/", scenario*"_forcings.csv")

    m = Model()
    set_dimension!(m, :time, start_year:1:start_year + nsteps - 1)

    # ---------------------------------------------
    # Add components to model
    # ---------------------------------------------
    add_comp!(m, co2_cycle, :co2_cycle)
    add_comp!(m, radiativeforcing, :radiativeforcing)
    add_comp!(m, temperature, :temperature)

    # ---------------------------------------------
    # Read in data
    # ---------------------------------------------
    emissions_data  = load(emissions_datafile, commentchar='#') |> DataFrame
    forcing_data    = load(forcing_datafile, commentchar='#') |> DataFrame

    # Find index for start year and subset data
    start_index     = findall(emissions_data.Year .== start_year)[1]
    emissions_data  = emissions_data[start_index:(start_index + nsteps-1), :]
    forcing_data    = forcing_data[start_index:(start_index + nsteps-1), :]

    # Create CO2 emissions variable and non-CO2 radiative forcing variable
    E   = (emissions_data.FossilCO2 + emissions_data.OtherCO2)
    Fext= forcing_data.SOLAR_RF + forcing_data.VOLCANIC_ANNUAL_RF + forcing_data.TOTAL_ANTHRO_RF - forcing_data.CO2_RF


    # ---------------------------------------------
    # Set component parameters
    # ---------------------------------------------

    #  CARBON CYCLE
    set_param!(m, :co2_cycle, :CO2_0, 278.0)
    set_param!(m, :co2_cycle, :r0, 35.0)
    set_param!(m, :co2_cycle, :rC, 0.019)
    set_param!(m, :co2_cycle, :rT, 4.165)
    set_param!(m, :co2_cycle, :iIRF_max, 97.0)
    set_param!(m, :co2_cycle, :a, [0.2173, 0.2240, 0.2824, 0.2763])
    set_param!(m, :co2_cycle, :τ, [10.0^6, 394.4, 36.54, 4.304])
    set_param!(m, :co2_cycle, :E, E)
    set_param!(m, :co2_cycle, :gtc2ppm, 2.1289)

    # RADIATIVE FORCING
    set_param!(m, :radiativeforcing, :C0, 278.0)
    set_param!(m, :radiativeforcing, :F2x, 3.74)
    set_param!(m, :radiativeforcing, :Fext, Fext)

    # TEMPERATURE
    set_param!(m, :temperature, :d, [239.0, 4.1])
    set_param!(m, :temperature, :q, [0.33, 0.41])
    set_param!(m, :temperature, :F2x, 3.74)

    # -----------------------------------------------
    # Create necessary connections between components
    # -----------------------------------------------
    connect_param!(m, :radiativeforcing, :C, :co2_cycle, :C)
    connect_param!(m, :temperature, :F, :radiativeforcing, :F)
    # Note: offset=1 => dependence is on on prior timestep, i.e., not a cycle
    connect_param!(m, :co2_cycle, :T, :temperature, :T)

    return m
end

end #module
