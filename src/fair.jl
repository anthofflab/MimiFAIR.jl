module fair

using Mimi
using DataFrames

include(joinpath(dirname(@__FILE__), "carboncycle.jl"))
include(joinpath(dirname(@__FILE__), "radiativeforcing.jl"))
include(joinpath(dirname(@__FILE__), "temperature.jl"))

export fair

#
# N.B. See fair-defmodel.jl for the @defmodel version of the following
#

const global scenario = "rcp8.5"
const global nsteps = 736
const global start_year = 1765

const global emissions_datafile = joinpath(dirname(@__FILE__),"../data/rcp_scenarios/", scenario*"_emissions.csv")
const global forcing_datafile = joinpath(dirname(@__FILE__),"../data/rcp_scenarios/", scenario*"_forcings.csv")

fair = Model()
set_dimension!(RICE, :time, start_year:1:start_year + nsteps - 1)

# ---------------------------------------------
# Add components to model
# ---------------------------------------------
addcomponent(fair, carboncycle, :carboncycle)
addcomponent(fair, radiativeforcing, :radiativeforcing)
addcomponent(fair, temperature :temperature)

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


# ---------------------------------------------
# Set component parameters
# ---------------------------------------------

#  CARBON CYCLE 
set_parameter!(fair, :carboncycle, :C0, 278.0)
set_parameter!(fair, :carboncycle, :r0, 32.4)
set_parameter!(fair, :carboncycle, :rC, 0.019)
set_parameter!(fair, :carboncycle, :rT, 4.165)
set_parameter!(fair, :carboncycle, :a, [0.2173, 0.2240, 0.2824, 0.2763])
set_parameter!(fair, :carboncycle, :τ, [10.0^6, 394.4, 36.54, 4.304])
set_parameter!(fair, :carboncycle, :E, E)

# RADIATIVE FORCING
set_parameter!(fair, :radiativeforcing, :C0, 278.0)
set_parameter!(fair, :radiativeforcing, :F2x, 3.74)
set_parameter!(fair, :radiativeforcing, :Fext, Fext)

# TEMPERATURE
set_parameter!(fair, :temperature, :d, [239.0, 4.1])
set_parameter!(fair, :temperature, :q, [0.33, 0.41])
set_parameter!(fair, :temperature, :F2x, 3.74)

# -----------------------------------------------
# Create necessary connections between components
# -----------------------------------------------
connect_parameter(fair, :radiativeforcing, :C, :carboncycle, :C, offset = 0)
connect_parameter(fair, :temperature, :F, :radiativeforcing, :F, offset = 0)
# Note: offset=1 => dependence is on on prior timestep, i.e., not a cycle
connect_parameter(fair, :carboncycle, :T, :temperature, :T, offset = 1)

add_connector_comps(fair)

end #module
