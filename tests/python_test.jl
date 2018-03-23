using PyCall
using Base.Test
using DataFrames

#= This file compares annual time series of temperature and CO2 produced by the julia and
   Python implementations of the FAIR model.  Current settings use the CO2 emissions (fossil
   + other sources) from the RCP8.5 scenario and set the non-CO2 radiative forcing to zero
   for all time periods.  Julia runs with scaling-factors saved from the Python version to
   allow for easier comparison. =#

#---------------------------------------------------------------------------------------------------
# Variable settings to change for different model runs
#---------------------------------------------------------------------------------------------------

start_year  = 1765      #First year to read data from emissions/forcing scenario
nsteps      = 736       #Number of time steps to run model for (one timestep = one year)
tolerance   = 1.0e-6   #Acceptable tolerance for differences between julia and Python results


#---------------------------------------------------------------------------------------------------
# Get Python version of Fair
#---------------------------------------------------------------------------------------------------

push!(pyimport("sys")["path"], dirname(@__FILE__))
@pyimport fair_python_version as FairPy


#---------------------------------------------------------------------------------------------------
# Read in emissions data and state-dependent scaling factors
#---------------------------------------------------------------------------------------------------

#Read data
emissions_data  = readtable(joinpath(dirname(@__FILE__), "../data/rcp_scenarios/rcp8.5_emissions.csv"), allowcomments=true)

# Find index for start year and subset data
start_index     = find(emissions_data[:Year] .== start_year)[1]
emissions_data  = emissions_data[start_index:(start_index + nsteps-1), :]

# Create CO2 emissions variable and non-CO2 radiative forcing variable
E = (emissions_data[:FossilCO2] + emissions_data[:OtherCO2])


#---------------------------------------------------------------------------------------------------
# Run julia and Python models and compare results
#---------------------------------------------------------------------------------------------------

#Calculate temperature and CO2 concentrations for Python FAIR
py_co2, py_temp = FairPy.fair_scm(in_driver = convert(Array,E))

#remove old code that ran Julia model
#=
#Construct base version of model
include(joinpath(dirname(@__FILE__), "../src/fair.jl"))

#Set number of time steps
julia_fair = constructfair(nsteps = nsteps, start_year = start_year)

#Set scaling factor, emissions, and non-co2 radiative forcing parameters
setparameter(julia_fair, :carboncycle, :E, E)
#setparameter(julia_fair, :carboncycle, :α, scale_factors)
setparameter(julia_fair, :radiativeforcing, :Fext, zeros(nsteps))

#Run model
run(julia_fair)
=#

#new code to run Julia model
## TODO:  do we need to set Fext to zeros here?
using Mimi
include(joinpath(dirname(@__FILE__), "../src/fair.jl"))
using fair
julia_fair = run(fair)

#---------------------------------------------------------------------------------------------------
# Run tests to compare differences in temperature and CO2 timeseries for Python and julia versions
#---------------------------------------------------------------------------------------------------

# Will throw an error if difference between two versions is greater than tolerance at any timestep
@test_approx_eq_eps maxabs(py_temp .- julia_fair[:temperature, :T]) 0. tolerance
@test_approx_eq_eps maxabs(py_co2 .- julia_fair[:carboncycle, :C]) 0. tolerance
