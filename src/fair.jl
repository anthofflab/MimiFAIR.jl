using Mimi

include(joinpath(dirname(@__FILE__), "carboncycle.jl"))
include(joinpath(dirname(@__FILE__), "radiativeforcing.jl"))
include(joinpath(dirname(@__FILE__), "temperature.jl"))

function constructfair(;nsteps=200)

    m = Model()
    setindex(m, :time, nsteps)

    # TODO Replace this with a proper emission and forcing loading routine
    E       = ones(200) .* 8.
    Fext    = zeros(200)

    # ---------------------------------------------
    # Add components to model
    # ---------------------------------------------
    addcomponent(m, carboncycle)
    addcomponent(m, radiativeforcing)
    addcomponent(m, temperature)

    # ---------------------------------------------
    # Read in data
    # ---------------------------------------------

    # TODO Add emissions & non-CO2 forcing data

    # ---------------------------------------------
    # Set component parameters
    # ---------------------------------------------
    setparameter(m, :carboncycle, :C0, 278.05158)
    setparameter(m, :carboncycle, :r0, 35.0)
    setparameter(m, :carboncycle, :rC, 0.02)
    setparameter(m, :carboncycle, :rT, 8.5)
    setparameter(m, :carboncycle, :a, [0.2173, 0.2240, 0.2824, 0.2763])
    setparameter(m, :carboncycle, :τ, [10.0^8, 394.4, 36.54, 4.304])
    setparameter(m, :carboncycle, :E, E)

    setparameter(m, :radiativeforcing, :C0, 278.05158)
    setparameter(m, :radiativeforcing, :F2x, 3.74/log(2))
    setparameter(m, :radiativeforcing, :Fext, Fext)

    setparameter(m, :temperature, :q, [0.46, 0.27])
    setparameter(m, :temperature, :d, [8.4, 409.5])

    # -----------------------------------------------
    # Create necessary connections between components
    # -----------------------------------------------
    connectparameter(m, :radiativeforcing, :C, :carboncycle, :C)
    connectparameter(m, :temperature, :F, :radiativeforcing, :F)
    connectparameter(m, :carboncycle, :T, :temperature, :T)

    # Return constructed model
    return m
end
