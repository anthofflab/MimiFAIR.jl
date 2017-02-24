include("../src/fair.jl")

# TODO run this with all timesteps
m = constructfair(nsteps=100)
run(m)
