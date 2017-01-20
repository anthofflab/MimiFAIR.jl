include(joinpath(dirname(@__FILE__), "fair.jl"))

#Construct a base version of the FAIR model
base = constructfair()

#Run model
run(base)
