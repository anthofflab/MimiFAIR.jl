include(joinpath(dirname(@__FILE__), "fair.jl"))

#Construct and run a base version of the FAIR model (default is 250 timesteps)
try
    base = constructfair()
    run(base)
    if isa(base, Mimi.Model)
        print("The base version worked")
    end
catch error
    print(error)
end

#Construct and run version with 100 timesteps
try
    less_steps = constructfair(nsteps=100)
    run(less_steps)
    if isa(less_steps, Mimi.Model)
        print("The version with less timesteps worked")
    end
catch error2
    print(error2)
end