using Mimi

include("fair.jl")
using .Fair

m = getfair()
run(m)

explore(m)
