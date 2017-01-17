using Mimi
using Plots

include("fair.jl")

m = Model()

setindex(m, :time, 200)
addcomponent(m, fair)

setparameter(m, :fair, :a, [0.2173, 0.2240, 0.2824, 0.2763])
setparameter(m, :fair, :τ, [10.0^8, 394.4, 36.54, 4.304])
setparameter(m, :fair, :C0, 278.05158)
setparameter(m, :fair, :r0, 35.0)
setparameter(m, :fair, :rC, 0.02)
setparameter(m, :fair, :rT, 8.5)
setparameter(m, :fair, :F2x, 3.74/log(2))
setparameter(m, :fair, :c, [0.46, 0.27])
setparameter(m, :fair, :d, [8.4, 409.5])

# Replace this with a proper emission and forcing loading routine
E = ones(200) .* 8.
Fext = zeros(200)

setparameter(m, :fair, :E, E)
setparameter(m, :fair, :Fext, Fext)

run(m)

# Plot temperature
plot(m[:fair,:T])
