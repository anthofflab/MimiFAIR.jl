using Mimi
using NLsolve

@defcomp carboncycle begin

    C0  = Parameter()               # Pre-industrial CO2 concentrations
    r0  = Parameter()               # Pre-industrial iIRF100
    rC  = Parameter()               # Increase in iIRF100 with cumulative carbon uptake (yr/GtC)
    rT  = Parameter()               # Increase in iIRF100 with warming (yr/K)
    a   = Parameter(index=[4])      # Fraction of emissions entering carbon pool (geological reabsorption[1], deep ocean invasion/equilibration[2], biospheric uptake/ocean thermocline invasion[3], rapid biospheric uptake/ocean mixed-layer invasion[4])
    τ   = Parameter(index=[4])      # Decay time constant for each carbon pool in 'a'
    E   = Parameter(index=[time])   # Annual CO2 emissions (in units of ppm/year with 1ppm = 2.12GtC)
    T   = Parameter(index=[time])   # Global mean surface temperature anomaly

    α   = Variable(index=[time])    # State-dependent carbon uptake scaling factor
    C   = Variable(index=[time])    # Total atmospheric CO2 concentrations
    Cacc= Variable(index=[time])    # Accumulated perturbation carbon stock (amount of emitted carbon that no longer resides in the atmosphere)
    R   = Variable(index=[time,4])  # CO2 concentration in each carbon pool
end

function run_timestep(s::carboncycle, t::Int)
    p = s.Parameters
    v = s.Variables

    if t==1

        # TODO Is this correct?
        for i=1:4
            v.R[t,i] = p.a[i] * p.E[t]
        end

        # Initical CO2 concentration
        v.C[t] = p.C0

        # Initial carbon stock perturbation
        v.Cacc[t] = sum(p.E[1:t]) - (v.C[t] - p.C0)

        # Initial state-dependent scaling factor
        v.α[1] = 0.1
    else

        # Use iIRF100 equation to solve for α with nlsolve
        iIRFT100 = p.r0 + p.rC * v.Cacc[t-1] + p.rT * p.T[t-1]

        # Function to pass to nlsolve
        function f!(x, fvec)
            α = x[1]
            fvec[1] = -iIRFT100 + sum(α .* p.a .* p.τ .* (1 .- exp.(-100 ./ α ./ p.τ)))
        end

        # Calculate α
        res = nlsolve(f!, [v.α[t-1]], autodiff=true)
        if !converged(res)
            error("Couldn't find a solution for α.")
        end
        v.α[t] = res.zero[1]

        # Calculate CO2 concentrations in each pool
        for i=1:4
            v.R[t,i] = v.R[t-1,i] + p.a[i] * p.E[t] - v.R[t-1,i] / v.α[t] / p.τ[i]
        end

        # Calculate total CO2 concentrations
        v.C[t] = p.C0 + sum(v.R[t,:])

        # Calculate accumulated perturbation of carbon stock
        # TODO Change to a recursive equation, should be a bit faster
        v.Cacc[t] = sum(p.E[1:t]) - (v.C[t] - p.C0)
    end
end
