using Mimi
using NLsolve

const gtc2ppm = 2.123                   # Conversion factor between GtC and ppm

@defcomp carboncycle begin

    C0      = Parameter()               # Pre-industrial CO2 concentrations.
    r0      = Parameter()               # Pre-industrial iIRF100.
    rC      = Parameter()               # Increase in iIRF100 with cumulative carbon uptake (yr/GtC).
    rT      = Parameter()               # Increase in iIRF100 with warming (yr/K).
    a       = Parameter(index=[4])      # Fraction of emissions entering carbon pool (geological reabsorption[1], deep ocean invasion/equilibration[2], biospheric uptake/ocean thermocline invasion[3], rapid biospheric uptake/ocean mixed-layer invasion[4]).
    τ       = Parameter(index=[4])      # Decay time constant for each carbon pool in 'a'.
    E       = Parameter(index=[time])   # Annual CO2 emissions (in units of ppm/year with 1ppm = 2.12GtC).
    T       = Parameter(index=[time])   # Global mean surface temperature anomaly.

    #α      = Variable(index=[time])    # State-dependent carbon uptake scaling factor.
    ΔCO2    = Variable(index=[time])    # Change in atmospheric CO2 concentrations.
    C       = Variable(index=[time])    # Total atmospheric CO2 concentrations.
    Cacc    = Variable(index=[time])    # Accumulated perturbation carbon stock (amount of emitted carbon that no longer resides in the atmosphere).
    R       = Variable(index=[time,4])  # CO2 concentration in each carbon pool.

    # TODO Change scaling factor parameter back into variable and unmute relevant equations. Just setting as a parameter to compare with Python version.
    α       = Parameter(index=[time])
end

function run_timestep(s::carboncycle, t::Int)
    p = s.Parameters
    v = s.Variables

    if t==1

        # "Initialise the carbon pools to be correct for first timestep in numerical method" -Python Version of FAIR
        for i=1:4
            v.R[t,i] = p.a[i] * p.E[t] / gtc2ppm * 0.5
        end

        # Initical change and total CO2 concentration
        v.ΔCO2[t] = 0.0
        v.C[t] = p.C0

        # Initial carbon stock perturbation
        v.Cacc[t] = p.E[t]

    else

        # Use iIRF100 equation to solve for α with nlsolve
        iIRFT100 = p.r0 + p.rC * v.Cacc[t-1] + p.rT * p.T[t-1]

# TODO Uncomment once Python version matches

        # Function to pass to nlsolve
        #function f!(x, fvec)
        #    α = x[1]
        #    fvec[1] = -iIRFT100 + sum(α .* p.a .* p.τ .* (1 .- exp.(-100 ./ α ./ p.τ)))
        #end

        # Calculate α
        #res = nlsolve(f!, [v.α[t-1]], autodiff=true)
        #if !converged(res)
        #    error("Couldn't find a solution for α.")
        #end
        #v.α[t] = res.zero[1]

        #Calculate updated carbon cycle time constants and CO2 concentrations in each pool
        for i=1:4
            b_new = p.τ[i] * p.α[t]
            v.R[t,i] = v.R[t-1,i]  * exp((-1.0/b_new)) + 0.5 * p.a[i] * (p.E[t] + p.E[t-1]) / gtc2ppm
        end

        # Calculate change in and total CO2 concentrations
        v.ΔCO2[t] = sum(v.R[t,:])
        v.C[t] = v.ΔCO2[t] + p.C0

        # Calculate accumulated perturbation of carbon stock
        v.Cacc[t] = v.Cacc[t-1] + p.E[t] - (v.C[t] -v.C[t-1]) * gtc2ppm
    end
end
