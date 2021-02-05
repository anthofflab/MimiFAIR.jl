# --------------------------------------------------
# Carbon cycle
# --------------------------------------------------

@defcomp co2_cycle begin

    CO2_0     = Parameter()              # Pre-industrial carbon dioxide concentrations (ppm).
    r0        = Parameter()              # Pre-industrial 100-year integrated impulse response function (iIRF100).
    rC        = Parameter()              # Increase in iIRF100 with cumulative carbon uptake (yr GtC⁻¹).
    rT        = Parameter()              # Increase in iIRF100 with warming (yr K⁻¹).
    gtc2ppm   = Parameter()              # Conversion factor between GtC and ppm.
    iIRF_max  = Parameter()              # Maximum value for iIRF100.
    a         = Parameter(index=[4])     # Fraction of emissions entering carbon pool (geological reabsorption[1], deep ocean invasion/equilibration[2], biospheric uptake/ocean thermocline invasion[3], rapid biospheric uptake/ocean mixed-layer invasion[4]).
    τ_CO₂     = Parameter(index=[4])     # Decay time constant for each carbon pool in 'a'.
    E_ox_CH₄  = Parameter(index=[time])  # Emissions from methane that has been oxidized to carbon dioxide (GtC yr⁻¹).
    E_CO₂     = Parameter(index=[time])  # Annual carbon dioxide emissions (GtC yr⁻¹).
    T         = Parameter(index=[time])  # Global mean surface temperature anomaly (K).

    E         = Variable(index=[time])   # Total carbon dioxide emissions that account for oxidized methane (GtC yr⁻¹).
    iIRFT100  = Variable(index=[time])   # 100-year integrated impulse response function.
    α         = Variable(index=[time])   # State-dependent carbon uptake scaling factor.
    C         = Variable(index=[time])   # Total atmospheric carbon dioxide concentrations (ppm).
    Cacc      = Variable(index=[time])   # Accumulated perturbation carbon stock (amount of emitted carbon that no longer resides in the atmosphere).
    R         = Variable(index=[time,4]) # Carbon dioxide concentration in each of the four carbon pools.

    function run_timestep(p, v, d, t)
        if is_first(t)

            # Calculate initial sum of CO₂ emissions and emissions from oxidized methane.
            v.E[t] = p.E_CO₂[t] + p.E_ox_CH₄[t]

            # "Initialise the carbon pools to be correct for first timestep in numerical method" -Python Version of FAIR
            for i=1:4
                v.R[t,i] = p.a[i] * v.E[t] / p.gtc2ppm
            end

            # Initial guess for state-dependet scaling factor.
            v.α[t] = 1.0

            # Initical change and total CO₂ concentration.
            v.C[t] = sum(v.R[t,:]) + p.CO2_0

            # Initial carbon stock perturbation.
            v.Cacc[t] = 0.0

        else

            # Calculate sum of CO₂ emissions and emissions from oxidized methane.
            v.E[t] = p.E_CO₂[t] + p.E_ox_CH₄[t]

            # Use iIRF100 equation to solve for α (bounded above by max value).
            v.iIRFT100[t] = min(p.r0 + p.rC * v.Cacc[t-1] + p.rT * p.T[t-1], p.iIRF_max)

            # Create function to pass to nlsolve to find α.
            function f!(fvec, x)
                fvec[1] = -v.iIRFT100[t] + sum(x[1] .* p.a .* p.τ_CO₂ .* (1 .- exp.(-100 ./ x[1] ./ p.τ_CO₂)))
            end

            # Solve for α.
            res = nlsolve(f!, [v.α[t-1]], autodiff=:forward)
            if !converged(res)
                error("Couldn't find a solution for α.")
            end
            v.α[t] = res.zero[1]

            #Calculate updated carbon cycle time constants and carbon concentrations in each pool.
            for i=1:4
                b_new = p.τ_CO₂[i] * v.α[t]
                v.R[t,i] = v.R[t-1,i]  * exp((-1.0/b_new)) + p.a[i] * v.E[t] / p.gtc2ppm
            end

            # Calculate atmospheric carbon dioxide concentrations.
            v.C[t] = sum(v.R[t,:]) + p.CO2_0

            # Calculate accumulated perturbation of carbon stock.
            v.Cacc[t] = v.Cacc[t-1] + 0.5*(v.E[t] + v.E[t-1]) - (v.C[t] -v.C[t-1]) * p.gtc2ppm
        end
    end
end
