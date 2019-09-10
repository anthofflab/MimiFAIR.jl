# -----------------------------------------------------------
# Global temperature change from a given radiative forcing
# -----------------------------------------------------------

@defcomp temperature begin

    F2x = Parameter()               # Radiative forcing due to doubling carbon dioxide (Wm⁻²).
    d   = Parameter(index=[2])      # Thermal response timescales: [1] Thermal equilibration of deep ocean & [2] Thermal admustment of upper ocean (years).
    q   = Parameter(index=[2])      # Raditive forcing coefficient: [1] Thermal equilibration of deep ocean & [2] Thermal admustment of upper ocean (K W⁻¹m²).
    F   = Parameter(index=[time])   # Total radiative forcing (Wm⁻²).

    T   = Variable(index=[time])    # Global mean surface temperature anomaly (K).
    Tj  = Variable(index=[time,2])  # Temperature change for two thermal response times (K).

    function run_timestep(p, v, d, t)

        #Calculate Temperatures.
        if is_first(t)

            #Set initial temperature change for two response times.
            for j=1:2
                v.Tj[t,j] = (p.q[j] / p.d[j]) * p.F[t]
            end

            #Sum thermal response boxes to get initial temperature anomaly.
            v.T[t] = sum(v.Tj[t,:])

        else

            #Calculate temperature change for the two different thermal response times.
            for j=1:2
                v.Tj[t,j] = v.Tj[t-1,j] * exp((-1.0)/p.d[j]) + p.q[j] * (1 - exp((-1.0)/p.d[j])) * p.F[t]
            end

            #Calculate global mean surface temperature anomaly.
            v.T[t] = sum(v.Tj[t,:])
        end
    end
end
