using Mimi

@defcomp temperature begin

    F2x = Parameter()               # Forcing due to CO2 doubling (w/m2)
    d   = Parameter(index=[thermresponse])      # Thermal response timescales (Thermal equilibration of deep ocean[1] & Thermal admustment of upper ocean[2]) - years.
    F   = Parameter(index=[time])   # Total radiative forcing (w/m2).

    q   = Parameter(index=[thermresponse])      # Forcing coefficient (Thermal equilibration of deep ocean[1] & Thermal admustment of upper ocean[2]) - K/(w/m^2).
    Tj  = Variable(index=[time,thermresponse])  # Temperature change for two thermal response times.
    T   = Variable(index=[time])    # Global mean surface temperature anomaly.

    function init(p, v, d)
        #Calculate Temperatures.
        t = 1
        #Set initial temperature change for two response times to zero.
        for j=1:2
            v.Tj[t,j] = 0.0
        end

        #Set initial surface temperature anomaly to zero.
        v.T[t] = sum(v.Tj[t,:])
    end

    function run_timestep(p, v, d, t)

        #Calculate Temperatures.
        if t > 1
            #Calculate temperature change for the two different thermal response times.
            for j=1:2
                v.Tj[t,j] = v.Tj[t-1,j] * exp((-1.0)/p.d[j]) + 0.5 *p.q[j]*(p.F[t-1] + p.F[t]) * (1 - exp((-1.0)/p.d[j]))
            end

            #Calculate global mean surface temperature anomaly
            v.T[t] = sum(v.Tj[t,:])
        end
    end
end
