using Mimi

@defcomp temperature begin

    F2x = Parameter()               # Forcing due to CO2 doubling (w/m2)
    d   = Parameter(index=[2])      # Thermal response timescales (Thermal equilibration of deep ocean[1] & Thermal admustment of upper ocean[2]) - years.
    F   = Parameter(index=[time])   # Total radiative forcing (w/m2).

    q   = Variable(index=[2])       # Forcing coefficient (Thermal equilibration of deep ocean[1] & Thermal admustment of upper ocean[2]) - K/(w/m^2).
    Tj  = Variable(index=[time,2])  # Temperature change for two thermal response times.
    T   = Variable(index=[time])    # Global mean surface temperature anomaly.
end


function run_timestep(s::temperature, t::Int)
    p = s.Parameters
    v = s.Variables

    #Calculate values needed for specific transient climate response (TCR) and Equilibrium Climate Sensitivity (ECS).
    # TODO: What does ecstcr represent? (taken from Python version of FAIR).
    ecstcr = [1.75, 2.5]

    #Calculate k (intermediate step used for calculating q).
    k = zeros(2)
    for j=1:2
        k[j] = 1.0 - (p.d[j]/70.0) * (1.0 - exp(-70.0/p.d[j]))
    end

    #Calculate forcing coefficients.
    v.q[1] = (1.0/p.F2x) * (1.0 / (k[1] - k[2])) * (ecstcr[1] - ecstcr[2] * k[2])
    v.q[2] = (1.0/p.F2x) * (1.0 / (k[1] - k[2])) * (ecstcr[2] * k[1] - ecstcr[1])

    #----------------------
    #Calculate Temperatures.
    #----------------------
    if t==1

        #Set initial temperature change for two response times to zero.
        for j=1:2
            v.Tj[t,j] = 0.0
        end

        #Set initial surface temperature anomaly to zero.
        v.T[t] = sum(v.Tj[t,:])
    else

        #Calculate temperature change for the two different thermal response times.
        for j=1:2
            v.Tj[t,j] = v.Tj[t-1,j] * exp((-1.0)/p.d[j]) + 0.5 *v.q[j]*(p.F[t-1] + p.F[t]) * (1 - exp((-1.0)/p.d[j]))
        end

        #Calculate global mean surface temperature anomaly
        v.T[t] = sum(v.Tj[t,:])
    end
end
