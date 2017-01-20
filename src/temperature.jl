using Mimi

@defcomp temperature begin

    q   = Parameter(index=[2])      # Forcing coefficient (Thermal equilibration of deep ocean[1] & Thermal admustment of upper ocean[2]) - K/(w/m^2)
    d   = Parameter(index=[2])      # Thermal response timescales (Thermal equilibration of deep ocean[1] & Thermal admustment of upper ocean[2]) - years
    F   = Parameter(index=[time])   # Total radiative forcing (w/m2)

    Tj  = Variable(index=[time,2])  # Temperature change for two thermal response times
    T   = Variable(index=[time])    # Global mean surface temperature anomaly
end


function run_timestep(s::temperature, t::Int)
    p = s.Parameters
    v = s.Variables

    if t==1

        #Set initial temperature change for two response times to 0
        for j=1:2
            v.Tj[t,j] = 0.
        end

        #Set initial surface temperature anomaly to 0
        v.T[t] = sum(v.Tj[t,:])
    else

        #Calculate temperature change for the two different thermal response times
        for j=1:2
            v.Tj[t,j] = v.Tj[t-1,j] + (p.q[j] * p.F[t] - v.T[t-1]) / p.d[j]
        end

        #Calculate global mean surface temperature anomaly
        v.T[t] = sum(v.Tj[t,:])
    end
end
