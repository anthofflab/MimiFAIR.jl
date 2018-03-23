using Mimi

@defcomp radiativeforcing begin

    C0  = Parameter()               # Pre-industrial CO2 concentrations
    F2x = Parameter()               # Forcing due to CO2 doubling (w/m2)
    Fext= Parameter(index=[time])   # Non-CO2 forcing (w/m2)
    C   = Parameter(index=[time])   # Total atmospheric CO2 concentrations

    F   = Variable(index=[time])    # Total radiative forcing (w/m2)

    function init(p, v, d)
        t = 1
        # Set initial radiative forcing to zero
        v.F[t] = 0.0
    end

    function run_timestep(p, v, d, t)
        if t>1
            # Calculate total radiative forcing
            v.F[t] = (p.F2x / log(2)) * log(p.C[t] / p.C0) + p.Fext[t]
        end
    end
end
