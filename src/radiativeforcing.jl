using Mimi

@defcomp radiativeforcing begin

    C0  = Parameter()               # Pre-industrial CO2 concentrations
    F2x = Parameter()               # Forcing due to CO2 doubling (w/m2)
    Fext= Parameter(index=[time])   # Non-CO2 forcing (w/m2)
    C   = Parameter(index=[time])   # Total atmospheric CO2 concentrations

    F   = Variable(index=[time])    # Total radiative forcing (w/m2)
end


function run_timestep(s::radiativeforcing, t::Int)
    p = s.Parameters
    v = s.Variables

    # Calculate total radiative forcing
    v.F[t] = p.F2x / log(2) * log(p.C[t] / p.C0) + p.Fext[t]
end
