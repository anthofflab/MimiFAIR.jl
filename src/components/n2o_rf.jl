# --------------------------------------------------------------
# Radiative forcing from nitrous oxide (Etminan 2016 equations)
# --------------------------------------------------------------

@defcomp n2o_rf begin
    a₂          = Parameter()             # Nitrous oxide forcing constant.
    b₂          = Parameter()             # Nitrous oxide forcing constant.
    c₂          = Parameter()             # Nitrous oxide forcing constant.
    N₂O_0       = Parameter()             # Initial (pre-industrial) nitrous oxide concentration (ppb).
    CO₂_0       = Parameter()             # Initial (pre-industrial) carbon dioxide concentration (ppm).
    CH₄_0       = Parameter()             # Initial (pre-industrial) methane concentration (ppb).
    N₂O         = Parameter(index=[time]) # Atmospheric nitrous oxide concentration (ppb).
    CO₂         = Parameter(index=[time]) # Atmospheric carbon dioxide concentration (ppb).
    CH₄         = Parameter(index=[time]) # Atmospheric methane concentration (ppb).

    forcing_N₂O = Variable(index=[time])  # Direct forcing from nitrous oxide concentrations (Wm⁻²).


    function run_timestep(p, v, d, t)

        # Create temporary averaging variables.
        M_bar = 0.5 * (p.CH₄[t] + p.CH₄_0)
        N_bar = 0.5 * (p.N₂O[t] + p.N₂O_0)
        C_bar = 0.5 * (p.CO₂[t] + p.CO₂_0)

        # Direct nitrous oxide radiative forcing.
        v.forcing_N₂O[t] = (p.a₂ * C_bar + p.b₂ * N_bar + p.c₂ * M_bar + 0.117) * (sqrt(p.N₂O[t]) - sqrt(p.N₂O_0))
    end
end
