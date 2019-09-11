# -------------------------------------------------------
# Radiative forcing from methane (Etminan 2016 equations)
# -------------------------------------------------------

@defcomp ch4_rf begin
    a₃              = Parameter()             # Methane forcing constant.
    b₃              = Parameter()             # Methane forcing constant.
    N₂O_0           = Parameter()             # Initial (pre-industrial) nitrous oxide concentration (ppb).
    scale_CH₄       = Parameter()             # Scaling factor for uncertainty in methane radiative forcing.
    H₂O_share       = Parameter()             # Share of direct methane forcing used to estimate stratospheric water vapor forcing due to methane oxidation.
    CH₄_0           = Parameter()             # Initial (pre-industrial) methane concentration (ppb).
    N₂O             = Parameter(index=[time]) # Atmospheric nitrous oxide concentration (ppb).
    CH₄             = Parameter(index=[time]) # Atmospheric methane concentration (ppb).

    forcing_CH₄_H₂O = Variable(index=[time])  # Additional indirect forcing due to production of stratospheric water vapor from methane oxidation (Wm⁻²).
    forcing_CH₄     = Variable(index=[time])  # Direct forcing from atmospheric methane concentrations (Wm⁻²).


    function run_timestep(p, v, d, t)

        # Create temporary averaging variables.
        M_bar = 0.5 * (p.CH₄[t] + p.CH₄_0)
        N_bar = 0.5 * (p.N₂O[t] + p.N₂O_0)

        # Direct methane radiative forcing.
        v.forcing_CH₄[t] = ((p.a₃ * M_bar + p.b₃ * N_bar + 0.043) * (sqrt(p.CH₄[t]) - sqrt(p.CH₄_0))) * p.scale_CH₄

        # Calculate indirect CH₄ forcing due to production of stratospheric H20 from CH4 oxidation.
        # NOTE: Removing forcing uncertainty term from above equation (water vapor forcing is based on direct methane forcing without scale_CH₄ term).
        v.forcing_CH₄_H₂O[t] = p.H₂O_share *  (v.forcing_CH₄[t] / p.scale_CH₄)
    end
end
