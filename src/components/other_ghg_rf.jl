# ----------------------------------------------------------------------------------------------
# Radiative forcing from of other well-mixed greenhouse gases (Kyoto and ozone-depleting gases)
# ----------------------------------------------------------------------------------------------

@defcomp other_ghg_rf begin

    other_ghg            = Index()                            # Index for other well-mixed greenhouse gases.

    other_ghg_0          = Parameter(index=[other_ghg])       # Initial (pre-industrial) concentration for other well-mixed greenhouse gases (ppt).
    radiative_efficiency = Parameter(index=[other_ghg])       # Radiative efficiencies for other well-mixed greenhouse gases - from IPCC AR5 WG1 Ch8 SM (Wm⁻²ppb⁻¹).
    conc_other_ghg       = Parameter(index=[time, other_ghg]) # Atmospheric concentrations for other well-mixed greenhouse gases (ppt).

    other_ghg_rf_net     = Variable(index=[time])             # Net radiative forcing for all well-mixed greenhouse gases (Wm⁻²).
    other_ghg_rf         = Variable(index=[time, other_ghg])  # Individual radiative forcings for each well-mixed greenhouse gas (Wm⁻²).


    function run_timestep(p, v, d, t)

        for g in d.other_ghg
            # Caluclate radiative forcing for individual gases.
            # Note: the factor of 0.001 here is because radiative efficiencies are given in Wm⁻²ppb⁻¹ and concentrations of minor gases are in ppt.
            v.other_ghg_rf[t,g] = (p.conc_other_ghg[t,g] - p.other_ghg_0[g]) * p.radiative_efficiency[g] * 0.001
        end

        # Calculate total radiative forcing as sum of individual gas forcings.
        v.other_ghg_rf_net[t] = sum(v.other_ghg_rf[t,:])
    end
end
