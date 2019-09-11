# --------------------------------------------------------------------------------------
# Concentrations of other well-mixed greenhouse gases (Kyoto and ozone-depleting gases)
# --------------------------------------------------------------------------------------

@defcomp other_ghg_cycles begin

    other_ghg                  = Index()                                            # Index for other well-mixed greenhouse gases
    ozone_depleting_substances = Index()                                            # Index for ozone-depleting substances.

    τ                          = Parameter(index=[other_ghg])                       # Atmospheric (e-folding) lifetime for each gas (years).
    other_ghg_0                = Parameter(index=[other_ghg])                       # Initial (pre-industrial) concentration for each gas (ppt).
    emiss2conc_other_ghg       = Parameter(index=[other_ghg])                       # Conversion between ppt concentrations and kt emissions.
    emiss_other_ghg            = Parameter(index=[time, other_ghg])                 # Emissions for other well-mixed greenhouse gases (kt yr⁻¹).

    conc_other_ghg             = Variable(index=[time, other_ghg])                  # Atmospheric concentrations for other well-mixed greenhouse gases (ppt).
    conc_ods                   = Variable(index=[time, ozone_depleting_substances]) # Concentrations for ozone-depleting substances (a subset of 'conc_other_ghg' - used in stratospheric O₃ component).


    function run_timestep(p, v, d, t)

        for g in d.other_ghg

            # Set initial concentration values.
            if isfirst(t)
                v.conc_other_ghg[t,g] = p.other_ghg_0[g]
            else
                # Calculate concentrations based on simple one-box exponential decay model.
                v.conc_other_ghg[t,g] = v.conc_other_ghg[t-1,g] - v.conc_other_ghg[t-1,g] * (1.0 - exp(-1/p.τ[g])) + 0.5 * (p.emiss_other_ghg[t-1,g] + p.emiss_other_ghg[t,g]) * (1.0/p.emiss2conc_other_ghg[g])
            end
        end

        # For convenience reasons, create a subset variable for ozone-depleting substances (indices 13-28) used in other model components.
        # TODO: Remove this type of indexing (leftover from orifinal FAIR code).
        v.conc_ods[t,:] = v.conc_other_ghg[t,13:28]
    end
end
