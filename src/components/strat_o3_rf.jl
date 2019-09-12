# --------------------------------------------
# Radiative forcing from stratospheric ozone
# --------------------------------------------

@defcomp strat_o3_rf begin

    ozone_depleting_substances = Index() # Index for ozone-depleting substances.

    δ1                         = Parameter()                                         # Fitting parameter to the stratospheric ozone ERF time series from AR5.
    δ2                         = Parameter()                                         # Fitting parameter to the stratospheric ozone ERF time series from AR5.
    δ3                         = Parameter()                                         # Fitting parameter to the stratospheric ozone ERF time series from AR5.
    Br                         = Parameter(index=[ozone_depleting_substances])       # Number of bromine atoms in each ozone-depleting species.
    Cl                         = Parameter(index=[ozone_depleting_substances])       # Number of chlorine atoms in each ozone-depleting species.
    FC                         = Parameter(index=[ozone_depleting_substances])       # Fractional stratospheric release for ozone-depleting substances (Reference: Daniel, J. and Velders, G.: A focus on information and options for policymakers, in: Scientific Assessment of Ozone Depletion, WMO, 2011).
    ODS₀                       = Parameter(index=[ozone_depleting_substances])       # Initial concentrations for ozone_depleting_substances (ppt).
    conc_ODS                   = Parameter(index=[time, ozone_depleting_substances]) # Atmospheric concentrations of ozone_depleting_substances (ppt).

    int_eecs                   = Variable(index=[time, ozone_depleting_substances])  # Intermediate terms to caluclate equivalent effective stratospheric chlorine (EESC) for ozone-depleting substances (just for convenience).
    eecs                       = Variable(index=[time])                              # Equivalent effective stratospheric chlorine (EESC) for ozone-depleting substances.
    forcing_strat_O₃           = Variable(index=[time])                              # Radiative forcing from stratospheric ozone (Wm⁻²).


    function run_timestep(p, v, d, t)

        for g in d.ozone_depleting_substances
            # Calculate intermediate variables for EECS (calculations done for each individual ozone-depleting substance).
            # Note from FAIR Code: "EESC takes ODS concentrations in ppb, we provide ppt."
            v.int_eecs[t,g] = p.Cl[g] * 1000.0 * (p.conc_ODS[t,g]-p.ODS₀[g]) * (p.FC[g]/p.FC[1]) + 45 * p.Br[g] * 1000 * (p.conc_ODS[t,g]-p.ODS₀[g]) * (p.FC[g]/p.FC[1])
        end

        # Calcualte equivalent effective stratospheric chlorine.
        v.eecs[t] = max((p.FC[1]*sum(v.int_eecs[t,:])), 0.0)

        # Calcualte stratospheric O₃ forcing.
        v.forcing_strat_O₃[t] = p.δ1 * (p.δ2 * v.eecs[t]) ^ p.δ3
    end
end
