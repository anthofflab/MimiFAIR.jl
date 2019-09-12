#--------------------------
# Total raditiave forcing
#--------------------------

@defcomp total_rf begin

    other_ghg = Index()                                            # Index for other well-mixed greenhouse gases.

                                                                   # NOTE: For a given gas, efficacy = global temperature response per unit forcing relative to the response to CO₂ forcing.
    efficacy_CO₂              = Parameter()                        # Carbon dioxide efficacy.
    efficacy_CH₄              = Parameter()                        # Methane efficacy.
    efficacy_CH₄_H₂O          = Parameter()                        # Stratospheric water vapor from methane oxidation efficay.
    efficacy_N₂O              = Parameter()                        # Nitrous oxide efficay.
    efficacy_trop_O₃          = Parameter()                        # Tropospheric ozone efficay.
    efficacy_strat_O₃         = Parameter()                        # Stratospheric ozone efficay.
    efficacy_aerosol_direct   = Parameter()                        # Direct aerosol effect efficay.
    efficacy_aerosol_indirect = Parameter()                        # Indirect aerosol effect efficay.
    efficacy_bcsnow           = Parameter()                        # Black carbon on snow efficay.
    efficacy_landuse          = Parameter()                        # Land-use surface albedo change efficay.
    efficacy_contrails        = Parameter()                        # Contrails efficay.
    efficacy_other_ghg        = Parameter(index=[other_ghg])       # Radiative forcing from other well-mixed greenhouse gases (Wm⁻²)
    F_CO₂                     = Parameter(index=[time])            # Radiative forcing from carbon dioxide (Wm⁻²).
    F_CH₄                     = Parameter(index=[time])            # Radiative forcing from methane (Wm⁻²).
    F_CH₄_H₂O                 = Parameter(index=[time])            # Radiative forcing from stratospheric water varpor from methane oxidiation (Wm⁻²).
    F_N₂O                     = Parameter(index=[time])            # Radiative forcing from nitrous oxide (Wm⁻²).
    F_trop_O₃                 = Parameter(index=[time])            # Radiative forcing from tropospheric ozone (Wm⁻²).
    F_strat_O₃                = Parameter(index=[time])            # Radiative forcing from stratospheric ozone (Wm⁻²).
    F_aerosol_direct          = Parameter(index=[time])            # Radiative forcing from direct aerosol effect (Wm⁻²).
    F_aerosol_indirect        = Parameter(index=[time])            # Radiative forcing from indirect aerosol effect (Wm⁻²).
    F_bcsnow                  = Parameter(index=[time])            # Radiative forcing from black carbon on snow (Wm⁻²).
    F_landuse                 = Parameter(index=[time])            # Radiative forcing from land-use surface albedo change (Wm⁻²).
    F_contrails               = Parameter(index=[time])            # Radiative forcing from contrails (Wm⁻²).
    F_volcanic                = Parameter(index=[time])            # Radiative forcing from volcanic eruptions (Wm⁻²).
    F_solar                   = Parameter(index=[time])            # Radiative forcing from solar irradience (Wm⁻²).
    F_exogenous               = Parameter(index=[time])            # Radiative forcing from exogenous sources not included in this module (Wm⁻²).
    F_other_ghg               = Parameter(index=[time, other_ghg]) # Radiative forcing from other well-mixed greenhouse gases (Wm⁻²).

    total_forcing             = Variable(index = [time])           # Total radiative forcing, with individual components scaled by their respective efficacy (Wm⁻²).


    function run_timestep(p, v, d, t)

        # Calculate total radiative forcing as the sum of individual radiative forcings scaled by efficacy factor.
        v.total_forcing[t] = p.F_CO₂[t]              * p.efficacy_CO₂ +
                             p.F_CH₄[t]              * p.efficacy_CH₄ +
                             p.F_CH₄_H₂O[t]          * p.efficacy_CH₄_H₂O +
                             p.F_N₂O[t]              * p.efficacy_N₂O +
                             sum(p.F_other_ghg[t,:] .* p.efficacy_other_ghg[:]) +
                             p.F_trop_O₃[t]          * p.efficacy_trop_O₃ +
                             p.F_strat_O₃[t]         * p.efficacy_strat_O₃ +
                             p.F_aerosol_direct[t]   * p.efficacy_aerosol_direct +
                             p.F_aerosol_indirect[t] * p.efficacy_aerosol_indirect +
                             p.F_bcsnow[t]           * p.efficacy_bcsnow +
                             p.F_landuse[t]          * p.efficacy_landuse +
                             p.F_contrails[t]        * p.efficacy_contrails +
                             p.F_solar[t] +
                             p.F_volcanic[t]
    end
end
