# --------------------------------------------
# Radiative forcing from tropospheric ozone
# --------------------------------------------

@defcomp trop_o3_rf begin

    mol_weight_N          = Parameter()             # Nolecular mass of N.
    mol_weight_NO         = Parameter()             # Molecular mass of NO.
    CH₄_0                 = Parameter()             # Initial (pre-industrial) atmospheric methane concentration (ppb).
    T0                    = Parameter()             # Initial global mean surface temperature anomaly (K).
    fix_pre1850_RCP       = Parameter{Bool}()             # Switch to use different relationship for 1750/65 to 1850 based on anthropogenic emissions from Skeie et al. (2011) (atmos-chem-phys.net/11/11827/2011).
    CH₄                   = Parameter(index=[time]) # Atmospheric methane concentration (ppb).
    NOx_emissions         = Parameter(index=[time]) # Nitrogen oxides emissions (MtN yr⁻¹).
    CO_emissions          = Parameter(index=[time]) # Carbon monoxide emissions (MtCO yr⁻¹).
    NMVOC_emissions       = Parameter(index=[time]) # Global non-methane volatile organic compounds emissions (Mt yr⁻¹).
    temperature           = Parameter(index=[time]) # Global mean surface temperature anomaly (K).

    F_CH₄                 = Variable(index=[time])  # Radiative forcing contribution from methane (Wm⁻²).
    F_CO                  = Variable(index=[time])  # Radiative forcing contribution from carbon monoxide (Wm⁻²).
    F_NMVOC               = Variable(index=[time])  # Radiative forcing contribution from non-methane olatile organic compounds (Wm⁻²).
    F_NOx                 = Variable(index=[time])  # Radiative forcing contribution from nitrogen oxides (Wm⁻²).
    forcing_trop_O₃       = Variable(index=[time])  # Radiative forcing from tropospheric ozone (Wm⁻²).


    function run_timestep(p, v, d, t)

        # Calculate individual forcing components for tropospheric ozone.
            # Note from original FAIR code: "The RCP scenarios give a negative forcing prior to ~1780. This is
                # because the anthropogenic emissions are given to be zero in RCPs but not zero in the Skeie
                # numbers which are used here. This can be fixed to give a more linear behaviour."
        if gettime(t)  >= 1850 || p.fix_pre1850_RCP == false
            v.F_CH₄[t] = 0.166/960  * (p.CH₄[t]-722.0)
            v.F_CO[t] = 0.058/681.8 * (p.CO_emissions[t]-170.0)
            v.F_NMVOC[t] = 0.035/155.84 * (p.NMVOC_emissions[t]-10.0)
            v.F_NOx[t] =  0.119/61.16  * (p.NOx_emissions[t] * p.mol_weight_NO / p.mol_weight_N - 4.29)
        else
            v.F_CH₄[t] = 0.166/960  * (p.CH₄[t]-722.0)
            v.F_CO[t] = 0.058/681.8 * 215.59  * p.CO_emissions[t] / 385.59
            v.F_NMVOC[t] = 0.035/155.84 * 51.97 * p.NMVOC_emissions[t] / 61.97
            v.F_NOx[t] =  0.119/61.16  * 7.31 * (p.NOx_emissions[t] * p.mol_weight_NO / p.mol_weight_N) / 11.6
        end

        # Calculate total tropospheric ozone radiative forcing, accounting for a temperature feedback.
            # Note from original FAIR code: "We fit a curve to the 2000, 2030 and 2100 best estimates of feedback based on middle-of-the-road temperature projections."
        if is_first(t)
            # Initial O₃ forcing temperature sensitivity assumed to be zero (since model has not estimated temperature yet).
            v.forcing_trop_O₃[t] = v.F_CH₄[t] + v.F_CO[t] + v.F_NMVOC[t] + v.F_NOx[t] + temperature_feedback(p.T0)
        else
            v.forcing_trop_O₃[t] = v.F_CH₄[t] + v.F_CO[t] + v.F_NMVOC[t] + v.F_NOx[t] + temperature_feedback(p.temperature[t-1])
        end
    end
end
