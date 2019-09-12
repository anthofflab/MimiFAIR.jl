# --------------------------------------------------
# Direct radiative forcing effect from aerosols.
# --------------------------------------------------

@defcomp aerosol_direct_rf begin

    β_SOx                = Parameter()             # Sulfur oxides radiative efficiency: Wm⁻²(Mt yr⁻¹)⁻¹.
    β_CO                 = Parameter()             # Carbon monoxide radiative efficiency: Wm⁻²(Mt yr⁻¹)⁻¹.
    β_NMVOC              = Parameter()             # Non-methane volatile organic compounds radiative efficiency: Wm⁻²(Mt yr⁻¹)⁻¹.
    β_NOx                = Parameter()             # Nitrogon oxides radiative efficiency: Wm⁻²(Mt yr⁻¹)⁻¹.
    β_BC                 = Parameter()             # Black carbon radiative efficiency: Wm⁻²(Mt yr⁻¹)⁻¹.
    β_OC                 = Parameter()             # Organic carbon radiative efficiency: Wm⁻²(Mt yr⁻¹)⁻¹.
    β_NH3                = Parameter()             # Ammonia radiative efficiency: Wm⁻²(Mt yr⁻¹)⁻¹.
    rf_scale_aero_direct = Parameter()             # Scaling factor to capture effective radiative forcing uncertainty.
    SOx_emiss            = Parameter(index=[time]) # Emissions (MtS yr⁻¹).
    CO_emiss             = Parameter(index=[time]) # Emissions (MtCO yr⁻¹).
    NMVOC_emiss          = Parameter(index=[time]) # Emissions (Mt yr⁻¹).
    NOx_emiss            = Parameter(index=[time]) # Emissions (MtN yr⁻¹).
    BC_emiss             = Parameter(index=[time]) # Emissions (Mt yr⁻¹).
    OC_emiss             = Parameter(index=[time]) # Emissions (Mt yr⁻¹).
    NH3_emiss            = Parameter(index=[time]) # Emissions (MtN yr⁻¹).

    F_SOx                = Variable(index=[time])  # Sulfur oxides forcing contribution.
    F_CO                 = Variable(index=[time])  # Carbon monoxide forcing contribution.
    F_NMVOC              = Variable(index=[time])  # Non-methane volatile organic compounds forcing contribution.
    F_NOx                = Variable(index=[time])  # Nitrogon oxides forcing contribution.
    F_BC                 = Variable(index=[time])  # Sulfur black carbon forcing contribution.
    F_OC                 = Variable(index=[time])  # Organic carbon forcing contribution.
    F_NH3                = Variable(index=[time])  # Ammonia forcing contribution.
    F_aerosol_direct     = Variable(index=[time])  # Direct radiative forcing from aerosols (Wm⁻²).


    function run_timestep(p, v, d, t)

        # Calculate forcing contributions from individual emissions.
        v.F_SOx[t]   = p.β_SOx   * p.SOx_emiss[t]
        v.F_CO[t]    = p.β_CO    * p.CO_emiss[t]
        v.F_NMVOC[t] = p.β_NMVOC * p.NMVOC_emiss[t]
        v.F_NOx[t]   = p.β_NOx   * p.NOx_emiss[t]
        v.F_BC[t]    = p.β_BC    * p.BC_emiss[t]
        v.F_OC[t]    = p.β_OC    * p.OC_emiss[t]
        v.F_NH3[t]   = p.β_NH3   * p.NH3_emiss[t]

        # Total direct aerosol forcing based on linear relationships between emissions and forcing in Aerocom models.
        # Reference: Myhre et al., 2013: https://www.atmos-chem-phys.net/13/1853/2013
        v.F_aerosol_direct[t] = (v.F_SOx[t] + v.F_CO[t] +  v.F_NMVOC[t] + v.F_NOx[t] + v.F_BC[t] + v.F_OC[t] +  v.F_NH3[t]) * (1.0 + p.rf_scale_aero_direct)
    end
end
