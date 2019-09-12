# -------------------------------------------------------
# Radiative forcing from contrails
# -------------------------------------------------------

@defcomp contrails_rf begin

    E_ref             = Parameter()             # Reference-year emissions of aviation nitrogen oxides  (Mt yr⁻¹).
    F_ref             = Parameter()             # Forcing from linear persistent contrails + contrail induced cirrus. The default for 2005 is from Lee et al, 2009 (Atmos. Environ., doi:10.1016/j.atmosenv.2009.04.024).
    ref_is_NO2::Bool  = Parameter()             # True if 'E_ref' is in units of NO₂ rather than N.
    mol_weight_NO₂    = Parameter()             # Molecular mass of nitrogen dioxide.
    mol_weight_N      = Parameter()             # Molecular mass of nitrogen.
    NOx_emiss         = Parameter(index=[time]) # Nitrogen oxides emissions.
    frac              = Parameter(index=[time]) # Fraction of total nitrogen oxides emissions due to aviation.

    forcing_contrails = Variable(index=[time])  # Forcing approximation from aviation contrails (Wm⁻²).


    function run_timestep(p, v, d, t)

        # Caluclate forcing from contrail.
        # Note: If 'E_ref' in units of NO₂ rather than N, use ratio of molecular weights for conversion.
        if p.ref_is_NO2
            v.forcing_contrails[t] = p.NOx_emiss[t] * p.frac[t] * (p.F_ref/p.E_ref) * (p.mol_weight_NO₂ / p.mol_weight_N)
        else
            v.forcing_contrails[t] = p.NOx_emiss[t] * p.frac[t] * (p.F_ref/p.E_ref)
        end
    end
end
