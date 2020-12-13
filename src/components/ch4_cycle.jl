# --------------------------------------------------
# Methane cycle.
# --------------------------------------------------

@defcomp ch4_cycle begin

    emiss2conc_ch4    = Parameter()             # Conversion between ppb/ppt concentrations and Mt/kt emissions.
    CH₄_0             = Parameter()             # Initial (pre-industrial) atmospheric methane concentration (ppb).
    τ_CH₄             = Parameter()             # Atmospheric (e-folding) lifetime of methane.
    gtc2ppm           = Parameter()             # Conversion factor between GtC and ppm.
    oxidation_frac    = Parameter()             # Fraction of methane lost through reaction with hydroxyl radical that is converted to carbon dioxide.
    mol_weight_CH₄    = Parameter()             # Molecular mass of methane.
    mol_weight_C      = Parameter()             # Molecular mass of carbon.
    fossil_emiss_CH₄  = Parameter(index=[time]) # Fossil-fuel methane emissions (Mt CH₄ yr⁻¹).
    natural_emiss_CH₄ = Parameter(index=[time]) # Natural methane emissions (Mt CH₄ yr⁻¹).
    fossil_frac       = Parameter(index=[time]) # Fraciton of anthropogenic methane attributable to fossil sources.

    CH₄              = Variable(index=[time])  # Atmospheric methane concentration (ppb).
    oxidised_CH₄     = Variable(index=[time])  # Methane that has been oxidized to carbon dioxide (ppb).
    oxidised_CH₄_GtC = Variable(index=[time])  # Methane that has been oxidized to carbon dioxide (GtC)

    function run_timestep(p, v, d, t)

        # Set initial methane concentration values.
        if is_first(t)
            v.CH₄[t] = p.CH₄_0
            v.oxidised_CH₄[t] = 0.0
            v.oxidised_CH₄_GtC[t] = 0.0
        else
            # Calculate atmospheric methane concentration.
            emiss_prev = p.fossil_emiss_CH₄[t-1] + p.natural_emiss_CH₄[t]
            emiss_curr = p.fossil_emiss_CH₄[t] + p.natural_emiss_CH₄[t]
            v.CH₄[t] = v.CH₄[t-1] - v.CH₄[t-1] * (1.0 - exp(-1/p.τ_CH₄)) + 0.5 * (emiss_prev + emiss_curr) * (1.0/p.emiss2conc_ch4)

            # Calculate carbon dioxide from oxidized methane (bounded below at 0.0).
            v.oxidised_CH₄[t] = max(((v.CH₄[t-1]-p.CH₄_0) * (1.0 - exp(-1.0/p.τ_CH₄)) * (p.mol_weight_C/p.mol_weight_CH₄ * 0.001 * p.oxidation_frac * p.fossil_frac[t])), 0.0)

            # Also provide oxidised CH₄ in units of GtC.
            v.oxidised_CH₄_GtC[t] = v.oxidised_CH₄[t] * p.gtc2ppm
        end
    end
end
