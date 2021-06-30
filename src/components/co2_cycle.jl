# --------------------------------------------------
# Carbon cycle
# --------------------------------------------------

@defcomp co2_cycle begin

    co2_0           = Parameter()   # Carbon dioxide concentration in initial model period (ppm).
    co2_pi          = Parameter()   # Pre-industrial carbon dioxide concentrations (ppm).
    g0_co2          = Parameter()   # Constant to set approximation of value for α equal to Millar et al. (2017) numerical solution for iIRF100 carbon cycle parameterization at α=1.
    g1_co2          = Parameter()   # Constant to set approximation of gradient for α equal to Millar et al. (2017) numerical solution for iIRF100 carbon cycle parameterization at α=1.
    emiss2conc_co2  = Parameter()   # Conversion factor between emissions (GtC) and concentrations (ppm).
    GU_co2_0        = Parameter()   # Initial model period value for cumulative uptake of agent (unit of E⁻¹).
    r0_co2          = Parameter()   # Strength of pre-industrial uptake from atmosphere.
    rA_co2          = Parameter()   # Sensitivity of uptake from atmosphere to current atmospheric burden of agent (unit of E⁻¹).
    rT_co2          = Parameter()   # Sensitivity of uptake from atmosphere to model temperature change since initialization (K⁻¹).
    rU_co2          = Parameter()   # Sensitivity of uptake from atmosphere to cumulative uptake of agent since model initialization (unit of E⁻¹).
    τ_co2           = Parameter(index=[4])  # Atmospheric lifetime of gas in iᵗʰ reservior (years).
    a_co2           = Parameter(index=[4])  # Fraction of emissions entering iᵗʰ carbon reservior.
    R0_co2          = Parameter(index=[4])  # Initial model period value for quantity of agent in iᵗʰ atmospheric reservior (unit of E).
    E_co2           = Parameter(index=[time])   # Annual carbon dioxide emissions (GtC yr⁻¹).
    Tj              = Parameter(index=[time,3]) # Temperature change for three thermal pools (K).

    α_co2           = Variable(index=[time])    # State-dependent multiplicative adjustment coefficient of reservior lifetimes.
    co2             = Variable(index=[time])    # Total atmospheric carbon dioxide concentrations (ppm).
    GA_co2          = Variable(index=[time])    # Atmospheric burden of agent above pre-industrial levels (unit of E).
    GU_co2          = Variable(index=[time])    # Cumulative uptake of agent since model initialization (unit of E⁻¹).
    iIRFT100_co2    = Variable(index=[time])    # 100-year integrated impulse response function (the average airborne fraction over a 100-year period).
    decay_rate_co2  = Variable(index=[time,4])  # Decay rates
    R_co2           = Variable(index=[time,4])  # Quantity of agent in iᵗʰ atmospheric reservior (unit of E).


    function run_timestep(p, v, d, t)

        if is_first(t)

            # Set initial values.
            v.GU_co2[t]  = p.GU_co2_0
            v.R_co2[t,:] = p.R0_co2
            v.co2[t]     = p.co2_0

            # Calculate initial burden above pre-industrial values.
            v.GA_co2[t] = sum(v.R_co2[t,:])

        else

            # Calculate iIRF100.
            v.iIRFT100_co2[t] = abs(p.r0_co2 + p.rU_co2 * (v.GU_co2[t-1] - v.GA_co2[t-1]) + p.rT_co2 * sum(p.Tj[t-1,:]) + p.rA_co2 * v.GA_co2[t-1])

            # Calculate state-dependent lifetime adjustment term based on iIRF100 value.
            v.α_co2[t] = p.g0_co2 * exp(v.iIRFT100_co2[t] / p.g1_co2)

            # Calculate concentrations in each carbon reservoir.
            for i = 1:4

                # Calculate decay rate for iᵗʰ reservoir.
                v.decay_rate_co2[t,i] = 1.0 / (v.α_co2[t] * p.τ_co2[i])

                # Calculate amount of carbon in iᵗʰ reservoir.
                v.R_co2[t,i] = p.E_co2[t] * p.a_co2[i] / v.decay_rate_co2[t,i] * (1.0 - exp(-v.decay_rate_co2[t,i])) + v.R_co2[t-1,i] * exp(-v.decay_rate_co2[t,i])
            end

            # Calcualte atmospheric burden above pre-industiral levels.
            v.GA_co2[t] = sum(v.R_co2[t,:])

            # Calculate cumulative emissions burden.
            v.GU_co2[t] = v.GU_co2[t-1] + p.E_co2[t]

            # Calculate atmospheric CO₂ concentration.
            v.co2[t] = p.co2_pi + p.emiss2conc_co2 * (v.GA_co2[t-1] + v.GA_co2[t]) / 2.0
        end
    end
end
