# --------------------------------------------------
# Methane Cycle
# --------------------------------------------------

@defcomp ch4_cycle begin

    ch4_0    		= Parameter()   # Methane concentration in initial model period (ppb).
    ch4_pi   		= Parameter()   # Pre-industrial methane concentrations (ppb).
    g0_ch4   		= Parameter()   # Constant to set approximation of value for α equal to Millar et al. (2017) numerical solution for iIRF100 methane cycle parameterization at α=1.
    g1_ch4   		= Parameter()   # Constant to set approximation of gradient for α equal to Millar et al. (2017) numerical solution for iIRF100 methane cycle parameterization at α=1.
    emiss2conc_ch4  = Parameter()   # Conversion factor between emissions and concentrations.
    GU_ch4_0        = Parameter()   # Initial model period value for cumulative uptake of agent (unit of E⁻¹).
    r0_ch4   		= Parameter()   # Strength of pre-industrial uptake from atmosphere.
    rA_ch4   		= Parameter()   # Sensitivity of uptake from atmosphere to current atmospheric burden of agent (unit of E⁻¹).
    rT_ch4   		= Parameter()   # Sensitivity of uptake from atmosphere to model temperature change since initialization (K⁻¹).
    rU_ch4   		= Parameter()   # Sensitivity of uptake from atmosphere to cumulative uptake of agent since model initialization (unit of E⁻¹).
    τ_ch4    		= Parameter(index=[4])  # Atmospheric lifetime of gas in iᵗʰ reservior (years).
    a_ch4    		= Parameter(index=[4])  # Fraction of emissions entering iᵗʰ methane reservior.
    R0_ch4    		= Parameter(index=[4])  # Initial model period value for quantity of agent in iᵗʰ atmospheric reservior (unit of E).
    E_ch4    		= Parameter(index=[time])   # Annual methane emissions (TgCH₄ yr⁻¹).
    Tj       		= Parameter(index=[time,3]) # Temperature change for three thermal pools (K).

    α_ch4           = Variable(index=[time])    # State-dependent multiplicative adjustment coefficient of reservior lifetimes.
    ch4             = Variable(index=[time])    # Total atmospheric methane concentrations (ppb).
    GA_ch4          = Variable(index=[time])    # Atmospheric burden of agent above pre-industrial levels (unit of E).
    GU_ch4          = Variable(index=[time])    # Cumulative uptake of agent since model initialization (unit of E⁻¹).
    iIRFT100_ch4    = Variable(index=[time])    # 100-year integrated impulse response function (the average airborne fraction over a 100-year period).
    decay_rate_ch4  = Variable(index=[time,4])  # Decay rates
    R_ch4           = Variable(index=[time,4])  # Quantity of agent in iᵗʰ atmospheric reservior (unit of E).


    function run_timestep(p, v, d, t)

        if is_first(t)

            # Set initial values.
            v.GU_ch4[t]  = p.GU_ch4_0
            v.R_ch4[t,:] = p.R0_ch4
            v.ch4[t]     = p.ch4_0

            # Calculate initial burden above pre-industrial values.
            v.GA_ch4[t] = sum(v.R_ch4[t,:])

        else

            # Calculate iIRF100.
            v.iIRFT100_ch4[t] = abs(p.r0_ch4 + p.rU_ch4 * (v.GU_ch4[t-1] - v.GA_ch4[t-1]) + p.rT_ch4 * sum(p.Tj[t-1,:]) + p.rA_ch4 * v.GA_ch4[t-1])

            # Calculate state-dependent lifetime adjustment term based on iIRF100 value.
            v.α_ch4[t] = p.g0_ch4 * exp(v.iIRFT100_ch4[t] / p.g1_ch4)

            # Calculate concentrations in each methane reservoir.
            for i = 1:4

                # Calculate decay rate for iᵗʰ reservoir.
                v.decay_rate_ch4[t,i] = 1.0 / (v.α_ch4[t] * p.τ_ch4[i])

                # Calculate amount of methane in iᵗʰ reservoir.
                v.R_ch4[t,i] = p.E_ch4[t] * p.a_ch4[i] / v.decay_rate_ch4[t,i] * (1.0 - exp(-v.decay_rate_ch4[t,i])) + v.R_ch4[t-1,i] * exp(-v.decay_rate_ch4[t,i])
            end

            # Calcualte atmospheric burden above pre-industiral levels.
            v.GA_ch4[t] = sum(v.R_ch4[t,:])

            # Calculate cumulative emissions burden.
            v.GU_ch4[t] = v.GU_ch4[t-1] + p.E_ch4[t]

            # Calculate atmospheric CH₄ concentration.
            v.ch4[t] = p.ch4_pi + p.emiss2conc_ch4 * (v.GA_ch4[t-1] + v.GA_ch4[t]) / 2.0
        end
    end
end
