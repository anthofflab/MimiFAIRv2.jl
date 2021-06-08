# --------------------------------------------------
# Nitrous Oxide Cycle
# --------------------------------------------------

@defcomp n2o_cycle begin

    n2o_0    		= Parameter() # Nitrous oxide concentration in initial model period (ppb).
    n2o_pi   		= Parameter()              # Pre-industrial nitrous oxide concentrations (ppb).
    g0_n2o   		= Parameter() # Constant to set approximation of value for α equal to Millar et al. (2017) numerical solution for iIRF100 nitrous oxide cycle parameterization at α=1.
    g1_n2o   		= Parameter() # Constant to set approximation of gradient for α equal to Millar et al. (2017) numerical solution for iIRF100 nitrous oxide cycle parameterization at α=1.
    emiss2conc_n2o  = Parameter()              # Conversion factor between emissions and concentrations.
    GU_n2o_0        = Parameter() # Initial model period value for cumulative uptake of agent (unit of E⁻¹).
    r0_n2o   		= Parameter()              # Strength of pre-industrial uptake from atmosphere.
    rA_n2o   		= Parameter() # Sensitivity of uptake from atmosphere to current atmospheric burden of agent (unit of E⁻¹).
    rT_n2o   		= Parameter() # Sensitivity of uptake from atmosphere to model temperature change since initialization (K⁻¹).
    rU_n2o   		= Parameter() # Sensitivity of uptake from atmosphere to cumulative uptake of agent since model initialization (unit of E⁻¹).
    τ_n2o    		= Parameter(index=[4])     # Atmospheric lifetime of gas in iᵗʰ reservior (years).
    a_n2o    		= Parameter(index=[4])     # Fraction of emissions entering iᵗʰ nitrous oxide reservior.
    R0_n2o    		= Parameter(index=[4]) # Initial model period value for quantity of agent in iᵗʰ atmospheric reservior (unit of E).
    E_n2o    		= Parameter(index=[time])  # Annual nitrous oxide emissions (TgN yr⁻¹).
    Tj       		= Parameter(index=[time,3])  # Temperature change for three thermal pools (K).

    α_n2o           = Variable(index=[time])   # State-dependent multiplicative adjustment coefficient of reservior lifetimes.
    n2o             = Variable(index=[time])   # Total atmospheric nitrous oxide concentrations (ppb).
    GA_n2o          = Variable(index=[time])   # Atmospheric burden of agent above pre-industrial levels (unit of E).
    GU_n2o          = Variable(index=[time])   # Cumulative uptake of agent since model initialization (unit of E⁻¹).
    iIRFT100_n2o    = Variable(index=[time])   # 100-year integrated impulse response function (the average airborne fraction over a 100-year period).
    decay_rate_n2o  = Variable(index=[time,4])
    R_n2o           = Variable(index=[time,4]) # Quantity of agent in iᵗʰ atmospheric reservior (unit of E).


    function run_timestep(p, v, d, t)

        if is_first(t)

            # Set initial values.
            v.GU_n2o[t]  = p.GU_n2o_0
            v.R_n2o[t,:] = p.R0_n2o
            v.n2o[t]     = p.n2o_0

            # Calculate initial burden above pre-industrial values.
            v.GA_n2o[t] = sum(v.R_n2o[t,:])

        else

            # Calculate iIRF100.
            v.iIRFT100_n2o[t] = abs(p.r0_n2o + p.rU_n2o * (v.GU_n2o[t-1] - v.GA_n2o[t-1]) + p.rT_n2o * sum(p.Tj[t-1,:]) + p.rA_n2o * v.GA_n2o[t-1])

            # Calculate state-dependent lifetime adjustment term based on iIRF100 value.
            v.α_n2o[t] = p.g0_n2o * exp(v.iIRFT100_n2o[t] / p.g1_n2o)

            # Calculate concentrations in each nitrous oxide reservoir.
            for i = 1:4

                # Calculate decay rate for iᵗʰ reservoir.
                v.decay_rate_n2o[t,i] = 1.0 / (v.α_n2o[t] * p.τ_n2o[i])

                # Calculate amount of N₂O in iᵗʰ reservoir.
                v.R_n2o[t,i] = p.E_n2o[t] * p.a_n2o[i] / v.decay_rate_n2o[t,i] * (1.0 - exp(-v.decay_rate_n2o[t,i])) + v.R_n2o[t-1,i] * exp(-v.decay_rate_n2o[t,i])
            end

            # Calcualte atmospheric burden above pre-industiral levels.
            v.GA_n2o[t] = sum(v.R_n2o[t,:])

            # Calculate cumulative emissions burden.
            v.GU_n2o[t] = v.GU_n2o[t-1] + p.E_n2o[t]

            # Calculate atmospheric N₂O concentration.
            v.n2o[t] = p.n2o_pi + p.emiss2conc_n2o * (v.GA_n2o[t-1] + v.GA_n2o[t]) / 2.0
        end
    end
end
