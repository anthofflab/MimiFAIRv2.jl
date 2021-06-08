# --------------------------------------------------
# Montreal Protocol Gas Cycles
# --------------------------------------------------

@defcomp montreal_cycles begin

    montreal_gases  = Index()                                            # Index for Ozone Depleting Substances controlled under the Montreal Protocol.

    montreal_0    		= Parameter(index=[montreal_gases]) # Montreal gas concentration in initial model period (ppb).
    montreal_pi   		= Parameter(index=[montreal_gases])              # Pre-industrial Montreal gas concentrations (ppb).
    g0_montreal   		= Parameter(index=[montreal_gases]) # Constant to set approximation of value for α equal to Millar et al. (2017) numerical solution for iIRF100 Montreal gas cycle parameterization at α=1.
    g1_montreal   		= Parameter(index=[montreal_gases]) # Constant to set approximation of gradient for α equal to Millar et al. (2017) numerical solution for iIRF100 Montreal gas cycle parameterization at α=1.
    emiss2conc_montreal = Parameter(index=[montreal_gases])              # Conversion factor between emissions and concentrations.
    GU_montreal_0       = Parameter(index=[montreal_gases]) # Initial model period value for cumulative uptake of agent (unit of E⁻¹).
    r0_montreal   		= Parameter(index=[montreal_gases])              # Strength of pre-industrial uptake from atmosphere.
    rA_montreal   		= Parameter(index=[montreal_gases]) # Sensitivity of uptake from atmosphere to current atmospheric burden of agent (unit of E⁻¹).
    rT_montreal   		= Parameter(index=[montreal_gases]) # Sensitivity of uptake from atmosphere to model temperature change since initialization (K⁻¹).
    rU_montreal   		= Parameter(index=[montreal_gases]) # Sensitivity of uptake from atmosphere to cumulative uptake of agent since model initialization (unit of E⁻¹).
    τ_montreal    		= Parameter(index=[montreal_gases, 4])     # Atmospheric lifetime of gas in iᵗʰ reservior (years).
    a_montreal    		= Parameter(index=[montreal_gases, 4])     # Fraction of emissions entering iᵗʰ Montreal gas reservior.
    R0_montreal    		= Parameter(index=[montreal_gases, 4]) # Initial model period value for quantity of agent in iᵗʰ atmospheric reservior (unit of E).
    E_montreal    		= Parameter(index=[time, montreal_gases])  # Annual Montreal gas emissions (Tg yr⁻¹).
    Tj       		    = Parameter(index=[time, 3])  # Temperature change for three thermal pools (K).

    α_montreal           = Variable(index=[time, montreal_gases])   # State-dependent multiplicative adjustment coefficient of reservior lifetimes.
    montreal_conc        = Variable(index=[time, montreal_gases])   # Total atmospheric Montreal gas concentrations (ppb).
    GA_montreal          = Variable(index=[time, montreal_gases])   # Atmospheric burden of agent above pre-industrial levels (unit of E).
    GU_montreal          = Variable(index=[time, montreal_gases])   # Cumulative uptake of agent since model initialization (unit of E⁻¹).
    iIRFT100_montreal    = Variable(index=[time, montreal_gases])   # 100-year integrated impulse response function (the average airborne fraction over a 100-year period).
    decay_rate_montreal  = Variable(index=[time, montreal_gases, 4])
    R_montreal           = Variable(index=[time, montreal_gases, 4]) # Quantity of agent in iᵗʰ atmospheric reservior (unit of E).


    function run_timestep(p, v, d, t)

    	for g in d.montreal_gases

        	if is_first(t)

            	# Set initial values.
            	v.GU_montreal[t,g]   = p.GU_montreal_0[g]
            	v.R_montreal[t,g,:]  = p.R0_montreal[g,:]
            	v.montreal_conc[t,g] = p.montreal_0[g]

            	# Calculate initial burden above pre-industrial values.
            	v.GA_montreal[t,g] = sum(v.R_montreal[t,g,:])

        	else

            	# Calculate iIRF100.
            	v.iIRFT100_montreal[t,g] = abs(p.r0_montreal[g] + p.rU_montreal[g] * (v.GU_montreal[t-1,g] - v.GA_montreal[t-1,g]) + p.rT_montreal[g] * sum(p.Tj[t-1,:]) + p.rA_montreal[g] * v.GA_montreal[t-1,g])

            	# Calculate state-dependent lifetime adjustment term based on iIRF100 value.
            	v.α_montreal[t,g] = p.g0_montreal[g] * exp(v.iIRFT100_montreal[t,g] / p.g1_montreal[g])

            	# Calculate concentrations in each Montreal gas reservoir.
            	for i = 1:4

                	# Calculate decay rate for iᵗʰ reservoir.
                	v.decay_rate_montreal[t,g,i] = 1.0 / (v.α_montreal[t,g] * p.τ_montreal[g,i])

                	# Calculate amount of Montreal gas in iᵗʰ reservoir.
                	v.R_montreal[t,g,i] = p.E_montreal[t,g] * p.a_montreal[g,i] / v.decay_rate_montreal[t,g,i] * (1.0 - exp(-v.decay_rate_montreal[t,g,i])) + v.R_montreal[t-1,g,i] * exp(-v.decay_rate_montreal[t,g,i])
	            end

    	        # Calcualte atmospheric burden above pre-industiral levels.
            	v.GA_montreal[t,g] = sum(v.R_montreal[t,g,:])

            	# Calculate cumulative emissions burden.
            	v.GU_montreal[t,g] = v.GU_montreal[t-1,g] + p.E_montreal[t,g]

            	# Calculate atmospheric CO₂ concentration.
            	v.montreal_conc[t,g] = p.montreal_pi[g] + p.emiss2conc_montreal[g] * (v.GA_montreal[t-1,g] + v.GA_montreal[t,g]) / 2.0
        	end
    	end
	end
end
