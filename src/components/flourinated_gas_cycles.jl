# --------------------------------------------------
# Flourinated Gas Cycles
# --------------------------------------------------

@defcomp flourinated_cycles begin

    flourinated_gases  = Index()                                            # Index for flourinated gases controlled under the Kyoto Protocol.

    flourinated_0    		= Parameter(index=[flourinated_gases]) # flourinated gas concentration in initial model period (ppb).
    flourinated_pi   		= Parameter(index=[flourinated_gases])              # Pre-industrial flourinated gas concentrations (ppb).
    g0_flourinated   		= Parameter(index=[flourinated_gases]) # Constant to set approximation of value for α equal to Millar et al. (2017) numerical solution for iIRF100 flourinated gas cycle parameterization at α=1.
    g1_flourinated   		= Parameter(index=[flourinated_gases]) # Constant to set approximation of gradient for α equal to Millar et al. (2017) numerical solution for iIRF100 flourinated gas cycle parameterization at α=1.
    emiss2conc_flourinated  = Parameter(index=[flourinated_gases])              # Conversion factor between emissions and concentrations.
    GU_flourinated_0        = Parameter(index=[flourinated_gases]) # Initial model period value for cumulative uptake of agent (unit of E⁻¹).
    r0_flourinated   		= Parameter(index=[flourinated_gases])              # Strength of pre-industrial uptake from atmosphere.
    rA_flourinated   		= Parameter(index=[flourinated_gases]) # Sensitivity of uptake from atmosphere to current atmospheric burden of agent (unit of E⁻¹).
    rT_flourinated   		= Parameter(index=[flourinated_gases]) # Sensitivity of uptake from atmosphere to model temperature change since initialization (K⁻¹).
    rU_flourinated   		= Parameter(index=[flourinated_gases]) # Sensitivity of uptake from atmosphere to cumulative uptake of agent since model initialization (unit of E⁻¹).
    τ_flourinated    		= Parameter(index=[flourinated_gases, 4])     # Atmospheric lifetime of gas in iᵗʰ reservior (years).
    a_flourinated    		= Parameter(index=[flourinated_gases, 4])     # Fraction of emissions entering iᵗʰ flourinated gas reservior.
    R0_flourinated    		= Parameter(index=[flourinated_gases, 4]) # Initial model period value for quantity of agent in iᵗʰ atmospheric reservior (unit of E).
    E_flourinated    		= Parameter(index=[time, flourinated_gases])  # Annual flourinated gas emissions (Tg yr⁻¹).
    Tj       		    = Parameter(index=[time, 3])  # Temperature change for three thermal pools (K).

    α_flourinated           = Variable(index=[time, flourinated_gases])   # State-dependent multiplicative adjustment coefficient of reservior lifetimes.
    flourinated_conc        = Variable(index=[time, flourinated_gases])   # Total atmospheric flourinated gas concentrations (ppb).
    GA_flourinated          = Variable(index=[time, flourinated_gases])   # Atmospheric burden of agent above pre-industrial levels (unit of E).
    GU_flourinated          = Variable(index=[time, flourinated_gases])   # Cumulative uptake of agent since model initialization (unit of E⁻¹).
    iIRFT100_flourinated    = Variable(index=[time, flourinated_gases])   # 100-year integrated impulse response function (the average airborne fraction over a 100-year period).
    decay_rate_flourinated  = Variable(index=[time, flourinated_gases, 4])
    R_flourinated           = Variable(index=[time, flourinated_gases, 4]) # Quantity of agent in iᵗʰ atmospheric reservior (unit of E).


    function run_timestep(p, v, d, t)

    	for g in d.flourinated_gases

        	if is_first(t)

            	# Set initial values.
            	v.GU_flourinated[t,g]   = p.GU_flourinated_0[g]
            	v.R_flourinated[t,g,:]  = p.R0_flourinated[g,:]
            	v.flourinated_conc[t,g] = p.flourinated_0[g]

            	# Calculate initial burden above pre-industrial values.
            	v.GA_flourinated[t,g] = sum(v.R_flourinated[t,g,:])

        	else

            	# Calculate iIRF100.
            	v.iIRFT100_flourinated[t,g] = abs(p.r0_flourinated[g] + p.rU_flourinated[g] * (v.GU_flourinated[t-1,g] - v.GA_flourinated[t-1,g]) + p.rT_flourinated[g] * sum(p.Tj[t-1,:]) + p.rA_flourinated[g] * v.GA_flourinated[t-1,g])

            	# Calculate state-dependent lifetime adjustment term based on iIRF100 value.
            	v.α_flourinated[t,g] = p.g0_flourinated[g] * exp(v.iIRFT100_flourinated[t,g] / p.g1_flourinated[g])

            	# Calculate concentrations in each flourinated gas reservoir.
            	for i = 1:4

                	# Calculate decay rate for iᵗʰ reservoir.
                	v.decay_rate_flourinated[t,g,i] = 1.0 / (v.α_flourinated[t,g] * p.τ_flourinated[g,i])

                	# Calculate amount of flourinated gas in iᵗʰ reservoir.
                	v.R_flourinated[t,g,i] = p.E_flourinated[t,g] * p.a_flourinated[g,i] / v.decay_rate_flourinated[t,g,i] * (1.0 - exp(-v.decay_rate_flourinated[t,g,i])) + v.R_flourinated[t-1,g,i] * exp(-v.decay_rate_flourinated[t,g,i])
	            end

    	        # Calcualte atmospheric burden above pre-industiral levels.
            	v.GA_flourinated[t,g] = sum(v.R_flourinated[t,g,:])

            	# Calculate cumulative emissions burden.
            	v.GU_flourinated[t,g] = v.GU_flourinated[t-1,g] + p.E_flourinated[t,g]

            	# Calculate atmospheric CO₂ concentration.
            	v.flourinated_conc[t,g] = p.flourinated_pi[g] + p.emiss2conc_flourinated[g] * (v.GA_flourinated[t-1,g] + v.GA_flourinated[t,g]) / 2.0
        	end
    	end
	end
end
