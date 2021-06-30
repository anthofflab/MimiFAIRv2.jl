# ---------------------------------------------------------------------------
# Tropospheric Ozone Precursors, Aerosols & Reactive Gas Cycles (Aerosol+)
# ---------------------------------------------------------------------------

@defcomp aerosol_plus_cycles begin

    aerosol_plus_gases  = Index()   # Index for tropospheric ozone precursors, aerosols and reactive gas

    aerosol_plus_0    		= Parameter(index=[aerosol_plus_gases])     # Aerosol+ gas concentration in initial model period (ppb).
    aerosol_plus_pi   		= Parameter(index=[aerosol_plus_gases])     # Pre-industrial Aerosol+ concentrations (ppb).
    g0_aerosol_plus   		= Parameter(index=[aerosol_plus_gases])     # Constant to set approximation of value for α equal to Millar et al. (2017) numerical solution for iIRF100 Aerosol+ gas cycle parameterization at α=1.
    g1_aerosol_plus   		= Parameter(index=[aerosol_plus_gases])     # Constant to set approximation of gradient for α equal to Millar et al. (2017) numerical solution for iIRF100 Aerosol+ gas cycle parameterization at α=1.
    emiss2conc_aerosol_plus = Parameter(index=[aerosol_plus_gases])     # Conversion factor between emissions and concentrations.
    GU_aerosol_plus_0       = Parameter(index=[aerosol_plus_gases])     # Initial model period value for cumulative uptake of agent (unit of E⁻¹).
    r0_aerosol_plus   		= Parameter(index=[aerosol_plus_gases])     # Strength of pre-industrial uptake from atmosphere.
    rA_aerosol_plus   		= Parameter(index=[aerosol_plus_gases])     # Sensitivity of uptake from atmosphere to current atmospheric burden of agent (unit of E⁻¹).
    rT_aerosol_plus   		= Parameter(index=[aerosol_plus_gases])     # Sensitivity of uptake from atmosphere to model temperature change since initialization (K⁻¹).
    rU_aerosol_plus   		= Parameter(index=[aerosol_plus_gases])     # Sensitivity of uptake from atmosphere to cumulative uptake of agent since model initialization (unit of E⁻¹).
    τ_aerosol_plus    		= Parameter(index=[aerosol_plus_gases, 4])  # Atmospheric lifetime of gas in iᵗʰ reservior (years).
    a_aerosol_plus    		= Parameter(index=[aerosol_plus_gases, 4])  # Fraction of emissions entering iᵗʰ Aerosol+ gas reservior.
    R0_aerosol_plus    		= Parameter(index=[aerosol_plus_gases, 4])  # Initial model period value for quantity of agent in iᵗʰ atmospheric reservior (unit of E).
    E_aerosol_plus    		= Parameter(index=[time, aerosol_plus_gases])  # Annual Aerosol+ emissions (Tg yr⁻¹).
    Tj       		        = Parameter(index=[time, 3])                # Temperature change for three thermal pools (K).

    α_aerosol_plus           = Variable(index=[time, aerosol_plus_gases])   # State-dependent multiplicative adjustment coefficient of reservior lifetimes.
    aerosol_plus_conc        = Variable(index=[time, aerosol_plus_gases])   # Total atmospheric Aerosol+ concentrations (ppb).
    GA_aerosol_plus          = Variable(index=[time, aerosol_plus_gases])   # Atmospheric burden of agent above pre-industrial levels (unit of E).
    GU_aerosol_plus          = Variable(index=[time, aerosol_plus_gases])   # Cumulative uptake of agent since model initialization (unit of E⁻¹).
    iIRFT100_aerosol_plus    = Variable(index=[time, aerosol_plus_gases])   # 100-year integrated impulse response function (the average airborne fraction over a 100-year period).
    decay_rate_aerosol_plus  = Variable(index=[time, aerosol_plus_gases, 4])    # Decay rates
    R_aerosol_plus           = Variable(index=[time, aerosol_plus_gases, 4])    # Quantity of agent in iᵗʰ atmospheric reservior (unit of E).


    function run_timestep(p, v, d, t)

    	for g in d.aerosol_plus_gases

        	if is_first(t)

            	# Set initial values.
            	v.GU_aerosol_plus[t,g]   = p.GU_aerosol_plus_0[g]
            	v.R_aerosol_plus[t,g,:]  = p.R0_aerosol_plus[g,:]
            	v.aerosol_plus_conc[t,g] = p.aerosol_plus_0[g]

            	# Calculate initial burden above pre-industrial values.
            	v.GA_aerosol_plus[t,g] = sum(v.R_aerosol_plus[t,g,:])

        	else

            	# Calculate iIRF100.
            	v.iIRFT100_aerosol_plus[t,g] = abs(p.r0_aerosol_plus[g] + p.rU_aerosol_plus[g] * (v.GU_aerosol_plus[t-1,g] - v.GA_aerosol_plus[t-1,g]) + p.rT_aerosol_plus[g] * sum(p.Tj[t-1,:]) + p.rA_aerosol_plus[g] * v.GA_aerosol_plus[t-1,g])

            	# Calculate state-dependent lifetime adjustment term based on iIRF100 value.
            	v.α_aerosol_plus[t,g] = p.g0_aerosol_plus[g] * exp(v.iIRFT100_aerosol_plus[t,g] / p.g1_aerosol_plus[g])

            	# Calculate concentrations in each Aerosol+ gas reservoir.
            	for i = 1:4

                	# Calculate decay rate for iᵗʰ reservoir.
                	v.decay_rate_aerosol_plus[t,g,i] = 1.0 / (v.α_aerosol_plus[t,g] * p.τ_aerosol_plus[g,i])

                	# Calculate amount of Aerosol+ gas in iᵗʰ reservoir.
                	v.R_aerosol_plus[t,g,i] = p.E_aerosol_plus[t,g] * p.a_aerosol_plus[g,i] / v.decay_rate_aerosol_plus[t,g,i] * (1.0 - exp(-v.decay_rate_aerosol_plus[t,g,i])) + v.R_aerosol_plus[t-1,g,i] * exp(-v.decay_rate_aerosol_plus[t,g,i])
	            end

    	        # Calcualte atmospheric burden above pre-industiral levels.
            	v.GA_aerosol_plus[t,g] = sum(v.R_aerosol_plus[t,g,:])

            	# Calculate cumulative emissions burden.
            	v.GU_aerosol_plus[t,g] = v.GU_aerosol_plus[t-1,g] + p.E_aerosol_plus[t,g]

            	# Calculate atmospheric CO₂ concentration.
            	v.aerosol_plus_conc[t,g] = p.aerosol_plus_pi[g] + p.emiss2conc_aerosol_plus[g] * (v.GA_aerosol_plus[t-1,g] + v.GA_aerosol_plus[t,g]) / 2.0
        	end
    	end
	end
end
