# --------------------------------------------------
# Radiative Forcing From All Sources
# --------------------------------------------------

@defcomp radiative_forcing begin

    montreal_gases  = Index()                                            # Index for Ozone Depleting Substances controlled under the Montreal Protocol.
    aerosol_plus_gases  = Index()                                            # Index for tropospheric ozone precursors, aerosols and reactive gas
    flourinated_gases  = Index()                                            # Index for flourinated gases controlled under the Kyoto Protocol.

    # CO2, CH4, and N2O parametes
    co2_conc = Parameter(index=[time])
    ch4_conc = Parameter(index=[time])
    n2o_conc = Parameter(index=[time])
    co2_pi = Parameter()
    ch4_pi = Parameter()
    n2o_pi = Parameter()
    co2_f    = Parameter(index=[3]) # f coefficient for
    ch4_f    = Parameter(index=[3]) # f coefficient for
    ch4_o3_f    = Parameter(index=[3]) # f coefficient for
    ch4_h2o_f    = Parameter(index=[3]) # f coefficient for
    n2o_f    = Parameter(index=[3]) # f coefficient for
    n2o_o3_f    = Parameter(index=[3]) # f coefficient for

    # Aerosol+ parameters
    aerosol_plus_conc = Parameter(index=[time, aerosol_plus_gases])
    aerosol_plus_pi = Parameter(index=[aerosol_plus_gases])
    bc_f      = Parameter(index=[3]) # f coefficient for
    bc_snow_f = Parameter(index=[3]) # f coefficient for
    bc_aci_f  = Parameter(index=[3]) # f coefficient for
    co_f    = Parameter(index=[3]) # f coefficient for
    co_o3_f = Parameter(index=[3]) # f coefficient for
    nh3_f               = Parameter(index=[3]) # f coefficient for
    nmvoc_f             = Parameter(index=[3]) # f coefficient for
    nmvoc_o3_f          = Parameter(index=[3]) # f coefficient for
    nox_f               = Parameter(index=[3]) # f coefficient for
    nox_o3_f            = Parameter(index=[3]) # f coefficient for
    nox_avi_f           = Parameter(index=[3]) # f coefficient for
    nox_avi_contrails_f = Parameter(index=[3]) # f coefficient for
    oc_f                = Parameter(index=[3]) # f coefficient for
    oc_aci_f            = Parameter(index=[3]) # f coefficient for
    so2_f               = Parameter(index=[3]) # f coefficient for
    so2_aci_f           = Parameter(index=[3]) # f coefficient for

    # Montreal gas parameters
    montreal_conc = Parameter(index=[time, montreal_gases])
    montreal_pi = Parameter(index=[montreal_gases])
    montreal_f = Parameter(index=[montreal_gases, 3]) # f coefficient
    montreal_ind_f = Parameter(index=[montreal_gases, 3]) # f coefficient for indirect forcing effect on O3

    # Flourinated gas parameters
    flourinated_conc = Parameter(index=[time, flourinated_gases])
    flourinated_pi = Parameter(index=[flourinated_gases])
    flourinated_f = Parameter(index=[flourinated_gases, 3]) # f coefficient

    exogenous_RF = Parameter(index=[time])

    bc_RF = Variable(index=[time])
    bc_snow_RF = Variable(index=[time])
    bc_aci_RF = Variable(index=[time])

    co_RF = Variable(index=[time])
    co_o3_RF = Variable(index=[time])

    nh3_RF = Variable(index=[time])
    nmvoc_RF = Variable(index=[time])
    nmvoc_o3_RF = Variable(index=[time])
    nox_RF = Variable(index=[time])
    nox_o3_RF = Variable(index=[time])
    nox_avi_RF = Variable(index=[time])
    nox_avi_contrails_RF = Variable(index=[time])
    oc_RF = Variable(index=[time])
    oc_aci_RF = Variable(index=[time])
    so2_RF = Variable(index=[time])
    so2_aci_RF = Variable(index=[time])

    co2_RF = Variable(index=[time])

    ch4_RF = Variable(index=[time])
    ch4_o3_RF = Variable(index=[time])
    ch4_h2o_RF = Variable(index=[time])

    n2o_RF = Variable(index=[time])
    n2o_o3_RF = Variable(index=[time])

    montreal_RF = Variable(index=[time, montreal_gases])
    montreal_ind_RF = Variable(index=[time, montreal_gases])

    flourinated_RF = Variable(index=[time, flourinated_gases])

    total_RF = Variable(index=[time])

    function run_timestep(p, v, d, t)


        v.co2_RF[t] = calculate_RF(p.co2_f, p.co2_conc[t], p.co2_pi)

        v.ch4_RF[t] = calculate_RF(p.ch4_f, p.ch4_conc[t], p.ch4_pi)
        v.ch4_o3_RF[t] = calculate_RF(p.ch4_o3_f, p.ch4_conc[t], p.ch4_pi)
        v.ch4_h2o_RF[t] = calculate_RF(p.ch4_h2o_f, p.ch4_conc[t], p.ch4_pi)

        v.n2o_RF[t] = calculate_RF(p.n2o_f, p.n2o_conc[t], p.n2o_pi)
        v.n2o_o3_RF[t] = calculate_RF(p.n2o_o3_f, p.n2o_conc[t], p.n2o_pi)

        # Radiative forcing for tropospheric ozone precursors, aerosols & reactive gas cycles (Aerosol+).
        # Note: gases in 'aerosol_plus_conc' ordered alphabetically [BC, CO, NH₃, NMVOCs, NOₓ, NOₓ(avi), OC, SO₂].

        # Radiative forcing for black carbon and indirect effects on 
        v.bc_RF[t]      = calculate_RF(p.bc_f,      p.aerosol_plus_conc[t,1], p.aerosol_plus_pi[1])
        v.bc_snow_RF[t] = calculate_RF(p.bc_snow_f, p.aerosol_plus_conc[t,1], p.aerosol_plus_pi[1])
        v.bc_aci_RF[t]  = calculate_RF(p.bc_aci_f,  p.aerosol_plus_conc[t,1], p.aerosol_plus_pi[1])

        # Radiative forcing for oc and indirect effects on 
        v.co_RF[t]    = calculate_RF(p.co_f,    p.aerosol_plus_conc[t,2], p.aerosol_plus_pi[2])
        v.co_o3_RF[t] = calculate_RF(p.co_o3_f, p.aerosol_plus_conc[t,2], p.aerosol_plus_pi[2])

        # Radiative forcing for nh3
        v.nh3_RF[t]    = calculate_RF(p.nh3_f,    p.aerosol_plus_conc[t,3], p.aerosol_plus_pi[3])

        # Radiative forcing for nmvocs and indirect effects on 
        v.nmvoc_RF[t]    = calculate_RF(p.nmvoc_f,    p.aerosol_plus_conc[t,4], p.aerosol_plus_pi[4])
        v.nmvoc_o3_RF[t]    = calculate_RF(p.nmvoc_o3_f,    p.aerosol_plus_conc[t,4], p.aerosol_plus_pi[4])

        # Radiative forcing for nox and indirect effects on 
        v.nox_RF[t]    = calculate_RF(p.nox_f,    p.aerosol_plus_conc[t,5], p.aerosol_plus_pi[5])
        v.nox_o3_RF[t]    = calculate_RF(p.nox_o3_f,    p.aerosol_plus_conc[t,5], p.aerosol_plus_pi[5])

        # Radiative forcing for nox-avi and indirect effects on 
        v.nox_avi_RF[t]    = calculate_RF(p.nox_avi_f,    p.aerosol_plus_conc[t,6], p.aerosol_plus_pi[6])
        v.nox_avi_contrails_RF[t]    = calculate_RF(p.nox_avi_contrails_f,    p.aerosol_plus_conc[t,6], p.aerosol_plus_pi[6])

        # Radiative forcing for oc and indirect effects on 
        v.oc_RF[t]    = calculate_RF(p.oc_f,    p.aerosol_plus_conc[t,7], p.aerosol_plus_pi[7])
        v.oc_aci_RF[t] = calculate_RF(p.oc_aci_f, p.aerosol_plus_conc[t,7], p.aerosol_plus_pi[7])

        v.so2_RF[t]    = calculate_RF(p.so2_f,    p.aerosol_plus_conc[t,8], p.aerosol_plus_pi[8])
        v.so2_aci_RF[t] = calculate_RF(p.so2_aci_f, p.aerosol_plus_conc[t,8], p.aerosol_plus_pi[8])

	    for g in d.montreal_gases
    		v.montreal_RF[t,g]     = calculate_RF(p.montreal_f[g,:],     p.montreal_conc[t,g], p.montreal_pi[g])
            v.montreal_ind_RF[t,g] = calculate_RF(p.montreal_ind_f[g,:], p.montreal_conc[t,g], p.montreal_pi[g])
    	end

        # Flourinated gases (note, no indirect forcing effects for flourinated gases.)
        for g in d.flourinated_gases
            v.flourinated_RF[t,g]     = calculate_RF(p.flourinated_f[g,:],     p.flourinated_conc[t,g], p.flourinated_pi[g])
        end

    	# Calculate total radiative forcing from all sources.
    	v.total_RF[t] = v.co2_RF[t] + v.ch4_RF[t] + v.n2o_RF[t] +
                        v.ch4_o3_RF[t] + v.ch4_h2o_RF[t] + v.n2o_o3_RF[t] +
                        v.bc_RF[t] + v.co_RF[t] + v.nh3_RF[t] + v.nmvoc_RF[t] + v.nox_RF[t] + v.nox_avi_RF[t] + v.oc_RF[t] + v.so2_RF[t] +
                        v.bc_snow_RF[t] + v.bc_aci_RF[t] + v.co_o3_RF[t] + v.nmvoc_o3_RF[t] + v.nox_o3_RF[t] + v.nox_avi_contrails_RF[t] + v.oc_aci_RF[t] + v.so2_aci_RF[t] +
                        sum(v.montreal_RF[t,:]) + sum(v.montreal_ind_RF[t,:]) + sum(v.flourinated_RF[t,:]) +
                        p.exogenous_RF[t]
    end
end
