#######################################################################################################################
# CREATE A FUNCTION TO RUN A MONTE CARLO WITH FAIRv2.0
#######################################################################################################################
# Description: This file contains functions to run a Monte Carlo with MimiFAIRv2 using the constrained parameters from
#              Leach et al. (2021). The first function, "create_fair_monte_carlo" loads the constrained parameter
#              samples from Leach et al. (2021) and cleans them up so they can be easily passed into Mimi-FAIRv2.0. It
#              then creates a second function, "fair_monte_carlo" that performs the actual analysis.
#
# Function Arguments:
#
#       n_samples:          Number of samples to randomly draw from original constrained parameter values (max = 9,995).
#       emissions_scenario: Current options are: "ssp119", "ssp126", "ssp245", "ssp370", "ssp585"
#       start_year:         First year to run the model (note, Mimi-FAIR requires user-supplied initial conditions if not starting in 1765).
#       end_year:           Final year to run the model.
#----------------------------------------------------------------------------------------------------------------------

# Load packages.
using CSVFiles, DataFrames, Mimi, MimiFAIRv2, StatsBase


function create_fair_monte_carlo_1000(;n_samples::Int=1000, emissions_scenario::String="ssp585", start_year::Int=2000, end_year::Int=2300)

    #Load constrained FAIR parameters from Leach et al. (2021).
    thermal_params           = DataFrame(load(joinpath(@__DIR__, "..", "data", "constrained_parameters", "julia_constrained_thermal_parameters_average_probs_1000.csv")))
    gas_params               = DataFrame(load(joinpath(@__DIR__, "..", "data", "constrained_parameters", "julia_constrained_gas_cycle_parameters_average_probs_1000.csv")))
    indirect_forcing_params  = DataFrame(load(joinpath(@__DIR__, "..", "data", "constrained_parameters", "julia_constrained_indirect_forcing_parameters_average_probs_1000.csv")))
    exogenous_forcing_params = DataFrame(load(joinpath(@__DIR__, "..", "data", "constrained_parameters", "julia_constrained_exogenous_forcing_parameters_average_probs_1000.csv")))

    # Load sample-specific initial conditions (defaults for the year 2000 as initial year).
    thermal_initial_conditions   = DataFrame(load(joinpath(@__DIR__, "..", "data", "constrained_parameters", "initial_thermal_conditions_constrained_average_probs_1000.csv")))
    gas_cycle_initial_conditions = DataFrame(load(joinpath(@__DIR__, "..", "data", "constrained_parameters", "initial_gas_cycle_conditions_constrained_average_probs_1000.csv")))

    # Get FAIR sample id values from subsetted data (for indexing data).
    sample_id_subset = sort(thermal_params.sample_id)

    # Isolate individual and group gas parameters for convenience.
    co2_p = filter(:gas_name  => ==("carbon_dioxide"), gas_params)
    ch4_p = filter(:gas_name  => ==("methane"), gas_params)
    n2o_p = filter(:gas_name  => ==("nitrous_oxide"), gas_params)

    montreal_gas_p     = filter(:gas_group => ==("montreal"), gas_params)
    montreal_indirect  = filter(:gas_group => ==("montreal"), indirect_forcing_params)
    flourinated_gas_p  = filter(:gas_group => ==("flourinated"), gas_params)
    aerosol_plus_gas_p = filter(:gas_group => ==("aerosol_plus"), gas_params)

    # Get number of gases in each gas group.
    n_montreal     = length(unique(montreal_gas_p.gas_name))
    n_flourinated  = length(unique(flourinated_gas_p.gas_name))
    n_aerosol_plus = length(unique(aerosol_plus_gas_p.gas_name))

    # Isolate thermal pool initial conditions.
    thermal_pool_init = Array(thermal_initial_conditions[:, [:thermal_box_1, :thermal_box_2, :thermal_box_3]])

    # Isolate initial conditions for primary gases and gas groupings.
    co2_init = filter(:gas_name  => ==("carbon_dioxide"), gas_cycle_initial_conditions)
    ch4_init = filter(:gas_name  => ==("methane"), gas_cycle_initial_conditions)
    n2o_init = filter(:gas_name  => ==("nitrous_oxide"), gas_cycle_initial_conditions)

    montreal_init     = filter(:gas_group  => ==("montreal"), gas_cycle_initial_conditions)
    flourinated_init  = filter(:gas_group  => ==("flourinated"), gas_cycle_initial_conditions)
    aerosol_plus_init = filter(:gas_group  => ==("aerosol_plus"), gas_cycle_initial_conditions)

    # Isolate primary gas cycle pool initial conditions.
    co2_pool_init = Array(co2_init[:, [:pool1, :pool2, :pool3, :pool4]])
    ch4_pool_init = Array(ch4_init[:, [:pool1, :pool2, :pool3, :pool4]])
    n2o_pool_init = Array(n2o_init[:, [:pool1, :pool2, :pool3, :pool4]])

    # Create arrays for 'a' and 'τ; parameter groups in primary gases.
    co2_a = Array(co2_p[:,[:a1,:a2,:a3,:a4]])
    ch4_a = Array(ch4_p[:,[:a1,:a2,:a3,:a4]])
    n2o_a = Array(n2o_p[:,[:a1,:a2,:a3,:a4]])
    co2_τ = Array(co2_p[:,[:tau1,:tau2,:tau3,:tau4]])
    ch4_τ = Array(ch4_p[:,[:tau1,:tau2,:tau3,:tau4]])
    n2o_τ = Array(n2o_p[:,[:tau1,:tau2,:tau3,:tau4]])

    # Create arrays for direct and indirect forcing 'f' parameters in primarty gases.
    co2_f     = Array(co2_p[:, [:f1, :f2, :f3]])
    ch4_f     = Array(ch4_p[:, [:f1, :f2, :f3]])
    n2o_f     = Array(n2o_p[:, [:f1, :f2, :f3]])
    ch4_o3_f  = indirect_forcing_params[(indirect_forcing_params.indirect_forcing_effect .== "methane|o3"), [:f1, :f2, :f3]]
    ch4_h2o_f = indirect_forcing_params[(indirect_forcing_params.indirect_forcing_effect .== "methane|strat_h2o"), [:f1, :f2, :f3]]
    n2o_o3_f  = indirect_forcing_params[(indirect_forcing_params.indirect_forcing_effect .== "nitrous_oxide|o3"), [:f1, :f2, :f3]]

    # Calculate constants for primary gases to approximate numerical solution for state-dependent timescale adjustment factor from Millar et al. (2017).
    g0_co2, g1_co2 =  MimiFAIRv2.calculate_g0_g1(co2_a, co2_τ)
    g0_ch4, g1_ch4 =  MimiFAIRv2.calculate_g0_g1(ch4_a, ch4_τ)
    g0_n2o, g1_n2o =  MimiFAIRv2.calculate_g0_g1(n2o_a, n2o_τ)

    # Iniitialize parameter arrays for different non-primary gas groups.
    montreal_a           = zeros(n_montreal, 4, n_samples)
    montreal_τ           = zeros(n_montreal, 4, n_samples)
    montreal_init_pools  = zeros(n_montreal, 4, n_samples)
    montreal_f           = zeros(n_montreal, 3, n_samples)
    montreal_ind_f       = zeros(n_montreal, 3, n_samples)
    montreal_r0          = zeros(n_montreal, n_samples)
    montreal_rC          = zeros(n_montreal, n_samples)
    montreal_rT          = zeros(n_montreal, n_samples)
    montreal_rA          = zeros(n_montreal, n_samples)
    montreal_pi_conc     = zeros(n_montreal, n_samples)
    montreal_init_conc   = zeros(n_montreal, n_samples)
    montreal_init_uptake = zeros(n_montreal, n_samples)


    flourinated_a        = zeros(n_flourinated, 4, n_samples)
    flourinated_τ        = zeros(n_flourinated, 4, n_samples)
    flourinated_init_pools  = zeros(n_flourinated, 4, n_samples)
    flourinated_f        = zeros(n_flourinated, 3, n_samples)
    flourinated_r0       = zeros(n_flourinated, n_samples)
    flourinated_rC       = zeros(n_flourinated, n_samples)
    flourinated_rT       = zeros(n_flourinated, n_samples)
    flourinated_rA       = zeros(n_flourinated, n_samples)
    flourinated_pi_conc  = zeros(n_flourinated, n_samples)
    flourinated_init_conc   = zeros(n_flourinated, n_samples)
    flourinated_init_uptake   = zeros(n_flourinated, n_samples)

    aerosol_plus_a       = zeros(n_aerosol_plus, 4, n_samples)
    aerosol_plus_τ       = zeros(n_aerosol_plus, 4, n_samples)
    aerosol_plus_init_pools  = zeros(n_aerosol_plus, 4, n_samples)
    aerosol_plus_r0      = zeros(n_aerosol_plus, n_samples)
    aerosol_plus_rC      = zeros(n_aerosol_plus, n_samples)
    aerosol_plus_rT      = zeros(n_aerosol_plus, n_samples)
    aerosol_plus_rA      = zeros(n_aerosol_plus, n_samples)
    aerosol_plus_pi_conc = zeros(n_aerosol_plus, n_samples)
    aerosol_plus_init_conc   = zeros(n_aerosol_plus, n_samples)
    aerosol_plus_init_uptake   = zeros(n_aerosol_plus, n_samples)

    # Iniitialize direct and indirect forcing parameter 'f' arrays for individual aerosol+ gaeses.
    bc_f                = zeros(n_samples, 3)
    bc_snow_f           = zeros(n_samples, 3)
    bc_aci_f            = zeros(n_samples, 3)
    co_f                = zeros(n_samples, 3)
    co_o3_f             = zeros(n_samples, 3)
    nh3_f               = zeros(n_samples, 3)
    nmvoc_f             = zeros(n_samples, 3)
    nmvoc_o3_f          = zeros(n_samples, 3)
    nox_f               = zeros(n_samples, 3)
    nox_o3_f            = zeros(n_samples, 3)
    nox_avi_f           = zeros(n_samples, 3)
    nox_avi_contrails_f = zeros(n_samples, 3)
    oc_f                = zeros(n_samples, 3)
    oc_aci_f            = zeros(n_samples, 3)
    so2_f               = zeros(n_samples, 3)
    so2_aci_f           = zeros(n_samples, 3)

    # Loop through constrained parameter samples and reshape into arrays that can be passed into Mimi-FAIRv2.0
    for i = 1:n_samples
        # Montreal gases.
        montreal_a[:,:,i]         .= montreal_gas_p[findall(x->x==sample_id_subset[i], montreal_gas_p.sample_id), [:a1,:a2,:a3,:a4]]
        montreal_τ[:,:,i]         .= montreal_gas_p[findall(x->x==sample_id_subset[i], montreal_gas_p.sample_id), [:tau1,:tau2,:tau3,:tau4]]
        montreal_f[:,:,i]         .= montreal_gas_p[findall(x->x==sample_id_subset[i], montreal_gas_p.sample_id), [:f1,:f2,:f3]]
        montreal_ind_f[:,:,i]     .= montreal_indirect[findall(x->x==sample_id_subset[i], montreal_indirect.sample_id), [:f1,:f2,:f3]]
        montreal_r0[:,i]          .= montreal_gas_p[findall(x->x==sample_id_subset[i], montreal_gas_p.sample_id), :r0]
        montreal_rC[:,i]          .= montreal_gas_p[findall(x->x==sample_id_subset[i], montreal_gas_p.sample_id), :rC]
        montreal_rT[:,i]          .= montreal_gas_p[findall(x->x==sample_id_subset[i], montreal_gas_p.sample_id), :rT]
        montreal_rA[:,i]          .= montreal_gas_p[findall(x->x==sample_id_subset[i], montreal_gas_p.sample_id), :rA]
        montreal_pi_conc[:,i]     .= montreal_gas_p[findall(x->x==sample_id_subset[i], montreal_gas_p.sample_id), :PI_conc]
        montreal_init_conc[:,i]   .= montreal_init[findall(x->x==sample_id_subset[i], montreal_init.sample_id), :concentration]
        montreal_init_uptake[:,i]   .= montreal_init[findall(x->x==sample_id_subset[i], montreal_init.sample_id), :cumulative_uptake]
        montreal_init_pools[:,:,i]   .= montreal_init[findall(x->x==sample_id_subset[i], montreal_init.sample_id), [:pool1, :pool2, :pool3, :pool4]]

        # Flourinated gases.
        flourinated_a[:,:,i]      .= flourinated_gas_p[findall(x->x==sample_id_subset[i], flourinated_gas_p.sample_id), [:a1,:a2,:a3,:a4]]
        flourinated_τ[:,:,i]      .= flourinated_gas_p[findall(x->x==sample_id_subset[i], flourinated_gas_p.sample_id), [:tau1,:tau2,:tau3,:tau4]]
        flourinated_f[:,:,i]      .= flourinated_gas_p[findall(x->x==sample_id_subset[i], flourinated_gas_p.sample_id), [:f1,:f2,:f3]]
        flourinated_r0[:,i]       .= flourinated_gas_p[findall(x->x==sample_id_subset[i], flourinated_gas_p.sample_id), :r0]
        flourinated_rC[:,i]       .= flourinated_gas_p[findall(x->x==sample_id_subset[i], flourinated_gas_p.sample_id), :rC]
        flourinated_rT[:,i]       .= flourinated_gas_p[findall(x->x==sample_id_subset[i], flourinated_gas_p.sample_id), :rT]
        flourinated_rA[:,i]       .= flourinated_gas_p[findall(x->x==sample_id_subset[i], flourinated_gas_p.sample_id), :rA]
        flourinated_pi_conc[:,i]  .= flourinated_gas_p[findall(x->x==sample_id_subset[i], flourinated_gas_p.sample_id), :PI_conc]
        flourinated_init_conc[:,i] .= flourinated_init[findall(x->x==sample_id_subset[i], flourinated_init.sample_id), :concentration]
        flourinated_init_uptake[:,i] .= flourinated_init[findall(x->x==sample_id_subset[i], flourinated_init.sample_id), :cumulative_uptake]
        flourinated_init_pools[:,:,i]   .= flourinated_init[findall(x->x==sample_id_subset[i], flourinated_init.sample_id), [:pool1, :pool2, :pool3, :pool4]]

        # Aerosol+ gases.
        aerosol_plus_a[:,:,i]     .= aerosol_plus_gas_p[findall(x->x==sample_id_subset[i], aerosol_plus_gas_p.sample_id), [:a1,:a2,:a3,:a4]]
        aerosol_plus_τ[:,:,i]     .= aerosol_plus_gas_p[findall(x->x==sample_id_subset[i], aerosol_plus_gas_p.sample_id), [:tau1,:tau2,:tau3,:tau4]]
        aerosol_plus_r0[:,i]      .= aerosol_plus_gas_p[findall(x->x==sample_id_subset[i], aerosol_plus_gas_p.sample_id), :r0]
        aerosol_plus_rC[:,i]      .= aerosol_plus_gas_p[findall(x->x==sample_id_subset[i], aerosol_plus_gas_p.sample_id), :rC]
        aerosol_plus_rT[:,i]      .= aerosol_plus_gas_p[findall(x->x==sample_id_subset[i], aerosol_plus_gas_p.sample_id), :rT]
        aerosol_plus_rA[:,i]      .= aerosol_plus_gas_p[findall(x->x==sample_id_subset[i], aerosol_plus_gas_p.sample_id), :rA]
        aerosol_plus_pi_conc[:,i] .= aerosol_plus_gas_p[findall(x->x==sample_id_subset[i], aerosol_plus_gas_p.sample_id), :PI_conc]
        aerosol_plus_init_conc[:,i] .= aerosol_plus_init[findall(x->x==sample_id_subset[i], aerosol_plus_init.sample_id), :concentration]
        aerosol_plus_init_uptake[:,i] .= aerosol_plus_init[findall(x->x==sample_id_subset[i], aerosol_plus_init.sample_id), :cumulative_uptake]
        aerosol_plus_init_pools[:,:,i]   .= aerosol_plus_init[findall(x->x==sample_id_subset[i], aerosol_plus_init.sample_id), [:pool1, :pool2, :pool3, :pool4]]

        # Aerosol+ forcing parameters.
        bc_f[i,:]                 = Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name .== "bc") .& (aerosol_plus_gas_p.sample_id .== sample_id_subset[i]), [:f1,:f2,:f3]])
        bc_snow_f[i,:]            = Array(indirect_forcing_params[(indirect_forcing_params.indirect_forcing_effect .== "bc|bc_on_snow") .& (indirect_forcing_params.sample_id .== sample_id_subset[i]), [:f1,:f2,:f3]])
        bc_aci_f[i,:]             = Array(indirect_forcing_params[(indirect_forcing_params.indirect_forcing_effect .== "bc|aci") .& (indirect_forcing_params.sample_id .== sample_id_subset[i]), [:f1,:f2,:f3]])
        co_f[i,:]                 = Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name .== "co") .& (aerosol_plus_gas_p.sample_id .== sample_id_subset[i]), [:f1,:f2,:f3]])
        co_o3_f[i,:]              = Array(indirect_forcing_params[(indirect_forcing_params.indirect_forcing_effect .== "co|o3") .& (indirect_forcing_params.sample_id .== sample_id_subset[i]), [:f1,:f2,:f3]])
        nh3_f[i,:]                = Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name .== "nh3") .& (aerosol_plus_gas_p.sample_id .== sample_id_subset[i]), [:f1,:f2,:f3]])
        nmvoc_f[i,:]              = Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name .== "nmvoc") .& (aerosol_plus_gas_p.sample_id .== sample_id_subset[i]), [:f1,:f2,:f3]])
        nmvoc_o3_f[i,:]           = Array(indirect_forcing_params[(indirect_forcing_params.indirect_forcing_effect .== "nmvoc|o3") .& (indirect_forcing_params.sample_id .== sample_id_subset[i]), [:f1,:f2,:f3]])
        nox_f[i,:]                = Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name .== "nox") .& (aerosol_plus_gas_p.sample_id .== sample_id_subset[i]), [:f1,:f2,:f3]])
        nox_o3_f[i,:]             = Array(indirect_forcing_params[(indirect_forcing_params.indirect_forcing_effect .== "nox|o3") .& (indirect_forcing_params.sample_id .== sample_id_subset[i]), [:f1,:f2,:f3]])
        nox_avi_f[i,:]            = Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name .== "nox_avi") .& (aerosol_plus_gas_p.sample_id .== sample_id_subset[i]), [:f1,:f2,:f3]])
        nox_avi_contrails_f[i,:]  = Array(indirect_forcing_params[(indirect_forcing_params.indirect_forcing_effect .== "nox_avi|contrails") .& (indirect_forcing_params.sample_id .== sample_id_subset[i]), [:f1,:f2,:f3]])
        oc_f[i,:]                 = Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name .== "oc") .& (aerosol_plus_gas_p.sample_id .== sample_id_subset[i]), [:f1,:f2,:f3]])
        oc_aci_f[i,:]             = Array(indirect_forcing_params[(indirect_forcing_params.indirect_forcing_effect .== "oc|aci") .& (indirect_forcing_params.sample_id .== sample_id_subset[i]), [:f1,:f2,:f3]])
        so2_f[i,:]                = Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name .== "so2") .& (aerosol_plus_gas_p.sample_id .== sample_id_subset[i]), [:f1,:f2,:f3]])
        so2_aci_f[i,:]            = Array(indirect_forcing_params[(indirect_forcing_params.indirect_forcing_effect .== "so2|aci") .& (indirect_forcing_params.sample_id .== sample_id_subset[i]), [:f1,:f2,:f3]])
    end

    # Iniitlaize arrays for g1 and g0 parameters.
    montreal_g0     = zeros(n_montreal, n_samples)
    montreal_g1     = zeros(n_montreal, n_samples)
    flourinated_g0  = zeros(n_flourinated, n_samples)
    flourinated_g1  = zeros(n_flourinated, n_samples)
    aerosol_plus_g0 = zeros(n_aerosol_plus, n_samples)
    aerosol_plus_g1 = zeros(n_aerosol_plus, n_samples)

    # Calculate g1 and g0 parameters for non-primary gas groups.
    for i = 1:n_samples
        montreal_g0[:,i],     montreal_g1[:,i]     = MimiFAIRv2.calculate_g0_g1(montreal_a[:,:,i], montreal_τ[:,:,i])
        flourinated_g0[:,i],  flourinated_g1[:,i]  = MimiFAIRv2.calculate_g0_g1(flourinated_a[:,:,i], flourinated_τ[:,:,i])
        aerosol_plus_g0[:,i], aerosol_plus_g1[:,i] = MimiFAIRv2.calculate_g0_g1(aerosol_plus_a[:,:,i], aerosol_plus_τ[:,:,i])
    end

    #Calculate thermal decay factors, defined as exp(-1/d).
    thermal_decay_factors = exp.(-1.0 ./ Array(thermal_params[:,[:d1,:d2,:d3]]))

    # Initialize an array to store FAIR projections.
    temperatures = zeros(length(start_year:end_year), n_samples)  # Global mean surface temperature anomaly (K)
    rf = zeros(length(start_year:end_year), n_samples)            # Total radiative forcing, with individual components scaled by their respective efficacy (Wm⁻²)
    co2 = zeros(length(start_year:end_year), n_samples)           # Total atmospheric carbon dioxide concentrations (ppm).
    ch4 = zeros(length(start_year:end_year), n_samples)           # Total atmospheric methane concentrations (ppb)
    n2o = zeros(length(start_year:end_year), n_samples)           # Total atmospheric nitrous oxide concentrations (ppb).

    # Load an instance of FAIR with user-specificed settings.
    fair_raw = MimiFAIRv2.get_model(emissions_forcing_scenario=emissions_scenario, start_year=start_year, end_year=end_year)

    # Create a model instance to speed things up.
    fair = Mimi.build(fair_raw)


    # Create a function to carry out the actual Monte Carlo analysis (passing in sampled constrained parameter values).
    function fair_monte_carlo_1000( ;co2_em_vals::Union{Nothing, Vector{Vector{Float64}}} = nothing,
                                    n2o_em_vals::Union{Nothing, Vector{Vector{Float64}}} = nothing,
                                    ch4_em_vals::Union{Nothing, Vector{Vector{Float64}}} = nothing)

        for i = 1:n_samples

            # ---- Emissions Trajectories ---- #
            if !(isnothing(co2_em_vals))
                update_param!(fair, :co2_cycle, :E_co2, co2_em_vals[i])
            end
            if !(isnothing(n2o_em_vals))
                update_param!(fair, :n2o_cycle, :E_n2o, n2o_em_vals[i])
            end
            if !(isnothing(ch4_em_vals))
                update_param!(fair, :ch4_cycle, :E_ch4, ch4_em_vals[i])
            end

            # ---- Global Temperature Anomaly ---- #
            update_param!(fair, :temperature, :decay_factor, thermal_decay_factors[i,:])
            update_param!(fair, :temperature, :q, Array(thermal_params[i,[:q1,:q2,:q3]]))
            update_param!(fair, :temperature, :Tj_0, thermal_pool_init[i,:])
            update_param!(fair, :temperature, :T_0, thermal_initial_conditions.global_temp_anomaly[i])

            # ---- Carbon Cycle ---- #
            update_param!(fair, :co2_cycle, :r0_co2, co2_p.r0[i])
            update_param!(fair, :co2_cycle, :rU_co2, co2_p.rC[i])
            update_param!(fair, :co2_cycle, :rT_co2, co2_p.rT[i])
            update_param!(fair, :co2_cycle, :rA_co2, co2_p.rA[i])
            update_param!(fair, :co2_cycle, :g0_co2, g0_co2[i])
            update_param!(fair, :co2_cycle, :g1_co2, g1_co2[i])
            update_param!(fair, :co2_cycle, :a_co2,  co2_a[i,:])
            update_param!(fair, :co2_cycle, :τ_co2,  co2_τ[i,:])
            update_param!(fair, :co2_cycle, :co2_0, co2_init.concentration[i])
            update_param!(fair, :co2_cycle, :GU_co2_0, co2_init.cumulative_uptake[i])
            update_param!(fair, :co2_cycle, :R0_co2, co2_pool_init[i,:])

            # ---- Methane Cycle ---- #
            update_param!(fair, :ch4_cycle, :r0_ch4, ch4_p.r0[i])
            update_param!(fair, :ch4_cycle, :rU_ch4, ch4_p.rC[i])
            update_param!(fair, :ch4_cycle, :rT_ch4, ch4_p.rT[i])
            update_param!(fair, :ch4_cycle, :rA_ch4, ch4_p.rA[i])
            update_param!(fair, :ch4_cycle, :g0_ch4, g0_ch4[i])
            update_param!(fair, :ch4_cycle, :g1_ch4, g1_ch4[i])
            update_param!(fair, :ch4_cycle, :a_ch4,  ch4_a[i,:])
            update_param!(fair, :ch4_cycle, :τ_ch4,  ch4_τ[i,:])
            update_param!(fair, :ch4_cycle, :ch4_0, ch4_init.concentration[i])
            update_param!(fair, :ch4_cycle, :GU_ch4_0, ch4_init.cumulative_uptake[i])
            update_param!(fair, :ch4_cycle, :R0_ch4, ch4_pool_init[i,:])

            # ---- Nitrous Oxide Cycle ---- #
            update_param!(fair, :n2o_cycle, :r0_n2o, n2o_p.r0[i])
            update_param!(fair, :n2o_cycle, :rU_n2o, n2o_p.rC[i])
            update_param!(fair, :n2o_cycle, :rT_n2o, n2o_p.rT[i])
            update_param!(fair, :n2o_cycle, :rA_n2o, n2o_p.rA[i])
            update_param!(fair, :n2o_cycle, :g0_n2o, g0_n2o[i])
            update_param!(fair, :n2o_cycle, :g1_n2o, g1_n2o[i])
            update_param!(fair, :n2o_cycle, :a_n2o,  n2o_a[i,:])
            update_param!(fair, :n2o_cycle, :τ_n2o,  n2o_τ[i,:])
            update_param!(fair, :n2o_cycle, :n2o_0, n2o_init.concentration[i])
            update_param!(fair, :n2o_cycle, :GU_n2o_0, n2o_init.cumulative_uptake[i])
            update_param!(fair, :n2o_cycle, :R0_n2o, n2o_pool_init[i,:])

            # ---- Flourinated Gas Cycles ---- #
            update_param!(fair, :flourinated_cycles, :r0_flourinated, flourinated_r0[:,i])
            update_param!(fair, :flourinated_cycles, :rU_flourinated, flourinated_rC[:,i])
            update_param!(fair, :flourinated_cycles, :rT_flourinated, flourinated_rT[:,i])
            update_param!(fair, :flourinated_cycles, :rA_flourinated, flourinated_rA[:,i])
            update_param!(fair, :flourinated_cycles, :g0_flourinated, flourinated_g0[:,i])
            update_param!(fair, :flourinated_cycles, :g1_flourinated, flourinated_g1[:,i])
            update_param!(fair, :flourinated_cycles, :a_flourinated,  flourinated_a[:,:,i])
            update_param!(fair, :flourinated_cycles, :τ_flourinated,  flourinated_τ[:,:,i])
            update_param!(fair, :flourinated_cycles, :flourinated_0, flourinated_init_conc[:,i])
            update_param!(fair, :flourinated_cycles, :GU_flourinated_0, flourinated_init_uptake[:,i])
            update_param!(fair, :flourinated_cycles, :R0_flourinated, flourinated_init_pools[:,:,i])

            # ---- Montreal Protocol Gas Cycles ---- #
            update_param!(fair, :montreal_cycles, :r0_montreal, montreal_r0[:,i])
            update_param!(fair, :montreal_cycles, :rU_montreal, montreal_rC[:,i])
            update_param!(fair, :montreal_cycles, :rT_montreal, montreal_rT[:,i])
            update_param!(fair, :montreal_cycles, :rA_montreal, montreal_rA[:,i])
            update_param!(fair, :montreal_cycles, :g0_montreal, montreal_g0[:,i])
            update_param!(fair, :montreal_cycles, :g1_montreal, montreal_g1[:,i])
            update_param!(fair, :montreal_cycles, :a_montreal,  montreal_a[:,:,i])
            update_param!(fair, :montreal_cycles, :τ_montreal,  montreal_τ[:,:,i])
            update_param!(fair, :montreal_cycles, :montreal_0, montreal_init_conc[:,i])
            update_param!(fair, :montreal_cycles, :GU_montreal_0, montreal_init_uptake[:,i])
            update_param!(fair, :montreal_cycles, :R0_montreal, montreal_init_pools[:,:,i])

            # ---- Tropospheric Ozone Precursors, Aerosols, & Reactive Gas Cycles (Aerosol+) ---- #
            update_param!(fair, :aerosol_plus_cycles, :r0_aerosol_plus, aerosol_plus_r0[:,i])
            update_param!(fair, :aerosol_plus_cycles, :rU_aerosol_plus, aerosol_plus_rC[:,i])
            update_param!(fair, :aerosol_plus_cycles, :rT_aerosol_plus, aerosol_plus_rT[:,i])
            update_param!(fair, :aerosol_plus_cycles, :rA_aerosol_plus, aerosol_plus_rA[:,i])
            update_param!(fair, :aerosol_plus_cycles, :g0_aerosol_plus, aerosol_plus_g0[:,i])
            update_param!(fair, :aerosol_plus_cycles, :g1_aerosol_plus, aerosol_plus_g1[:,i])
            update_param!(fair, :aerosol_plus_cycles, :a_aerosol_plus,  aerosol_plus_a[:,:,i])
            update_param!(fair, :aerosol_plus_cycles, :τ_aerosol_plus,  aerosol_plus_τ[:,:,i])
            update_param!(fair, :aerosol_plus_cycles, :aerosol_plus_0, aerosol_plus_init_conc[:,i])
            update_param!(fair, :aerosol_plus_cycles, :GU_aerosol_plus_0, aerosol_plus_init_uptake[:,i])
            update_param!(fair, :aerosol_plus_cycles, :R0_aerosol_plus, aerosol_plus_init_pools[:,:,i])

            # ---- Radiative Forcing ---- #
            update_param!(fair, :radiative_forcing, :co2_f, co2_f[i,:])
            update_param!(fair, :radiative_forcing, :ch4_f, ch4_f[i,:])
            update_param!(fair, :radiative_forcing, :ch4_o3_f, ch4_o3_f[i,:])
            update_param!(fair, :radiative_forcing, :ch4_h2o_f, ch4_h2o_f[i,:])
            update_param!(fair, :radiative_forcing, :n2o_f, n2o_f[i,:])
            update_param!(fair, :radiative_forcing, :n2o_o3_f, n2o_o3_f[i,:])
            update_param!(fair, :radiative_forcing, :montreal_f, montreal_f[:,:,i])
            update_param!(fair, :radiative_forcing, :montreal_ind_f, montreal_ind_f[:,:,i])
            update_param!(fair, :radiative_forcing, :flourinated_f,flourinated_f[:,:,i])
            update_param!(fair, :radiative_forcing, :bc_f, bc_f[i,:])
            update_param!(fair, :radiative_forcing, :bc_snow_f, bc_snow_f[i,:])
            update_param!(fair, :radiative_forcing, :bc_aci_f, bc_aci_f[i,:])
            update_param!(fair, :radiative_forcing, :co_f, co_f[i,:])
            update_param!(fair, :radiative_forcing, :co_o3_f, co_o3_f[i,:])
            update_param!(fair, :radiative_forcing, :nh3_f, nh3_f[i,:])
            update_param!(fair, :radiative_forcing, :nmvoc_f, nmvoc_f[i,:])
            update_param!(fair, :radiative_forcing, :nmvoc_o3_f, nmvoc_o3_f[i,:])
            update_param!(fair, :radiative_forcing, :nox_f, nox_f[i,:])
            update_param!(fair, :radiative_forcing, :nox_o3_f, nox_o3_f[i,:])
            update_param!(fair, :radiative_forcing, :nox_avi_f, nox_avi_f[i,:])
            update_param!(fair, :radiative_forcing, :nox_avi_contrails_f, nox_avi_contrails_f[i,:])
            update_param!(fair, :radiative_forcing, :oc_f, oc_f[i,:])
            update_param!(fair, :radiative_forcing, :oc_aci_f, oc_aci_f[i,:])
            update_param!(fair, :radiative_forcing, :so2_f, so2_f[i,:])
            update_param!(fair, :radiative_forcing, :so2_aci_f, so2_aci_f[i,:])
            update_param!(fair, :radiative_forcing, :solar_f, exogenous_forcing_params[i, :solar])
            update_param!(fair, :radiative_forcing, :landuse_f, exogenous_forcing_params[i, :land_use])
            update_param!(fair, :radiative_forcing, :volcanic_f, exogenous_forcing_params[i, :volcanic])

            # --- Set Parameters Common to Multiple Components ---- #
            update_param!(fair, :shared_co2_pi, co2_p[i, :PI_conc])
            update_param!(fair, :shared_ch4_pi, ch4_p[i, :PI_conc])
            update_param!(fair, :shared_n2o_pi, n2o_p[i, :PI_conc])
            update_param!(fair, :shared_montreal_pi, montreal_pi_conc[:,i])
            update_param!(fair, :shared_flourinated_pi, flourinated_pi_conc[:,i])
            update_param!(fair, :shared_aerosol_plus_pi, aerosol_plus_pi_conc[:,i])

            # Run model.
            run(fair)

            # Store projections.
            temperatures[:,i] = fair[:temperature, :T]  # Global mean surface temperature anomaly (K)
            rf[:, i] = fair[:radiative_forcing, :total_RF] # Total radiative forcing, with individual components scaled by their respective efficacy (Wm⁻²)
            co2[:, i] = fair[:co2_cycle, :co2]  # Total atmospheric carbon dioxide concentrations (ppm)
            ch4[:, i] = fair[:ch4_cycle, :ch4]  # Total atmospheric methane concentrations (ppb)
            n2o[:, i] = fair[:n2o_cycle, :n2o]  # Total atmospheric nitrous oxide concentrations (ppb)

        end

        # Return temperature and radiative forcing 
        return Dict(:temperatures => temperatures, :rf => rf, :co2 => co2, :ch4 => ch4, :n2o => n2o)
    end

    # Return 'fair_monte_carlo' function.
    return fair_monte_carlo_1000
end
