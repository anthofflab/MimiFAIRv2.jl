module MimiFAIRv2

# Load required packages.
using CSVFiles, DataFrames, Mimi

# Load helper functions and MimiFAIRv2 model commponent files.
include("helper_functions.jl")
include("components/co2_cycle.jl")
include("components/ch4_cycle.jl")
include("components/n2o_cycle.jl")
include("components/montreal_gas_cycles.jl")
include("components/flourinated_gas_cycles.jl")
include("components/aerosol_plus_gas_cycles.jl")
include("components/radiative_forcing.jl")
include("components/temperature.jl")

# load monte Carlo
include("monte_carlo.jl")
include("monte_carlo_1000.jl")
include("monte_carlo_10k.jl")

"""
    get_model(;emissions_forcing_scenario::String="ssp585", start_year::Int=1750, 
                end_year::Int=2500, TCR::Float64=1.79, RWF::Float64=0.552, 
                F2x::Float64=3.759)

Return a constructed model with default settings for emissions forcing scenario 
(default to ssp585), start year (default to 1750), end year (default to 2500),
transient climate response (default to 1.79), realized warming fraction (default 
to 5.22), and forcing from a doubling of CO₂ (default to 3.759).
"""
function get_model(;emissions_forcing_scenario::String="ssp585", start_year::Int=1750, end_year::Int=2500, TCR::Float64=1.79, RWF::Float64=0.552, F2x::Float64=3.759)

    # TODO turning this off for now so it doesn't drive us nuts with 10k runs while we're generalizing this function
    # if start_year !== 1750
    #     error("FAIRv2 model monte carlo simulation should not be set to start with a year differing from 1750 as initial conditions are not calibrated for a different start year!")
    # end 

 	# ---------------------------------------------
	# ---------------------------------------------
	# Set Up Data and Parameter Values
	# ---------------------------------------------
 	# ---------------------------------------------

	# Load emissions and forcing data and crop to appropriate model years (scenarios span 1750-2500 by default).
	scenario_indices = indexin(start_year:end_year, 1750:2500)
	forcing_data     = DataFrame(load(joinpath(@__DIR__, "..", "data", "rcmip_"*emissions_forcing_scenario*"_effective_radiative_forcing_1750_to_2500.csv"), skiplines_begin=6))[scenario_indices,:]
	emissions_data   = DataFrame(load(joinpath(@__DIR__, "..", "data", "rcmip_"*emissions_forcing_scenario*"_emissions_1750_to_2500.csv"), skiplines_begin=6))[scenario_indices,:]

	# Load FAIR default gas cycle (gas) and indirect radiative forcing (irf_p) parameters.
	gas_p = DataFrame(load(joinpath(@__DIR__, "..", "data", "default_gas_cycle_parameters.csv"), skiplines_begin=6))
	irf_p = DataFrame(load(joinpath(@__DIR__, "..", "data", "default_indirect_radiative_forcing_parameters.csv"), skiplines_begin=7))

	# Load initial conditions for 1750 under the RCMIP emissions & forcing scenario.
	init_gas_vals     = DataFrame(load(joinpath(@__DIR__, "..", "data", "fair_initial_gas_cycle_conditions_1750.csv"), skiplines_begin=7))
	init_thermal_vals = DataFrame(load(joinpath(@__DIR__, "..", "data", "fair_initial_thermal_conditions_1750.csv"), skiplines_begin=7))

	# Isolate specific gas parameters for convenience.
	co2_p              = filter(:gas_name  => ==("carbon_dioxide"), gas_p)
	ch4_p              = filter(:gas_name  => ==("methane"), gas_p)
	n2o_p              = filter(:gas_name  => ==("nitrous_oxide"), gas_p)
	montreal_gas_p     = filter(:gas_group => ==("montreal"), gas_p)
	flourinated_gas_p  = filter(:gas_group => ==("flourinated"), gas_p)
	aerosol_plus_gas_p = filter(:gas_group => ==("aerosol_plus"), gas_p)

	# Sort arrays of parameters for multiple gases so they are listed alphabetically.
	sort!(flourinated_gas_p,  :gas_name)
	sort!(montreal_gas_p,     :gas_name)
	sort!(aerosol_plus_gas_p, :gas_name)

	#Isolate Montreal indirect forcing parameters and sort alphabetically.
	montreal_irf_p = filter(:gas_group => ==("montreal"), irf_p)
	sort!(montreal_irf_p, :gas_name)

	# Isolate initial model conditions (initial conditions currently default to year 1750).
	co2_init          = filter(:gas_name  => ==("carbon_dioxide"), init_gas_vals)
	ch4_init          = filter(:gas_name  => ==("methane"), init_gas_vals)
	n2o_init          = filter(:gas_name  => ==("nitrous_oxide"), init_gas_vals)
	montreal_init     = filter(:gas_group => ==("montreal"), init_gas_vals)
	flourinated_init  = filter(:gas_group => ==("flourinated"), init_gas_vals)
	aerosol_plus_init = filter(:gas_group => ==("aerosol_plus"), init_gas_vals)

	# Sort arrays of initial conditions for multiple gases so they are listed alphabetically.
	sort!(montreal_init,     :gas_name)
	sort!(flourinated_init,  :gas_name)
	sort!(aerosol_plus_init, :gas_name)

	# Extract emissions arrays for multi-gas groupings.
	montreal_emissions     = emissions_data[:, Symbol.(montreal_init.gas_name)]
	flourinated_emissions  = emissions_data[:, Symbol.(flourinated_init.gas_name)]
	aerosol_plus_emissions = emissions_data[:, Symbol.(aerosol_plus_init.gas_name)]

    # Create helper arrays of indices to use for pulling parameter groups.
    a_idxs = [:a1,:a2,:a3,:a4]
    τ_idxs = [:tau1,:tau2,:tau3,:tau4]

	# Create arrays for 'a' parameter groups. 
	co2_a          = vec(Array(co2_p[:, a_idxs]))
	ch4_a          = vec(Array(ch4_p[:, a_idxs]))
	n2o_a          = vec(Array(n2o_p[:, a_idxs]))
	montreal_a     = Array(montreal_gas_p[:, a_idxs])
	flourinated_a  = Array(flourinated_gas_p[:, a_idxs])
	aerosol_plus_a = Array(aerosol_plus_gas_p[:, a_idxs])

	# Create arrays for 'τ' parameter groups.
	co2_τ 		   = vec(Array(co2_p[:, τ_idxs]))
	ch4_τ 		   = vec(Array(ch4_p[:, τ_idxs]))
	n2o_τ 		   = vec(Array(n2o_p[:, τ_idxs]))
	montreal_τ 	   = Array(montreal_gas_p[:, τ_idxs])
	flourinated_τ  = Array(flourinated_gas_p[:, τ_idxs])
	aerosol_plus_τ = Array(aerosol_plus_gas_p[:, τ_idxs])

	# Calculate constants to approximate numerical solution for state-dependent timescale adjustment factor from Millar et al. (2017).
	g0_co2, g1_co2 					 =  calculate_g0_g1(co2_a, co2_τ)
	g0_ch4, g1_ch4 				     =  calculate_g0_g1(ch4_a, ch4_τ)
	g0_n2o, g1_n2o 					 =  calculate_g0_g1(n2o_a, n2o_τ)
	g0_montreal, g1_montreal 		 =  calculate_g0_g1(montreal_a, montreal_τ)
	g0_flourinated, g1_flourinated   =  calculate_g0_g1(flourinated_a, flourinated_τ)
	g0_aerosol_plus, g1_aerosol_plus =  calculate_g0_g1(aerosol_plus_a, aerosol_plus_τ)

	# Calculate default thermal parameter values 
    # Defaults pulled from get_model function defaults are are as follows:
    #   transient climate response = 1.79
    #   realized warming fraction = 0.552
    #   forcing from a doubling of CO₂ = 3.759
	thermal_p = get_thermal_parameter_defaults(TCR = TCR, RWF= RWF, F2x = F2x)

	# Calculate thermal decay factors, defined as exp(-1/d).
	thermal_decay_factors = exp.(-1.0 ./ vec(Array(thermal_p[:,[:d1,:d2,:d3]])))


 	# ---------------------------------------------
 	# ---------------------------------------------
    # Initialize Mimi model.
    # ---------------------------------------------
 	# ---------------------------------------------

 	# Create a Mimi model.
    m = Model()

    # Set time and gas-grouping indices.
    set_dimension!(m, :time, start_year:end_year)
    set_dimension!(m, :montreal_gases, montreal_gas_p[:,:gas_name])
    set_dimension!(m, :flourinated_gases, flourinated_gas_p[:,:gas_name])
    set_dimension!(m, :aerosol_plus_gases, aerosol_plus_gas_p[:,:gas_name])

    # ---------------------------------------------
    # Add components to model
    # ---------------------------------------------

    add_comp!(m, co2_cycle)
    add_comp!(m, ch4_cycle)
    add_comp!(m, n2o_cycle)
    add_comp!(m, montreal_cycles)
    add_comp!(m, flourinated_cycles)
    add_comp!(m, aerosol_plus_cycles)
    add_comp!(m, radiative_forcing)
    add_comp!(m, temperature)

	# ---------------------------------------------
    # Set component-specific parameters
    # ---------------------------------------------

 	# ---- Carbon Cycle ---- #
	update_param!(m, :co2_cycle, :co2_0, co2_init.concentration[1])
	update_param!(m, :co2_cycle, :r0_co2, co2_p.r0[1])
	update_param!(m, :co2_cycle, :rU_co2, co2_p.rC[1])
	update_param!(m, :co2_cycle, :rT_co2, co2_p.rT[1])
	update_param!(m, :co2_cycle, :rA_co2, co2_p.rA[1])
	update_param!(m, :co2_cycle, :GU_co2_0, co2_init.cumulative_uptake[1])
	update_param!(m, :co2_cycle, :g0_co2, g0_co2)
	update_param!(m, :co2_cycle, :g1_co2, g1_co2)
	update_param!(m, :co2_cycle, :R0_co2, vec(Array(co2_init[:,[:pool1,:pool2,:pool3,:pool4]])))
	update_param!(m, :co2_cycle, :emiss2conc_co2, co2_p.emis2conc[1])
	update_param!(m, :co2_cycle, :a_co2, co2_a)
	update_param!(m, :co2_cycle, :τ_co2, co2_τ)
	update_param!(m, :co2_cycle, :E_co2, emissions_data.carbon_dioxide)

	# ---- Methane Cycle ---- #
	update_param!(m, :ch4_cycle, :ch4_0, ch4_init.concentration[1])
	update_param!(m, :ch4_cycle, :r0_ch4, ch4_p.r0[1])
	update_param!(m, :ch4_cycle, :rU_ch4, ch4_p.rC[1])
	update_param!(m, :ch4_cycle, :rT_ch4, ch4_p.rT[1])
	update_param!(m, :ch4_cycle, :rA_ch4, ch4_p.rA[1])
	update_param!(m, :ch4_cycle, :GU_ch4_0, ch4_init.cumulative_uptake[1])
	update_param!(m, :ch4_cycle, :g0_ch4, g0_ch4)
	update_param!(m, :ch4_cycle, :g1_ch4, g1_ch4)
	update_param!(m, :ch4_cycle, :R0_ch4, vec(Array(ch4_init[:,[:pool1,:pool2,:pool3,:pool4]])))
	update_param!(m, :ch4_cycle, :emiss2conc_ch4, ch4_p.emis2conc[1])
	update_param!(m, :ch4_cycle, :a_ch4, ch4_a)
	update_param!(m, :ch4_cycle, :τ_ch4, ch4_τ)
	update_param!(m, :ch4_cycle, :E_ch4, emissions_data.methane)

	# ---- Nitrous Oxide Cycle ---- #
	update_param!(m, :n2o_cycle, :n2o_0, n2o_init.concentration[1])
	update_param!(m, :n2o_cycle, :r0_n2o, n2o_p.r0[1])
	update_param!(m, :n2o_cycle, :rU_n2o, n2o_p.rC[1])
	update_param!(m, :n2o_cycle, :rT_n2o, n2o_p.rT[1])
	update_param!(m, :n2o_cycle, :rA_n2o, n2o_p.rA[1])
	update_param!(m, :n2o_cycle, :GU_n2o_0, n2o_init.cumulative_uptake[1])
	update_param!(m, :n2o_cycle, :g0_n2o, g0_n2o)
	update_param!(m, :n2o_cycle, :g1_n2o, g1_n2o)
	update_param!(m, :n2o_cycle, :R0_n2o, vec(Array(n2o_init[:,[:pool1,:pool2,:pool3,:pool4]])))
	update_param!(m, :n2o_cycle, :emiss2conc_n2o, n2o_p.emis2conc[1])
	update_param!(m, :n2o_cycle, :a_n2o, n2o_a)
	update_param!(m, :n2o_cycle, :τ_n2o, n2o_τ)
	update_param!(m, :n2o_cycle, :E_n2o, emissions_data.nitrous_oxide)

	# ---- Flourinated Gas Cycles ---- #
	update_param!(m, :flourinated_cycles, :flourinated_0, flourinated_init[:, :concentration])
	update_param!(m, :flourinated_cycles, :r0_flourinated, flourinated_gas_p[:,:r0])
	update_param!(m, :flourinated_cycles, :rU_flourinated, flourinated_gas_p[:,:rC])
	update_param!(m, :flourinated_cycles, :rT_flourinated, flourinated_gas_p[:,:rT])
	update_param!(m, :flourinated_cycles, :rA_flourinated, flourinated_gas_p[:,:rA])
	update_param!(m, :flourinated_cycles, :GU_flourinated_0, flourinated_init[:, :cumulative_uptake])
	update_param!(m, :flourinated_cycles, :g0_flourinated, g0_flourinated)
	update_param!(m, :flourinated_cycles, :g1_flourinated, g1_flourinated)
	update_param!(m, :flourinated_cycles, :R0_flourinated, Array(flourinated_init[:,[:pool1,:pool2,:pool3,:pool4]]))
	update_param!(m, :flourinated_cycles, :emiss2conc_flourinated, flourinated_gas_p[:,:emis2conc])
	update_param!(m, :flourinated_cycles, :a_flourinated, flourinated_a)
	update_param!(m, :flourinated_cycles, :τ_flourinated, flourinated_τ)
	update_param!(m, :flourinated_cycles, :E_flourinated, Array(flourinated_emissions))

	# ---- Montreal Protocol Gas Cycles ---- #
	update_param!(m, :montreal_cycles, :montreal_0, montreal_init[:, :concentration])
	update_param!(m, :montreal_cycles, :r0_montreal, montreal_gas_p[:,:r0])
	update_param!(m, :montreal_cycles, :rU_montreal, montreal_gas_p[:,:rC])
	update_param!(m, :montreal_cycles, :rT_montreal, montreal_gas_p[:,:rT])
	update_param!(m, :montreal_cycles, :rA_montreal, montreal_gas_p[:,:rA])
	update_param!(m, :montreal_cycles, :GU_montreal_0, montreal_init[:, :cumulative_uptake])
	update_param!(m, :montreal_cycles, :g0_montreal, g0_montreal)
	update_param!(m, :montreal_cycles, :g1_montreal, g1_montreal)
	update_param!(m, :montreal_cycles, :R0_montreal, Array(montreal_init[:,[:pool1,:pool2,:pool3,:pool4]]))
	update_param!(m, :montreal_cycles, :emiss2conc_montreal, montreal_gas_p[:,:emis2conc])
	update_param!(m, :montreal_cycles, :a_montreal, montreal_a)
	update_param!(m, :montreal_cycles, :τ_montreal, montreal_τ)
	update_param!(m, :montreal_cycles, :E_montreal, Array(montreal_emissions))

	# ---- Tropospheric Ozone Precursors, Aerosols, & Reactive Gas Cycles (Aerosol+) ---- #
	update_param!(m, :aerosol_plus_cycles, :aerosol_plus_0, aerosol_plus_init[:, :concentration])
	update_param!(m, :aerosol_plus_cycles, :r0_aerosol_plus, aerosol_plus_gas_p[:,:r0])
	update_param!(m, :aerosol_plus_cycles, :rU_aerosol_plus, aerosol_plus_gas_p[:,:rC])
	update_param!(m, :aerosol_plus_cycles, :rT_aerosol_plus, aerosol_plus_gas_p[:,:rT])
	update_param!(m, :aerosol_plus_cycles, :rA_aerosol_plus, aerosol_plus_gas_p[:,:rA])
	update_param!(m, :aerosol_plus_cycles, :GU_aerosol_plus_0, aerosol_plus_init[:, :cumulative_uptake])
	update_param!(m, :aerosol_plus_cycles, :g0_aerosol_plus, g0_aerosol_plus)
	update_param!(m, :aerosol_plus_cycles, :g1_aerosol_plus, g1_aerosol_plus)
	update_param!(m, :aerosol_plus_cycles, :R0_aerosol_plus, Array(aerosol_plus_init[:,[:pool1,:pool2,:pool3,:pool4]]))
	update_param!(m, :aerosol_plus_cycles, :emiss2conc_aerosol_plus, aerosol_plus_gas_p[:,:emis2conc])
	update_param!(m, :aerosol_plus_cycles, :a_aerosol_plus, aerosol_plus_a)
	update_param!(m, :aerosol_plus_cycles, :τ_aerosol_plus, aerosol_plus_τ)
	update_param!(m, :aerosol_plus_cycles, :E_aerosol_plus, Array(aerosol_plus_emissions))

	# ---- Radiative Forcing ---- #
	update_param!(m, :radiative_forcing, :co2_f, vec(Array(co2_p[:, [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :ch4_f, vec(Array(ch4_p[:, [:f1, :f2, :f3]])))
   	update_param!(m, :radiative_forcing, :ch4_o3_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="methane|o3"), [:f1, :f2, :f3]])) )
   	update_param!(m, :radiative_forcing, :ch4_h2o_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="methane|strat_h2o"), [:f1, :f2, :f3]])) )
  	update_param!(m, :radiative_forcing, :n2o_f, vec(Array(n2o_p[:, [:f1, :f2, :f3]])))
   	update_param!(m, :radiative_forcing, :n2o_o3_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="nitrous_oxide|o3"), [:f1, :f2, :f3]])) )
	update_param!(m, :radiative_forcing, :montreal_f, Array(montreal_gas_p[:,[:f1, :f2, :f3]]))
	update_param!(m, :radiative_forcing, :montreal_ind_f, Array(montreal_irf_p[:,[:f1, :f2, :f3]]))
	update_param!(m, :radiative_forcing, :flourinated_f, Array(flourinated_gas_p[:,[:f1, :f2, :f3]]))
	update_param!(m, :radiative_forcing, :bc_f, vec(Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name.=="bc"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :bc_snow_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="bc|bc_on_snow"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :bc_aci_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="bc|aci"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :co_f, vec(Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name.=="co"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :co_o3_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="co|o3"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :nh3_f, vec(Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name.=="nh3"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :nmvoc_f, vec(Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name.=="nmvoc"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :nmvoc_o3_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="nmvoc|o3"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :nox_f, vec(Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name.=="nox"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :nox_o3_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="nox|o3"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :nox_avi_f, vec(Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name.=="nox_avi"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :nox_avi_contrails_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="nox_avi|contrails"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :oc_f, vec(Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name.=="oc"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :oc_aci_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="oc|aci"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :so2_f, vec(Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name.=="so2"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :so2_aci_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="so2|aci"), [:f1, :f2, :f3]])))
  	update_param!(m, :radiative_forcing, :solar_f, 1.0)
  	update_param!(m, :radiative_forcing, :landuse_f, 1.0)
  	update_param!(m, :radiative_forcing, :volcanic_f, 1.0)
  	update_param!(m, :radiative_forcing, :solar_RF, forcing_data[:, :solar])
  	update_param!(m, :radiative_forcing, :landuse_RF, forcing_data[:, :land_use])
  	update_param!(m, :radiative_forcing, :volcanic_RF, forcing_data[:, :volcanic])
  	update_param!(m, :radiative_forcing, :other_RF, zeros(length(start_year:end_year)))

    # ---- Global Temperature Anomaly ---- #
    update_param!(m, :temperature, :Tj_0, vec(Array(init_thermal_vals[:,[:thermal_box_1,:thermal_box_2,:thermal_box_3]])))
    update_param!(m, :temperature, :T_0, init_thermal_vals.global_temp_anomaly[1])
    update_param!(m, :temperature, :q, vec(Array(thermal_p[:,[:q1,:q2,:q3]])))
    update_param!(m, :temperature, :decay_factor, thermal_decay_factors)

    # --- Set Parameters Common to Multiple Components ---- #
	add_shared_param!(m, :shared_co2_pi, co2_p.PI_conc[1])
    connect_param!(m, :co2_cycle, :co2_pi, :shared_co2_pi)
    connect_param!(m, :radiative_forcing, :co2_pi, :shared_co2_pi)

	add_shared_param!(m, :shared_ch4_pi, ch4_p.PI_conc[1])
    connect_param!(m, :ch4_cycle, :ch4_pi, :shared_ch4_pi)
    connect_param!(m, :radiative_forcing, :ch4_pi, :shared_ch4_pi)

	add_shared_param!(m, :shared_n2o_pi, n2o_p.PI_conc[1])
    connect_param!(m, :n2o_cycle, :n2o_pi, :shared_n2o_pi)
    connect_param!(m, :radiative_forcing, :n2o_pi, :shared_n2o_pi)

	add_shared_param!(m, :shared_montreal_pi, montreal_gas_p[:,:PI_conc], dims = [:montreal_gases])
    connect_param!(m, :montreal_cycles, :montreal_pi, :shared_montreal_pi)
    connect_param!(m, :radiative_forcing, :montreal_pi, :shared_montreal_pi)

	add_shared_param!(m, :shared_flourinated_pi, flourinated_gas_p[:,:PI_conc], dims = [:flourinated_gases])
    connect_param!(m, :flourinated_cycles, :flourinated_pi, :shared_flourinated_pi)
    connect_param!(m, :radiative_forcing, :flourinated_pi, :shared_flourinated_pi)

	add_shared_param!(m, :shared_aerosol_plus_pi, aerosol_plus_gas_p[:,:PI_conc], dims = [:aerosol_plus_gases])
    connect_param!(m, :aerosol_plus_cycles, :aerosol_plus_pi, :shared_aerosol_plus_pi)
    connect_param!(m, :radiative_forcing, :aerosol_plus_pi, :shared_aerosol_plus_pi)

 	# ---------------------------------------------
    # Create Connections Between Mimi Components
    # ---------------------------------------------

    # Syntax is :component_needing_a_parameter_input => :name_of_that_parameter, :component_calculating_required_values => :name_of_variable_output
    connect_param!(m, :montreal_cycles     => :Tj,                :temperature         => :Tj)
    connect_param!(m, :flourinated_cycles  => :Tj,                :temperature         => :Tj)
    connect_param!(m, :aerosol_plus_cycles => :Tj,                :temperature         => :Tj)
    connect_param!(m, :co2_cycle           => :Tj,                :temperature         => :Tj)
    connect_param!(m, :ch4_cycle           => :Tj,                :temperature         => :Tj)
    connect_param!(m, :n2o_cycle           => :Tj,                :temperature         => :Tj)
    connect_param!(m, :radiative_forcing   => :co2_conc,          :co2_cycle     	   => :co2)
    connect_param!(m, :radiative_forcing   => :ch4_conc,          :ch4_cycle     	   => :ch4)
    connect_param!(m, :radiative_forcing   => :n2o_conc,          :n2o_cycle     	   => :n2o)
    connect_param!(m, :radiative_forcing   => :montreal_conc,     :montreal_cycles     => :montreal_conc)
    connect_param!(m, :radiative_forcing   => :flourinated_conc,  :flourinated_cycles  => :flourinated_conc)
    connect_param!(m, :radiative_forcing   => :aerosol_plus_conc, :aerosol_plus_cycles => :aerosol_plus_conc)
    connect_param!(m, :temperature   	   => :F,                 :radiative_forcing   => :total_RF)

    # Return model.
    return m
end

end # module
