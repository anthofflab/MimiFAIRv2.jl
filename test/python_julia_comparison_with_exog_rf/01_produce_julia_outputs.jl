using DataFrames
using CSVFiles
using Mimi

# This script was written with the MimiFAIRv2.0.0 implementation which uses 
# exogenous radiative forcing in the radiative_forcing.jl component, and compared
# against the Python FAIR v2.0.0:
#   - paper: FaIRv2.0.0: A Generalized Impulse Response Model for Climate Uncertainty and Future Scenario Exploration.
#   - Github: https://github.com/njleach/FAIR/tree/47c6eec031d2edcf09424394dbb86581a1b246ba
#   - paper replication code and other information see paper above, section Code and data availability

# As of this PR: https://github.com/FrankErrickson/MimiFAIRv2.jl/pull/9, 
# Running FAIR with constrained parameters requires individually scaling volcanic, 
# land use, and solar exogenous forcing by an uncertain term. These forcings were
# previously aggregated into a single parameter time series, "exogenous_RF". This
# commit splits them up and also adds an additional parameter, "other_RF", so that 
# a user could provide additional exogenous forcing sources if they wanted to 
# (these default to 0.0 in all periods). As such this code will not run with the current 
# MimiFAIRv2 and will be adapted, but serves as a record of testing done.

include("../../src/MimiFAIRv2.0.jl")

output_dir = joinpath(@__DIR__, "MimiFAIRv2_data")
mkpath(output_dir)

#Load initial conditions (just need parameter names) and extract gas names.
init_gas_vals     = DataFrame(load(joinpath(@__DIR__, "..", "..", "data", "fair_initial_gas_cycle_conditions_1750.csv"), skiplines_begin=7))
montreal_init     = filter(:gas_group => ==("montreal"), init_gas_vals)
flourinated_init  = filter(:gas_group => ==("flourinated"), init_gas_vals)
aerosol_plus_init = filter(:gas_group => ==("aerosol_plus"), init_gas_vals)

 # Sort arrays of initial conditions for multiple gases so they are listed alphabetically.
 sort!(montreal_init,     :gas_name)
 sort!(flourinated_init,  :gas_name)
 sort!(aerosol_plus_init, :gas_name)

 # Set SSP (options = ['ssp119','ssp126','ssp245','ssp370','ssp585'])
 for ssp in ["ssp119", "ssp126", "ssp245", "ssp370", "ssp585"]

    # Load replication emissions and forcing from Python.
    python_emiss   = DataFrame(load(joinpath(@__DIR__, "..", "..", "data", "python_replication", ssp*"_emissions.csv")))
    python_exog_RF = DataFrame(load(joinpath(@__DIR__, "..", "..", "data", "python_replication", ssp*"_forcing.csv")))[:,:External]

    # Extract emissions arrays for multi-gas groupings.
    montreal_emissions     = python_emiss[:, Symbol.(montreal_init.gas_name)]
    flourinated_emissions  = python_emiss[:, Symbol.(flourinated_init.gas_name)]
    aerosol_plus_emissions = python_emiss[:, Symbol.(aerosol_plus_init.gas_name)]

    # Get a model instance for that SSP.
    m = get_model(emissions_forcing_scenario=ssp, start_year=1750, end_year=2100)

    # Update exogenous forcing and emissions data.
    update_param!(m, :exogenous_RF, python_exog_RF)
    update_param!(m, :E_co2, python_emiss.carbon_dioxide)
    update_param!(m, :E_ch4, python_emiss.methane)
    update_param!(m, :E_n2o, python_emiss.nitrous_oxide)
    update_param!(m, :E_flourinated, Array(flourinated_emissions))
    update_param!(m, :E_montreal, Array(montreal_emissions))
    update_param!(m, :E_aerosol_plus, Array(aerosol_plus_emissions))

    run(m)

    # Put CO₂ concentration, total radiative forcing, and global temperature anomaly into a dataframe.
    results = DataFrame(co2 = m[:co2_cycle, :co2], rf = m[:radiative_forcing, :total_RF], temp = m[:temperature, :T])

    # Save results.
    results |> save(joinpath(output_dir, "Mimi_FAIR_replication_"*ssp*".csv"))
end

# Call the R Script 02_julia_python_comparison_plot.R to create comparison plots