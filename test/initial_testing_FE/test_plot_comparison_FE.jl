using DataFrames
using CSVFiles
using Mimi
using MimiFAIRv2

# include(joinpath(@__DIR__, "..", "src/MimiFAIRv2.jl"))

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
    m = MimiFAIRv2.get_model(emissions_forcing_scenario=ssp, start_year=1750, end_year=2100)

    # Update exogenous forcing and emissions data.

    # LFR CHANGE 2023
    # https://github.com/FrankErrickson/MimiFAIRv2.jl/pull/9 - after these tests were run we
    # split up exogenous forcing to be able to run the constrained parameters so no longer
    # can include this in the testing graphs.

    # update_param!(m, :radiative_forcing, :exogenous_RF, python_exog_RF)
    update_param!(m, :co2_cycle, :E_co2, python_emiss.carbon_dioxide)
    update_param!(m, :ch4_cycle, :E_ch4, python_emiss.methane)
    update_param!(m, :n2o_cycle, :E_n2o, python_emiss.nitrous_oxide)
    update_param!(m, :flourinated_cycles, :E_flourinated, Array(flourinated_emissions))
    update_param!(m, :montreal_cycles, :E_montreal, Array(montreal_emissions))
    update_param!(m, :aerosol_plus_cycles, :E_aerosol_plus, Array(aerosol_plus_emissions))

    run(m)

    # Put CO₂ concentration, total radiative forcing, and global temperature anomaly into a dataframe.
    results = DataFrame(co2 = m[:co2_cycle, :co2], rf = m[:radiative_forcing, :total_RF], temp = m[:temperature, :T])

    # Save results.
    results |> save(joinpath(@__DIR__, "output", "Mimi_FAIR_replication_"*ssp*".csv"))
end
