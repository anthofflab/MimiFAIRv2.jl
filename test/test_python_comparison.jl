module testPythonComparison

using DataFrames
using CSVFiles
using Test
using MimiFAIRv2
using Mimi

using MimiFAIRv2: get_model, create_fair_monte_carlo # load `get_model` function to avoid need for `MimiFAIRv2.` prefix

# run model for each possible SSP and compare variable values between julia
# and python versons

#Load initial conditions (just need parameter names) and extract gas names.
init_gas_vals     = DataFrame(load(joinpath(@__DIR__, "..", "data", "fair_initial_gas_cycle_conditions_1750.csv"), skiplines_begin=7))
montreal_init     = filter(:gas_group => ==("montreal"), init_gas_vals)
flourinated_init  = filter(:gas_group => ==("flourinated"), init_gas_vals)
aerosol_plus_init = filter(:gas_group => ==("aerosol_plus"), init_gas_vals)

# Sort arrays of initial conditions for multiple gases so they are listed alphabetically.
sort!(montreal_init,     :gas_name)
sort!(flourinated_init,  :gas_name)
sort!(aerosol_plus_init, :gas_name)

for ssp in ["ssp119", "ssp126", "ssp245", "ssp370", "ssp585"]

    # Load replication emissions and forcing from Python.
    python_emiss   = load(joinpath(@__DIR__, "..", "data", "python_replication", ssp*"_emissions.csv")) |> DataFrame

    # Extract emissions arrays for multi-gas groupings.
    montreal_emissions     = python_emiss[:, Symbol.(montreal_init.gas_name)]
    flourinated_emissions  = python_emiss[:, Symbol.(flourinated_init.gas_name)]
    aerosol_plus_emissions = python_emiss[:, Symbol.(aerosol_plus_init.gas_name)]

    # Get a model instance for that SSP.
    m = MimiFAIRv2.get_model(emissions_forcing_scenario=ssp, start_year=1750, end_year=2100)

    # Update dmissions data.
    update_param!(m, :co2_cycle, :E_co2, python_emiss.carbon_dioxide)
    update_param!(m, :ch4_cycle, :E_ch4, python_emiss.methane)
    update_param!(m, :n2o_cycle, :E_n2o, python_emiss.nitrous_oxide)
    update_param!(m, :flourinated_cycles, :E_flourinated, Array(flourinated_emissions))
    update_param!(m, :montreal_cycles, :E_montreal, Array(montreal_emissions))
    update_param!(m, :aerosol_plus_cycles, :E_aerosol_plus, Array(aerosol_plus_emissions))
    
    run(m)

    # Put CO₂ concentration, total radiative forcing, and global temperature anomaly into a dataframe.
    julia_results = DataFrame(co2 = m[:co2_cycle, :co2], rf = m[:radiative_forcing, :total_RF], temp = m[:temperature, :T])

    # compare temperatures
    python_results = load(joinpath(@__DIR__, "..", "data", "python_replication", string(ssp, "_temperature.csv"))) |> DataFrame
    @test maximum(abs, julia_results[!, :temp] .- python_results[!, :default]) ≈ 0. atol = 1e-3
    
    # compare radiative forcing
    python_results = load(joinpath(@__DIR__, "..", "data", "python_replication", string(ssp, "_forcing.csv"))) |> DataFrame
    select!(python_results, "Forcing component", "Total")
    @test maximum(abs, julia_results[!, :rf] .- python_results[!, :Total]) ≈ 0. atol = 1e-3

    # compare emissions
    conv = select!(DataFrame(load(joinpath(@__DIR__, "..", "data", "default_gas_cycle_parameters.csv"), skiplines_begin=6)), :gas_name, :emis2conc)
    filter(:gas_name  => ==("carbon_dioxide"), conv)

    python_results = load(joinpath(@__DIR__, "..", "data", "python_replication", string(ssp, "_concentrations.csv"))) |> DataFrame
    select!(python_results, "Gas name", "carbon_dioxide")
    @test maximum(abs, julia_results[!, :co2] .- python_results[!, :carbon_dioxide]) ≈ 0. atol = 1e-3

end

end