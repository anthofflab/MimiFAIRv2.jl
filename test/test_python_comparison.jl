module testPythonComparison

using DataFrames
using CSVFiles
using Test
using MimiFAIRv2
using Mimi

using MimiFAIRv2: get_model, create_fair_monte_carlo 

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
    atol = ssp == "ssp370" ? 1e-1 : 1e-6 # ssp370 outputs from python should be double checked
    python_results = load(joinpath(@__DIR__, "..", "data", "python_replication", string(ssp, "_temperature.csv"))) |> DataFrame   
    @test julia_results[1:end-1, :temp] ≈ python_results[1:end-1, :default] atol = atol # leaving off year 2100, pyton data has a jump in this year

    # df1 = DataFrame(:time => collect(1750:2100), :temp => julia_results[:, :temp], :implementation => "Julia-Mimi")
    # df2 = DataFrame(:time => collect(1750:2100), :temp => python_results[:, :default], :implementation => "Python")
    # vcat(df1, df2) |> @vlplot(:line, x="time:t", y=:temp, color=:implementation, width = 700, height = 500)
    # maximum(abs.(df1.temp .- df2.temp))

    # compare radiative forcing
    atol = ssp == "ssp370" ? 5e-1 : 1e-6 # ssp370 outputs from python should be double checked
    python_results = load(joinpath(@__DIR__, "..", "data", "python_replication", string(ssp, "_forcing.csv"))) |>  DataFrame
    @test julia_results[1:end-1, :rf] ≈ python_results[1:end-1, :Total] atol = atol # leaving off year 2100, pyton data has a jump in this year

    # df1 = DataFrame(:time => collect(1750:2100), :rf => julia_results[:, :rf], :implementation => "Julia-Mimi")
    # df2 = DataFrame(:time => collect(1750:2100), :rf => python_results[:, :Total], :implementation => "Python")
    # vcat(df1, df2) |> @vlplot(:line, x="time:t", y=:rf, color=:implementation, width = 700, height = 500)
    # maximum(abs.(df1.rf .- df2.rf))

    # compare co2 concentrations
    atol = ssp == "ssp370" ? 5e-1 : 1e-5 # ssp370 outputs from python should be double checked
    python_results = load(joinpath(@__DIR__, "..", "data", "python_replication", "$(ssp)_concentrations.csv")) |> DataFrame
    @test julia_results[:, :co2] ≈ python_results[:, :carbon_dioxide] atol = atol
    
    # df1 = DataFrame(:time => collect(1750:2100), :co2 => julia_results[:, :co2], :implementation => "Julia-Mimi")
    # df2 = DataFrame(:time => collect(1750:2100), :co2 => python_results[:, :carbon_dioxide], :implementation => "Python")
    # vcat(df1, df2) |> @vlplot(:line, x="time:t", y=:co2, color=:implementation, width = 700, height = 500)
    # maximum(abs.(df1.co2 .- df2.co2))

end

end