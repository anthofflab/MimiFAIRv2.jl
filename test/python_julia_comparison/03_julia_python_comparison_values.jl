module testJuliaPythonComparison

using DataFrames
using CSVFiles
using Test
using MimiFAIRv2
using Mimi

# run model for each possible SSP and compare variable values between julia
# and python versons
for ssp in ["ssp119", "ssp126", "ssp245", "ssp370", "ssp585"]

    # Get a model instance for that SSP.
    m = MimiFAIRv2.get_model(emissions_forcing_scenario=ssp, start_year=1750, end_year=2100)

    # From PR https://github.com/FrankErrickson/MimiFAIRv2.jl/pull/9 the Python 
    # data we have for comparison uses one exogenous forcing file, so we will use
    # that here to match
    update_param!(m, :radiative_forcing, :solar_RF, fill(0., length(1750:2100)))
    update_param!(m, :radiative_forcing, :volcanic_RF, fill(0., length(1750:2100)))
    update_param!(m, :radiative_forcing, :landuse_RF, fill(0., length(1750:2100)))

    python_exog_RF = DataFrame(load(joinpath(@__DIR__, "..", "..", "data", "python_replication", ssp*"_forcing.csv")))[:,:External]
    update_param!(m, :radiative_forcing, :other_RF, python_exog_RF)

    run(m)

    # Put CO₂ concentration, total radiative forcing, and global temperature anomaly into a dataframe.
    julia_results = DataFrame(co2 = m[:co2_cycle, :co2], rf = m[:radiative_forcing, :total_RF], temp = m[:temperature, :T])

    # compare temperatures
    atol = 1e-6
    python_results = load(joinpath(@__DIR__, "..", "..", "data", "python_replication", string(ssp, "_temperature.csv"))) |> DataFrame   
    @test julia_results[:, :temp] ≈ python_results[:, :default] atol = atol

    # df1 = DataFrame(:time => collect(1750:2100), :temp => julia_results[:, :temp], :implementation => "Julia-Mimi")
    # df2 = DataFrame(:time => collect(1750:2100), :temp => python_results[:, :default], :implementation => "Python")
    # vcat(df1, df2) |> @vlplot(:line, x="time:t", y=:temp, color=:implementation, width = 700, height = 500)
    # maximum(abs.(df1.temp .- df2.temp))

    # compare radiative forcing
    atol = 1e-6
    python_results = load(joinpath(@__DIR__, "..", "..", "data", "python_replication", string(ssp, "_forcing.csv"))) |>  DataFrame
    @test julia_results[:, :rf] ≈ python_results[:, :Total] atol = atol

    # df1 = DataFrame(:time => collect(1750:2100), :rf => julia_results[:, :rf], :implementation => "Julia-Mimi")
    # df2 = DataFrame(:time => collect(1750:2100), :rf => python_results[:, :Total], :implementation => "Python")
    # vcat(df1, df2) |> @vlplot(:line, x="time:t", y=:rf, color=:implementation, width = 700, height = 500)
    # maximum(abs.(df1.rf .- df2.rf))

    # compare co2 concentrations
    atol = 1e-5 
    python_results = load(joinpath(@__DIR__, "..", "..", "data", "python_replication", "$(ssp)_concentrations.csv")) |> DataFrame
    @test julia_results[:, :co2] ≈ python_results[:, :carbon_dioxide] atol = atol
    
    # df1 = DataFrame(:time => collect(1750:2100), :co2 => julia_results[:, :co2], :implementation => "Julia-Mimi")
    # df2 = DataFrame(:time => collect(1750:2100), :co2 => python_results[:, :carbon_dioxide], :implementation => "Python")
    # vcat(df1, df2) |> @vlplot(:line, x="time:t", y=:co2, color=:implementation, width = 700, height = 500)
    # maximum(abs.(df1.co2 .- df2.co2))

end

end
