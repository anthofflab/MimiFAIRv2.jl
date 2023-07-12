using DataFrames
using CSVFiles
using Mimi
using MimiFAIRv2
using Query
using VegaLite

for ssp in ["ssp119", "ssp126", "ssp245", "ssp370", "ssp585"]

    # Get a model instance for that SSP.
    m = MimiFAIRv2.get_model(emissions_forcing_scenario=ssp, start_year=1750, end_year=2100)
    run(m)

    # Put CO₂ concentration, total radiative forcing, and global temperature anomaly into a dataframe.
    # and Save results.
    mkpath(joinpath(@__DIR__, "Mimi_FAIR_replication"))
    DataFrame(co2 = m[:co2_cycle, :co2], rf = m[:radiative_forcing, :total_RF], temp = m[:temperature, :T]) |>
        save(joinpath(@__DIR__, "Mimi_FAIR_replication", "Mimi_FAIR_replication_"*ssp*".csv"))
end 

# Plotting Temperature
df = DataFrame()
for ssp in ["ssp119", "ssp126", "ssp245", "ssp370", "ssp585"]
    ssp_df = load(joinpath(@__DIR__, "Mimi_FAIR_replication", "Mimi_FAIR_replication_"*ssp*".csv")) |> DataFrame
    insertcols!(ssp_df, :ssp => ssp, :model => :MimiFAIRv2, :time => string.(1750:2100))
    append!(df, ssp_df)
end

df |> @vlplot(
    :line,
    x = {"time:t", title = "Year"},
    y = {:temp, title = "Global Surface Temperature Anomaly (K)"},
    color = {:ssp, title = "SSP"},
    width = 600,
    height = 400,
    title = ["Python vs. Julia-Mimi versions of FaIR v2.0";"RCMIP Emission & Forcing Scenarios (1750-2100)"]
) |> save(joinpath(@__DIR__, "Mimi_FAIR_replication", "PythonReplication_Test_Temperature.png"))

df |> @vlplot(
    :line,
    x = {"time:t", title = "Year"},
    y = {:rf, title = "Total Radiative Forcing (W/m2)"},
    color = {:ssp, title = "SSP"},
    width = 600,
    height = 400,
    title = ["Python vs. Julia-Mimi versions of FaIR v2.0";"RCMIP Emission & Forcing Scenarios (1750-2100)"]
) |> save(joinpath(@__DIR__, "Mimi_FAIR_replication", "PythonReplication_Test_RF.png"))

df |> @vlplot(
    :line,
    x = {"time:t", title = "Year"},
    y = {:co2, title = "CO₂ Concentrations (ppm)"},
    color = {:ssp, title = "SSP"},
    width = 600,
    height = 400,
    title = ["Python vs. Julia-Mimi versions of FaIR v2.0";"RCMIP Emission & Forcing Scenarios (1750-2100)"]
) |> save(joinpath(@__DIR__, "Mimi_FAIR_replication", "PythonReplication_Test_co2.png"))
