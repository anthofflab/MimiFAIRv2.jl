using DataFrames
using CSVFiles
using Mimi
using MimiFAIRv2
using Query
using VegaLite

outdir = joinpath(@__DIR__, "output/test_plot_comparison")
mkpath(outdir)

for ssp in ["ssp119", "ssp126", "ssp245", "ssp370", "ssp585"]

    # Get a model instance for that SSP.
    m = MimiFAIRv2.get_model(emissions_forcing_scenario=ssp, start_year=1750, end_year=2100)
    run(m)

    # Put CO₂ concentration, total radiative forcing, and global temperature anomaly into a dataframe.
    # and Save results.
    mkpath(joinpath(outdir))
    DataFrame(co2 = m[:co2_cycle, :co2], rf_total = m[:radiative_forcing, :total_RF], rf_co2 = m[:radiative_forcing, :co2_RF], temp = m[:temperature, :T]) |>
        save(joinpath(outdir, "Mimi_FAIR_replication_"*ssp*".csv"))
end 

# Plotting Temperature
df = DataFrame()
for ssp in ["ssp119", "ssp126", "ssp245", "ssp370", "ssp585"]

    # Add MimiFAIRv2 replication variables
    ssp_df = load(joinpath(outdir, "Mimi_FAIR_replication_"*ssp*".csv")) |> DataFrame
    insertcols!(ssp_df, :ssp => ssp, :implementation => "FAIRv2.0 (Julia-Mimi version)", :time => string.(1750:2100))
    append!(df, ssp_df)

    # add Python replication variables
    temp = load(joinpath(@__DIR__, "..", "data", "python_replication", "$(ssp)_temperature.csv")) |> DataFrame
    rf = load(joinpath(@__DIR__, "..", "data", "python_replication", "$(ssp)_forcing.csv")) |> DataFrame
    co2 = load(joinpath(@__DIR__, "..", "data", "python_replication", "$(ssp)_concentrations.csv")) |> DataFrame
    append!(df, DataFrame(:co2 => co2.carbon_dioxide, :rf_total => rf.Total, :rf_co2 => rf.carbon_dioxide, :temp => temp.default, :ssp => ssp, :implementation => "FAIRv2.0 (original Python version)", :time => string.(1750:2100))) 
end

# TEMPERATURE ANOMALY 

df |> @vlplot(
        x = {"time:t", title = "Year"},
        y = {:temp, title = "Global Surface Temperature Anomaly (K)"},
        color = {:ssp, title = "SSP", legend = {symbolOpacity = 1.}},
        width = 700,
        height = 500,
        title = ["Python vs. Julia-Mimi versions of FaIR v2.0";"RCMIP Emission & Forcing Scenarios (1750-2100)"; "points Julia-Mimi; solid line original Python"]
    ) +
    @vlplot(
        data = df |> @filter(_.implementation == "FAIRv2.0 (original Python version)"),
        :line
    )  +
    @vlplot(
        data = df |> @filter(_.implementation == "FAIRv2.0 (Julia-Mimi version)" && parse(Int64, _.time) in collect(1750:5:2100)),
        {:point, filled = true, stroke = "black", strokeWidth = 1.5}
    ) |> save(joinpath(outdir, "PythonReplication_Test_Temperature.svg"))

df |>
    @vlplot(
        x = {"time:t", title = "Year"},
        y = {:temp, title = "Global Surface Temperature Anomaly (K)"},
        color = {:ssp, title = "SSP", legend = {symbolOpacity = 1.}},
        width = 700,
        height = 500,
        title = ["Python vs. Julia-Mimi versions of FaIR v2.0";"RCMIP Emission & Forcing Scenarios (1750-2100)"; "points Julia-Mimi; solid line original Python"]
    ) +
    @vlplot(
        data = df |> @filter(_.implementation == "FAIRv2.0 (original Python version)" && parse(Int64, _.time) >= 2000),
        :line
    )  +
    @vlplot(
        data = df |> @filter(_.implementation == "FAIRv2.0 (Julia-Mimi version)" && parse(Int64, _.time) in collect(2000:2:2100)),
        {:point, filled = true, stroke = "black", strokeWidth = 1.5}
    ) |> save(joinpath(outdir, "PythonReplication_Test_Temperature_ZoomTime.svg"))

# TOTAL RADIATVE FORCING 

df |> @vlplot(
        x = {"time:t", title = "Year"},
        y = {:rf_total, title = "Total Radiative Forcing (W/m2)"},
        color = {:ssp, title = "SSP", legend = {symbolOpacity = 1.}},
        width = 700,
        height = 500,
        title = ["Python vs. Julia-Mimi versions of FaIR v2.0";"RCMIP Emission & Forcing Scenarios (1750-2100)"; "points Julia-Mimi; solid line original Python"]
    ) +
    @vlplot(
        data = df |> @filter(_.implementation == "FAIRv2.0 (original Python version)"),
        :line
    )  +
    @vlplot(
        data = df |> @filter(_.implementation == "FAIRv2.0 (Julia-Mimi version)" && parse(Int64, _.time) in collect(1750:5:2100)),
        {:point, filled = true, stroke = "black", strokeWidth = 1.5}
    ) |> save(joinpath(outdir, "PythonReplication_Test_Total_RF.svg"))

df |> @vlplot(
        x = {"time:t", title = "Year"},
        y = {:rf_total, title = "Total Radiative Forcing (W/m2)"},
        color = {:ssp, title = "SSP", legend = {symbolOpacity = 1.}},
        width = 700,
        height = 500,
        title = ["Python vs. Julia-Mimi versions of FaIR v2.0";"RCMIP Emission & Forcing Scenarios (1750-2100)"; "points Julia-Mimi; solid line original Python"]
    ) +
    @vlplot(
        data = df |> @filter(_.implementation == "FAIRv2.0 (original Python version)" && parse(Int64, _.time) >= 2000),
        :line
    )  +
    @vlplot(
        data = df |> @filter(_.implementation == "FAIRv2.0 (Julia-Mimi version)" && parse(Int64, _.time) in collect(2000:2:2100)),
        {:point, filled = true, stroke = "black", strokeWidth = 1.5}
    ) |> save(joinpath(outdir, "PythonReplication_Test_Total_RF_ZoomTime.svg"))

# CO2 RADIATIVE FORCING

df |> @vlplot(
        x = {"time:t", title = "Year"},
        y = {:rf_co2, title = "CO₂ Radiative Forcing (W/m2)"},
        color = {:ssp, title = "SSP", legend = {symbolOpacity = 1.}},
        width = 700,
        height = 500,
        title = ["Python vs. Julia-Mimi versions of FaIR v2.0";"RCMIP Emission & Forcing Scenarios (1750-2100)"; "points Julia-Mimi; solid line original Python"]
    ) +
    @vlplot(
        data = df |> @filter(_.implementation == "FAIRv2.0 (original Python version)"),
        :line
    )  +
    @vlplot(
        data = df |> @filter(_.implementation == "FAIRv2.0 (Julia-Mimi version)" && parse(Int64, _.time) in collect(1750:5:2100)),
        {:point, filled = true, stroke = "black", strokeWidth = 1.5}
    ) |> save(joinpath(outdir, "PythonReplication_Test_CO2_RF.svg"))

df |> @vlplot(
        x = {"time:t", title = "Year"},
        y = {:rf_co2, title = "CO₂ Radiative Forcing (W/m2)"},
        color = {:ssp, title = "SSP", legend = {symbolOpacity = 1.}},
        width = 700,
        height = 500,
        title = ["Python vs. Julia-Mimi versions of FaIR v2.0";"RCMIP Emission & Forcing Scenarios (1750-2100)"; "points Julia-Mimi; solid line original Python"]
    ) +
    @vlplot(
        data = df |> @filter(_.implementation == "FAIRv2.0 (original Python version)" && parse(Int64, _.time) >= 2000),
        :line
    )  +
    @vlplot(
        data = df |> @filter(_.implementation == "FAIRv2.0 (Julia-Mimi version)" && parse(Int64, _.time) in collect(2000:2:2100)),
        {:point, filled = true, stroke = "black", strokeWidth = 1.5}
    ) |> save(joinpath(outdir, "PythonReplication_Test_CO2_RF_ZoomTime.svg"))

# CO2 CONCENTRATIONS

df |> @vlplot(
        x = {"time:t", title = "Year"},
        y = {:co2, title = "CO₂ Concentrations (ppm)"},
        color = {:ssp, title = "SSP", legend = {symbolOpacity = 1.}},
        width = 700,
        height = 500,
        title = ["Python vs. Julia-Mimi versions of FaIR v2.0";"RCMIP Emission & Forcing Scenarios (1750-2100)"; "points Julia-Mimi; solid line original Python"]
    ) +
    @vlplot(
        data = df |> @filter(_.implementation == "FAIRv2.0 (original Python version)"),
        :line
    )  +
    @vlplot(
        data = df |> @filter(_.implementation == "FAIRv2.0 (Julia-Mimi version)" && parse(Int64, _.time) in collect(1750:5:2100)),
        {:point, filled = true, stroke = "black", strokeWidth = 1.5}
    ) |> save(joinpath(outdir, "PythonReplication_Test_CO2_Concentrations_ZoomTime.svg"))


df |> @vlplot(
        x = {"time:t", title = "Year"},
        y = {:co2, title = "CO₂ Concentrations (ppm)"},
        color = {:ssp, title = "SSP", legend = {symbolOpacity = 1.}},
        width = 700,
        height = 500,
        title = ["Python vs. Julia-Mimi versions of FaIR v2.0";"RCMIP Emission & Forcing Scenarios (1750-2100)"; "points Julia-Mimi; solid line original Python"]
    ) +
    @vlplot(
        data = df |> @filter(_.implementation == "FAIRv2.0 (original Python version)" && parse(Int64, _.time) >= 2000),
        :line
    )  +
    @vlplot(
        data = df |> @filter(_.implementation == "FAIRv2.0 (Julia-Mimi version)" && parse(Int64, _.time) in collect(2000:2:2100)),
        {:point, filled = true, stroke = "black", strokeWidth = 1.5}
    ) |> save(joinpath(outdir, "PythonReplication_Test_CO2_Concentrations_ZoomTime.svg"))
    
