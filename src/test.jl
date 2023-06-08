using CSVFiles, DataFrames, Interpolations, Query, MimiFAIRv2, Mimi, Statistics, VegaLite

# EMISSIONS
benchmark_df = load("/Users/lisarennels/JuliaProjects/IAM-RDM/data/benchmark_ssps/benchmark_ssps_allgases.csv") |> DataFrame
benchmark_df = gas_df |> @filter(_.Gas == "N2O") |> DataFrame
benchmark_df = select(benchmark_df, Not(:Gas))
benchmark_df = stack(benchmark_df, Not(:Year))
benchmark_df.value = benchmark_df.value .* 28/44 ./ 1e3 

benchmark_df |> @vlplot(:line, x= "Year:q", y=:value, color = :variable, width = 500)

fair_df = DataFrame()
for ssp in ["ssp119", "ssp126", "ssp245", "ssp370", "ssp585"]
    df = load("/Users/lisarennels/.julia/dev/MimiFAIRv2/data/rcmip_$(ssp)_emissions_1750_to_2500.csv", skiplines_begin=6) |> DataFrame
    titles = names(df)
    titles[1] = "year"
    rename!(df, titles)
    df = select(df, :year, :nitrous_oxide)
    insertcols!(df, :SSP => ssp)
    append!(fair_df, df)
end
rename!(fair_df, [:Year, :value, :variable])

fair_df |> @vlplot(:line, x= "Year:q", y=:value, color = :variable, width = 500)

insertcols!(fair_df, :model => "fair")
insertcols!(benchmark_df, :model => "bench")

df = vcat(fair_df, benchmark_df)

df |> @filter(_.Year in collect(1900:2100)) |> @vlplot(
    # row = :model,
    :line, 
    x = "Year:q", 
    y = "value:q",
    color = "variable:n",
    width = 500
)

# TEMPERATURE 

n = 1000
my_fair_monte_carlo = MimiFAIRv2.create_fair_monte_carlo(n, delete_downloaded_data = false)

# default emissions
df = my_fair_monte_carlo()

T = DataFrame(df[:temperatures], Symbol.(1:1000))
insertcols!(T, :year => 1750:2300)
T = stack(T, Not(:year))
T.year = string.(T.year)

T = T |>
    @groupby({_.year}) |>
    @map({key(_)..., q05 = quantile(_.value, 0.05), q25 = quantile(_.value, 0.25), median = quantile(_.value, 0.5), q75 = quantile(_.value, 0.75), q95 = quantile(_.value, 0.95), mean = mean(_.value)}) |>
    DataFrame

T |> @filter(parse.(Int64, _.year) in (1990:2010)) |>
        @vlplot(
            width = 500,
            height = 500,
            x = {"year:t",timeUnit = :year,title = "year"},
            layer = [
                @vlfrag(
                    mark = {:errorband, opacity = 0.1},
                    y = {"q05:q",axis = {grid = true}, title = "T Anomaly (deg C)"},
                    y2 = "q95:q"
                ),
                @vlfrag(
                    mark = {:errorband, opacity = 0.2},
                    y = {"q25:q",axis = {grid = true}, title = "T Anomaly (deg C)"},
                    y2 = "q75:q"
                ),
                @vlfrag(
                    mark = {:line, strokeDash = "4,4"},
                    y = {"median:q",title = "T Anomaly (deg C)"}
                ),
                @vlfrag(
                    mark = :line,
                    y = {"mean:q",title = "T Anomaly (deg C)"}
                ) 
            ]
        )

# with other emissions
default_emissions = load("/Users/lisarennels/.julia/dev/MimiFAIRv2/data/rcmip_ssp585_emissions_1750_to_2500.csv", skiplines_begin=6) |> DataFrame
titles = names(default_emissions)
titles[1] = "year"
rename!(default_emissions, titles)

gas_df = load("/Users/lisarennels/JuliaProjects/IAM-RDM/data/benchmark_ssps/benchmark_ssps_allgases.csv") |> DataFrame
gas_df = gas_df |>
        @select(:Year, :Gas, :SSP2) |>
        DataFrame

idxs = indexin(1950:2300, 1750:2500)

co2_em_vals = gas_df |> @filter(_.Gas == "CO2") |> DataFrame
itp =  LinearInterpolation(co2_em_vals.Year, co2_em_vals.SSP2)
co2_em_vals = itp[collect(1950:2300)] .* 12/44 ./ 1e3 

co2_emissions = default_emissions.carbon_dioxide
co2_emissions[idxs] = co2_em_vals

n2o_em_vals = gas_df |> @filter(_.Gas == "N2O") |> DataFrame
itp =  LinearInterpolation(n2o_em_vals.Year, n2o_em_vals.SSP2)
n2o_em_vals = itp[collect(1950:2300)] .* 28/44 ./ 1e3 

n2o_emissions = default_emissions.nitrous_oxide
n2o_emissions[idxs] = n2o_em_vals

ch4_em_vals = gas_df |> @filter(_.Gas == "CH4") |> DataFrame
itp =  LinearInterpolation(ch4_em_vals.Year, ch4_em_vals.SSP2)
ch4_em_vals = itp[collect(1950:2300)]

ch4_emissions = default_emissions.methane
ch4_emissions[idxs] = ch4_em_vals

m = MimiFAIRv2.get_model()

update_param!(m, :co2_cycle, :E_co2, co2_emissions)
update_param!(m, :n2o_cycle, :E_n2o, n2o_emissions)
update_param!(m, :ch4_cycle, :E_ch4, ch4_emissions)

run(m)

df = my_fair_monte_carlo(co2_em_vals = fill(co2_emissions[1:551], n), n2o_em_vals = fill(n2o_emissions[1:551], n), ch4_em_vals = fill(ch4_emissions[1:551], n))

T = DataFrame(df[:temperatures], Symbol.(1:1000))
insertcols!(T, :year => 1750:2300)
T = stack(T, Not(:year))
T.year = string.(T.year)

T = T |>
    @groupby({_.year}) |>
    @map({key(_)..., q05 = quantile(_.value, 0.05), q25 = quantile(_.value, 0.25), median = quantile(_.value, 0.5), q75 = quantile(_.value, 0.75), q95 = quantile(_.value, 0.95), mean = mean(_.value)}) |>
    DataFrame

T |> @filter(parse.(Int64, _.year) > 2000) |>
        @vlplot(
            width = 800,
            height = 500,
            x = {"year:t",timeUnit = :year,title = "year"},
            layer = [
                @vlfrag(
                    mark = {:errorband, opacity = 0.1},
                    y = {"q05:q",axis = {grid = true}, title = "T Anomaly (deg C)"},
                    y2 = "q95:q"
                ),
                @vlfrag(
                    mark = {:errorband, opacity = 0.2},
                    y = {"q25:q",axis = {grid = true}, title = "T Anomaly (deg C)"},
                    y2 = "q75:q"
                ),
                @vlfrag(
                    mark = {:line, strokeDash = "4,4"},
                    y = {"median:q",title = "T Anomaly (deg C)"}
                ),
                @vlfrag(
                    mark = :line,
                    y = {"mean:q",title = "T Anomaly (deg C)"}
                ) 
            ]
        )