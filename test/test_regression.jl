module testRegression

using DataFrames
using CSVFiles
using Test
using MimiFAIRv2
using Mimi

using MimiFAIRv2: get_model, create_fair_monte_carlo 
 
outdir = joinpath(@__DIR__, "regression_testing")

##
## INITIAL CREATION OF BASELINE DATA TO COMPARE AGAINST (run with a new version name if expect changes)
##

# version = "saved_XXX"

# # build up dataframe of results
# df = DataFrame()
# for ssp in ["ssp119", "ssp126", "ssp245", "ssp370", "ssp585"]
#     m = MimiFAIRv2.get_model(emissions_forcing_scenario=ssp, start_year=1750, end_year=2100)
#     run(m)
#     append!(df, DataFrame(co2 = m[:co2_cycle, :co2], rf = m[:radiative_forcing, :total_RF], temp = m[:temperature, :T], ssp = ssp))
# end
# df |> save(joinpath(outdir, "output_data_$version.csv"))

##
## COMPARISON
##

regression_baseline_version = "saved_07122023_rebased"
baseline_data = load(joinpath(outdir, "output_data_$(regression_baseline_version).csv")) |> DataFrame

df = DataFrame()
for ssp in ["ssp119", "ssp126", "ssp245", "ssp370", "ssp585"]
    m = MimiFAIRv2.get_model(emissions_forcing_scenario=ssp, start_year=1750, end_year=2100)
    run(m)
    append!(df, DataFrame(co2 = m[:co2_cycle, :co2], rf = m[:radiative_forcing, :total_RF], temp = m[:temperature, :T], ssp = ssp))
end

for var in [:co2, :rf, :temp]
    @test baseline_data[:, var] ≈ df[:, var] atol = 1e-9
end

end