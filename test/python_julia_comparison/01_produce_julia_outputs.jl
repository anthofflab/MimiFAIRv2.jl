using DataFrames
using CSVFiles
using Mimi
using MimiFAIRv2
using Query

outdir = joinpath(@__DIR__, "julia_data")
mkpath(outdir)

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
    # and Save results.
    mkpath(joinpath(outdir))
    DataFrame(co2 = m[:co2_cycle, :co2], rf_total = m[:radiative_forcing, :total_RF], rf_co2 = m[:radiative_forcing, :co2_RF], temp = m[:temperature, :T]) |>
        save(joinpath(outdir, "Mimi_FAIR_replication_"*ssp*".csv"))
end 
