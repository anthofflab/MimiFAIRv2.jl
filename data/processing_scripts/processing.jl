using DataFrames, CSVFiles, Query, Interpolations

for ssp in ["ssp119", "ssp126", "ssp245", "ssp370", "ssp585"]
    df = load(joinpath(@__DIR__, "rcmip_$(ssp)_emissions_1750_to_2500.csv")) |> DataFrame
    titles = names(df)
    titles[1] = "year"
    rename!(df, titles)

    for col in names(df)[2:end]
        idxs = (!ismissing).(df[:,col])
        itp = LinearInterpolation(df[:, :year][idxs], df[:,col][idxs])
        df[:,col] .= itp[df.year]
    end

    df |> save(joinpath(@__DIR__, "../rcmip_$(ssp)_emissions_1750_to_2500.csv"))
end

# : `ch2cl, chcl3, methyl_bromine, methyl_chlorine, so2, nox, co, nnvoc, bc, nh3, oc` 

main_df = DataFrame()
for ssp in ["ssp119", "ssp126", "ssp245", "ssp370", "ssp585"]
    df = load(joinpath(@__DIR__, "rcmip_$(ssp)_emissions_1750_to_2500_itp.csv"), skiplines_begin = 6) |> @select(:year, :nox) |> DataFrame
    insertcols!(df, :ssp => ssp)
    insertcols!(df, :type => "new")
    append!(main_df, df)
end

for ssp in ["ssp119", "ssp126", "ssp245", "ssp370", "ssp585"]
    df = load(joinpath(@__DIR__, "..", "rcmip_$(ssp)_emissions_1750_to_2500.csv"), skiplines_begin = 6) |> @select(:nox) |> DataFrame
    insertcols!(df, :ssp => ssp)
    insertcols!(df, :type => "old")
    insertcols!(df, 1, :year => collect(1750:2500))
    append!(main_df, df)
end

main_df |> @vlplot(row = :ssp, :line, x = :year, y = :nox, color = :type, height = 100, width = 600)
