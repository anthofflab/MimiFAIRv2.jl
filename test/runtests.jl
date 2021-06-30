using DataFrames
using CSVFiles
using Test

include(joinpath(@__DIR__, "..", "src", "MimiFAIRv2.jl"))
using Main.MimiFAIRv2

@testset "UI" begin
    
    # run basic steps presented in README
    m = get_model()
    run(m)
    fair_temps = m[:temperature, :T]

    # Run model with altered keyword args and make sure results reflect the
    # appropriate changes, or at least a change at all
    # - start and end years
    # - emissions scenario
    # - TCR, RWF, and F2x

    # end year
    m = get_model(end_year = 2100)
    run(m)
    fair_temps = m[:temperature, :T]
    @test length(fair_temps) == length(1750:2100)

    # start year - should throw a warning if use a start year other than 1750
    @test_logs (:warn, "Model should not be set to start with a year differing from 1750.") get_model(start_year = 1760)
    m = get_model(start_year = 1760)
    run(m)
    fair_temps = m[:temperature, :T]
    @test length(fair_temps) == length(1760:2500)

    # emissions scenario
    m1 = get_model(emissions_forcing_scenario = "ssp126")
    m2 = get_model(emissions_forcing_scenario = "ssp585")
    run(m1)
    run(m2)
    m1_fair_temps = m1[:temperature, :T]
    m2_fair_temps = m2[:temperature, :T]
    # Mimi.plot(m1, :temperature, :T)
    # Mimi.plot(m2, :temperature, :T)
    @test m1_fair_temps !== m2_fair_temps
    @test maximum(m1_fair_temps) < maximum(m2_fair_temps)

    # TCR - transient climate response (default to 1.79)
    m1 = get_model()
    m2 = get_model(TCR = 2.0)
    run(m1)
    run(m2)
    m1_fair_temps = m1[:temperature, :T]
    m2_fair_temps = m2[:temperature, :T]
    # Mimi.plot(m1, :temperature, :T)
    # Mimi.plot(m2, :temperature, :T)
    @test m1_fair_temps !== m2_fair_temps
    @test maximum(m1_fair_temps) < maximum(m2_fair_temps)

    # RWF - realized warming fraction (default to 0.552)
    m1 = get_model()
    m2 = get_model(RWF = 0.5)
    run(m1)
    run(m2)
    m1_fair_temps = m1[:temperature, :T]
    m2_fair_temps = m2[:temperature, :T]
    # Mimi.plot(m1, :temperature, :T)
    # Mimi.plot(m2, :temperature, :T)
    @test m1_fair_temps !== m2_fair_temps
    @test maximum(m1_fair_temps) < maximum(m2_fair_temps)

    # F2x - forcing from a doubling of CO₂ (default to 3.759)
    m1 = get_model()
    m2 = get_model(F2x = 3.5)
    run(m1)
    run(m2)
    m1_fair_temps = m1[:temperature, :T]
    m2_fair_temps = m2[:temperature, :T]
    # Mimi.plot(m1, :temperature, :T)
    # Mimi.plot(m2, :temperature, :T)
    @test m1_fair_temps !== m2_fair_temps
    @test maximum(m1_fair_temps) < maximum(m2_fair_temps)

end

@testset "Python Comparison" begin 

    # run model for each possible SSP and compare variable values between julia
    # and python versons

    for SSP in ["ssp119", "ssp126", "ssp245", "ssp370", "ssp585"]

        m = get_model(;emissions_forcing_scenario = SSP, end_year = 2100)
        run(m)

        # compare temperatures (TODO - this is failing)
        julia_temps = m[:temperature, :T]

        data = load(joinpath(@__DIR__, "..", "data", "python_replication_data", string(SSP, "_temperature.csv"))) |> DataFrame
        rename!(data, [:year, :T])
        python_temps = data[!, :T]

        @test maximum(abs, julia_temps .- python_temps) ≈ 0. atol = 1e-3

        # TODO - what other variables do we want to compare?
        
    end
end
