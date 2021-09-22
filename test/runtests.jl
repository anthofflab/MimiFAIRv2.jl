using DataFrames
using CSVFiles
using Test
using MimiFAIRv2
using Mimi

using MimiFAIRv2: get_model, create_fair_monte_carlo # load `get_model` function to avoid need for `MimiFAIRv2.` prefix

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

    # start year - should error if use a start year other than 1750
    @test_throws ErrorException get_model(start_year = 1760)

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

@testset "Python Comparison - With Python Replication Data" begin 

    # run model for each possible SSP and compare variable values between julia
    # and python versons

    for ssp in ["ssp119", "ssp126", "ssp245", "ssp370", "ssp585"]

        #Load initial conditions (just need parameter names) and extract gas names.
        init_gas_vals     = DataFrame(load(joinpath(@__DIR__, "..", "data", "fair_initial_gas_cycle_conditions_1750.csv"), skiplines_begin=7))
        montreal_init     = filter(:gas_group => ==("montreal"), init_gas_vals)
        flourinated_init  = filter(:gas_group => ==("flourinated"), init_gas_vals)
        aerosol_plus_init = filter(:gas_group => ==("aerosol_plus"), init_gas_vals)

        # Sort arrays of initial conditions for multiple gases so they are listed alphabetically.
        sort!(montreal_init,     :gas_name)
        sort!(flourinated_init,  :gas_name)
        sort!(aerosol_plus_init, :gas_name)

        # Load replication emissions and forcing from Python.
        python_emiss   = DataFrame(load(joinpath(@__DIR__, "..", "data", "python_replication_data", ssp*"_emissions.csv")))
        python_exog_RF = DataFrame(load(joinpath(@__DIR__, "..", "data", "python_replication_data", ssp*"_forcing.csv")))[:,:External]

        # Extract emissions arrays for multi-gas groupings.
        montreal_emissions     = python_emiss[:, Symbol.(montreal_init.gas_name)]
        flourinated_emissions  = python_emiss[:, Symbol.(flourinated_init.gas_name)]
        aerosol_plus_emissions = python_emiss[:, Symbol.(aerosol_plus_init.gas_name)]

        # Get a model instance for that SSP.
        m = MimiFAIRv2.get_model(emissions_forcing_scenario = ssp, start_year=1750, end_year=2100)

        # Update exogenous forcing and emissions data.
        update_param!(m, :radiative_forcing, :exogenous_RF, python_exog_RF)
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
        python_results = load(joinpath(@__DIR__, "..", "data", "python_replication_data", string(ssp, "_temperature.csv"))) |> DataFrame
        @test maximum(abs, julia_results[!, :temp] .- python_results[!, :default]) ≈ 0. atol = 1e-3
       
        # compare radiative forcing
        python_results = load(joinpath(@__DIR__, "..", "data", "python_replication_data", string(ssp, "_forcing.csv"))) |> DataFrame
        select!(python_results, "Forcing component", "Total")
        @test maximum(abs, julia_results[!, :rf] .- python_results[!, :Total]) ≈ 0. atol = 1e-3

        # compare emissions
        conv = select!(DataFrame(load(joinpath(@__DIR__, "..", "data", "default_gas_cycle_parameters.csv"), skiplines_begin=6)), :gas_name, :emis2conc)
        filter(:gas_name  => ==("carbon_dioxide"), conv)

        python_results = load(joinpath(@__DIR__, "..", "data", "python_replication_data", string(ssp, "_concentrations.csv"))) |> DataFrame
        select!(python_results, "Gas name", "carbon_dioxide")
        @test maximum(abs, julia_results[!, :co2] .- python_results[!, :carbon_dioxide]) ≈ 0. atol = 1e-3

    end
end

@testset "Python Comparison - With Default Model" begin 

    # run model for each possible SSP and compare variable values between julia
    # and python versons

    for ssp in ["ssp119", "ssp126", "ssp245", "ssp370", "ssp585"]

        # Get a model instance for that SSP.
        m = MimiFAIRv2.get_model(emissions_forcing_scenario = ssp, start_year=1750, end_year=2100)
        run(m)

        # Put CO₂ concentration, total radiative forcing, and global temperature anomaly into a dataframe.
        julia_results = DataFrame(co2 = m[:co2_cycle, :co2], rf = m[:radiative_forcing, :total_RF], temp = m[:temperature, :T])
        
        python_results = DataFrame(:year => collect(1750:2100))

        # add temperature
        res = load(joinpath(@__DIR__, "..", "data", "python_replication_data", string(ssp, "_temperature.csv"))) |> DataFrame
        insertcols!(python_results, 2, :temp => res[!, :default])

        # compare radiative forcing
        res = load(joinpath(@__DIR__, "..", "data", "python_replication_data", string(ssp, "_forcing.csv"))) |> DataFrame
        insertcols!(python_results, 2, :rf => res[!, :Total])

        # compare emissions
        res = load(joinpath(@__DIR__, "..", "data", "python_replication_data", string(ssp, "_concentrations.csv"))) |> DataFrame
        insertcols!(python_results, 2, :co2 => res[!, :carbon_dioxide])

        for name in [:co2, :rf, :temp]
            maxdiff = maximum(abs, julia_results[!, name] .- python_results[!, name])
            if maxdiff > 1e-3
                @warn("Maximum absolute difference for $(ssp)'s $name is $maxdiff")
            end
        end
    end
end

@testset "Monte Carlo" begin
    
    n_samples = 5

    # test the create_fair_monte_carlo function
    mcs_function_ssp585 = create_fair_monte_carlo(n_samples; delete_downloaded_data = false) # make delete_downloaded_data false for speed, will delete them at the end
    @test_throws ErrorException create_fair_monte_carlo(n_samples; start_year = 2000) # should error with a different start year

    # test the returned function
    results_ssp585 = mcs_function_ssp585()
    results_ssp585_zero_emissions= mcs_function_ssp585(co2_em_vals = fill(zeros(551), n_samples))
    @test results_ssp585[:temperatures][end] > results_ssp585_zero_emissions[:temperatures][end] # zero emissiosn should have lower end temperature than ssp585

    @test_throws ErrorException mcs_function_ssp585(co2_em_vals = fill(zeros(551), n_samples + 1)) # wrong number of samples
    @test_throws ErrorException mcs_function_ssp585(co2_em_vals = fill(zeros(250), n_samples)) # wrong vector length

    mcs_function_ssp245 = create_fair_monte_carlo(n_samples; emissions_scenario = "ssp245", delete_downloaded_data = false) # make delete_downloaded_data false for speed, will delete them at the end
    results_ssp245 = mcs_function_ssp245()

    @test results_ssp585[:temperatures][end] > results_ssp245[:temperatures][end] # ssp585 has higher end temperature than ssp245

    mcs_function_deterministic = create_fair_monte_carlo(n_samples; delete_downloaded_data = false, sample_id_subset = ["mem100057","mem10010","mem100210","mem100229", "mem100232"])
    results_deterministic_run1 = mcs_function_deterministic()
    results_deterministic_run2 = mcs_function_deterministic()
    @test results_deterministic_run1[:temperatures] == results_deterministic_run2[:temperatures] # should always pull the same samples

    # deletes the data because delete_downloaded_data should default to true
    MimiFAIRv2.create_fair_monte_carlo(n_samples)

end