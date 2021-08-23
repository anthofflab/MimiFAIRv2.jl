using DataFrames, CSVFiles

# Helper script to randomly sample n indices of constrained parameters and save
# the resulting files

function _get_constrained_parameters_subset(n_samples::Int, data_dir::String, outdir::String)

    #Load constrained FAIR parameters from Leach et al. (2021).
    thermal_params           = DataFrame(load(joinpath(data_dir, "julia_constrained_thermal_parameters_average_probs.csv")))
    gas_params               = DataFrame(load(joinpath(data_dir, "constrained_parameters", "julia_constrained_gas_cycle_parameters_average_probs.csv")))
    indirect_forcing_params  = DataFrame(load(joinpath(data_dir, "data", "constrained_parameters", "julia_constrained_indirect_forcing_parameters_average_probs.csv")))
    exogenous_forcing_params = DataFrame(load(joinpath(data_dir, "data", "constrained_parameters", "julia_constrained_exogenous_forcing_parameters_average_probs.csv")))

    # Load sample-specific initial conditions (defaults for the year 2000 as initial year).
    thermal_initial_conditions   = DataFrame(load(joinpath(@__DIR__, "..", "data", "constrained_parameters", "initial_thermal_conditions_constrained_average_probs.csv")))
    gas_cycle_initial_conditions = DataFrame(load(joinpath(@__DIR__, "..", "data", "constrained_parameters", "initial_gas_cycle_conditions_constrained_average_probs.csv")))
    
    # Create random indices and create a subset of FAIR sample id values (for indexing data).
    rand_indices     = sort(sample(1:93995, n_samples, replace=false))
    sample_id_subset = thermal_params[rand_indices, :sample_id]

    # Subset constrained parameters using random sample id values.
    thermal_params           = thermal_params[findall((in)(sample_id_subset), thermal_params.sample_id), :]
    gas_params               = gas_params[findall((in)(sample_id_subset), gas_params.sample_id), :]
    indirect_forcing_params  = indirect_forcing_params[findall((in)(sample_id_subset), indirect_forcing_params.sample_id), :]
    exogenous_forcing_params = exogenous_forcing_params[findall((in)(sample_id_subset), exogenous_forcing_params.sample_id), :]

    thermal_initial_conditions = thermal_initial_conditions[findall((in)(sample_id_subset), thermal_initial_conditions.sample_id), :]
    gas_cycle_initial_conditions = gas_cycle_initial_conditions[findall((in)(sample_id_subset), gas_cycle_initial_conditions.sample_id), :]

    # save
    thermal_params           |> save(joinpath(outdir, "julia_constrained_thermal_parameters_average_probs_$(n).csv"))
    gas_params               |> save(joinpath(outdir, "constrained_parameters", "julia_constrained_gas_cycle_parameters_average_probs_$(n).csv"))
    indirect_forcing_params  |> save(joinpath(outdir, "data", "constrained_parameters", "julia_constrained_indirect_forcing_parameters_average_probs_$(n).csv"))
    exogenous_forcing_params |> save(joinpath(data_dir, "data", "constrained_parameters", "julia_constrained_exogenous_forcing_parameters_average_probs_$(n).csv"))
    
    thermal_initial_conditions  |> save(joinpath(@__DIR__, "..", "data", "constrained_parameters", "initial_thermal_conditions_constrained_average_probs_$(n).csv"))
    gas_cycle_initial_conditions|> save(joinpath(@__DIR__, "..", "data", "constrained_parameters", "initial_gas_cycle_conditions_constrained_average_probs_$(n).csv"))
    
end
