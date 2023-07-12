module testPythonComparison

using DataFrames
using CSVFiles
using Test
using MimiFAIRv2
using Mimi

using MimiFAIRv2: get_model, create_fair_monte_carlo # load `get_model` function to avoid need for `MimiFAIRv2.` prefix
 
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