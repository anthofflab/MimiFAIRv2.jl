module testUI

using DataFrames
using CSVFiles
using Test
using MimiFAIRv2
using Mimi

using MimiFAIRv2: get_model, create_fair_monte_carlo 

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
@test_logs (:warn, "FAIRv2 model monte carlo simulation should not be set to start with a year differing from 1750 as initial conditions are not calibrated for a different start year!") get_model(start_year = 1760)

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