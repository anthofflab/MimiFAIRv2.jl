# MimiFAIRv2.jl

This is a work-in-progress repository for a Julia-Mimi implementation of the FAIRv2.0 simple climate model. The model description paper can be found at [FaIRv2.0.0: A Generalized Impulse Response Model for Climate Uncertainty and Future Scenario Exploration](https://gmd.copernicus.org/articles/14/3007/2021/gmd-14-3007-2021.html). 

## Preparing the Software Environment

To add the package to your current environment, run the following command at the julia package REPL:

```julia
pkg> add https://github.com/FrankErrickson/MimiFAIRv2.jl.git
```

```julia
include("src/MimiFAIRv2.jl") # load the MimiFAIRv2 module 
```

You probably also want to install the Mimi package into your julia environment, so that you can use some of the tools in there:

```julia
pkg> add Mimi
```

## Running the Model

The model uses the Mimi framework and it is highly recommended to read the Mimi  documentation first to understand the code structure. The basic way to access a copy of the default MimiFAIRv2 model and explore the resuts is the following:

```julia
using Mimi 
using MimiFAIRv2

# Create an instance of MimiFAIRv2.
m = MimiFAIRv2.get_model() 

# Run the model.
run(m)

# Access the temperature output.
fair_temps = m[:temperature, :T]

# Explore interactive plots of all the model output.
explore(m)
```

The `get_model()` function currently has the following keyword arguments:  

* `emissions_forcing_scenario`: One of the RCMIP scenarios from the original FAIRv2.0 paper. Current options include "ssp119", "ssp126", "ssp245", "ssp370", and"ssp585". The default is "ssp585".  
* `start_year`: The model has an option to be initialized at different time periods, however this is only currently set up to start in 1750.
* `end_year`: The model can be run out to 2500 (the default final year).  
* `TCR`: The transient climate response (default = 1.79).  
* `RWF`: The realized warming fraction, defined as the TCR/ECS ratio (default = 0.552).  
* `F2x`: The forcing from a doubling of CO2 (default = 3.759).  

\
\
![Python vs. Julia temperature comparison](https://github.com/FrankErrickson/MimiFAIRv2.jl/blob/main/data/python_replication_data/Python_Mimi_FAIR2_temperature_comparison.png)

## Running a Monte Carlo Simulation

See `monte_carlo.jl` for details on running a Monte Carlo Simulation, and don't hesitate to contact the developers with any questions, or post an Issue on Github.
