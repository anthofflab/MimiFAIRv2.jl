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

### Overview

See `monte_carlo.jl` for the script and details on running a Monte Carlo Simulation, and don't hesitate to contact the developers with any questions, or post an Issue on Github.

This file contains functions to run a Monte Carlo with MimiFAIRv2 using the constrained parameters from Leach et al. (2021). The first function, `create_fair_monte_carlo`, loads the constrained parameter samples from Leach et al. (2021) and cleans them up so they can be easily passed into Mimi-FAIRv2.0. It then creates a second function, in the script called `fair_monte_carlo` that performs the actual analysis. This function returns a dictionary with temperature, atmospheric co2, ch4, and n2o, and radiative forcing. Adding other outputs is doable please add an Issue on Github if you would like this to be done. This returned function can optionally take
a vector of vectors to force co2, ch4, and n2o emissions trajectories for the `n_samples` runs, providing these in the native MimiFAIRv2 emissions units, using optional arguments `co2_em_vals`, `ch4_em_vals`, and `n2o_em_vals` respectively. 

The required constrained parameter files are uploaded to Zenodo with the citation Frank Errickson, & Lisa Rennels. (2021). MimiFAIRV2 Large Data File Storage (0.1.0-DEV) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5513221

### Function Details

The `create_fair_monte_carlo` function is the primary user-facing function provided for the monte carlo simulation and has the signature and function arguments as follows:

```julia
function create_fair_monte_carlo(
    n_samples::Int; 
    emissions_scenario::String="ssp585", 
    start_year::Int=1750, 
    end_year::Int=2300, 
    data_dir::String = joinpath(@__DIR__, "..", "data", "large_constrained_parameter_files"),
    sample_id_subset::Union{Vector, Nothing} = nothing,
    delete_downloaded_data::Bool = true
)
```

Function Arguments:

- n_samples: Number of samples to randomly draw from original constrained parameter values (max = 93,995).
- emissions_scenario: Current options are: "ssp119", "ssp126", "ssp245", "ssp370", "ssp585"
- start_year: First year to run the model (note, Mimi-FAIR requires user-supplied initial conditions if not starting in 1750 so this is not yet supported for a year other than 1750).
- end_year: Final year to run the model.
- data_dir: Location of the data files
- sample_id_subset: IDs of the subset of samples to use from original constrained parameter values (each ID between 1 and 93,995  inclusive). If this argument is set, it will be used instead of the step to create n random indices. The length  should match n_samples, and if it is longer we will choose the first 1:n.
- delete_downloaded_data: Boolean (defaulting to false) to recursively delete all downloaded large data files of constrained parameters at the end of the script.  Should be set to `false` if one will be running this several times and want to avoid re-downloading each time.

Calling `create_fair_monte_carlo` as follows:

```julia
my_fair_monte_carlo = create_fair_monte_carlo(1_000)
```

will return another function, in this case named `my_fair_monte_carlo`, that performs the actual analysis. This function returns a dictionary with keys temperature, atmospheric co2, ch4, and n2o, and radiative forcing. Adding other outputs is doable please add an Issue on Github if you would like this to be done. This returned function can optionally take a vector of vectors to force co2, ch4, and n2o emissions trajectories for the `n_samples` runs, providing these in the native MimiFAIRv2 emissions units, using optional arguments `co2_em_vals`, `ch4_em_vals`, and `n2o_em_vals` respectively. 

```julia
function fair_monte_carlo( ; 
    co2_em_vals::Union{Nothing, Vector{Vector{T1}}}  = nothing,
    n2o_em_vals::Union{Nothing, Vector{Vector{T2}}} = nothing,
    ch4_em_vals::Union{Nothing, Vector{Vector{T3}}} = nothing
) where {T1, T2, T3}
```

Function Arguments: 

- co2_em_vals: A vector with n_samples elements, each of which is a vector spanning the full time of the model being run and holding an exogenous stream of co2 emissions in the native FAIRv2 units, GtC (convert to GtCO2 with multiplier 44/12)
- n2o_em_vals: A vector with n_samples elements, each of which is a vector spanning the full time of the model being run and holding an exogenous stream of n2o emissions in the native FAIRv2 units, GtN (convert to GtN2O with multiplier 44/28)
- ch4_em_vals: A vector with n_samples elements, each of which is a vector spanning the full time of the model being run and holding an exogenous stream of ch4 emissions in the native FAIRv2 units, GtCH4.
