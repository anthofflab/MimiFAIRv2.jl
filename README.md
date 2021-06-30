# MimiFAIRv2.jl

This is a work-in-progress respository for a Julia-Mimi implementation of the FAIRv2.0 simple climate model. The model description paper can be found at [FaIRv2.0.0: A Generalized Impulse Response Model for Climate Uncertainty and Future Scenario Exploration](https://gmd.copernicus.org/articles/14/3007/2021/gmd-14-3007-2021.html). 

## Running Mimi-FAIRv2
To run the model, execute the following code:

```julia

# Load the model code.
include("src/MimiFAIRv2.jl") # load the MimiFAIRv2 module
using Main.MimiFAIRv2 # bring the module into your namespace

# Create an instance of Mimi-FAIRv2.
m = get_model() 

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