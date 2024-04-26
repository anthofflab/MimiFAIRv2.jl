using DataFrames
using CSVFiles
using Test
using MimiFAIRv2
using Mimi

@testset "MimiFAIRv2" begin

    @info("test_UI.jl")
    @time include("test_UI.jl")

    @info("test_mcs.jl")
    @time include("test_mcs.jl")

    @info("python_julia_comparison/03_julia_python_comparison_values.jl")
    @time(include("python_julia_comparison/03_julia_python_comparison_values.jl"))

    # Can locally run code to produce plots for comparison using 
    # python_julia_comparison/01 and 02 files
    
    @info("test_regression.jl")
    @time include("test_regression.jl")

end
