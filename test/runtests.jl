using DataFrames
using CSVFiles
using Test
using MimiFAIRv2
using Mimi

using MimiFAIRv2: get_model, create_fair_monte_carlo # load `get_model` function to avoid need for `MimiFAIRv2.` prefix

# @testset "MimiFAIRv2" begin

    @info("test_UI.jl")
    @time include("test_UI.jl")

    @info("test_python_comparison.jl")
    @time include("test_python_comparison.jl")

    @info("test_mcs.jl")
    @time include("test_mcs.jl")

    @info("test_plot_comparison.jl")
    @time include("test_plot_comparison.jl") # View results in test/output/test_plot_comparison

    @info("test_regression.jl")
    @time include("test_regression.jl")

# end
