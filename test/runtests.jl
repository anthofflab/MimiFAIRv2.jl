using DataFrames
using CSVFiles
using Test
using MimiFAIRv2
using Mimi

using MimiFAIRv2: get_model, create_fair_monte_carlo 

@testset "MimiFAIRv2" begin

    @info("test_UI.jl")
    @time include("test_UI.jl")

    @info("test_python_comparison.jl")
    @time include("test_python_comparison.jl")

    @info("test_mcs.jl")
    @time include("test_mcs.jl")

    # run locally and see test/output/test_plot_comparison for graphs
    # @info("test_plot_comparison.jl")
    # @time include("test_plot_comparison.jl")

    @info("test_regression.jl")
    @time include("test_regression.jl")

end
