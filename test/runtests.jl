using RIXSCalculator
using Test

@testset "RIXSCalculator.jl" begin

    # include subfile for abstract type tests
    include("tests_abstract_types.jl")

    # include subfile for concrete type tests
    include("tests_concrete_types.jl")

    # include subfile for eigensysten tests
    include("tests_eygensystem.jl")

end