
using Test
using SystemDynamicsBuildR 

@testset "SystemDynamicsBuildR tests" begin
    println("About to include plusTwo.jl")
    include("plusTwo.jl")
    println("About to include add_numbers.jl")
    include("add_numbers.jl")
    println("Done including files")

    include("test_sdbuildR_units.jl")
    include("test_custom_func.jl")
    include("test_clean.jl")
    include("test_ensemble.jl")
    include("test_unit_func.jl")
end