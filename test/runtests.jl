using InvariantManifolds, StaticArrays, OrdinaryDiffEq, DataInterpolations
using Test

@testset "InvariantManifolds.jl" begin
    include("nsstate.jl")
    include("ns_vector_field_constructors.jl")
    include("smoothone.jl")
    include("nonsmoothone.jl")
end
