using InvariantManifolds, StaticArrays, OrdinaryDiffEq, DataInterpolations, LinearAlgebra
using Test

@testset "InvariantManifolds.jl" begin
    include("nsstate.jl")
    include("ns_vector_field_constructors.jl")
    include("gen_seg_disk.jl")
    include("probs.jl")
end
