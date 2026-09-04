module InvariantManifolds

using LinearAlgebra, StaticArrays, OrdinaryDiffEq, DataInterpolations, NearestNeighbors, FiniteDiff

import Base: -, +, *, /, length, eltype, getindex, size, zero, show

export NSState, NSSetUp

export PiecewiseV, BilliardV, PiecewiseImpactV

export OneDManifoldProblem, NSOneDManifoldProblem, VTwoDManifoldProblem

export TwoDManifoldProblem, NSVTwoDManifoldProblem

export setmap, ns_solver

export gen_segment, gen_disk, initialize

export grow!, growmanifold, Saddle, findsaddle, iscontact

export thin_curve, ns_monodromy


include("nsstate_construction.jl")
include("basic_constructors.jl")
include("curve_thinning.jl")

# internal helpers for counting map/time-T-map calls
@inline function _count_map_call!(counter)
    counter === nothing || Threads.atomic_add!(counter, 1)
    nothing
end

function _new_map_counter()
    Threads.Atomic{Int}(0)
end

function _counting_map(f, counter)
    (x, p) -> (Threads.atomic_add!(counter, 1); f(x, p))
end

# the algorithm for smooth nonlinear maps' one dimensional manifolds
include("smooth-one.jl")

# simple algorithm for smooth nonlinear maps' two dimensional manifolds
include("smooth_maps_two.jl")

# simple algorithm for autonomous vector fields ' two dimensional manifolds
include("smooth_vectorfield_two.jl")

# construct the time-T-map for different non-smooth ODE systems
include("piecewise.jl")
include("impact.jl")
include("piecewise_impact.jl")
include("monodromy.jl")

# non-smooth manifold algorithm
include("nonsmooth-one.jl")
include("nonsmooth_vectorfield_two.jl")

# newton method to locate the saddles of smooth vector fields' time-T-map
include("newton.jl")
end
