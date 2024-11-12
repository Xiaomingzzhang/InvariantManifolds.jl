module InvariantManifolds

using LinearAlgebra, StaticArrays, OrdinaryDiffEq, DataInterpolations, NearestNeighbors

import Base: -, +, *, /, length, eltype, getindex, size, show

export NSState, NSSolution, NSSetUp

export PiecewiseV, BilliardV, SFilippovV

export OneDManifoldProblem, OneDManifold, NSOneDManifoldProblem, NSOneDManifold, FlawPoint

export TwoDManifoldProblem, TwoDManifold, VTwoDManifoldProblem, VTwoDManifold

export setmap, timetmap, ns_solver

export gen_segment, gen_disk, initialize

export grow!, growmanifold, Saddle, findsaddle

include("basic_constructors.jl")
include("nsstate_construction.jl")

# the algorithm for smooth nonlinear maps' one dimensional manifolds
include("smooth-one.jl")

# simple algorithm for smooth nonlinear maps' two dimensional manifolds
include("smooth_maps_two.jl")

# simple algorithm for autonomous vector fields ' two dimensional manifolds
include("smooth_vectorfield_two.jl")

# construct the time-T-map for different non-smooth ODE systems
include("piecewise.jl")
include("impact.jl")
include("simple_filippov.jl")

# the one dimensional non-smooth manifold algorithm
include("nonsmooth-one.jl")


# newton method to locate the saddles of smooth vector fields' time-T-map
include("newton.jl")

end