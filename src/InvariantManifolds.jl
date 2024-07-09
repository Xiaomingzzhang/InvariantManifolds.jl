module InvariantManifolds

using LinearAlgebra, StaticArrays, OrdinaryDiffEq, DataInterpolations, NearestNeighbors

import Base: -, +, *, /, length, eltype, getindex, size, show, insert!

import DataInterpolations: AbstractInterpolation, show

export NSState, NSSolution, NSSetUp, AnnulusBoundaries

export PiecewiseV, BilliardV, SFilippovV

export setmap, timetmap, ns_solver

export inintialise_mesh, initialise_curve, ns_initialise_curve

export grow_line!, grow_surface!, generate_curves, generate_surface

include("basic_constructors.jl")
include("nsstate_construction.jl")

# the algorithm for smooth nonlinear maps' one dimensional manifolds
include("smooth-one.jl")

# the algorithm for autonomous smooth vector fields' two dimensional manifolds
include("smooth_vectorfield_two.jl")

# simple algorithm for smooth nonlinear maps' two dimensional manifolds
include("smooth_maps_two.jl")

# construct the time-T-map for different non-smooth ODE systems
include("piecewise.jl")
include("impact.jl")
include("simple_filippov.jl")

# the one dimensional non-smooth manifold algorithm
include("nonsmooth-one.jl")

include("nonsmooth_maps_two.jl")

# newton method to locate the saddles of smooth vector fields' time-T-map
include("newton.jl")

end