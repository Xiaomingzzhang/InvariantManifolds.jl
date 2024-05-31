module InvariantManifolds

using LinearAlgebra, StaticArrays, OrdinaryDiffEq, DataInterpolations

import Base: -, show, insert!

import DataInterpolations.LinearInterpolation

export State, NSState, IterationCurve, NSSolution, NSSetUp, AnnulusBoundaries

export PiecewiseV, BilliardV, SFilippovV

export setmap, timetmap, ns_solver

export segment, generate_curves, generate_surface

include("basic_constrators.jl")

# the one dimensional smooth manifold algorithm
include("smooth-one.jl")

# construct the time-T-map for different nonsmooth ODE systems
include("piecewise.jl")
include("impact.jl")
include("simple_filippov.jl")

# construct the time-T-map for mixed types of nonsmooth ODE systems, 
#e.g., imapct&&piecewise smooth, filippov&&impact.
# include("mixed.jl")

# the one dimensional non-smooth manifold algorithm
include("nonsmooth-one.jl")

# vector fields's two dimesional manifolds
include("smooth-vectorfield-two.jl")

end