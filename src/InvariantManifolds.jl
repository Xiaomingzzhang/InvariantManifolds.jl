module InvariantManifolds

using LinearAlgebra, StaticArrays, OrdinaryDiffEq, Plots

import Interpolations.linear_interpolation
import Base.-, Base.show
import Plots.plot

export State, IterationCurve, show, segment, generate_curves, plot

export PiecewiseV, BilliardV, SFillippoV, NSState, NSSetUp, setmap, gen_prob

# the one dimensional smooth manifold algorithm
include("smooth-one.jl")

# construct the time-T-map for different nonsmooth ODE systems
include("piecewise.jl")
include("impact.jl")
include("simple_filippov.jl")

# construct the time-T-map for mixed types of nonsmooth ODE systems, 
#e.g., imapct&&piecewise smooth, filippov&&impact.
include("mixed.jl")

# the one dimensional non-smooth manifold algorithm
include("nonsmooth-one.jl")

# simple Newton method to find the saddle
include("newton.jl")

# floquet theory to compute the jacobian matrix of time-T-map
include("floq.jl")

end