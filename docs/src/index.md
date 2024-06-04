# Home

[![](https://img.shields.io/badge/docs-online-blue.svg)](https://Xiaomingzzhang.github.io/InvariantManifolds.jl/dev/)
[![Build Status](https://github.com/Xiaomingzzhang/InvariantManifolds.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Xiaomingzzhang/InvariantManifolds.jl/actions/workflows/CI.yml?query=branch%3Amaster)

## InvariantManifolds.jl

<img src="/assets/henon.png" alt="henon" width="600"/>

<!-- ```@example fast_code
using InvariantManifolds, LinearAlgebra, StaticArrays, OrdinaryDiffEq, GLMakie

# First define the lorenz vector field.
# Note that we have to rescale the system to ensure the uniform extention of the flow.
function lorenz(x, p, t)
    σ, ρ, β = p
    v = SA[σ*(x[2]-x[1]),
        ρ*x[1]-x[2]-x[1]*x[3],
        x[1]*x[2]-β*x[3]
    ]
    v / sqrt(1 + norm(v)^2)
end

function lorenz_map(x, p)
    prob = ODEProblem{false}(lorenz, x, (0.0, -2.0), p)
    sol = solve(prob, Tsit5())
    sol[end]
end

function eigenv(p)
    σ, ρ, β = p
    (SA[0.0, 0.0, 1.0], SA[-(-1 + σ + sqrt(1 - 2 * σ + 4 * ρ * σ + σ^2))/(2*ρ), 1, 0])
end

second(x) = x[2]

function myplot(annulus)
    fig = Figure()
    axes = LScene(fig[1, 1], show_axis=false,scenekw = (backgroundcolor=:white, clear=true))
    for i in eachindex(annulus)
        points = annulus[i].outer.pcurve.u
        scatter!(axes, first.(points), second.(points), last.(points), fxaa=true)
    end
    fig
end

para = [10, 28, 9 / 3]
lorenz_manifold = generate_surface(lorenz_map, para, SA[0.0, 0.0, 0.0], eigenv(para)..., 70, 2, 2)
myplot(lorenz_manifold)
``` -->