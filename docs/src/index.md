# Home

[![](https://img.shields.io/badge/docs-online-blue.svg)](https://Xiaomingzzhang.github.io/InvariantManifolds.jl/dev/)
[![Build Status](https://github.com/Xiaomingzzhang/InvariantManifolds.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Xiaomingzzhang/InvariantManifolds.jl/actions/workflows/CI.yml?query=branch%3Amaster)

This is a package to compute the invariant manifolds of a dynamical system. Currently, I am still work on the one-dimensional manifolds of saddle points.

Main features:

- Compute saddles' one-dimensional manifolds of smooth mapping;
- Compute saddles' one-dimensional manifolds of non-smooth mapping: these mapping are the time-T-map of a non-smooth ODE systems, e.g., impact systems, piecewise smooth systems, and simple Fillippov systems;
- Compute saddles' two-dimensional manifolds of autonomous vector field.


# Basic example: Unstable manifold of Henon map
Consider the Henon map:

```math
\begin{aligned}
x'&=1-\alpha x^2+y,\\
y'&=\beta x,
\end{aligned}
```

where $\alpha=1.4,\beta=0.3$. This map has a saddle located at $(0.6313544770895048, 0.18940634312685142)$, and its unstable eigenvector is $(-6.412462860511356, 1.0)$. 

To compute the unstable manifolds of the saddle numerically, InvariantManifolds.jl just needs a segment of unstable manifold, whose start point is the saddle.
It's resonable to choose a short unstable eigenvector as the segment. You don't have to shorten the eigenvector started at the saddle yourself. We provide a function `segment` to do this automatically. The `segment` can generate equal distributed `n` points at one point, with given length and direction.
```julia
using InvaraiantManifolds, StaticArrays, Plots

function henonmap(x, p)
    SA[1 - p[1] * x[1]^2 + x[2], p[2] * x[1]]
end

seg = segment(SA[0.6313544770895048, 0.18940634312685142], SA[-6.412462860511356, 1.0], 150, 0.01)
result = generate_curves(henonmap, SA[1.4, 0.3], seg, 0.001, 8)
```

This package does not provide any function to plot the manifolds. However, it's simple to plot it using the stardard julia ploting library.

To do this, just define
```julia
using GLMakie

function myplot(v::Vector{IterationCurve})
    figure = Figure()
    axes = Axis(figure[1,1])
    for k in eachindex(v)
        data = v[k].pcurve.u
        lines!(axes, first.(data), last.(data))
    end
    figure
end

myplot(result)
```
![henon](/docs/src/assets/henon.png)

# An advanced example: Unstable manifold of the periodic pertubed system:

```julia
using GLMakie

function myplot(v::Vector{IterationCurve})
    figure = Figure()
    axes = Axis(figure[1,1])
    for k in eachindex(v)
        data = v[k].pcurve.u
        lines!(axes, first.(data), last.(data))
    end
    figure
end

myplot(result)
```