# InvariantManifolds

[![](https://img.shields.io/badge/docs-online-blue.svg)](https://Xiaomingzzhang.github.io/InvariantManifolds.jl/dev/)
[![Build Status](https://github.com/Xiaomingzzhang/InvariantManifolds.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Xiaomingzzhang/InvariantManifolds.jl/actions/workflows/CI.yml?query=branch%3Amaster)

This is a package to compute the one-dimensional invariant manifolds of saddle points for nonlinear mappings.

Main features:

- Compute saddles' one-dimensional manifolds of smooth mapping;
- Compute these manifolds in any precisions;
- Compute saddles' one-dimensional manifolds of non-smooth mapping: these mapping are the time-T-map of a non-smooth ODE systems, e.g., impact systems, piecewise smooth systems, and simple Filippov systems.

# Basic example: Unstable manifold of Henon map
Consider the Henon map:

$$
x'=1-\alpha x^2+y,\\
y'=\beta x,
$$

where $\alpha=1.4,\beta=0.3$. This map has a saddle located at $(0.6313544770895048, 0.18940634312685142)$, and its unstable eigenvector is $(-0.9880577559947047, 0.15408397327012555)$. 

To compute the unstable manifolds of the saddle numerically, InvariantManifolds.jl just needs a segment of unstable manifold, whose start point is the saddle.
It's reasonable to choose a short unstable eigenvector as the segment. You don't have to shorten the eigenvector started at the saddle yourself. We provide a function `segment` to do this automatically. The `segment` can generate equal distributed `n` points at one point, with given length and direction.
```julia
using InvariantManifolds, StaticArrays, Plots

function henonmap(x, p)
    y1 = 1 - p[1] * x[1]^2 + x[2]
    y2 = p[2] * x[1]
    SA[y1, y2]
end

seg = segment(SA[0.6313544770895048, 0.18940634312685142], SA[-0.9880577559947047, 0.15408397327012555], 150, 0.01)
result = generate_curves(henonmap, SA[1.4, 0.3], seg, 0.005, 13)
plot(result)
```
![henon](/docs/src/assets/henon.svg)

You can use `Plotly` backends to see the details of manifolds. To do this, just define
```julia
function myplot(v::Vector{IterationCurve{N,T}}) where {N,T}
    plotlyjs()
    figure = plot(;legend=false)
    for i in eachindex(v)
        nn = length(v[i].states)
        id = range(0, 1, length=2 * nn)
        data = v[i].pcurve.(id)
        plot!(first.(data), last.(data))
    end
    figure
end
```
Run `myplot(result)` to get the figure that can be zoom in.

See more examples in document.

These features may be added in the future:

- Compute saddles' high-dimensional manifolds of smooth mapping;
- Compute invariant tori.