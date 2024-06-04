# InvariantManifolds.jl

[![](https://img.shields.io/badge/docs-online-blue.svg)](https://Xiaomingzzhang.github.io/InvariantManifolds.jl/dev/)
[![Build Status](https://github.com/Xiaomingzzhang/InvariantManifolds.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Xiaomingzzhang/InvariantManifolds.jl/actions/workflows/CI.yml?query=branch%3Amaster)

This is a package to compute the invariant manifolds of a dynamical system.

Main features:

- Compute saddles' one-dimensional manifolds of smooth mapping;
- Compute saddles' two-dimensional manifolds of autonomous vector field;
- Compute saddles' one-dimensional manifolds of non-smooth mapping: these mapping are the time-T-map of a non-smooth ODE systems, e.g., impact systems, piecewise smooth systems, and simple Fillippov systems.

To use this package, first install [julia](https://julialang.org/), then run 
```julia
using Pkg;
Pkg.add("https://github.com/Xiaomingzzhang/InvariantManifolds.jl")
```

## Basic example: Unstable manifold of Henon map
Consider the Henon map:

$$
\begin{aligned}
x'&=1-\alpha x^2+y,\\
y'&=\beta x,
\end{aligned}
$$

where $\alpha=1.4,\beta=0.3$. This map has a saddle located at $(0.6313544770895048, 0.18940634312685142)$, and its unstable eigenvector is $(-0.9880577559947047, 0.15408397327012555)$. 

To compute the unstable manifolds of the saddle numerically, InvariantManifolds.jl just needs a segment of unstable manifold, whose start point is the saddle.
It's resonable to choose a short unstable eigenvector as the segment. You don't have to shorten the eigenvector started at the saddle yourself. We provide a function `segment` to do this automatically. The `segment` can generate equal distributed `n` points at one point, with given length and direction.
```julia
using InvaraiantManifolds, StaticArrays

function henonmap(x, p)
    y1 = 1 - p[1] * x[1]^2 + x[2]
    y2 = p[2] * x[1]
    SA[y1, y2]
end

henonmap2(x, p)=henonmap(henonmap(x, p), p)

# For technical reason, we have to use the double iteration of the map, since the eigenvalue is less than -1.

seg = segment(SA[0.6313544770895048, 0.18940634312685142], SA[-0.9880577559947047, 0.15408397327012555], 150, 0.01)
result = generate_curves(henonmap2, SA[1.4, 0.3], seg, 0.001, 8)
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

## An advanced example: Unstable manifold of the periodic pertubed system:

```julia
using InvariantManifolds, LinearAlgebra, StaticArrays, OrdinaryDiffEq, GLMakie
# Duffing Equation

f(x, p, t) = SA[x[2], x[1] - x[1]^3 + p[1]*cos(t)]

# First find the fixed point and its unstable direction by using the Newton iteration.

function timemap(x,p)
    prob = ODEProblem{false}(f, x, (0.0, 2pi), p)
    solve(prob, Vern9())[end]
end

function jac(x, p)
    prob = ODEProblem{false}(f, x, (0.0, 2pi), p)
    sol = solve(prob, Vern9())
    function df(x, p, t)
        SA[0 1; 1-3*(sol(t)[1])^2 0] * x
    end
    ii = SA[1.0 0.0; 0.0 1.0]
    nprob = ODEProblem{false}(df, ii, (0.0, 2pi), p)
    solve(nprob, Vern9())[end]
end

function newton(x, p; n=100, atol=1e-8)
    xn = x - inv(jac(x, p) - I) * (timemap(x, p) - x)
    data = typeof(x)[x, xn]
    i = 1
    while norm(data[2] - data[1]) > atol && i <= n
        data[1] = data[2]
        data[2] = data[1] - inv(jac(data[1], p) - I) * (timemap(data[1], p) - data[1])
        i = i + 1
    end
    if norm(data[2] - data[1]) < atol
        println("Fixed point found successfully:")
        data[2]
    else
        println("Failed to find a fixed point after $n times iterations. The last point is:")
        data[2]
    end
end

function myplot(v)
    figure = Figure()
    axes = Axis(figure[1,1])
    for k in eachindex(v)
        data = v[k].pcurve.u
        lines!(axes, first.(data), last.(data))
    end
    figure
end

para = [0.1]
fixedpoint = newton(SA[-0.05, 0.0], para)
unstable_direction = eigen(jac(fixedpoint, para)).vectors[:, 2]
seg = segment(fixedpoint, unstable_direction, 150, 0.01)
result = generate_curves(timemap, para, seg, 0.002, 3)
myplot(result)
```
![duffing](/docs/src/assets/duffing.png)

## Lorenz manifold:

```julia
using InvariantManifolds, LinearAlgebra, StaticArrays, OrdinaryDiffEq, GLMakie

# First define the lorenz vector field. Note that we have to resized the system to ensure the uniform extention.
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
lorenz_manifold = generate_surface(lorenz_map, para, SA[0.0, 0.0, 0.0], eigenv(para)..., 120, 1, 1; n=150)

myplot(lorenz_manifold)
```
![lorenz](/docs/src/assets/lorenz.png)

See more examples in the docs.