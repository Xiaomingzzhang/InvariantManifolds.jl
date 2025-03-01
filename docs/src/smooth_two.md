# Two-Dimensional Smooth Manifolds

Essentially, we haven't introduced new algorithms; the core functions for computing two-dimensional manifolds are the same as those for one-dimensional manifolds. We simply represent the two-dimensional manifold as circles of one-dimensional manifolds that are close enough to each other.

```@setup lorenz
using InvariantManifolds, LinearAlgebra, StaticArrays, OrdinaryDiffEq, CairoMakie
```
## Autonomous Vector Field: Lorenz Manifold
First, let's load the required packages and define the Lorenz vector field:
```@repl lorenz
using InvariantManifolds, LinearAlgebra, StaticArrays, OrdinaryDiffEq, CairoMakie
function lorenz(x, p, t)
    σ, ρ, β = p
    v = SA[σ*(x[2]-x[1]),
        ρ*x[1]-x[2]-x[1]*x[3],
        x[1]*x[2]-β*x[3]
    ]
    v / sqrt(0.1 + norm(v)^2)
end
```

It's worth noting that we performed an approximate normalization of the vector field to keep its magnitude within a small range. This ensures uniform expansion of the manifold. Under the classical parameters:
```@repl lorenz
para = [10.0, 28.0, 8/3]
```
the Jacobian matrix at the origin equilibrium point has two stable directions, which are:
```@example lorenz
function eigenv(p)
    σ, ρ, β = p
    [SA[0.0, 0.0, 1.0], SA[-(-1 + σ + sqrt(1 - 2 * σ + 4 * ρ * σ + σ^2))/(2*ρ), 1, 0]]
end
eigenv(para)
```
Then we can create a [`Saddle`](@ref) structure to store this saddle point:
```@repl lorenz
saddle = Saddle(SA[0, 0, 0.0], eigenv(para), [1.0, 1.0])
```
The magnitude of the eigenvalues can be specified arbitrarily, as it won't affect the computation results. Since we're computing the stable manifold, we need to evolve the flow backward. Let's define the following mapping:
```@repl lorenz
function lorenz_map(x, p)
    prob = ODEProblem{false}(lorenz, x, (0.0, -1.0), p)
    sol = solve(prob, Vern9(), abstol = 1e-10)
    sol[end]
end
```

Now we can create the problem:
```@repl lorenz
prob = VTwoDManifoldProblem(lorenz_map, para, d=1.0, amax=1.0, dsmin=1e-3)
```
For the meaning of these keyword arguments, please refer to [`VTwoDManifoldProblem`](@ref).

Similar to one-dimensional manifolds, we need a local manifold to start the extension. The corresponding function for creating a local manifold is [`gen_disk`](@ref):
```@repl lorenz
disk = gen_disk(saddle, r=1.0)
```
For detailed information about the `gen_disk` function, please refer to [`gen_disk`](@ref). Now we can proceed with the extension:
```@repl lorenz
manifold = growmanifold(prob, disk, 200)
```

We can also define a plotting function to visualize the results:
```@example lorenz
using CairoMakie
function manifold_plot(annulus)
    fig = Figure()
    axes = LScene(fig[1, 1], show_axis=false, scenekw=(backgroundcolor=:white, clear=true))
    second(x) = x[2]
    for i in eachindex(annulus)
        points = annulus[i].u
        lines!(axes, first.(points), second.(points), last.(points), fxaa=true)
    end
    fig
end
manifold_plot(manifold.data)
```


## Nonlinear Mapping

Consider the following nonlinear mapping:

```math
f(X)=\varphi\circ\Lambda\circ\varphi^{-1}(X)
```
where $\varphi(x,y,z)=(x,y,z-\alpha x^2-\beta y^2)$ is a nonlinear mapping, and $\Lambda$ is a diagonal matrix whose diagonal elements can be used to control the Jacobian matrix of mapping $f$ near the origin. Here's the code for computing its invariant manifold:
```@example nonlinearmap
using InvariantManifolds, LinearAlgebra, StaticArrays, OrdinaryDiffEq, CairoMakie
Λ = SDiagonal(SA[2.1, 6.3, 0.6])
φ(x, p)= SA[x[1],x[2],x[3]-p[1]*x[1]^2-p[2]*x[2]^2]
iφ(x, p)= SA[x[1],x[2],x[3]+p[1]*x[1]^2+p[2]*x[2]^2]
f(x,p) = φ(Λ*iφ(x, p),p)

para = [1.2,-1.2]
saddle = Saddle(SA[0.0, 0.0, 0.0], [SA[1.0, 0.0, 0.0], SA[0.0, 1.0, 0.0]], [2.1, 6.3])
prob = TwoDManifoldProblem(f, para, dcircle=0.05, d = 0.02, dsmin=1e-3)

disk = gen_disk(saddle, times=4, r= 0.05)
manifold = growmanifold(prob, disk, 3)
function manifold_plot(data)
    fig = Figure()
    axes = Axis3(fig[1,1])
    second(x) = x[2]
    for k in eachindex(data)
        for j in eachindex(data[k])
            points=data[k][j].u
            scatter!(axes,first.(points),second.(points),last.(points),fxaa=true)
        end
    end
    fig
end
manifold_plot(manifold.data)
```