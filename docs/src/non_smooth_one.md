# One-Dimensional Non-Smooth Manifolds
Perhaps the most notable feature of this package is its ability to compute non-smooth manifolds. Currently, it supports the computation of non-smooth manifolds for two types of systems:
- One-dimensional manifold computation for time-periodic non-smooth differential equations
- Two-dimensional manifold computation for non-smooth autonomous systems

The manifolds in both cases are invariant manifolds of saddle points. The former requires taking time-periodic mappings, while the latter requires taking fixed-step time-T-mappings. The non-smooth factors in these two types of systems can be diverse, including piecewise smooth functions, collisions, and combinations thereof. Users don't need to solve these three types of non-smooth systems themselves, as we provide three encapsulated structures:
- [`PiecewiseV`](@ref)
- [`BilliardV`](@ref)
- [`PiecewiseImpactV`](@ref)

We use the Callback functionality from [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) to compute time mappings. The following three examples will demonstrate the methods for computing invariant manifolds of these three types of non-smooth systems.

!!! warning 
    The computation of non-smooth manifolds heavily depends on the ODE solving algorithm and precision. When the solution fails or performs poorly, try changing the algorithm, increasing the ODE solving precision, or reducing the values of parameters `ϵ`, `d`, `amax` in [`NSOneDManifoldProblem`](@ref).


## Piecewise Smooth Systems
Consider a simple piecewise smooth system:

```math
\begin{aligned}
\dot{x}&=y,\\
\dot{y}&=f(x) + \epsilon \sin(2\pi t),
\end{aligned}
```
where

```math
f(x) =
\begin{cases}
-k_1 x& \text{if } x < -d,\\
k_2 x & \text{if } -d<x<d,\\
-k_3 x& \text{if } x > d.
\end{cases}
```

$k_1,k_2,k_3,d>0$ are all positive constants. We will compute the invariant manifold of the time-periodic mapping. Note that when the periodic perturbation is small, the saddle point should be close to the origin.

First, let's load the required packages
```@setup piecewise
using InvariantManifolds, LinearAlgebra, StaticArrays, OrdinaryDiffEq, CairoMakie
```
```@repl piecewise
using InvariantManifolds, LinearAlgebra, StaticArrays, OrdinaryDiffEq, CairoMakie
```
Next, define the piecewise smooth vector field:
```@example piecewise
f1(x, p, t) = SA[x[2], p[1]*x[1]+p[4]*sin(2pi * t)]
f2(x, p, t) = SA[x[2], -p[2]*x[1]+p[4]*sin(2pi * t)]
f3(x, p, t) = SA[x[2], -p[3]*x[1]+p[4]*sin(2pi * t)]

hyper1(x, p, t) = x[1] - p[5]
hyper2(x, p, t) = x[1] + p[5]

dom1(x, p, t) = -p[5] < x[1] < p[5]
dom2(x, p, t) = x[1] > p[5]
dom3(x, p, t) = x[1] < -p[5]

vectorfield = PiecewiseV((f1, f2, f3), (dom1, dom2, dom3), (hyper1, hyper2))
```
The parameters passed to the `PiecewiseV` structure are: vector fields, their domains, and the hyperplanes that divide these domains. For more details, refer to [`PiecewiseV`](@ref).

Next, we'll encapsulate the key information for solving the time-periodic mapping in another structure [`NSSetUp`](@ref):

```@repl piecewise
setup = setmap(vectorfield, (0.0, 1.0), Tsit5(), abstol=1e-8)
```
The function [`setmap`](@ref) is used to encapsulate the time mapping computation information. Now we have defined everything needed to solve the time-periodic mapping.

Next, to generate the local manifold, we also need to locate the saddle point and its unstable eigenvector. Let's set the parameters:
```@repl piecewise
para = para = [2, 5, 5, 0.6, 2]
```
Since the perturbation is small, the saddle-type periodic orbit should still be in `dom1`. Therefore, we can use `findsaddle` to calculate the position of the saddle point:
```@example piecewise
function df1(x, p, t)
    SA[0 1; p[1] 0]
end
initialguess = SA[0.0, 0.0]
saddle = findsaddle(f1, df1, (0.0,1.0), initialguess, para, abstol=1e-10)
```

Next, create the problem, generate the local manifold, and perform the extension
```@repl piecewise
prob = NSOneDManifoldProblem(setup, para, ϵ = 1e-3)
segment = gen_segment(saddle)
manifold = growmanifold(prob, segment, 8)
```

Note that the data type of `manifold.data` is `Vector{Vector{S}}`, where `S` is an interpolation function. So we need to use the following function to plot the results:
```@example piecewise
using CairoMakie
function manifold_plot(data)
    fig = Figure()
    axes = Axis(fig[1,1])
    for k in eachindex(data)
        for j in eachindex(data[k])
            points=data[k][j].u
            lines!(axes,first.(points),last.(points))
        end
    end
    fig
end
manifold_plot(manifold.data)
```


## Impact Systems

Consider the following forced inverted pendulum equation:

```math
\begin{aligned}
\dot{x}&= y,\\
\dot{y}&= \sin(x) - \epsilon \cos(2\pi t),
\end{aligned}
```
Assume there are walls on both sides of the inverted pendulum: when $x=\xi$ or $x=-\xi$, we have $y\rightarrow - ry$. Similarly, we need to first construct the non-smooth vector field
```@setup impact
using InvariantManifolds, LinearAlgebra, StaticArrays, OrdinaryDiffEq, CairoMakie
function manifold_plot(data)
    fig = Figure()
    axes = Axis(fig[1,1])
    for k in eachindex(data)
        for j in eachindex(data[k])
            points=data[k][j].u
            lines!(axes,first.(points),last.(points))
        end
    end
    fig
end
```
```@example impact
f(x, p, t) = SA[x[2], sin(x[1])-p[1]*cos(2 * pi * t)]

hyper1(x, p, t) = x[1] + p[2]
hyper2(x, p, t) = x[1] - p[2]

rule1(x, p, t) = SA[x[1], -p[3]*x[2]]
rule2(x, p, t) = SA[x[1], -p[3]*x[2]]

vectorfield = BilliardV(f, (hyper1, hyper2), (rule1, rule2))
```
Next, encapsulate the information for solving the time-periodic mapping:
```@example impact
setup = setmap(vectorfield, (0.0, 1.0), Vern9(), abstol=1e-10)
```
Find the saddle point:
```@example impact
function df(x, p, t)
    SA[0 1; cos(x[1]) 0]
end
para = [0.2, pi / 4, 0.98]
initialguess = SA[0.0, 0.0]
saddle = findsaddle(f, df, (0.0,1.0), initialguess, para, abstol=1e-10)
```
Next, create the problem, generate the local manifold, and perform the extension
```@repl impact
prob = NSOneDManifoldProblem(setup, para)
segment = gen_segment(saddle)
manifold = growmanifold(prob, segment, 11)
```

Finally, plot the results: 
```@example impact
manifold_plot(manifold.data)
```

## Combined Piecewise Smooth and Impact ODE Systems
```@setup piecewiseimpact
using InvariantManifolds, LinearAlgebra, StaticArrays, OrdinaryDiffEq, CairoMakie
function manifold_plot(data)
    fig = Figure()
    axes = Axis(fig[1,1])
    for k in eachindex(data)
        for j in eachindex(data[k])
            points=data[k][j].u
            lines!(axes,first.(points),last.(points))
        end
    end
    fig
end
```
Now consider an ODE system with both piecewise smooth functions and impacts:
```math
\begin{aligned}
\dot{x}&=y,\\
\dot{y}&=f(x) + \epsilon \sin(2\pi t),
\end{aligned}
```
where

```math
f(x) =
\begin{cases}
-k_1 x& \text{if } x < -d,\\
k_2 x & \text{if } -d<x<d
\end{cases}
```

$k_1,k_2,d>0$ are all positive constants. When $x=d$, we have $\dot{y}->-r\dot{y}$. We will compute the invariant manifold of the time-periodic mapping. Note that when the periodic perturbation is small, the saddle point should be close to the origin. First, let's load the required packages
```@setup piecewiseimpact
using InvariantManifolds, LinearAlgebra, StaticArrays, OrdinaryDiffEq, CairoMakie
```
```@repl piecewiseimpact
using InvariantManifolds, LinearAlgebra, StaticArrays, OrdinaryDiffEq, CairoMakie
```
Next, define the non-smooth vector field:
```@example piecewiseimpact
f1(x, p, t) = SA[x[2], p[1]*x[1]+p[3]*sin(2pi * t)]
f2(x, p, t) = SA[x[2], -p[2]*x[1]+p[3]*sin(2pi * t)]

hyper1(x, p, t) = x[1] - p[4]
hyper2(x, p, t) = x[1] + p[4]

dom1(x, p, t) = -p[4] < x[1] < p[4]
dom2(x, p, t) = x[1] < -p[4]

impact_rule(x, p, t) = SA[x[1], -p[5]*x[2]]
id(x,p,t) = x

vectorfield = PiecewiseImpactV((f1, f2), (dom1, dom2), (hyper1, hyper2), (impact_rule, id), [1])
```
The parameters passed to the `PiecewiseImpactV` structure are: vector fields, their domains, the hyperplanes that divide these domains, the rules that act on the hyperplanes, and a list of rules that have impact effects. For more details, refer to [`PiecewiseImpactV`](@ref).

Next, we'll encapsulate the key information for solving the time-periodic mapping in another structure [`NSSetUp`](@ref):

```@repl piecewiseimpact
setup = setmap(vectorfield, (0.0, 1.0), Tsit5(), abstol=1e-8, reltol=1e-8)
```
The function [`setmap`](@ref) is used to encapsulate the time mapping computation information. Now we have defined everything needed to solve the time-periodic mapping.

Next, to generate the local manifold, we also need to locate the saddle point and its unstable eigenvector. Let's set the parameters:
```@repl piecewiseimpact
para = [2, 5, 0.6, 2, 0.98]
```
Since the perturbation is small, the saddle-type periodic orbit should still be in `dom1`. Therefore, we can use `findsaddle` to calculate the position of the saddle point:
```@example piecewiseimpact
function df1(x, p, t)
    SA[0 1; p[1] 0]
end
initialguess = SA[0.0, 0.0]
saddle = findsaddle(f1, df1, (0.0,1.0), initialguess, para, abstol=1e-10)
```

Next, create the problem, generate the local manifold, and perform the extension
```@repl piecewiseimpact
prob = NSOneDManifoldProblem(setup, para)
segment = gen_segment(saddle)
manifold = growmanifold(prob, segment, 9)
```

Note that the data type of `manifold.data` is `Vector{Vector{S}}`, where `S` is an interpolation function. So we need to use the following function to plot the results:
```@example piecewiseimpact
using CairoMakie
function manifold_plot(data)
    fig = Figure()
    axes = Axis(fig[1,1])
    for k in eachindex(data)
        for j in eachindex(data[k])
            points=data[k][j].u
            lines!(axes,first.(points),last.(points))
        end
    end
    fig
end
manifold_plot(manifold.data)
```




