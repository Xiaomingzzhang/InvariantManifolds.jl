```@setup smooth_one
using InvariantManifolds, StaticArrays, OrdinaryDiffEq, CairoMakie
```
# Getting Started: One-Dimensional Smooth Manifolds

## Nonlinear Mapping
Consider the following Henon map:

```math
\begin{aligned}
x'&=1-\alpha x^2+y,\\
y'&=\beta x,
\end{aligned}
```

where $\alpha,\beta$ are parameters. This mapping has fixed points:

```math
\begin{aligned}
(x_1,y_1)&=(\frac{-\sqrt{4 \alpha +\beta ^2-2 \beta +1}+\beta -1}{2 \alpha },\frac{1}{2} \left(\frac{\beta ^2}{\alpha }-\frac{\beta  \sqrt{4 \alpha +\beta ^2-2 \beta +1}}{\alpha }-\frac{\beta }{\alpha }\right)),\\
(x_2,y_2)&=(\frac{\sqrt{4 \alpha +\beta ^2-2 \beta +1}+\beta -1}{2 \alpha },\frac{1}{2} \left(\frac{\beta ^2}{\alpha }+\frac{\beta  \sqrt{4 \alpha +\beta ^2-2 \beta +1}}{\alpha }-\frac{\beta }{\alpha }\right)),
\end{aligned}
```

Let's calculate the eigenvalues of these two fixed points under the classical parameters $\alpha=1.4,\beta=0.3$:

```@example smooth_one
using StaticArrays, LinearAlgebra
function fixedpoints(p)
    a , b = p
    x1 = (-sqrt(4 * a + b^2 - 2 * b + 1) + b - 1) / (2 * a)
    y1 = (1 / 2) * (b^2 / a - b * sqrt(4 * a + b^2 - 2 * b + 1) / a - b / a)
    x2 = (sqrt(4 * a + b^2 - 2 * b + 1) + b - 1) / (2 * a)
    y2 = (1 / 2) * (b^2 / a + b * sqrt(4 * a + b^2 - 2 * b + 1) / a - b / a)
    return SA[x1, y1], SA[x2, y2]
end

function jacobian(x, p)
    a, b = p
    J = @SMatrix [-2 * a * x[1] 1.0; b 0.0]
    return J
end
```

```@repl smooth_one
eigen(jacobian(fixedpoints([1.4, 0.3])[1], [1.4, 0.3]))
```

```@repl smooth_one
eigen(jacobian(fixedpoints([1.4, 0.3])[2], [1.4, 0.3]))
```

As we can see, under the classical parameters, both fixed points are unstable. Next, we'll use the InvariantManifolds.jl package to compute one branch of the unstable manifold of the second fixed point.

The InvariantManifolds.jl package has an interface similar to many Julia packages. First, we need to load the package in Julia and define the Henon map:

```@repl smooth_one
using InvariantManifolds
function henonmap(x, p)
    y1 = 1 - p[1] * x[1]^2 + x[2]
    y2 = p[2] * x[1]
    SA[y1, y2]
end
```

Since the unstable eigenvalue at the saddle point is:
```@repl smooth_one
eigen(jacobian(fixedpoints([1.4, 0.3])[2], [1.4, 0.3])).values[1]
```
We need to iterate this mapping twice to ensure the manifold doesn't reverse during extension:

```@repl smooth_one
henonmap2(x, p)=henonmap(henonmap(x, p), p)
```

Now let's define a problem for computing a one-dimensional manifold of the smooth mapping:

```@repl smooth_one
para = [1.4, 0.3]
prob = OneDManifoldProblem(henonmap2, para)
```

To compute the manifold, we need a small local manifold segment starting at the saddle point. Usually, an unstable eigenvector starting at the saddle point with a very small length will suffice. InvariantManifolds.jl provides a function [`gen_segment`](@ref) to generate such a local manifold:

```@example smooth_one
saddle = fixedpoints(para)[2]
unstable_direction = eigen(jacobian(fixedpoints([1.4, 0.3])[2], [1.4, 0.3])).vectors[:,1]
segment = gen_segment(saddle, unstable_direction)
```

Under default keyword arguments, this function will generate a local manifold starting at the saddle point, with a length of 150 units and a step size of 0.01. Now let's use this local manifold to compute the smooth manifold:

```@repl smooth_one
manifold = growmanifold(prob, segment, 8)
```

This package doesn't provide plotting functionality. However, since the computation results are stored in `manifold.data`, and `manifold.data` is actually a vector whose elements are interpolation functions from the [DataInterpolations.jl](https://github.com/SciML/DataInterpolations.jl) package:
```@repl smooth_one
manifold.data
```

Therefore, we can define the following function to plot the smooth manifold:

```@example smooth_one
using CairoMakie
function manifold_plot(data)
    figure = Figure()
    axes = Axis(figure[1,1])
    for k in eachindex(data)
        points = Point2f.(data[k].u)
        lines!(axes, points)
    end
    figure
end
manifold_plot(manifold.data)
```

Full codes without comments:
```julia
using StaticArrays, LinearAlgebra, InvariantManifolds, CairoMakie
function fixedpoints(p)
    a , b = p
    x1 = (-sqrt(4 * a + b^2 - 2 * b + 1) + b - 1) / (2 * a)
    y1 = (1 / 2) * (b^2 / a - b * sqrt(4 * a + b^2 - 2 * b + 1) / a - b / a)
    x2 = (sqrt(4 * a + b^2 - 2 * b + 1) + b - 1) / (2 * a)
    y2 = (1 / 2) * (b^2 / a + b * sqrt(4 * a + b^2 - 2 * b + 1) / a - b / a)
    return SA[x1, y1], SA[x2, y2]
end
function jacobian(x, p)
    a, b = p
    J = @SMatrix [-2 * a * x[1] 1.0; b 0.0]
    return J
end
function henonmap(x, p)
    y1 = 1 - p[1] * x[1]^2 + x[2]
    y2 = p[2] * x[1]
    SA[y1, y2]
end
function henonmap2(x, p)
    henonmap(henonmap(x, p), p)
end
para = [1.4, 0.3]
prob = OneDManifoldProblem(henonmap2, para)
saddle = fixedpoints(para)[2]
unstable_direction = eigen(jacobian(fixedpoints([1.4, 0.3])[2], [1.4, 0.3])).vectors[:,1]
segment = gen_segment(saddle, unstable_direction)
manifold = growmanifold(prob, segment, 8)
function manifold_plot(data)
    figure = Figure()
    axes = Axis(figure[1,1])
    for k in eachindex(data)
        points = Point2f.(data[k].u)
        lines!(axes, points)
    end
    figure
end
manifold_plot(manifold.data)
```

## Oscillator with Periodic Forcing

Now let's consider a higher-order example. Consider the following oscillator with periodic forcing:
```math
\begin{aligned}
\dot{x}&=y,\\
\dot{y}&=x-\delta x^3+\gamma \cos(\omega t).
\end{aligned}
```

When $\gamma=0$, the system has a saddle point at $(0,0)$. After a small periodic perturbation, this saddle point becomes a saddle periodic orbit, which is a saddle point of the mapping $T:X\mapsto \phi(X,2\pi/\omega,0)$, where $\phi(X,t,t_0)$ is the solution of the system under the initial condition $X(t_0)=X\in\mathbb{R}^2$. Fortunately, we can obtain the Jacobian matrix of the mapping $T$ using the solution of the variational equation. The saddle point position and unstable direction of the mapping $T$ can also be obtained through numerical methods.

InvariantManifolds.jl provides a function [`findsaddle`](@ref) to obtain the saddle point position and unstable direction of $T$. Here's the code demonstrating how to use this function:

```@example smooth_one
using InvariantManifolds, LinearAlgebra, StaticArrays, OrdinaryDiffEq, CairoMakie
f(x, p, t) = SA[x[2], x[1] - p[1]*(x[1]^3) + p[2]*cos(p[3]*t)]
df(x, p, t) = SA[0.0 1.0; 1-p[1]*3*(x[1]^2) 0.0]
initial_guess = SA[0.0, 0.0]
para = [1.0, 0.1, 2.2]
timespan = (0.0, 2pi/para[3])
saddle = findsaddle(f, df, timespan, initial_guess, para)
```

The `gen_segment` function can act directly on the [`Saddle`](@ref) structure. Therefore, we can use the following code to generate a local manifold:
```@repl smooth_one
segment = gen_segment(saddle)
```

Now we can define the nonlinear mapping:
```@repl smooth_one
function timeTmap(x, p)
    prob = ODEProblem{false}(f, x, (0.0, 2pi/p[3]), p)
    solve(prob, Vern9(), abstol=1e-10)[end]
end
```

Then create the problem and solve it:
```@repl smooth_one
prob = OneDManifoldProblem(timeTmap, para)
manifold = growmanifold(prob, segment, 7)
```
Finally, use the function defined in the previous section to plot the results:
```@example smooth_one
manifold_plot(manifold.data)
```

Full codes without comments:
```julia
using StaticArrays, LinearAlgebra, InvariantManifolds, CairoMakie
f(x, p, t) = SA[x[2], x[1] - p[1]*(x[1]^3) + p[2]*cos(p[3]*t)]
df(x, p, t) = SA[0.0 1.0; 1-p[1]*3*(x[1]^2) 0.0]
initial_guess = SA[0.0, 0.0]
para = [1.0, 0.1, 2.2]
timespan = (0.0, 2pi/para[3])
saddle = findsaddle(f, df, timespan, initial_guess, para)
segment = gen_segment(saddle)
function timeTmap(x, p)
    prob = ODEProblem{false}(f, x, (0.0, 2pi/p[3]), p)
    solve(prob, Vern9(), abstol=1e-10)[end]
end
prob = OneDManifoldProblem(timeTmap, para)
manifold = growmanifold(prob, segment, 7)
function manifold_plot(data)
    figure = Figure()
    axes = Axis(figure[1,1])
    for k in eachindex(data)
        points = Point2f.(data[k].u)
        lines!(axes, points)
    end
    figure
end
manifold_plot(manifold.data)
```