# Curve thinning and saltation (monodromy) matrices

This page documents two extensions that are used by the one-dimensional
continuation algorithms:

- curvature-aware point thinning after every continuation step;
- computation of the time-map derivative with jump (saltation) matrices for
  non-smooth systems.

## Thinning over-resolved curve points

During a continuation the manifold often becomes over-resolved: many almost
collinear points remain in regions that were refined because of an earlier sharp
turn. [`thin_curve`](@ref) removes such points after every continuation step
without evaluating the map. The keyword `thinning=true` is the default in
`grow!`/`growmanifold`; use `thinning=false` to keep all points.

```julia
manifold = growmanifold(prob, segment, N)        # thinning enabled by default
manifold = growmanifold(prob, segment, N; thinning=false)
```

A point is deleted only when it is sufficiently close to the chord connecting its
neighbours and when deleting it does not increase the local turning angle beyond
`amax`. Endpoints of smooth components and genuine sharp turns are kept.

## Scale-aware angle-times-length criterion

The default acceptance rule uses the Goodman–Wróbel style condition

```math
\alpha \le \alpha^*, \qquad \alpha\,\delta_1 \le \frac{d^*\alpha^*}{10}.
```

The second condition prevents a large raw angle on an extremely short segment from
being recorded as a numerical flaw. The bound can be changed with the `adbnd`
keyword of `grow!`/`growmanifold`:

```julia
manifold = growmanifold(prob, segment, N; adbnd=1e-4)
```

## Monodromy matrices and jump matrices

For a non-smooth system the derivative of the time-`T`-map cannot be obtained by
simply integrating the variational equation, because at every impact/switching
event the linearized flow has to be multiplied by a saltation (jump) matrix.

[`ns_monodromy`](@ref) solves the augmented system

```math
\dot{x} = f(x,t),\qquad \dot{X} = Df(x,t)X,
```

and applies the jump matrix at each event. For an impact system (BilliardV):

```julia
saddle = findsaddle(vectorfield, df, timespan, initialguess, para)
DF = ns_monodromy(vectorfield, df, timespan, initialguess, para)
```

For a piecewise smooth or combined system the derivative argument must be a tuple
containing the Jacobians of all smooth vector fields:

```julia
dvs = (df1, df2)
DF = ns_monodromy(vectorfield, dvs, timespan, x0, para)
saddle = findsaddle(vectorfield, dvs, timespan, x0, para)
```

The monodromy matrix is used by Newton's method in the `findsaddle` overloads so
that saddle periodic orbits which collide with impact surfaces are located
accurately.
