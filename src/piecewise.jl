function _region_detect(regions, x, p, t)
    i = 1
    n = length(regions)
    while i <= n && regions[i](x, p, t) == false
        i = i + 1
    end
    if i > n
        error("Cannot detect $x's location")
    end
    i
end


"""
    setmap(v::T, timespan, alg, N, T; region_detect=_region_detect, extra...)

The function `setmap` is to get a `NSSetUp`.

# Parameters
- `v` a `ContinuousVectorField` or `JumpVectorField` like `PiecewiseV` or `BilliardV`.
- `timespan` the time span of the time-T-map.
- `alg` algorithm in `OrdinaryDiffEq` to solve ODE.
- `N` the dimension of the vector field.
- `T` number type used in computation. Usually, `Float64` is enough.

To ensure type stable, the numbers in `timespan` should be type of `T`.
The last two parameters has to be specified, since we have to store the event data.

# Keyword arguments
For `ContinuousVectorField` such as `PiecewiseV`, we have keyword argument:
- `region_detect=_region_detect` the region detect function to determine which domain the state in.

You can also pass the keywords of `solve` of [OrdinaryDiffEq](https://github.com/SciML/OrdinaryDiffEq.jl) to this function, 
except the `callback` and saving related keywords.
"""
function setmap(v::PiecewiseV, timespan, alg, N, T; region_detect=_region_detect, extra...)
    nn = length(v.hypers)
    event_at = Int[]
    event_state = SVector{N,T}[]
    event_t = T[]
    function affect!(integrator, idx)
        t0 = integrator.t + 1 // 20
        p = integrator.p
        u0 = integrator.sol(t0)
        i = region_detect(v.regions, u0, p, t0)
        integrator.f.f == i
        append!(event_at, [idx])
        append!(event_state, [integrator.u])
        append!(event_t, [integrator.t])
    end
    function condition(out, u, t, integrator)
        for i in eachindex(v.hypers)
            out[i] = v.hypers[i](u, integrator.p, t)
        end
    end
    vcb = VectorContinuousCallback(condition, affect!, nn)
    function tmap(State::State{N,T}, para) where {N,T}
        x = State.state
        v.n = region_detect(v.regions, x, para, timespan[1])
        prob = ODEProblem{false}(v, x, timespan, para)
        sol = solve(prob, alg, callback=vcb; extra...)
        newv_event_at = copy(event_at)
        newv_event_t = copy(event_t)
        newv_event_state = copy(event_state)
        empty!(event_at)
        empty!(event_t)
        empty!(event_state)
        NSState(sol[end], newv_event_t, newv_event_state, newv_event_at, State.s)
    end
    NSSetUp(v, timespan, tmap)
end

"""
    timetmap(v::T, para, timespan, alg; region_detect=_region_detect, extra...)

The function `timetmap` is similar to `setmap`. The output of this function is a function
which maps a `SVector` and parameters of ODE to a `SVector`, i.e. the time-T-map.

# Parameters
- `v` a `ContinuousVectorField` or `JumpVectorField` like `PiecewiseV` or `BilliardV`.
- `timespan` the time span of the time-T-map.
To ensure type stable, the numbers `timespan` should be type of float numbers you are using.
- `alg` algorithm in `OrdinaryDiffEq` to solve ODE.

# Keyword arguments
For `ContinuousVectorField` such as `PiecewiseV`, we have keyword argument:
- `region_detect=_region_detect` the region detect function to determine which domain the state in.


You can also pass the keywords of `solve` of [OrdinaryDiffEq](https://github.com/SciML/OrdinaryDiffEq.jl) to this function, 
except the `callback` and saving related keywords.
"""
function timetmap(v::PiecewiseV, timespan, alg ;region_detect=_region_detect, extra...)
    nn = length(v.hypers)
    function affect!(integrator, idx)
        t0 = integrator.t + 1 // 20
        p = integrator.p
        u0 = integrator.sol(t0)
        i = region_detect(v.regions, u0, p, t0)
        integrator.f.f == i
    end
    function condition(out, u, t, integrator)
        for i in eachindex(v.hypers)
            out[i] = v.hypers[i](u, integrator.p, t)
        end
    end
    vcb = VectorContinuousCallback(condition, affect!, nn)
    function tmap(x, para)
        v.n = region_detect(v.regions, x, para, timespan[1])
        prob = ODEProblem{false}(v, x, timespan, para)
        sol = solve(prob, alg, callback=vcb; extra...)
        sol[end]
    end
end

"""
    ns_solver(v::T, para, timespan, alg, N, T)

The function `ns_solver` is similar to `timetmap`. The output of this function is a function
which maps a `SVector` to a `NSSolution`. This `NSSolution` contain all data of an non-smooth ODE solution.

# Parameters
- `v` a `ContinuousVectorField` or `JumpVectorField` like `PiecewiseV` or `BilliardV`.
- `para` the parameter of the vector field.
- `timespan` the time span of the time-T-map.
- `alg` algorithm in `OrdinaryDiffEq` to solve ODE.
- `N` the dimension of the vector field.
- `T` number type used in computation.

To ensure type stable, the numbers in `para` and `timespan` should be type of `T`.
The last two parameters have to be specified, since we need to store the event data.
You can also pass the keywords of `solve` of [OrdinaryDiffEq](https://github.com/SciML/OrdinaryDiffEq.jl) to this function, 
except the `callback` and saving related keywords.
"""
function ns_solver(v::PiecewiseV, timespan, alg, N, T;region_detect=_region_detect, extra...)
    nn = length(v.hypers)
    event_at = Int[]
    event_state = SVector{N,T}[]
    event_t = T[]
    function affect!(integrator, idx)
        t0 = integrator.t + 1 // 20
        p = integrator.p
        u0 = integrator.sol(t0)
        i = region_detect(v.regions, u0, p, t0)
        integrator.f.f == i
        append!(event_at, [idx])
        append!(event_state, [integrator.u])
        append!(event_t, [integrator.t])
    end
    function condition(out, u, t, integrator)
        for i in eachindex(v.hypers)
            out[i] = v.hypers[i](u, integrator.p, t)
        end
    end
    vcb = VectorContinuousCallback(condition, affect!, nn)
    function tmap(x, para)
        v.n = region_detect(v.regions, x, para, timespan[1])
        prob = ODEProblem{false}(v, x, timespan, para)
        sol = solve(prob, alg, callback=vcb; extra...)
        newv_event_at = copy(event_at)
        newv_event_t = copy(event_t)
        newv_event_state = copy(event_state)
        empty!(event_at)
        empty!(event_t)
        empty!(event_state)
        NSSolution(sol, newv_event_t, newv_event_state, newv_event_at)
    end
end