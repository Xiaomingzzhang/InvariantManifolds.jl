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
    setmap(v, timespan::Tuple{T,T}, alg) where {T}

The function `setmap` is to get a `NSSetUp`.

# Parameters
- `v` a nonsmooth vector field like [`PiecewiseV`](@ref) or [`BilliardV`](@ref).
- `timespan` the time span of the time-T-map.
- `alg` algorithm in `OrdinaryDiffEq` to solve ODE.

To ensure type stable, the numbers in `timespan` should be type of `T`.

# Keyword arguments
- `repeat_nudge=1//100` this is used to set the next testing point after a previously found zero.

For vector fields [`PiecewiseV`](@ref) and [`PiecewiseImpactV`](@ref), we have two special keyword arguments:
- `cross_time= 0.01` when the solution `sol` hits the hypersurface at time `t`, we need to know which domain it enters. We choose the state `sol(t+cross_time)` to determine which domain it enters.
- `region_detect=_region_detect` the region detect function to determine which domain the state in.


You can also pass the keywords of `solve` of [OrdinaryDiffEq](https://github.com/SciML/OrdinaryDiffEq.jl) to this function, 
except the `callback` and saving related keywords.
"""
function setmap(v::PiecewiseV, timespan::Tuple{T,T}, alg;
    cross_time=0.01, region_detect=_region_detect, repeat_nudge=1//100, extra...) where {T}
    nn = length(v.hypers)
    event_at = Int[]
    function affect!(integrator, idx)
        t0 = integrator.t + cross_time
        p = integrator.p
        u0 = integrator.sol(t0)
        i = region_detect(v.regions, u0, p, t0)
        integrator.f.f.n = i
        append!(event_at, [idx])
    end
    function condition(out, u, t, integrator)
        for i in eachindex(v.hypers)
            out[i] = v.hypers[i](u, integrator.p, t)
        end
    end
    vcb = VectorContinuousCallback(condition, affect!, nn; repeat_nudge=repeat_nudge)
    function tmap(X::NSState{N,T}, para) where {N,T}
        x = X.state
        event = copy(X.event_at)
        v.n = region_detect(v.regions, x, para, timespan[1])
        prob = ODEProblem{false}(v, x, timespan, para)
        sol = solve(prob, alg, callback=vcb; extra...)
        newv_event_at = copy(event_at)
        append!(event, newv_event_at)
        empty!(event_at)
        NSState(sol[end], newv_event_at)
    end
    NSSetUp(v, timespan, tmap)
end


"""
    ns_solver(v::T, para, timespan, alg, N, T)

The function `ns_solver` is similar to `timetmap`. The output of this function is a function
which maps a `SVector` to a [`NSSolution`](@ref). This `NSSolution` contain all data of an non-smooth ODE solution.

# Parameters
- `v` vector fields like [`PiecewiseV`](@ref) or [`BilliardV`](@ref).
- `para` the parameter of the vector field.
- `timespan` the time span of the time-T-map.
- `alg` algorithm in `OrdinaryDiffEq` to solve ODE.
- `N` the dimension of the vector field.
- `T` number type used in computation.

To ensure type stable, the numbers in `para` and `timespan` should be type of `T`.
The last two parameters have to be specified, since we need to store the event data.

# Keyword arguments
- `repeat_nudge=1//100` this is used to set the next testing point after a previously found zero.

For vector fields [`PiecewiseV`](@ref) and [`PiecewiseImpactV`](@ref), we have two special keyword arguments:
- `cross_time= 0.01` when the solution `sol` hits the hypersurface at time `t`, we need to know which domain it enters. We choose the state `sol(t+cross_time)` to determine which domain it enters.
- `region_detect=_region_detect` the region detect function to determine which domain the state in.

You can also pass the keywords of `solve` of [OrdinaryDiffEq](https://github.com/SciML/OrdinaryDiffEq.jl) to this function, 
except the `callback` and saving related keywords.
"""
function ns_solver(v::PiecewiseV, timespan, alg, N, T;
    cross_time=0.01, region_detect=_region_detect, repeat_nudge=1//100, extra...)
    nn = length(v.hypers)
    event_at = Int[]
    event_state = SVector{N,T}[]
    event_t = T[]
    function affect!(integrator, idx)
        t0 = integrator.t + cross_time
        p = integrator.p
        u0 = integrator.sol(t0)
        i = region_detect(v.regions, u0, p, t0)
        integrator.f.f.n = i
        append!(event_at, [idx])
        append!(event_state, [integrator.u])
        append!(event_t, [integrator.t])
    end
    function condition(out, u, t, integrator)
        for i in eachindex(v.hypers)
            out[i] = v.hypers[i](u, integrator.p, t)
        end
    end
    vcb = VectorContinuousCallback(condition, affect!, nn; repeat_nudge=repeat_nudge)
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

"""
    iscontact(setup, saddle, para)

Check if the saddle point of the time-T-map has contacts with the hypersurfaces.

# Parameters
- `setup` a [`NSSetUp`](@ref).
- `saddle` a [`Saddle`](@ref) or a `SVector`.
- `para` the parameter of the vector field.

# Returns
- `true` if the saddle point has contact with the hypersurfaces.
- `false` if the saddle point has no contact with the hypersurfaces.
"""
function iscontact(setup, saddle::Saddle{N,T,S}, para) where {N,T,S}
    event = setup.timetmap(NSState(saddle.saddle), para).event_at
    !(isempty(event))
end

function iscontact(setup, saddle::SVector{N,T}, para) where {N,T}
    event = setup.timetmap(NSState(saddle), para).event_at
    !(isempty(event))
end