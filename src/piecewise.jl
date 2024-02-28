"""
    PiecewiseV

Piecewise smooth vector field. Fields:
- `fs` is a vector of smooth vector fields in different regions.
- `regions` is a vector of the region functions: `[r1,r2,...]`, where `r1(x,p,t)` should return a Bool value to indicate that `x` is in this region or not.
- `hypers` is a vector of the hypersurfaces separating the regions.
"""
struct PiecewiseV
    fs::Vector{Function}
    regions::Vector{Function}
    hypers::Vector{Function}
end

function (v::PiecewiseV)(x, p, t)
    n = convert(Int, p[end])
    v.fs[n](x, p, t)
end


"""
    NSSetUp{T}

`NSSetUp` is a struct to contatin all the information needed in continuing the manifold. Fields:
- `f::T` the Non-smooth vector field, like `PiecewiseV`;
- `p` the parameter;
- `timespan` the time span of time-T-map;
- `timetmap` the time-t-map of nonsmooth ODE, which maps a `State` to a `NSState`;
- `alg` algorithm used to solve ODE.
"""
struct NSSetUp{T}
    f::T
    p
    timespan
    timetmap
    alg
end

function show(io::IO, NS::NSSetUp{T}) where {T}
    p = NS.p
    timespan = NS.timespan
    print(io,
        "NSSetUp{$T}
parameter:$p
timespan:$timespan")
end

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
    setmap(v::PiecewiseV, para, timespan, alg; region_detect=_region_detect, extra...)

Function to get a `NSSetUp`.
"""
function setmap(v::PiecewiseV, para, timespan, alg; region_detect=_region_detect, extra...)
    nn = length(v.hypers)
    ctimes = zeros(Int, nn)
    function affect!(integrator, idx)
        t0 = integrator.t + 1 // 20
        p = integrator.p
        u0 = integrator.sol(t0)
        i = region_detect(v.regions, u0, p, t0 + 1 // 20)
        p[end] = i
        ctimes[idx] = ctimes[idx] + 1
    end
    function condition(out, u, t, integrator)
        for i in eachindex(v.hypers)
            out[i] = v.hypers[i](u, integrator.p, t)
        end
    end
    vcb = VectorContinuousCallback(condition, affect!, nn)
    function tmap(State::State{N,T}) where {N,T}
        x = State.state
        ss = State.s
        para[end] = region_detect(v.regions, x, para, timespan[1])
        prob = ODEProblem(v, x, timespan, para)
        sol = solve(prob, alg, callback=vcb; extra...)
        newctimes = Vector{Int}(undef, nn)
        newctimes .= ctimes
        ctimes .= zeros(Int, nn)
        NSState(sol[end], newctimes, ss)
    end
    NSSetUp(v, para, timespan, tmap, alg)
end

function gen_prob(v::PiecewiseV, x, para, timespan; region_detect=_region_detect)
    nn = length(v.hypers)
    function affect!(integrator, idx)
        t0 = integrator.t + 1 // 20
        p = integrator.p
        u0 = integrator.sol(t0)
        i = region_detect(v.regions, u0, p, t0 + 1 // 20)
        p[end] = i
    end
    function condition(out, u, t, integrator)
        for i in eachindex(v.hypers)
            out[i] = v.hypers[i](u, integrator.p, t)
        end
    end
    vcb = VectorContinuousCallback(condition, affect!, nn)
    para[end] = region_detect(v.regions, x, para, timespan[1])
    ODEProblem(v, x, timespan, para, callback=vcb)
end