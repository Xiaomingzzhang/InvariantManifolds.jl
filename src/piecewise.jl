struct PiecewiseV
    fs::Vector{Function}
    regions::Vector{Function}
    hypers::Vector{Function}
end

function (v::PiecewiseV)(x, p, t)
    n = convert(Int, p[end])
    v.fs[n](x, p, t)
end

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
    timetmap(vecs,hypers,map,timespan)
This is the function to generate a time-T-map from a piecewise smooth ODE system.
The output of this function is a function, whose input is the state in state space,
output is the image of the time-T-map and the integer vector to represent the cross times with
hypersurfaces.
- `vecs`, vector of ODE functions, with type f(x,p,t);
- `hypers`, hypersurfaces of switch surface; if `vecs=[f1,f2,f3]` and `hypers=[h1,h2,h3]`,
it is assumed that `f1` is at h1<0, `f2` is at 0<h1 and h2<0, `f3` is at h2>0;
- `map`, vector of mapping function on hypersurfaces;
- `timespan`, the timespan of time-T-map.
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