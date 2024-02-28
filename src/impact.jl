"""
    BilliardV

The struct `BilliardV` has three fields:
- `f` is the vector field, of type `f(x,p,t)`, and its output is a SVector;
- `hypers` is vector of hypersurfaces:`[h1,h2,...]`, `h1(x,p,t)`;
- `irules` is vector of rules on hypersurfaces:`[r1,r2,r3,...]`; the rules must be independent of time `t`, i.e., `r1(x,p)`.
"""
struct BilliardV
    f
    hypers::Vector{Function}
    irules::Vector{Function}
end

function (v::BilliardV)(x, p, t)
    v.f(x, p, t)
end

function setmap(v::BilliardV, para, timespan, alg; extra...)
    nn = length(v.hypers)
    ctimes = zeros(Int, nn)
    function affect!(integrator, idx)
        p0 = integrator.p
        u0 = integrator.u
        u0 .= v.irules[idx](u0, p0)
        ctimes[idx] = ctimes[idx] + 1
    end
    function condition(out, u, t, integrator)
        for i in eachindex(v.hypers)
            out[i] = v.hypers[i](u, integrator.p, t)
        end
    end
    vcb = VectorContinuousCallback(condition, affect!, nn)
    function tmap(State)
        x = State.state
        ss = State.s
        prob = ODEProblem(v, x, timespan, para)
        sol = solve(prob, alg, callback=vcb; extra...)
        newctimes = Vector{Int}(undef, nn)
        newctimes .= ctimes
        ctimes .= zeros(Int, nn)
        NSState(sol[end], newctimes, ss)
    end
    NSSetUp(v, para, timespan, tmap, alg)
end

function gen_prob(v::BilliardV, x, para, timespan)
    nn = length(v.hypers)
    function affect!(integrator, idx)
        p0 = integrator.p
        u0 = integrator.u
        type = typeof(u0[1])
        n0 = length(u0)
        newu0 = Vector{type}(undef, n0)
        newu0 .= v.irules[idx](u0, p0)
        integrator.u = SVector{n0}(newu0)
    end
    function condition(out, u, t, integrator)
        for i in eachindex(v.hypers)
            out[i] = v.hypers[i](u, integrator.p, t)
        end
    end
    vcb = VectorContinuousCallback(condition, affect!, nn)
    ODEProblem(v, x, timespan, para, callback=vcb)
end