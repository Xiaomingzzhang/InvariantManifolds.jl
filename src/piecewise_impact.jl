function setmap(v::PiecewiseImpactV, timespan, alg;
    cross_time=0.01, region_detect=_region_detect, repeat_nudge=1 // 100, extra...)
    nn = length(v.hypers)
    event_at = Int[]
    function affect!(integrator, idx)
        if idx in v.idxs
            integrator.u = v.rules[idx](integrator.u, integrator.p, integrator.t)
        else
            t0 = integrator.t + cross_time
            u0 = integrator.sol(t0)
            p = integrator.p
            integrator.f.f.n = region_detect(v.regions, u0, p, t0)
        end
        push!(event_at, idx)
    end
    function condition(out, u, t, integrator)
        for i in eachindex(v.hypers)
            out[i] = v.hypers[i](u, integrator.p, t)
        end
    end
    vcb = VectorContinuousCallback(condition, affect!, nn, repeat_nudge=repeat_nudge)
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

function ns_solver(v::PiecewiseImpactV, timespan, alg, N, T;
    cross_time=0.01, region_detect=_region_detect, repeat_nudge=1 // 100, extra...)
    nn = length(v.hypers)
    event_at = Int[]
    event_state = SVector{N,T}[]
    event_t = T[]
    function affect!(integrator, idx)
        if idx in v.idxs
            integrator.u = v.rules[idx](integrator.u, integrator.p, integrator.t)
        else
            t0 = integrator.t + cross_time
            u0 = integrator.sol(t0)
            p = integrator.p
            integrator.f.f.n = region_detect(v.regions, u0, p, t0)
        end
        push!(event_at, idx)
        push!(event_state, integrator.u)
        push!(event_t, integrator.t)
    end
    function condition(out, u, t, integrator)
        for i in eachindex(v.hypers)
            out[i] = v.hypers[i](u, integrator.p, t)
        end
    end
    vcb = VectorContinuousCallback(condition, affect!, nn, repeat_nudge=repeat_nudge)
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