function setmap(v::BilliardV, timespan, alg; repeat_nudge=1 // 100, extra...)
    nn = length(v.hypers)
    event_at = Int[]
    function affect!(integrator, idx)
        p0 = integrator.p
        u0 = integrator.u
        t0 = integrator.t
        integrator.u = v.rules[idx](u0, p0, t0)
        append!(event_at, [idx])
    end
    function condition(out, u, t, integrator)
        for i in eachindex(v.hypers)
            out[i] = v.hypers[i](u, integrator.p, integrator.t)
        end
    end
    vcb = VectorContinuousCallback(condition, affect!, nn; repeat_nudge=repeat_nudge)
    function tmap(X::NSState{N,T}, para) where {N,T}
        x = X.state
        event = copy(X.event_at)
        prob = ODEProblem{false}(v, x, timespan, para)
        sol = solve(prob, alg, callback=vcb; extra...)
        newv_event_at = copy(event_at)
        append!(event, newv_event_at)
        empty!(event_at)
        NSState(sol[end], event)
    end
    NSSetUp(v, timespan, tmap)
end


function ns_solver(v::BilliardV, timespan, alg, N, T; repeat_nudge=1 // 100, extra...)
    nn = length(v.hypers)
    event_at = Int[]
    event_state = SVector{N,T}[]
    event_t = T[]
    function affect!(integrator, idx)
        p0 = integrator.p
        u0 = integrator.u
        integrator.u = v.rules[idx](u0, p0, integrator.t)
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