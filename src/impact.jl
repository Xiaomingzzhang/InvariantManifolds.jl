function setmap(v::BilliardV, timespan, alg, N, T; extra...)
    nn = length(v.hypers)
    event_at = Int[]
    function affect!(integrator, idx)
        p0 = integrator.p
        u0 = integrator.u
        integrator.u = v.rules[idx](u0, p0)
        append!(event_at, [idx])
    end
    function condition(out, u, t, integrator)
        for i in eachindex(v.hypers)
            out[i] = v.hypers[i](u, integrator.p)
        end
    end
    vcb = VectorContinuousCallback(condition, affect!, nn)
    function tmap(X::NSState{N,T}, para) where {N,T}
        x = X.state
        event = copy(X.event_at)
        prob = ODEProblem{false}(v, x, timespan, para)
        sol = solve(prob, alg, callback=vcb; extra...)
        newv_event_at = copy(event_at)
        append!(event, newv_event_at)
        empty!(event_at)
        NSState(sol[end], event, false, 0, 0)
    end
    NSSetUp(v, timespan, tmap)
end

function timetmap(v::BilliardV, timespan, alg;extra...)
    nn = length(v.hypers)
    function affect!(integrator, idx)
        p0 = integrator.p
        u0 = integrator.u
        integrator.u = v.rules[idx](u0, p0, integrator.t)
    end
    function condition(out, u, t, integrator)
        for i in eachindex(v.hypers)
            out[i] = v.hypers[i](u, integrator.p, t)
        end
    end
    vcb = VectorContinuousCallback(condition, affect!, nn)
    function tmap(x,para)
        prob = ODEProblem{false}(v, x, timespan, para)
        sol = solve(prob, alg, callback=vcb; extra...)
        sol[end]
    end
end

function ns_solver(v::BilliardV, timespan, alg, N, T;extra...)
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
    vcb = VectorContinuousCallback(condition, affect!, nn)
    function tmap(x,para)
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