function setmap(v::SFilippovV, timespan, alg, N, T; extra...)
    f1 = v.fs[1]
    f2 = v.fs[2]
    event_at = Int[]
    event_state = SVector{N,T}[]
    event_t = T[]
    function affect!(integrator)
        v1 = integrator.u
        t1 = integrator.t
        p0 = integrator.p
        append!(event_state, [v1])
        append!(event_t, [t1])
        n = Int(p0[end])
        if n <= 2
            if dot(v.dhyper(v1, p0, t1), f1(v1, p0, t1)) > 0 && dot(v.dhyper(v1, p0, t1), f2(v1, p0, t1)) < 0
                integrator.f.f.n = 3
                append!(event_at, [1])
            elseif dot(v.dhyper(v1, p0, t1), f1(v1, p0, t1)) < 0 && dot(v.dhyper(v1, p0, t1), f2(v1, p0, t1)) < 0
                integrator.f.f.n = 1
                append!(event_at, [2])
            else
                integrator.f.f.n = 2
                append!(event_at, [2])
            end
        else # exit slide surface
            append!(event_at, [3])
            if abs(dot(v.dhyper(v1, p0, t1), f1(v1, p0, t1))) < 1e-8
                integrator.f.f.n = 1
            else
                integrator.f.f.n = 2
            end
        end
    end
    function condition(u, t, integrator)
        n = integrator.f.f.n
        if n > 2
            v.exit(u, integrator.p, t)
        elseif n <= 2
            v.hyper(u, integrator.p, t)
        end
    end
    cb = ContinuousCallback(condition, affect!)
    function tmap(State::State{N,T}, para) where {N,T}
        x = State.state
        if v.hyper(x, para, timespan[1]) < 0
            v.n = 1
        elseif v.hyper(x, para, timespan[1]) > 0
            v.n = 2
        elseif dot(v.dhyper(x, para, timespan[1]), f1(x, para, timespan[1])) > 0 && dot(v.dhyper(x, para, timespan[1]), f2(x, para, timespan[1])) < 0
            v.n = 3
        elseif dot(v.dhyper(x, para, timespan[1]), f1(x, para, timespan[1])) < 0 && dot(v.dhyper(x, para, timespan[1]), f2(x, para, timespan[1])) > 0
            error("The solution is not defined at $x")
        end
        prob = ODEProblem{false}(v, x, timespan, para)
        sol = solve(prob, alg, callback=cb; extra...)
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


function timetmap(v::SFilippovV, timespan, alg; extra...)
    f1 = v.fs[1]
    f2 = v.fs[2]
    function affect!(integrator)
        v1 = integrator.u
        t1 = integrator.t
        p0 = integrator.p
        n = Int(p0[end])
        if n <= 2
            if dot(v.dhyper(v1, p0, t1), f1(v1, p0, t1)) > 0 && dot(v.dhyper(v1, p0, t1), f2(v1, p0, t1)) < 0
                integrator.f.f.n = 3
            elseif dot(v.dhyper(v1, p0, t1), f1(v1, p0, t1)) < 0 && dot(v.dhyper(v1, p0, t1), f2(v1, p0, t1)) < 0
                integrator.f.f.n = 1
            else
                integrator.f.f.n = 2
            end
        else # exit slide surface
            if abs(dot(v.dhyper(v1, p0, t1), f1(v1, p0, t1))) < 1e-8
                integrator.f.f.n = 1
            else
                integrator.f.f.n = 2
            end
        end
    end
    function condition(u, t, integrator)
        n = v.n
        if n > 2
            v.exit(u, integrator.p, t)
        elseif n <= 2
            v.hyper(u, integrator.p, t)
        end
    end
    cb = ContinuousCallback(condition, affect!)
    function tmap(x, para)
        if v.hyper(x, para, timespan[1]) < 0
            v.n = 1
        elseif v.hyper(x, para, timespan[1]) > 0
            v.n = 2
        elseif dot(v.dhyper(x, para, timespan[1]), f1(x, para, timespan[1])) > 0 && dot(v.dhyper(x, para, timespan[1]), f2(x, para, timespan[1])) < 0
            v.n = 3
        elseif dot(v.dhyper(x, para, timespan[1]), f1(x, para, timespan[1])) < 0 && dot(v.dhyper(x, para, timespan[1]), f2(x, para, timespan[1])) > 0
            error("The solution is not defined at $x")
        end
        prob = ODEProblem{false}(v, x, timespan, para)
        sol = solve(prob, alg, callback=cb; extra...)
        sol[end]
    end
end

function ns_solver(v::SFilippovV, timespan, alg, N, T; extra...)
    event_at = Int[]
    event_state = SVector{N,T}[]
    event_t = T[]
    function affect!(integrator)
        f1 = v.fs[1]
        f2 = v.fs[2]
        v1 = integrator.u
        t1 = integrator.t
        p0 = integrator.p
        append!(event_state, [v1])
        append!(event_t, [t1])
        n = Int(p0[end])
        if n <= 2
            if dot(v.dhyper(v1, p0, t1), f1(v1, p0, t1)) > 0 && dot(v.dhyper(v1, p0, t1), f2(v1, p0, t1)) < 0
                integrator.f.f.n = 3
                append!(event_at, [1])
            elseif dot(v.dhyper(v1, p0, t1), f1(v1, p0, t1)) < 0 && dot(v.dhyper(v1, p0, t1), f2(v1, p0, t1)) < 0
                integrator.f.f.n = 1
                append!(event_at, [2])
            else
                integrator.f.f.n = 2
                append!(event_at, [2])
            end
        else # exit slide surface
            append!(event_at, [3])
            if abs(dot(v.dhyper(v1, p0, t1), f1(v1, p0, t1))) < 1e-8
                integrator.f.f.n = 1
            else
                integrator.f.f.n = 1
            end
        end
    end
    function condition(u, t, integrator)
        n = v.n
        if n > 2
            v.exit(u, integrator.p, t)
        elseif n <= 2
            v.hyper(u, integrator.p, t)
        end
    end
    cb = ContinuousCallback(condition, affect!)
    function tmap(x, para)
        f1 = v.fs[1]
        f2 = v.fs[2]
        if v.hyper(x, para, timespan[1]) < 0
            v.n = 1
        elseif v.hyper(x, para, timespan[1]) > 0
            v.n = 2
        elseif dot(v.dhyper(x, para, timespan[1]), f1(x, para, timespan[1])) > 0 && dot(v.dhyper(x, para, timespan[1]), f2(x, para, timespan[1])) < 0
            v.n = 3
        elseif dot(v.dhyper(x, para, timespan[1]), f1(x, para, timespan[1])) < 0 && dot(v.dhyper(x, para, timespan[1]), f2(x, para, timespan[1])) > 0
            error("The solution is not defined at $x")
        end
        prob = ODEProblem{false}(v, x, timespan, para)
        sol = solve(prob, alg, callback=cb; extra...)
        newv_event_at = copy(event_at)
        newv_event_t = copy(event_t)
        newv_event_state = copy(event_state)
        empty!(event_at)
        empty!(event_t)
        empty!(event_state)
        NSSolution(sol, newv_event_t, newv_event_state, newv_event_at)
    end
end
