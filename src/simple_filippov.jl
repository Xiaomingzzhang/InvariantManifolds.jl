"""
    SFilippovV

`SFilippovV` means simple Fillippov vector fields, which means that there only exists one hypersurface to separate the phase space. Fields:
- `fs` vector fields in two sides of hypersurface. The slide vector field can be generated automatically.
- `hyper` hypersurface.
- `dhyper` grad of hypersurface. Warn!!! The grad must point to the second of `fs`.
- `exit` conditions to exit the hypersurface, which can also be generated automatically
"""
struct SFilippovV
    fs::SVector{3,Function}
    hyper
    dhyper
    exit
end

function SFilippovV(fs, h, ∇h)
    f1 = fs[1]
    f2 = fs[2]
    function sv(x, p, t)
        α = dot(∇h(x, p, t), f1(x, p, t)) / dot(∇h(x, p, t), f1(x, p, t) - f2(x, p, t))
        (1 - α) * f1(x, p, t) + α * f2(x, p, t)
    end
    function exit(x, p, t)
        dot(∇h(x, p, t), f2(x, p, t)) * dot(∇h(x, p, t), f1(x, p, t))
    end
    SFilippovV(SA[f1, f2, sv], h, ∇h, exit)
end

function (v::SFilippovV)(x, p, t)
    n = convert(Int, p[end])
    v.fs[n](x, p, t)
end

function gen_prob(v::SFilippovV, x, timespan, para)
    f1 = v.fs[1]
    f2 = v.fs[2]
    if v.hyper(x, para, timespan[1]) < 0
        para[end] = 1
    elseif v.hyper(x, para, timespan[1]) > 0
        para[end] = 2
    elseif dot(v.dhyper(x, para, timespan[1]), f1(x, para, timespan[1])) > 0 && dot(v.dhyper(v1, p0, t1), f2(x, para, timespan[1])) < 0
        para[end] = 3
    elseif dot(v.dhyper(x, para, timespan[1]), f1(x, para, timespan[1])) < 0 && dot(v.dhyper(v1, p0, t1), f2(x, para, timespan[1])) > 0
        error("The solution is not defined at $x")
    end
    function affect!(integrator)
        v1 = integrator.u
        t1 = integrator.t
        p0 = integrator.p
        n = convert(Int, p0[end])
        if n <= 2
            if dot(v.dhyper(v1, p0, t1), f1(v1, p0, t1)) > 0 && dot(v.dhyper(v1, p0, t1), f2(v1, p0, t1)) < 0
                p0[end] = 3
            elseif dot(v.dhyper(v1, p0, t1), f1(v1, p0, t1)) < 0 && dot(v.dhyper(v1, p0, t1), f2(v1, p0, t1)) < 0
                p0[end] = 1
            else
                p0[end] = 2
            end
        else # exit slide surface
            if abs(dot(v.dhyper(v1, p0, t1), f1(v1, p0, t1))) < 1e-8
                p0[end] = 1
            else
                p0[end] = 2
            end
        end
    end
    function condition(u, t, integrator)
        n = convert(Int, integrator.p[end])
        if n > 2
            v.exit(u, integrator.p, t)
        elseif n <= 2
            v.hyper(u, integrator.p, t)
        end
    end
    cb = ContinuousCallback(condition, affect!)
    ODEProblem(v, x, timespan, para, callback=cb)
end



function setmap(v::SFilippovV, para, timespan, alg; extra...)
    ctimes = zeros(Int, 3)
    f1 = v.fs[1]
    f2 = v.fs[2]
    function affect!(integrator)
        v1 = integrator.u
        t1 = integrator.t
        p0 = integrator.p
        n = convert(Int, p0[end])
        if n <= 2
            if dot(v.dhyper(v1, p0, t1), f1(v1, p0, t1)) > 0 && dot(v.dhyper(v1, p0, t1), f2(v1, p0, t1)) < 0
                p0[end] = 3
                ctimes[3] = ctimes[3] + 1
            elseif dot(v.dhyper(v1, p0, t1), f1(v1, p0, t1)) < 0 && dot(v.dhyper(v1, p0, t1), f2(v1, p0, t1)) < 0
                p0[end] = 1
                if n==1
                    nothing
                else
                    ctimes[1] = ctimes[1] + 1
                end
            else
                p0[end] = 2
                if n==2
                    nothing
                else
                    ctimes[2] = ctimes[2] + 1
                end
            end
        else # exit slide surface
            if abs(dot(v.dhyper(v1, p0, t1), f1(v1, p0, t1))) < 1e-8
                p0[end] = 1
                ctimes[1] = ctimes[1] + 1
            else
                p0[end] = 2
                ctimes[2] = ctimes[2] + 1
            end
        end
    end
    function condition(u, t, integrator)
        n = convert(Int,p[end])
        if n > 2
            v.exit(u, integrator.p, t)
        elseif n <= 2
            v.hyper(u, integrator.p, t)
        end
    end
    cb = ContinuousCallback(condition, affect!)
    function tmap(State::State{N,T}) where {N,T}
        x = State.state
        ss = State.s
        if v.hyper(x, para, timespan[1]) < 0
            para[end] = 1
        elseif v.hyper(x, para, timespan[1]) > 0
            para[end] = 2
        elseif dot(v.dhyper(x, para, timespan[1]), f1(x, para, timespan[1])) > 0 && dot(v.dhyper(v1, p0, t1), f2(x, para, timespan[1])) < 0
            para[end] = 3
        elseif dot(v.dhyper(x, para, timespan[1]), f1(x, para, timespan[1])) < 0 && dot(v.dhyper(v1, p0, t1), f2(x, para, timespan[1])) > 0
            error("The solution is not defined at $x")
        end
        prob = ODEProblem(v, x, timespan, para)
        sol = solve(prob, alg, callback=cb; extra...)
        newctimes = Vector{Int}(undef, 3)
        newctimes .= ctimes
        ctimes .= zeros(Int, 3)
        NSState(sol[end], newctimes, ss)
    end
    NSSetUp(v, para, timespan, tmap, alg)
end