
"""
    NSOneDManifoldProblem{F,T}

`NSOneDManifoldProblem` is a struct to contain the main information for continuing the non-smooth one-dimensional manifold of
the time-T-map of a non-smooth ODE.
# Fields
- `f` the struct [`NSSetUp`](@ref);
- `para` the parameters of the nonlinear map;
- `amax` the maximum angle between points when continuing the manifold;
- `d` the maximum distance between points when continuing the manifold;
- `ϵ` the max value of the following expression: ``\\max\\{|H(x_0,T)|,H(x_1,T)\\},`` where ``H(x,t)`` is the hypersurface the manifold cross, ``x_0`` and ``x_1`` are points before and after the cross, ``T`` is the end of the time-``T``-map (from 0 to ``T``).
- `dsmin` the minimum arc length allowing; note that if in a continuation point, this value is achieved and the angle as well as the distance values are not achieved, then we will record this point as a [`FlawPoint`](@ref).

Convenient consturctors are `NSOneDManifoldProblem(f)` and `NSOneDManifoldProblem(f,para)`
"""
struct NSOneDManifoldProblem{F,T}
    f::F
    para::Vector{T}
    amax::T
    d::T
    ϵ::T
    dsmin::T
end

function NSOneDManifoldProblem(f; amax=0.5, d=0.001, ϵ=0.00001, dsmin=0.0001)
    NSOneDManifoldProblem(f, Float64[], amax, d, ϵ, dsmin)
end


function NSOneDManifoldProblem(f, para::AbstractVector{T};
    amax=T(0.5), d=T(0.001), ϵ=T(0.00001), dsmin=T(0.0001)) where {T}
    NSOneDManifoldProblem(f, para, amax, d, ϵ, dsmin)
end

function show(io::IO, m::MIME"text/plain", A::NSOneDManifoldProblem)
    printstyled(io, "NSOneDManifoldProblem:"; color=:cyan)
    println(io)
    print(io, "f: NSSetUp generated by non-smooth vector field: ")
    show(io, A.f.f)
    println(io)
    print(io, "para: ")
    show(io, A.para)
    println(io)
    print(io, "amax: ")
    show(io, m, A.amax)
    println(io)
    print(io, "d: ")
    show(io, m, A.d)
    println(io)
    print(io, "ϵ: ")
    show(io, m, A.ϵ)
    println(io)
    print(io, "dsmin: ")
    show(io, m, A.dsmin)
end

"""
    NSOneDManifold{F,S,N,T}

`NSOneDManifold` is a struct contains all the information of the non-smooth one-dimensional numerical manifold.
# Fields
- `prob` the problem [`NSOneDManifoldProblem`](@ref);
- `data` the numerical data that should be `Vector{Vector{Vector{S}}}`, where `S` is the interpolation curve (we use `DataInterpolation` in this package);
- `flawpoints` the flaw points generated during continuation.
"""
mutable struct NSOneDManifold{F,S,N,T}
    prob::NSOneDManifoldProblem{F,T}
    data::Vector{Vector{S}}
    flawpoints::Vector{FlawPoint{N,T}}
end

function Base.show(io::IO, m::MIME"text/plain", A::NSOneDManifold)
    m = 0
    n = 0
    for i in eachindex(A.data)
        n = n + length(A.data[i])
        for j in eachindex(A.data[i])
            m = m + length(A.data[i][j].t)
        end
    end
    k = length(A.flawpoints)
    printstyled(io, "Non-smooth one-dimensional manifold"; bold=true, color=:cyan)
    println(io)
    printstyled(io, "Curves number: "; color=:cyan)
    println(io, "$n")
    printstyled(io, "Points number: "; color=:cyan)
    println(io, "$m")
    if  k == 0
        printstyled(io, "Flaw points number: "; color=:cyan)
        print(io, "0")
    else
        amax = A.prob.amax
        d = A.prob.d
        prend = findall(x -> x.d > d, A.flawpoints)
        nd = length(prend)
        prenc = findall(x -> x.α > amax, A.flawpoints)
        nc = length(prenc)
        printstyled(io, "Flaw points number: "; color=:cyan)
        println(io, "$k")
        printstyled(io, "Distance failed points number: "; color=:cyan)
        println(io, "$nd")
        printstyled(io, "Curvature failed points number: "; color=:cyan)
        println(io, "$nc")
        dα = maximum([x.α*x.d for x in A.flawpoints])
        printstyled(io, "Max dα in Flaw Points: "; color=:cyan)
        print(io, "$dα")
    end
end

function partition(v::Vector{NSState{N,T}}, s0::Vector{T}; interp=QuadraticInterpolation) where {N,T}
    ctime1 = v[1].event_at
    n = length(v)
    a0 = NSState{N,T}[]
    result = Vector{NSState{N,T}}[]
    result_s = Vector{T}[]
    s00 = T[]
    for j in 1:n
        if v[j].event_at == ctime1
            append!(a0, [v[j]])
            append!(s00, [s0[j]])
        else
            append!(result, [a0])
            append!(result_s, [s00])
            a0 = [v[j]]
            s00 = [s0[j]]
            ctime1 = v[j].event_at
        end
    end
    append!(result, [a0])
    append!(result_s, [s00])
    for i in eachindex(result_s)
        first_s = result_s[i][1]
        for j in eachindex(result_s[i])
            result_s[i][j] = result_s[i][j] - first_s
        end
    end
    [paramise(result[i], result_s[i], interp=interp) for i in eachindex(result)]
end

take_state(x) = x.state
take_u(x) = x.u
take_t(x) = x.t


function union_same_event(v::Vector{S}; interp=QuadraticInterpolation) where {S}
    if length(v) == 1
        v[1]
    else
        totals = v[1].t
        totalu = vcat(take_u.(v)...)
        l = length(v)
        T = typeof(totals[1])
        s0 = T(0)
        for i in 1:l-1
            last_s = v[i].t[end]
            add_s = norm(v[i].u[end] - v[i+1].u[1])
            s0 = s0 + last_s + add_s
            final_s = s0 .+ v[i+1].t
            append!(totals, final_s)
        end
        paramise(totalu, totals, interp=interp)
    end
end

# Union broken lines with the same event
function union_lines(v::Vector{S}; interp=QuadraticInterpolation) where {S}
    result = S[]
    event = v[1].u[1].event_at
    line = S[]
    for j in eachindex(v)
        if (v[j].u[1].event_at) == event
            append!(line, [v[j]])
        else
            append!(result, [union_same_event(line, interp=interp)])
            line = [v[j]]
            event = v[j].u[1].event_at
        end
    end
    append!(result, [union_same_event(line, interp=interp)])
    result
end

function ispartitioned(v::Vector{NSState{N,T}}) where {N,T}
    ctime1 = v[1].event_at
    n = length(v)
    result = true
    j = 1
    while j <= n && result
        result = v[j].event_at == ctime1
        j = j + 1
    end
    result
end


"""
    InvariantManifolds.ns_addpoints!(tmap, p, d, dsmin, oldcurve, newu, olds, αmax, tend, hypers, ϵ, flawpoints) -> Vector{T}

Add points to ensure proper spacing and accuracy when computing the non-smooth one-dimensional manifold.

# Arguments
- `tmap`: Time map function that evolves states forward
- `p`: Vector of parameters
- `d`: Maximum allowed distance between consecutive points
- `dsmin`: Minimum allowed arc length between points
- `oldcurve`: Previous curve data used for interpolation
- `newu`: Vector of new states to be processed
- `olds`: Vector of arc length parameters
- `αmax`: Maximum allowed angle between consecutive segments
- `tend`: End time of the time map
- `hypers`: Vector of hypersurface functions
- `ϵ`: Maximum allowed error in hypersurface intersection
- `flawpoints`: Vector to store problematic points encountered

# Returns
- Vector of arc length parameters for the processed points

# Details
The function adaptively adds points to maintain:
1. Maximum distance `d` between consecutive points
2. Maximum angle `αmax` between segments
3. Accuracy `ϵ` at hypersurface intersections

If constraints cannot be satisfied within `dsmin`, points are marked as flaws.
"""
@inline function ns_addpoints!(tmap, p, d, dsmin, oldcurve,
    newu::Vector{NSState{N,T}}, olds::Vector{T}, αmax, tend, hypers, ϵ, flawpoints) where {N,T}
    n = length(newu)
    i = 1
    newpara = T[0]
    oldcurve_event = oldcurve.u[1].event_at
    @inbounds while i + 1 <= n
        if i + 2 <= n
            u0 = newu[i]
            u1 = newu[i+1]
            u2 = newu[i+2]
            if u0.event_at == u1.event_at == u2.event_at
                δ1 = norm(u0 - u1)
                δ2 = norm(u1 - u2)
                if δ1 == 0 || (δ1 != 0 && δ2 == 0)
                    deleteat!(newu, i + 1)
                    deleteat!(olds, i + 1)
                    n = n - 1
                else
                    baru0 = u1 + (u1 - u2) * δ1 / δ2
                    α = norm(baru0 - u0) / δ1
                    if δ1 <= d && α <= αmax
                        i = i + 1
                        dd = newpara[end]
                        append!(newpara, [dd + δ1])
                    else
                        if olds[i+1] - olds[i] > dsmin
                            s0 = olds[i]
                            s1 = olds[i+1]
                            paras = (s0 + s1) / 2
                            addps = tmap(NSState(oldcurve(paras), copy(oldcurve_event)), p)
                            insert!(newu, i + 1, addps)
                            insert!(olds, i + 1, paras)
                            n = n + 1
                        else
                            i = i + 1
                            append!(flawpoints, [FlawPoint(u0.state, α, δ1)])
                            dd = newpara[end]
                            append!(newpara, [dd + δ1])
                        end
                    end
                end
            elseif u0.event_at == u1.event_at && u1.event_at != u2.event_at
                δ = norm(u0 - u1)
                if δ <= d
                    i = i + 1
                    dd = newpara[end]
                    append!(newpara, [dd + δ])
                else
                    if olds[i+1] - olds[i] > dsmin
                        s0 = olds[i]
                        s1 = olds[i+1]
                        paras = (s0 + s1) / 2
                        addps = tmap(NSState(oldcurve(paras), copy(oldcurve_event)), p)
                        insert!(newu, i + 1, addps)
                        insert!(olds, i + 1, paras)
                        n = n + 1
                    else
                        i = i + 1
                        append!(flawpoints, [FlawPoint(u0.state, T(0), δ)])
                        dd = newpara[end]
                        append!(newpara, [dd + δ])
                    end
                end
            elseif u0.event_at == u1.event_at[1:end-1]
                idx = u1.event_at[end]
                ϵ0 = max(abs(hypers[idx](u0, p, tend)), abs(hypers[idx](u1, p, tend)))
                if ϵ0 <= ϵ
                    i = i + 1
                    dd = newpara[end]
                    append!(newpara, [dd + ϵ0])
                else
                    if olds[i+1] - olds[i] > 1e-12
                        s0 = olds[i]
                        s1 = olds[i+1]
                        paras = (s0 + s1) / 2
                        addps = tmap(NSState(oldcurve(paras), copy(oldcurve_event)), p)
                        insert!(newu, i + 1, addps)
                        insert!(olds, i + 1, paras)
                        n = n + 1
                    else
                        error("Cannot locate the intersection points at desired accurcy. Please increase the value of ϵ")
                    end
                end
            elseif u0.event_at[1:end-1] == u1.event_at
                idx = u0.event_at[end]
                ϵ0 = max(abs(hypers[idx](u0, p, tend)), abs(hypers[idx](u1, p, tend)))
                if ϵ0 <= ϵ
                    i = i + 1
                    dd = newpara[end]
                    append!(newpara, [dd + ϵ0])
                else
                    if olds[i+1] - olds[i] > 1e-12
                        s0 = olds[i]
                        s1 = olds[i+1]
                        paras = (s0 + s1) / 2
                        addps = tmap(NSState(oldcurve(paras), copy(oldcurve_event)), p)
                        insert!(newu, i + 1, addps)
                        insert!(olds, i + 1, paras)
                        n = n + 1
                    else
                        error("Cannot locate the intersection points. Please increase the value of ϵ")
                    end
                end
            else
                if olds[i+1] - olds[i] > 10eps()
                    s0 = olds[i]
                    s1 = olds[i+1]
                    paras = (s0 + s1) / 2
                    addps = tmap(NSState(oldcurve(paras), copy(oldcurve_event)), p)
                    insert!(newu, i + 1, addps)
                    insert!(olds, i + 1, paras)
                    n = n + 1
                else
                    error("The manifold cross two events in an almost zero arc length")
                end
            end
        else
            u0 = newu[i]
            u1 = newu[i+1]
            if u0.event_at == u1.event_at
                δ = norm(u0 - u1)
                if δ <= d
                    i = i + 1
                    dd = newpara[end]
                    append!(newpara, [dd + δ])
                else
                    if olds[i+1] - olds[i] > dsmin
                        s0 = olds[i]
                        s1 = olds[i+1]
                        paras = (s0 + s1) / 2
                        addps = tmap(NSState(oldcurve(paras), copy(oldcurve_event)), p)
                        insert!(newu, i + 1, addps)
                        insert!(olds, i + 1, paras)
                        n = n + 1
                    else
                        i = i + 1
                        append!(flawpoints, [FlawPoint(u0.state, T(0), δ)])
                        dd = newpara[end]
                        append!(newpara, [dd + δ])
                    end
                end
            elseif u0.event_at == u1.event_at[1:end-1]
                idx = u1.event_at[end]
                ϵ0 = max(abs(hypers[idx](u0, p, tend)), abs(hypers[idx](u1, p, tend)))
                if ϵ0 <= ϵ
                    i = i + 1
                    dd = newpara[end]
                    append!(newpara, [dd + ϵ0])
                else
                    if olds[i+1] - olds[i] > 1e-12
                        s0 = olds[i]
                        s1 = olds[i+1]
                        paras = (s0 + s1) / 2
                        addps = tmap(NSState(oldcurve(paras), copy(oldcurve_event)), p)
                        insert!(newu, i + 1, addps)
                        insert!(olds, i + 1, paras)
                        n = n + 1
                    else
                        error("Cannot locate the intersection points. Please increase the value of ϵ")
                    end
                end
            elseif u0.event_at[1:end-1] == u1.event_at
                idx = u0.event_at[end]
                ϵ0 = max(abs(hypers[idx](u0, p, tend)), abs(hypers[idx](u1, p, tend)))
                if ϵ0 <= ϵ
                    i = i + 1
                    dd = newpara[end]
                    append!(newpara, [dd + ϵ0])
                else
                    if olds[i+1] - olds[i] > 1e-12
                        s0 = olds[i]
                        s1 = olds[i+1]
                        paras = (s0 + s1) / 2
                        addps = tmap(NSState(oldcurve(paras), copy(oldcurve_event)), p)
                        insert!(newu, i + 1, addps)
                        insert!(olds, i + 1, paras)
                        n = n + 1
                    else
                        error("Cannot locate the intersection points. Please increase the value of ϵ")
                    end
                end
            else
                if olds[i+1] - olds[i] > 1e-12
                    s0 = olds[i]
                    s1 = olds[i+1]
                    paras = (s0 + s1) / 2
                    addps = tmap(NSState(oldcurve(paras), copy(u0.event_at)), p)
                    insert!(newu, i + 1, addps)
                    insert!(olds, i + 1, paras)
                    n = n + 1
                else
                    error("The manifold cross two events in an almost zero arc length")
                end
            end
        end
    end
    newpara
end


function initialize(prob::NSOneDManifoldProblem, seg::Vector{SVector{N,T}}; interp=QuadraticInterpolation, event=Int[]) where {N,T}
    parameters = prob.para
    tmap = prob.f.timetmap
    d = prob.d
    αmax = prob.amax
    tend = prob.f.timespan[end]
    hypers = prob.f.f.hypers
    dsmin = prob.dsmin
    ϵ = prob.ϵ
    flawpoints = FlawPoint{N,T}[]
    n = length(seg)
    points = Vector{NSState{N,T}}(undef, n)
    for i in eachindex(points)
        points[i] = NSState(seg[i], event)
    end
    result = [[paramise(points, interp=interp)]]
    n = length(points)
    image = Vector{NSState{N,T}}(undef, n)
    for i in eachindex(image)
        image[i] = tmap(points[i], parameters)
    end
    oldcurve = result[1][1]
    olds = copy(oldcurve.t)
    ns_addpoints!(tmap, parameters, d, dsmin, oldcurve,
        image, olds, αmax, tend, hypers, ϵ, flawpoints)
    if ispartitioned(image) == false
        error("The initial curve has to be chosen more small")
    end
    p0 = first(points)
    pn = last(points)
    j = 1
    dd = norm(pn - p0)
    while norm(image[j] - p0) < dd
        j = j + 1
    end
    j = j + 1
    first_iteration_curve = image[j:end]
    event = copy(first_iteration_curve[1].event_at)
    prepend!(first_iteration_curve, [NSState(pn.state, event)])
    append!(result, [[paramise(first_iteration_curve, interp=interp)]])
    NSOneDManifold(prob, result, flawpoints)
end

function grow!(manifold::NSOneDManifold{F,S,N,T}; interp=QuadraticInterpolation) where {F,S,N,T}
    αmax = manifold.prob.amax
    d = manifold.prob.d
    v = manifold.prob.f
    ϵ = manifold.prob.ϵ
    dsmin = manifold.prob.dsmin
    para = manifold.prob.para
    tmap = v.timetmap
    tend = v.timespan[end]
    hypers = v.f.hypers
    data = manifold.data
    nic = data[end]
    flawpoints = manifold.flawpoints
    n = length(manifold.data) + 1
    result = empty(manifold.data[1])
    @inbounds for j in eachindex(nic)
        curve = nic[j]
        n = length(curve.u)
        ic2_states = Vector{NSState{N,T}}(undef, n)
        for k in eachindex(ic2_states)
            ic2_states[k] = tmap(curve.u[k], para)
        end
        olds = copy(curve.t)
        newpara = ns_addpoints!(tmap, para, d, dsmin, curve,
            ic2_states, olds, αmax, tend, hypers, ϵ, flawpoints)
        _result = partition(ic2_states, newpara, interp=interp)
        append!(result, _result)
    end
    append!(data, [union_lines(result, interp=interp)])
end

function growmanifold(prob::NSOneDManifoldProblem, segment, N; interp=QuadraticInterpolation)
    manifold = initialize(prob, segment, interp=interp)
    for i in 1:N
        grow!(manifold, interp=interp)
    end
    manifold
end