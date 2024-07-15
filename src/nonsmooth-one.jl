function partition(v::Vector{NSState{N,T}}, s0::Vector{T}; interp=LinearInterpolation) where {N,T}
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
    [interp(result[i], result_s[i]) for i in eachindex(result)]
end

take_state(x) = x.state
take_u(x) = x.u
take_t(x) = x.t


function union_same_event(v::Vector{S}; interp=LinearInterpolation) where {S<:AbstractInterpolation}
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
        @show length(totalu)
        @show length(totals)
        interp(totalu, totals)
    end
end

# Union broken lines with the same event
function union_lines(v::Vector{S}; interp=LinearInterpolation) where {S<:AbstractInterpolation}
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

# function myreplace!(a::Array{T,1}, i::Integer, b::Array{S,1}) where {T,S}
#     # Throw convert error before changing the shape of the array
#     _b = S == T ? b : convert(Array{T,1}, b)::Array{T,1}
#     n = length(b)
#     Base._growat!(a, i, n - 1)
#     @inbounds for j in 1:n
#         a[i+j-1] = _b[j]
#     end
#     return a
# end

@inline function nsnorm(a::NSState{N,T}, b::NSState{N,T}, hypers, p, time_end, dtimes) where {N,T}
    if a.event_at == b.event_at
        norm(a - b)
    elseif abs(length(a.event_at) - length(b.event_at)) == 1
        if a.event_at == b.event_at[1:end-1]
            idx = b.event_at[end]
            max(abs(hypers[idx](a.state, p, time_end)), abs(hypers[idx](b.state, p, time_end))) * dtimes
        elseif a.event_at[1:end-1] == b.event_at
            idx = a.event_at[end]
            max(abs(hypers[idx](a.state, p, time_end)), abs(hypers[idx](b.state, p, time_end))) * dtimes
        else
            T(2)
        end
    else
        T(2)
    end
end


@inline function ns_addpoints!(tmap, p, δ, oldcurve, newu::Vector{NSState{N,T}}, olds::Vector{T}, dtimes, tend, hypers) where {N,T}
    n = length(newu)
    i = 1
    newpara = T[0]
    event = oldcurve.u[1].event_at
    @inbounds while i + 1 <= n
        dist = nsnorm(newu[i], newu[i+1], hypers, p, tend, dtimes)
        if dist > δ
            # @show dist
            m = ceil(Int, dist / δ)
            if m == 1
                m = 2
            end
            s0 = olds[i]
            s1 = olds[i+1]
            plengh = (s1 - s0) / m
            paras = [s0 + plengh * i for i in 1:m-1]
            addps = Vector{NSState{N,T}}(undef, m - 1)
            for j in 1:m-1
                addps[j] = tmap(NSState(oldcurve(paras[j]), copy(event)), p)
            end
            insert!(newu, i + 1, addps)
            insert!(olds, i + 1, paras)
            n = n + m - 1
        else
            if newu[i].event_at == newu[i+1].event_at
                newdist = dist
            else
                newdist = dist / dtimes
            end
            dd = newpara[end]
            append!(newpara, [dd + newdist])
            i = i + 1
        end
    end
    newpara
end

function del_extra!(newpara, newu, min)
    k = 1
    M = length(newpara)
    @inbounds while k + 2 <= M
        s1 = newpara[k+1] - newpara[k]
        s2 = newpara[k+2] - newpara[k+1]
        s3 = norm(newu[k+2] - newu[k])
        if s1 < min && s2 < min && s3 < min
            deleteat!(newpara, k + 1)
            deleteat!(newu, k + 1)
            M = M - 1
        else
            k = k + 1
        end
    end
end





"""
    InvariantManifolds.ns_initialise_curve(points, tmap, para)

Initialise the curve of non-smooth ODE's time-T-map.
"""
function ns_initialise_curve(f::NSSetUp, para, saddle, direction, nn, d, δ; abstol=1e-4, interp=LinearInterpolation)
    points = NSState.(segment(saddle, direction, nn, d))
    result = [[paramise(points, interp=interp)]]
    tmap = f.timetmap
    tend = f.timespan[end]
    hypers = f.f.hypers
    dtimes = δ / abstol
    n = length(points)
    m = length(saddle)
    type = typeof(saddle[1][1])
    image = Vector{NSState{m,type}}(undef, n)
    for i in eachindex(image)
        image[i] = tmap(points[i], para)
    end
    oldcurve = result[1][1]
    olds = copy(result[1][1].t)
    ns_addpoints!(tmap, para, δ, oldcurve, image, olds, dtimes, tend, hypers)
    if ispartitioned(image) == false
        error("The initial curve has to be chosen more small")
    end
    p0 = first(points)
    pn = last(points)
    j = 1
    d = norm(pn - p0)
    while norm(image[j] - p0) < d
        j = j + 1
    end
    j = j + 1
    first_iteration_curve = image[j:end]
    prepend!(first_iteration_curve, [pn])
    append!(result, [[paramise(first_iteration_curve, interp=interp)]])
end

function grow_line!(v::NSSetUp, para, data::Vector{Vector{P}}, δ;
    abstol=1e-4, interp=LinearInterpolation) where {P}
    nic = data[end]
    result = P[]
    N = length(nic[1].u[1])
    T = typeof(nic[1].u[1][1])
    tmap = v.timetmap
    tend = v.timespan[end]
    hypers = v.f.hypers
    dtimes = δ / abstol
    @inbounds for j in eachindex(nic)
        curve = nic[j]
        n = length(curve.u)
        ic2_states = Vector{NSState{N,T}}(undef, n)
        for k in eachindex(ic2_states)
            ic2_states[k] = tmap(curve.u[k], para)
        end
        olds = copy(curve.t)
        newpara = ns_addpoints!(tmap, para, δ, curve, ic2_states, olds, dtimes, tend, hypers)
        _result = partition(ic2_states, newpara)
        # insert left zeros and right zeros
        for i in eachindex(_result)
            del_extra!(_result[i].t, _result[i].u, δ)
        end
        # delete extra points
        append!(result, _result)
    end
    append!(data, [union_lines(result, interp=interp)])
    # append!(data, [result])
end


function grow_linetest!(v::NSSetUp, para, data::Vector{Vector{P}}, δ;
    abstol=1e-4, interp=LinearInterpolation) where {P}
    nic = data[end]
    result = P[]
    N = length(nic[1].u[1])
    T = typeof(nic[1].u[1][1])
    tmap = v.timetmap
    tend = v.timespan[end]
    hypers = v.f.hypers
    dtimes = δ / abstol
    @inbounds for j in eachindex(nic)
        curve = nic[j]
        n = length(curve.u)
        ic2_states = Vector{NSState{N,T}}(undef, n)
        for k in eachindex(ic2_states)
            ic2_states[k] = tmap(curve.u[k], para)
        end
        olds = copy(curve.t)
        newpara = ns_addpoints!(tmap, para, δ, curve, ic2_states, olds, dtimes, tend, hypers)
        _result = partition(ic2_states, newpara)
        # insert left zeros and right zeros
        for i in eachindex(_result)
            del_extra!(_result[i].t, _result[i].u, δ)
        end
        # delete extra points
        append!(result, _result)
    end
    result
    # append!(data, [result])
end


function generate_curves(f::NSSetUp, p, saddle, direction, δ, N; interp=LinearInterpolation, n=150, initial_d=0.01, abstol=1e-4)
    curves = ns_initialise_curve(f, p, saddle, direction, n, initial_d, δ; abstol=abstol, interp=interp)
    for i in 1:N
        grow_line!(f, p, curves, δ; interp=interp, abstol=abstol)
    end
    curves
end