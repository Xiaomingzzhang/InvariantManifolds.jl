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
    [interp(result[i], result_s[i]) for i in eachindex(result)]
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

@inline function nsnorm(a::NSState{N,T}, b::NSState{N,T}, mirrors, p, time_end, dtimes) where {N,T}
    if a.event_at == b.event_at
        norm(a - b)
    elseif abs(length(a.event_at) - length(b.event_at)) == 1
        if a.event_at == b.event_at[1:end-1]
            norm(mirrors[b.event_at[end]](a.state, p) - b.state) * dtimes
        elseif a.event_at[1:end-1] == b.event_at
            norm(mirrors[a.event_at[end]](b.state, p) - a.state) * dtimes
        else
            T(2)
        end
    else
        T(2)
    end
end


@inline function ns_addpoints!(tmap, p, δ, oldcurve, newu::Vector{NSState{N,T}}, olds::Vector{T}, dtimes, tend, mirrors) where {N,T}
    n = length(newu)
    i = 1
    newpara = T[0]
    event = oldcurve.u[1].event_at
    @inbounds while i + 1 <= n
        dist = nsnorm(newu[i], newu[i+1], mirrors, p, tend, dtimes)
        if dist > δ
            @show dist
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


take_state(x) = x.state


"""
    InvariantManifolds.ns_initialise_curve(points, tmap, para)

Initialise the curve of non-smooth ODE's time-T-map.
"""
function ns_initialise_curve(f::NSSetUp, para, saddle, direction, nn, d, δ; dtimes=100, interp=LinearInterpolation)
    points = NSState.(segment(saddle, direction, nn, d))
    result = [[paramise(points, interp=interp)]]
    tmap = f.timetmap
    tend = f.timespan[end]
    mirrors = f.f.mirrors
    n = length(points)
    m = length(saddle)
    type = typeof(saddle[1][1])
    image = Vector{NSState{m,type}}(undef, n)
    for i in eachindex(image)
        image[i] = tmap(points[i], para)
    end
    oldcurve = result[1][1]
    olds = copy(result[1][1].t)
    ns_addpoints!(tmap, para, δ, oldcurve, image, olds, dtimes, tend, mirrors)
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
    dtimes=100, interp=LinearInterpolation) where {P}
    nic = data[end]
    result = P[]
    N = length(nic[1].u[1])
    T = typeof(nic[1].u[1][1])
    tmap = v.timetmap
    tend = v.timespan[end]
    mirrors = v.f.mirrors
    @inbounds for j in eachindex(nic)
        curve = nic[j]
        n = length(curve.u)
        ic2_states = Vector{NSState{N,T}}(undef, n)
        for k in eachindex(ic2_states)
            ic2_states[k] = tmap(curve.u[k], para)
        end
        olds = copy(curve.t)
        newpara = ns_addpoints!(tmap, para, δ, curve, ic2_states, olds, dtimes, tend, mirrors)
        _result = partition(ic2_states, newpara)
        # insert left zeros and right zeros

        # delete extra points
        append!(result, _result)
    end
    append!(data, [result])
end


function generate_curves(f::NSSetUp, p, saddle, direction, δ, N; interp=LinearInterpolation, n=150, initial_d=0.01, dtimes=100)
    curves = ns_initialise_curve(f, p, saddle, direction, n, initial_d, δ; dtimes=dtimes, interp=interp)
    for i in 1:N
        grow_line!(f, p, curves, δ; interp=interp, dtimes=dtimes)
    end
    curves
end