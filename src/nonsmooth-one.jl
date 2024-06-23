function partition(v::Vector{NSState{N,T}}) where {N,T}
    ctime1 = v[1].event_at
    n = length(v)
    a0 = NSState{N,T}[]
    result = Vector{NSState{N,T}}[]
    for j in 1:n
        if v[j].event_at == ctime1
            append!(a0, [v[j]])
        else
            append!(result, [a0])
            a0 = [v[j]]
            ctime1 = v[j].event_at
        end
    end
    append!(result, [a0])
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

@inline function nsnorm(a::NSState{N,T}, b::NSState{N,T}, rules, p, time_end, dtimes) where {N,T}
    if a.event_at == b.event_at
        norm(a - b)
    elseif abs(length(a.event_at) - length(b.event_at)) == 1
        if a.event_at == b.event_at[1:end-1]
            norm(rules[b.event_at[end]](a.state, p, time_end) - b.state) * dtimes
        elseif a.event_at[1:end-1] == b.event_at
            norm(rules[a.event_at[end]](b.state, p, time_end) - a.state) * dtimes
        else
            T(2)
        end
    else
        T(2)
    end
end

take_state(x) = x.state


"""
    InvariantManifolds.ns_initialise_curve(points, tmap, para)

Initialise the curve of non-smooth ODE's time-T-map.
"""
function ns_initialise_curve(f::NSSetUp, para, saddle, direction, nn, d; interp=LinearInterpolation)
    points = segment(saddle, direction, nn, d)
    result = [[paramise(points, interp=interp)]]
    tmap = f.timetmap
    n = length(points)
    m = length(saddle)
    type = typeof(saddle[1][1])
    _image = Vector{NSState{m,type}}(undef, n)
    for i in eachindex(_image)
        _image[i] = tmap(points[i], para)
    end
    if ispartitioned(_image) == false
        error("The initial curve has to be chosen more small")
    end
    image = take_state.(_image)
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

@inline function grow_line!(v::NSSetUp, para, data::Vector{Vector{P}}, δ;
    dtimes=100, interp=LinearInterpolation) where {P}
    nic = data[end]
    result = P[]
    N = length(nic[1].u[1])
    T = typeof(nic[1].u[1][1])
    rules = v.f.rules
    t_end = v.timespan[end]
    f = v.timetmap
    @inbounds for j in eachindex(nic)
        curve = nic[j]
        n = length(curve.u)
        ic2_states = Vector{NSState{N,T}}(undef, n)
        for k in eachindex(ic2_states)
            ic2_states[k] = f(curve.u[k], para)
        end
        i = 1
        news = T[0]
        olds = copy(curve.t)
        while i + 1 <= n
            dist = nsnorm(ic2_states[i], ic2_states[i+1], rules, para, t_end, dtimes)
            if dist > δ
                m = ceil(Int, dist / δ)
                if m == 1
                    m = 2
                end
                s0 = olds[i]
                s1 = olds[i+1]
                plengh = (s1 - s0) / m
                paras = [s0 + plengh * i for i in 1:m-1]
                addps = Vector{NSState{N,T}}(undef, m - 1)
                for kk in 1:m-1
                    addps[kk] = f(curve(paras[kk]), para)
                end
                insert!(ic2_states, i + 1, addps)
                insert!(olds, i + 1, paras)
                n = n + m - 1
            else
                if ic2_states[i].event_at == ic2_states[i+1].event_at
                    newdist = dist
                else
                    newdist = dist / dtimes
                end
                dd = news[end]
                append!(news, [dd + newdist])
                i = i + 1
            end
        end
        __result = partition(ic2_states)
        mm = length(__result)
        _result = Vector{P}(undef, mm)
        l0 = 1
        lend = length(__result[1])
        for k in 1:mm
            _result[k] = interp(take_state.(__result[k]), news[l0:lend])
            if k <= mm - 1
                l0 = lend + 1
                lend = lend + length(__result[k+1])
            end
        end
        append!(result, _result)
    end
    append!(data, [result])
end


function generate_curves(f::NSSetUp, p, saddle, direction, δ, N; interp=LinearInterpolation, n=150, initial_d=0.01, dtimes=100)
    curves = ns_initialise_curve(f, p, saddle, direction, n, initial_d; interp=interp)
    for i in 1:N
        grow_line!(f, p, curves, δ; interp=interp, dtimes=dtimes)
    end
    curves
end