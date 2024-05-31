function show(io::IO, v::Vector{IterationCurve{N,T}}) where {N,T}
    n = 0
    for i in eachindex(v)
        n = n + length(v[i].states)
    end
    print(io, "Vector of IterationCurve{$N,$T} with total $n points")
end

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

@inline function nsnorm(a::NSState{N,T}, b::NSState{N,T}, rules, p, time_end, ntimes) where {N,T}
    if a.event_at == b.event_at
        norm(a - b)
    elseif abs(length(a.event_at) - length(b.event_at)) == 1
        if a.event_at == b.event_at[1:end-1]
            norm(rules[b.event_at[end]](a.state, p, time_end) - b.state) * ntimes
        elseif a.event_at[1:end-1] == b.event_at
            norm(rules[a.event_at[end]](b.state, p, time_end) - a.state) * ntimes
        else
            T(2)
        end
    else
        T(2)
    end
end

@inline function addpoints(f, p, nic::Vector{IterationCurve{N,T}}, rules, time_end, min, ntimes) where {N,T}
    result = IterationCurve{N,T}[]
    @inbounds for j in eachindex(nic)
        n = length(nic[j].states)
        ic2_states = Vector{NSState{N,T}}(undef, n)
        for k in eachindex(ic2_states)
            ic2_states[k] = f(nic[j].states[k])
        end
        curve = nic[j].pcurve
        i = 1
        newpara = T[0]
        while i + 1 <= n
            dist = nsnorm(ic2_states[i], ic2_states[i+1], rules, time_end, p, ntimes)
            if dist > min
                m = ceil(Int, dist / min)
                if m == 1
                    m = 2
                end
                s0 = ic2_states[i].s
                s1 = ic2_states[i+1].s
                plengh = (s1 - s0) / m
                paras = [s0 + plengh * i for i in 1:m-1]
                addps = Vector{NSState{N,T}}(undef, m - 1)
                for kk in 1:m-1
                    addps[kk] = f(State(curve(paras[kk]), paras[kk]))
                end
                insert!(ic2_states, i + 1, addps)
                n = n + m - 1
            else
                i = i + 1
                dd = newpara[end]
                append!(newpara, [dd + dist])
            end
        end
        __result = partition(ic2_states)
        mm = length(__result)
        _result = Vector{IterationCurve{N,T}}(undef, mm)
        l0 = 1
        lend = length(__result[1])
        for k in 1:mm
            newstates = Vector{State{N,T}}(undef, length(__result[k]))
            for i in eachindex(newstates)
                newstates[i] = State(__result[k][i].state, newpara[l0+i-1])
            end
            _result[k] = IterationCurve(newstates, LinearInterpolation(take_state.(__result[k]), newpara[l0:lend]))
            if k <= mm - 1
                l0 = lend + 1
                lend = lend + length(__result[k+1])
            end
        end
        append!(result, _result)
    end
    result
end

"""
    InvariantManifolds.initialise_curve(points, tmap)

Initialise the curve of ODE's time-T-map.
"""
function initialise_curve(points, tmap)
    map = x -> tmap(State(x, 0))
    n = length(points)
    m = length(points[1])
    type = typeof(points[1][1])
    _image = Vector{NSState{m,type}}(undef, n)
    for i in eachindex(_image)
        _image[i] = map(points[i])
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
    result = image[j:end]
    prepend!(result, [pn])
    paramise(result)
end


ID(x, p, t) = x

function generate_curves(v::NSSetUp{S}, seg, d, n; ntimes=100) where {S<:ContinuousVectorField}
    N = length(seg[1])
    T = typeof(seg[1][1])
    m = length(v.f.hypers)
    result = Vector{Vector{IterationCurve{N,T}}}(undef, 2)
    _seg = initialise_curve(seg, v.timetmap)
    result[1] = [paramise(seg)]
    result[2] = [_seg]
    for i in 1:n
        seg2 = addpoints(v.timetmap, v.p, result[end], fill(ID, m), v.timespan[end], d, ntimes)
        append!(result, [seg2])
    end
    result
end

function generate_curves(v::NSSetUp{S}, seg, d, n; ntimes=100) where {S<:JumpVectorField}
    N = length(seg[1])
    T = typeof(seg[1][1])
    result = Vector{Vector{IterationCurve{N,T}}}(undef, 2)
    _seg = initialise_curve(seg, v.timetmap)
    result[1] = [paramise(seg)]
    result[2] = [_seg]
    for i in 1:n
        seg2 = addpoints(v.timetmap, v.p, result[end], v.f.rules, v.timespan[end], d, ntimes)
        append!(result, [seg2])
    end
    result
end