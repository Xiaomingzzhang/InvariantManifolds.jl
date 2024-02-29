struct State{N,T<:Number}
    state::SVector{N,T}
    s::T
end

function State(v::SVector{N,T}, t::M) where {M<:Number,N,T<:Number}
    datetype = promote_type(M, T)
    newv = convert(SVector{N,datetype}, v)
    newt = convert(datetype, t)
    State(newv, newt)
end

struct IterationCurve{N,T<:Number}
    states::Vector{State{N,T}}
    pcurve
end

function show(io::IO, v::IterationCurve{N,T}) where {N,T<:Number}
    n = length(v.states)
    print(io, "IterationCurve{$N,$T} with $n points")
end

function take_s(x)
    x.s
end

function take_state(x)
    x.state
end

function linear_interpolation(v::Vector{State{N,T}}) where {N,T<:Number} # Vector{T} 中的T必须是具体的
    linear_interpolation(take_s.(v), take_state.(v))
end

function -(a, b)
    a.state - b.state
end

"""
    segment(point, direction, n, d)

Generating `n` points at `point` in the `direction`, with length `d`
"""
function segment(point, direction, n, d)
    tangent = normalize(direction)
    datatype = typeof(point)
    data = Vector{datatype}(undef, n)
    for i in eachindex(data)
        data[i] = point + ((d * (i - 1)) / (n - 1)) * tangent
    end
    data
end

function paramise(image)
    datatype = typeof(image[1][1])
    l = length(image[1])
    m = length(image)
    if m == 1
        error("inintial points are too small")
    end
    final = Vector{State{l,datatype}}(undef, m)
    @inbounds for q in 1:m
        final[q] = State(image[q], (q - 1) // (m - 1))
    end
    IterationCurve(final, linear_interpolation(final))
end

"""
    InvariantManifolds.initialise_curve(points, map, parameters)

This is a function to get the iteration curve. From `[p0,...,pn]`, where `p0` is the
saddle point, to get `[pn,...,f(pn)]`.
- `points` stands for distributed points in a local manifolds, `points[1]` is the saddle;
- `map` is the map function;
- `parameters` is the parameters of the `map`.
"""
function initialise_curve(points, map, parameters)
    n = length(points)
    image = Vector{eltype(points)}(undef, n)
    for i in eachindex(image)
        image[i] = map(points[i], parameters)
    end
    p0 = first(points)
    pn = last(points)
    j = 1
    d = norm(pn - p0)
    while norm(image[j] - p0) < d
        j = j + 1
    end
    j = j + 1
    result = image[j:end]
    result=prepend!(result, [pn])
    paramise(result)
end


function myinsert!(a::Array{T,1}, i::Integer, b::Array{S,1}) where {T,S}
    # Throw convert error before changing the shape of the array
    _b = S == T ? b : convert(Array{T,1}, b)::Array{T,1}
    n = length(b)
    Base._growat!(a, i, n)
    # _growat! already did bound check
    @inbounds for j in 1:n
        a[i+j-1] = _b[j]
    end
    return a
end

"""
    InvariantManifolds.addpoints(f, p, ic1::IterationCurve{N,T}, min)
    
`addpoints` will add enough points from `ic1` so that its image of points are
dense, i.e., the distance of nearby points will less than `min`.
"""
@inline function addpoints(f, p, ic1::IterationCurve{N,T}, min) where{N, T}
    n = length(ic1.states)
    if n < 3
        error("The length of states must be greater than 2")
    end
    datatype = typeof(ic1.states[1])
    ic2_states = Vector{datatype}(undef, n)
    @inbounds @simd for i in eachindex(ic2_states)
        ic2_states[i] = State(f(ic1.states[i].state, p), ic1.states[i].s)
    end
    curve = ic1.pcurve
    i = 1
    @inbounds while i + 1 <= n
        dist = norm(ic2_states[i+1] - ic2_states[i])
        if dist > min
            m = ceil(Int, dist / min)
            if m == 1
                m = 2
            end
            s0 = ic2_states[i].s
            s1 = ic2_states[i+1].s
            plengh = (s1 - s0) / m
            paras = [s0 + plengh * i for i in 1:m-1]
            addps = Vector{datatype}(undef, m - 1)
            for j in 1:m-1
                addps[j] = State(f(curve(paras[j]), p), paras[j])
            end
            myinsert!(ic2_states, i + 1, addps)
            n = n + m - 1
        else
            i = i + 1
        end
    end
    IterationCurve(ic2_states, linear_interpolation(ic2_states))
end

function iterate!(f, p, result, min, n)
    for i in 1:n
        append!(result, [addpoints(f, p, result[end], min)])
    end
    result
end


function generate_curves(f, p, seg, d, n)
    N = length(seg[1])
    T = typeof(seg[1][1])
    result = Vector{IterationCurve{N,T}}(undef, 2)
    _seg = initialise_curve(seg, f, p)
    result[1] = paramise(seg)
    result[2] = _seg
    for i in 1:n
        seg2 = addpoints(f, p, result[end], d)
        append!(result, [seg2])
    end
    result
end

function plot(v::Vector{IterationCurve{N,T}}) where {N,T}
    figure = plot(; legend=false)
    for i in eachindex(v)
        nn = length(v[i].states)
        id = range(0, 1, length=2 * nn)
        data = v[i].pcurve.(id)
        plot!(first.(data), last.(data))
    end
    figure
end