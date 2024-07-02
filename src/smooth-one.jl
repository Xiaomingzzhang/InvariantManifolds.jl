"""
    segment(saddle, direction, n, d)

Generating `n` points at `saddle` in the `direction`, with length `d`
"""
function segment(saddle, direction, n, d)
    tangent = normalize(direction)
    datatype = typeof(saddle)
    data = Vector{datatype}(undef, n)
    for i in eachindex(data)
        data[i] = saddle + ((d * (i - 1)) / (n - 1)) * tangent
    end
    data
end

function paramise(data::Vector{S}; interp=LinearInterpolation) where {S}
    m = length(data)
    if m == 1
        error("there is just one point!")
    end
    T = typeof(data[1][1])
    s0 = Vector{T}(undef, m)
    s0[1] = 0
    for i in 2:m
        dd = norm(data[i] - data[i-1])
        s0[i] = s0[i-1] + dd
    end
    interp(data, s0)
end

# extend the function Base.insert! so that we can insert many elements at a certain position.
function Base.insert!(a::Array{T,1}, i::Integer, b::Array{T,1}) where {T}
    # Throw convert error before changing the shape of the array
    n = length(b)
    Base._growat!(a, i, n)
    # _growat! already did bound check
    @inbounds for j in 1:n
        a[i+j-1] = b[j]
    end
    return a
end

function Base.insert!(a::Vector{Vector{T}}, i::Integer, b::Vector{T}) where {T}
    # Throw convert error before changing the shape of the array
    Base._growat!(a, i, 1)
    # _growat! already did bound check
    a[i] = b
    return a
end

@inline function addpoints!(f, p, min, oldcurve, newu::Vector{SVector{N,T}}, olds::Vector{T}; del_extra=false) where {N,T}
    n = length(newu)
    i = 1
    newpara = T[0]
    # first add enough points
    @inbounds while i + 1 <= n
        dist = norm(newu[i+1] - newu[i])
        if dist > min
            m = ceil(Int, dist / min)
            if m == 1
                m = 2
            end
            s0 = olds[i]
            s1 = olds[i+1]
            plengh = (s1 - s0) / m
            paras = [s0 + plengh * i for i in 1:m-1]
            addps = Vector{SVector{N,T}}(undef, m - 1)
            for j in 1:m-1
                addps[j] = f(oldcurve(paras[j]), p)
            end
            insert!(newu, i + 1, addps)
            insert!(olds, i + 1, paras)
            n = n + m - 1
        else
            i = i + 1
            dd = newpara[end]
            append!(newpara, [dd + dist])
        end
    end
    # delete unnesscery points
    if del_extra == true
        k = 1
        M = length(newpara)
        @inbounds while k + 2 <= M
            s1 = newpara[k+1] - newpara[k]
            s2 = newpara[k+2] - newpara[k+1]
            s3 = norm(newu[k+2] - newu[k])
            if s1 < min && s2 < min && s3 < min
                deleteat!(newpara, k + 1)
                deleteat!(olds, k + 1)
                deleteat!(newu, k + 1)
                M = M - 1
            else
                k = k + 1
            end
        end
    end
    newpara
end


"""
    InvariantManifolds.initialise_curve(points, map, parameters)

This is a function to get the iteration curve. From `[p0,...,pn]`, where `p0` is the
saddle point, to get `[pn,...,f(pn)]`.
# Parameters
- `points` stands for distributed points in a local manifolds, `points[1]` is the saddle;
- `map` is the map function;
- `parameters` is the parameters of the `map`.
"""
function initialise_curve(map, parameters, saddle, direction, nn, d, min; interp=LinearInterpolation, del_extra=false)
    points = segment(saddle, direction, nn, d)
    result = [paramise(points, interp=interp)]
    n = length(points)
    data = Vector{eltype(points)}(undef, n)
    for i in eachindex(data)
        data[i] = map(points[i], parameters)
    end
    curve = paramise(points)
    olds = copy(curve.t)
    addpoints!(map, parameters, min, curve, data, olds, del_extra = del_extra)
    p0 = first(points)
    pn = last(points)
    j = 1
    dist = norm(pn - p0)
    while norm(data[j] - p0) < dist
        j = j + 1
    end
    j = j + 1
    data = data[j:end]
    prepend!(data, [pn])
    append!(result, [paramise(data, interp=interp)])
end



"""
    InvariantManifolds.addpoints(f, p, ic1::IterationCurve{N,T}, min)
    
`addpoints` will add enough points from `ic1` so that its data of points are
dense, i.e., the distance of nearby points will less than `min`. 
This function will also delete extra points so that the distance of nearby points aren't two small.
"""
function grow_line!(f, p, data, min; interp=LinearInterpolation, del_extra = false)
    ic1 = copy(data[end].u)
    olds = copy(data[end].t)
    n = length(ic1)
    if n < 3
        error("The length of states must be greater than 2")
    end
    datatype = typeof(ic1[1])
    T = typeof(ic1[1][1])
    ic2 = Vector{datatype}(undef, n)
    @inbounds @simd for i in eachindex(ic1)
        ic2[i] = f(ic1[i], p)
    end
    curve = data[end]
    newpara = addpoints!(f, p, min, curve, ic2, olds, del_extra = del_extra)
    newintep = interp(ic2, newpara)
    append!(data, [newintep])
end

"""
    generate_curves(f, p, seg, d, n)

This function is to generate one-dimensional manifold of smooth mapping or a time-T-map of a non-smooth ODE.
# Parameters
- `f` a discrete map of type `f(x,p)=SA[...]` or a `NSSetUp` of non-smooth ODE;
- `para` the parameters of systems;
- `seg` initial segment of manifolds, which can be generated by function `segment` if the unstable direction is known;
- `d` max distance between points;
- `n` iteration times.
"""
function generate_curves(f, p, saddle, direction, d, N; interp=LinearInterpolation, n=150, initial_d=0.01, del_extra=false)
    curves = initialise_curve(f, p, saddle, direction, n, initial_d, d; interp=interp, del_extra = del_extra)
    for i in 1:N
        grow_line!(f, p, curves, d; interp=interp, del_extra = del_extra)
    end
    curves
end