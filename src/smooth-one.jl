
"""
    OneDManifoldProblem{F,T}

`OneDManifoldProblem` is a struct to contain the main information for continuing the one-dimensional manifold of a nonlinear map.
# Fields
- `f` the nonlinear map, which should has the form `f(x,p)` and return a `SVector`;
- `para` the parameters of the nonlinear map;
- `amax` the maximum angle between points when continuing the manifold;
- `d` the maximum distance between points when continuing the manifold;
- `dsmin` the minimum arc length allowing; note that if in a continuation point, this value is achieved and the angle as well as the distance values are not achieved, then we will record this point as a [`FlawPoint`](@ref).

Convenient consturctors are `OneDManifoldProblem(f)` and `OneDManifoldProblem(f,para)`
"""
struct OneDManifoldProblem{F,T}
    f::F
    para::Vector{T}
    amax::T
    d::T
    dsmin::T
end


"""
    FlawPoint{N,T}

`FlawPoint` is a struct to record the points that don't satisfy the angle and distance request while the 
program has reached the minimum arc length.
# Fields
- `point` flaw point in the process of continuation;
- `α` the angle recorded;
- `d` the distance recorded.
"""
struct FlawPoint{N,T}
    point::SVector{N,T}
    α::T
    d::T
end


function OneDManifoldProblem(f; amax=0.5, d=0.001, dsmin=1e-5)
    OneDManifoldProblem(f, Float64[], amax, d, dsmin)
end


function OneDManifoldProblem(f, para::Vector{T};
    amax=T(0.5), d=T(0.001), dsmin=T(1e-5)) where {T}
    OneDManifoldProblem(f, para, amax, d, dsmin)
end


"""
    OneDManifold{F,S,N,T}

`OneDManifold` is a struct contains all the information of the one-dimensional numerical manifold.
# Fields
- `prob` the problem `OneDManifoldProblem`;
- `data` the numerical data that should be `Vector{Vector{S}}`, where `S` is the interpolation curve (we use `DataInterpolation` in this package);
- `flawpoints` the flaw points generated during continuation.
"""
mutable struct OneDManifold{F,S,N,T}
    prob::OneDManifoldProblem{F,T}
    data::Vector{S}
    flawpoints::Vector{FlawPoint{N,T}}
end

function Base.show(io::IO, m::MIME"text/plain", A::OneDManifold)
    m = 0
    n = length(A.data)
    for i in eachindex(A.data)
        m = m + length(A.data[i].t)
    end
    k = length(A.flawpoints)
    println(io, "Two-dimensional manifold")
    println(io, "Curves number: $n")
    println(io, "Points number: $m")
    amax = A.prob.amax
    d = A.prob.d
    prend = findall(x -> x.d > d, A.flawpoints)
    nd = length(prend)
    prenc = findall(x -> x.α > amax, A.flawpoints)
    nc = length(prenc)
    println(io, "Flaw points number: $k")
    println(io, "Distance failed  points number: $nd")
    println(io, "Curvature failed points number: $nc")
end

"""
    segment(saddle, direction)

Generating `n` points at `saddle` in the `direction`, with length `d`, with default `n=150` and `d=0.01`.
Another Convenient consturctor is `segment(p::Saddle)`.
"""
function gen_segment(saddle::SVector{N,T}, direction; n=150, d=0.01) where {N,T}
    tangent = normalize(direction)
    data = Vector{SVector{N,T}}(undef, n)
    for i in eachindex(data)
        data[i] = saddle + ((d * (i - 1)) / (n - 1)) * tangent
    end
    data
end

function gen_segment(p::Saddle{N,T,S}; n=150, d=0.01) where {N,T,S}
    saddle = p.saddle
    direction = p.unstable_directions[1]
    tangent = normalize(direction)
    data = Vector{SVector{N,T,S}}(undef, n)
    for i in eachindex(data)
        data[i] = saddle + ((d * (i - 1)) / (n - 1)) * tangent
    end
    data
end

@inline function addpoints!(f, p, d, oldcurve, newu::Vector{SVector{N,T}}, olds::Vector{T},
    dsmin, αmax, flawpoints) where {N,T}
    n = length(newu)
    i = 1
    newpara = T[0]
    # first add enough points
    @inbounds while i + 1 <= n
        if i + 2 <= n
            u0 = newu[i]
            u1 = newu[i+1]
            u2 = newu[i+2]
            δ = norm(u0 - u1)
            baru0 = u1 + (u1 - u2) * norm(u1 - u0) / norm(u1 - u2)
            α = norm(baru0 - u0) / (norm(u1 - u0))
            if δ <= d && α <= αmax
                i = i + 1
                dd = newpara[end]
                append!(newpara, [dd + δ])
            else
                if olds[i+1] - olds[i] > dsmin
                    s0 = olds[i]
                    s1 = olds[i+1]
                    paras = (s0 + s1) / 2
                    addps = f(oldcurve(paras), p)
                    insert!(newu, i + 1, addps)
                    insert!(olds, i + 1, paras)
                    n = n + 1
                else
                    i = i + 1
                    append!(flawpoints, [FlawPoint(u0, α, δ)])
                    dd = newpara[end]
                    append!(newpara, [dd + δ])
                end
            end
        else
            u0 = newu[i]
            u1 = newu[i+1]
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
                    addps = f(oldcurve(paras), p)
                    insert!(newu, i + 1, addps)
                    insert!(olds, i + 1, paras)
                    n = n + 1
                else
                    i = i + 1
                    append!(flawpoints, [FlawPoint(u0, T(0), δ)])
                    dd = newpara[end]
                    append!(newpara, [dd + δ])
                end
            end
        end
    end
    newpara
end

"""
    grow!(manifold)

One time iteration to grow the manifold.
# Parameters
- `manifold` the manifold struct.
# Keyword argument
- `interp` the interpolation method used, default to be `LinearInterpolation`.
"""
function grow!(manifold::OneDManifold; interp=LinearInterpolation)
    αmax = manifold.prob.amax
    d = manifold.prob.d
    f = manifold.prob.f
    para = manifold.prob.para
    Δsmin = manifold.prob.dsmin
    flawpoints = manifold.flawpoints
    data = manifold.data
    curve = data[end]
    ic1 = copy(curve.u)
    olds = copy(curve.t)
    n = length(ic1)
    datatype = typeof(ic1[1])
    ic2 = Vector{datatype}(undef, n)
    Threads.@threads for i in eachindex(ic1)
        ic2[i] = f(ic1[i], para)
    end
    newpara = addpoints!(f, para, d, curve, ic2, olds, Δsmin, αmax, flawpoints)
    newintep = interp(ic2, newpara)
    append!(data, [newintep])
end

function paramise(data::Vector{S}; interp=LinearInterpolation) where {S}
    m = length(data)
    if m == 1
        T = typeof(data[1][1])
        interp(data, [T(0)])
    else
        T = typeof(data[1][1])
        s0 = Vector{T}(undef, m)
        s0[1] = 0
        for i in 2:m
            dd = norm(data[i] - data[i-1])
            s0[i] = s0[i-1] + dd
        end
        interp(data, s0)
    end
end


"""
    initialize(prob, points)

This is a function to initialize the continuation process. Its output is a manifold struct.
# Parameters
- `prob` the problem such as `OneDManifoldProblem`.
- `points` the points in the local manifold. For one dimensional manifolds, these points should be a `Vector{SVector}` and the start point should be the saddle. For two dimensional manifolds, these points should be a `Vector{Vector{S}}` and its first element should like `[saddle, saddle, saddle]`. Note that in the both cases, the functions [`gen_segment`](@ref), [`gen_disk`](@ref), and [`gen_circle`](@ref) can generate these points easily.

# Keyword argument
- `interp` the interpolation method used, default to be `LinearInterpolation`.
"""
function initialize(prob::OneDManifoldProblem, points::Vector{SVector{N,T}}; interp=LinearInterpolation) where {N,T}
    parameters = prob.para
    map = prob.f
    αmax = prob.amax
    d = prob.d
    Δsmin = prob.dsmin
    flawpoints = FlawPoint{N,T}[]
    Δsmin = prob.dsmin
    result = [paramise(points, interp=interp)]
    n = length(points)
    data = Vector{SVector{N,T}}(undef, n)
    for i in eachindex(data)
        data[i] = map(points[i], parameters)
    end
    curve = paramise(points)
    olds = copy(curve.t)
    addpoints!(map, parameters, d, curve, data, olds, Δsmin, αmax, flawpoints)
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
    OneDManifold(prob, result, flawpoints)
end


"""
    growmanifold(prob, points, N)

This is the mani function to continuate the numerical manifolds. Its output is a manifold struct.
# Parameters
- `prob` the problem such as `OneDManifoldProblem`.
- `points` the points in the local manifold. For one dimensional manifolds, these points should be a `Vector{SVector}` and the start point should be the saddle. For two dimensional manifolds, these points should be a `Vector{Vector{S}}` and its first element should like `[saddle, saddle, saddle]`. Note that in the both cases, the functions [`segment`](@ref) and [`disk`](@ref) can generate these points easily.
- `N` the number of iterations.

# Keyword argument
- `interp` the interpolation method used, default to be `LinearInterpolation`.
"""
function growmanifold(prob::OneDManifoldProblem, points, N; interp=LinearInterpolation)
    manifold = initialize(prob, points, interp=interp)
    for i in 1:N
        grow!(manifold, interp=interp)
    end
    manifold
end