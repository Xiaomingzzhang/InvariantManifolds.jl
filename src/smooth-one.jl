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

function OneDManifoldProblem(f; amax=0.5, d=1e-3, dsmin=1e-6)
    OneDManifoldProblem(f, Float64[], amax, d, dsmin)
end


function OneDManifoldProblem(f, para::AbstractVector{T};
    amax=T(0.5), d=T(1e-3), dsmin=T(1e-6)) where {T}
    OneDManifoldProblem(f, para, amax, d, dsmin)
end


function show(io::IO, m::MIME"text/plain", A::OneDManifoldProblem)
    printstyled(io, "OneDManifoldProblem:"; color=:cyan)
    println(io)
    print(io, "f:")
    show(io, A.f)
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
    print(io, "dsmin: ")
    show(io, m, A.dsmin)
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


"""
    OneDManifold{F,S,N,T}

`OneDManifold` is a struct contains all the information of the one-dimensional numerical manifold.
# Fields
- `prob` the problem [`OneDManifoldProblem`](@ref);
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
    arclength = 0.0
    for i in eachindex(A.data)
        m = m + length(A.data[i].t)
        arclength = arclength + A.data[i].t[end]
    end
    k = length(A.flawpoints)
    printstyled(io, "One-dimensional manifold"; bold=true, color=:cyan)
    println(io)
    printstyled(io, "Curves number: "; color=:cyan)
    println(io, "$n")
    printstyled(io, "Points number: "; color=:cyan)
    println(io, "$m")
    printstyled(io, "Total arc length: "; color=:cyan)
    println(io, "$arclength")
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

"""
    gen_segment(saddle, direction)

Generating `n` points at `saddle` in the `direction`, with length `d`, with default `n=50` and `d=0.01`.
Another Convenient consturctor is `gen_segment(p::Saddle)`.
"""
function gen_segment(saddle::SVector{N,T}, direction; n=50, d=0.01) where {N,T}
    tangent = normalize(direction)
    data = Vector{SVector{N,T}}(undef, n)
    for i in eachindex(data)
        data[i] = saddle + ((d * (i - 1)) / (n - 1)) * tangent
    end
    data
end

function gen_segment(p::Saddle{N,T,S}; n=50, d=0.01) where {N,T,S}
    saddle = p.saddle
    direction = p.unstable_directions[1]
    tangent = normalize(direction)
    data = Vector{SVector{N,T}}(undef, n)
    for i in eachindex(data)
        data[i] = saddle + ((d * (i - 1)) / (n - 1)) * tangent
    end
    data
end

function newasin(x::T) where {T}
    if  abs(x) <= 1
        return asin(x)
    elseif x > 1
        return T(pi)/2
    else
        return -T(pi)/2
    end
end

"""
    InvariantManifolds.addpoints!(f, p, d, oldcurve, newu, olds, dsmin, αmax, flawpoints)

Add points to ensure proper spacing and curvature along a one-dimensional manifold curve.

# Arguments
- `f`: The nonlinear map function
- `p`: Parameters for the nonlinear map
- `d`: Maximum allowed distance between consecutive points
- `oldcurve`: Previous curve interpolation
- `newu::Vector{SVector{N,T}}`: New points to be processed
- `olds::Vector{T}`: Parameter values of the points
- `dsmin`: Minimum allowed arc length between points
- `αmax`: Maximum allowed angle between consecutive segments
- `flawpoints`: Vector to store problematic points that don't meet criteria

# Details
The function iteratively processes points, adding new ones when:
- The distance between consecutive points exceeds `d`
- The angle between consecutive segments exceeds `αmax`
- The arc length is greater than `dsmin`

If a point fails these criteria and the minimum arc length is reached, it's recorded as a flaw point.

# Returns
Vector of parameter values for the new points
"""
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
            δ1 = norm(u0 - u1)
            δ2 = norm(u1 - u2)
            if δ1 == 0 || (δ1 != 0 && δ2 == 0)
                deleteat!(newu, i + 1)
                deleteat!(olds, i + 1)
                n = n - 1
            else
                baru0 = u1 + (u1 - u2) * δ1 / δ2
                α = 2*newasin(norm(baru0 - u0) / (2*δ1))
                if δ1 <= d && α <= αmax
                    i = i + 1
                    dd = newpara[end]
                    append!(newpara, [dd + δ1])
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
                        append!(flawpoints, [FlawPoint(u0, α, δ1)])
                        dd = newpara[end]
                        append!(newpara, [dd + δ1])
                    end
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
- `interp` the interpolation method used, default to be `QuadraticInterpolation`.
"""
function grow!(manifold::OneDManifold; interp=QuadraticInterpolation)
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
    manifold
end

function paramise(data::Vector{S}; interp=QuadraticInterpolation) where {S}
    m = length(data)
    T = typeof(data[1][1])
    if m == 1
        interp([data[1], data[1], data[1]], [T(0), T(0), T(0)])
    elseif m == 2
        s1 = norm(data[2] - data[1])
        interp([data[1], data[2], data[2]], [0, s1, s1])
    else
        s0 = Vector{T}(undef, m)
        s0[1] = T(0)
        for i in 2:m
            dd = norm(data[i] - data[i-1])
            s0[i] = s0[i-1] + dd
        end
        interp(data, s0)
    end
end

function paramise(data::Vector{S}, s0::Vector{T}; interp=QuadraticInterpolation) where {S,T}
    m = length(data)
    if m == 1
        interp([data[1], data[1], data[1]], [T(0), T(0), T(0)])
    elseif m == 2
        s1 = norm(data[2] - data[1])
        interp([data[1], data[2], data[2]], [0, s1, s1])
    else
        interp(data, s0)
    end
end


"""
    initialize(prob, points)

This is a function to initialize the continuation process. Its output is a manifold struct.
# Parameters
- `prob` the problem such as [`OneDManifoldProblem`](@ref).
- `points` the points in the local manifold. For one dimensional manifolds, these points should be a `Vector{SVector}` and the start point should be the saddle. For two dimensional manifolds, these points should be a `Vector{Vector{S}}` and its first element should like `[saddle, saddle, saddle]`. Note that in the both cases, the functions [`gen_segment`](@ref) and [`gen_disk`](@ref) can generate these points easily.

# Keyword argument
- `interp` the interpolation method used, default to be `QuadraticInterpolation`.
"""
function initialize(prob::OneDManifoldProblem, points::Vector{SVector{N,T}}; interp=QuadraticInterpolation) where {N,T}
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
- `prob` the problem such as [`OneDManifoldProblem`](@ref).
- `points` the points in the local manifold. For one dimensional manifolds, these points should be a `Vector{SVector}` and the start point should be the saddle. For two dimensional manifolds, these points should be a `Vector{Vector{S}}` and its first element should like `[saddle, saddle, saddle]`. Note that in the both cases, the functions [`gen_segment`](@ref) and [`gen_disk`](@ref) can generate these points easily.
- `N` the number of iterations.

# Keyword argument
- `interp` the interpolation method used, default to be `QuadraticInterpolation`.
"""
function growmanifold(prob::OneDManifoldProblem, points, N; interp=QuadraticInterpolation)
    manifold = initialize(prob, points, interp=interp)
    for i in 1:N
        grow!(manifold, interp=interp)
    end
    manifold
end