"""
    VTwoDManifoldProblem{F,T}

`VTwoDManifoldProblem` is a struct to contain the main information for continuing the two-dimensional manifold of an autonomous vector field.
# Fields
- `f` the time flow map of the vector field, which should has the form `f(x,p)` and return a `SVector`; note that the vector field should be regularized, e.g., ``\\dot{x}=v(x)`` should be rewritten as ``\\dot{x}=v(x)/\\sqrt{1+|v(x)|}``;
- `para` the parameters of the time flow map;
- `amax` the maximum angle between points when continuing the manifold;
- `d` the maximum distance between points when continuing the manifold;
- `dsmin` the minimum arc length allowing; note that if in a continuation point, this value is achieved and the angle as well as the distance values are not achieved, then we will record this point as a [`FlawPoint`](@ref).

Convenient consturctors are `VTwoDManifoldProblem(f)` and `VTwoDManifoldProblem(f,para)`
"""
struct VTwoDManifoldProblem{F,T}
    f::F
    para::Vector{T}
    amax::T
    d::T
    dsmin::T
end

function VTwoDManifoldProblem(f; amax=0.5, d=0.001, dsmin=1e-5)
    VTwoDManifoldProblem(f, Float64[], amax, d, dsmin)
end


function VTwoDManifoldProblem(f, para::AbstractVector{T};
    amax=T(0.5), d=T(0.001), dsmin=T(1e-5)) where {T}
    VTwoDManifoldProblem(f, para, amax, d, dsmin)
end

function show(io::IO, m::MIME"text/plain", A::VTwoDManifoldProblem)
    printstyled(io, "VTwoDManifoldProblem:"; color=:cyan)
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
    VTwoDManifold{F,S,N,T}

`VTwoDManifold` is a struct contains all the information of the two-dimensional numerical manifold of an autonomous vector field.
# Fields
- `prob` the problem [`VTwoDManifoldProblem`](@ref);
- `data` the numerical data that should be `Vector{S}`, where `S` is the interpolation curve (we use `DataInterpolation` in this package);
- `flawpoints` the flaw points generated during continuation.
"""
mutable struct VTwoDManifold{F,S,N,T}
    prob::VTwoDManifoldProblem{F,T}
    data::Vector{S}
    flawpoints::Vector{FlawPoint{N,T}}
end


function Base.show(io::IO, A::VTwoDManifold)
    m = 0
    n = length(A.data)
    for i in eachindex(A.data)
        m = m + length(A.data[i].t)
    end
    k = length(A.flawpoints)
    printstyled(io, "Two-dimensional manifold"; bold=true, color=:cyan)
    println(io)
    printstyled(io, "Circles number: "; color=:cyan)
    println(io, "$n")
    printstyled(io, "Points number: "; color=:cyan)
    println(io, "$m")
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
    print(io, "$nc")
end

function initialize(prob::VTwoDManifoldProblem, disk::Vector{Vector{SVector{N,T}}}; interp=QuadraticInterpolation) where {N,T}
    para = prob.para
    f = prob.f
    αmax = prob.amax
    d = prob.d
    dsmin = prob.dsmin
    flawpoints = FlawPoint{N,T}[]
    circles = [paramise(disk[i], interp=interp) for i in eachindex(disk)]
    VTwoDManifold(prob, circles, flawpoints)
end

function grow!(manifold::VTwoDManifold{F,S,N,T}; interp=QuadraticInterpolation) where {F,S,N,T}
    prob = manifold.prob
    para = prob.para
    f = prob.f
    αmax = prob.amax
    d = prob.d
    dsmin = prob.dsmin
    flawpoints = manifold.flawpoints
    dsmin = prob.dsmin
    data = manifold.data
    circle = data[end]
    olds = copy(circle.t)
    newcircle = Vector{SVector{N,T}}(undef, length(olds))
    points = circle.u
    Threads.@threads for i in eachindex(newcircle)
        newcircle[i] = f(points[i], para)
    end
    newpara = addpoints!(f, para, d, circle, newcircle, olds, dsmin, αmax, flawpoints)
    finalnewcircle = interp(newcircle, newpara)
    append!(data, [finalnewcircle])
    manifold
end

function growmanifold(prob::VTwoDManifoldProblem, disk, N; interp=QuadraticInterpolation)
    manifold = initialize(prob, disk, interp=interp)
    for i in 1:N
        grow!(manifold; interp=interp)
    end
    manifold
end
