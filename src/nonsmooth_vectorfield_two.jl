"""
    NSVTwoDManifoldProblem{F,T}

`NSVTwoDManifoldProblemVTwoDManifoldProblem` is a struct to contain the main information for continuing the two-dimensional manifold of an autonomous vector field.
# Fields
- `f` the `NSSetUp` of a nonsmooth vector field;
- `para` the parameters of the time flow map;
- `amax` the maximum angle between points when continuing the manifold;
- `d` the maximum distance between points when continuing the manifold;
- `ϵ` the max value of the following expression: ``\\max\\{|H(x_0,T)|,H(x_1,T)\\},`` where ``H(x,t)`` is the hypersurface of the manifold cross, ``x_0`` and ``x_1`` are points before and after the cross, ``T`` is the end of the time-``T``-map (from 0 to ``T``).
- `dsmin` the minimum arc length allowing; note that if in a continuation point, this value is achieved and the angle as well as the distance values are not achieved, then we will record this point as a [`FlawPoint`](@ref).

Convenient consturctors are `NSVTwoDManifoldProblem(f)` and `NSVTwoDManifoldProblem(f,para)`
"""
struct NSVTwoDManifoldProblem{F,T}
    f::F
    para::Vector{T}
    amax::T
    d::T
    ϵ::T
    dsmin::T
end

function NSVTwoDManifoldProblem(f; amax=0.5, d=0.001, ϵ=0.00001, dsmin=1e-5)
    NSVTwoDManifoldProblem(f, Float64[], amax, d, ϵ, dsmin)
end


function NSVTwoDManifoldProblem(f, para::AbstractVector{T};
    amax=T(0.5), d=T(0.001), ϵ=T(0.00001), dsmin=T(1e-5)) where {T}
    NSVTwoDManifoldProblem(f, para, amax, d, ϵ, dsmin)
end

function show(io::IO, m::MIME"text/plain", A::NSVTwoDManifoldProblem)
    printstyled(io, "NSVTwoDManifoldProblem:"; color=:cyan)
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
    NSVTwoDManifold{F,S,N,T}

`NSVTwoDManifold` is a struct contains all the information of the non-smooth two-dimensional numerical manifold of an autonomous vector field.
# Fields
- `prob` the problem [`NSVTwoDManifoldProblem`](@ref);
- `data` the numerical data that should be `Vector{Vector{S}}`, where `S` is the interpolation curve (we use `DataInterpolation` in this package);
- `flawpoints` the flaw points generated during continuation.
"""
mutable struct NSVTwoDManifold{F,S,N,T}
    prob::NSVTwoDManifoldProblem{F,T}
    data::Vector{Vector{S}}
    flawpoints::Vector{FlawPoint{N,T}}
end


function Base.show(io::IO, m::MIME"text/plain", A::NSVTwoDManifold)
    m = 0
    n = length(A.data)
    for i in eachindex(A.data)
        for j in eachindex(A.data[i])
            m = m + length(A.data[i][j].t)
        end
    end
    k = length(A.flawpoints)
    printstyled(io, "Non-smooth two-dimensional manifold"; bold=true, color=:cyan)
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

function initialize(prob::NSVTwoDManifoldProblem, disk::Vector{Vector{SVector{N,T}}}; interp=LinearInterpolation) where {N,T}
    flawpoints = FlawPoint{N,T}[]
    newdisk = Vector{NSState{N,T}}[]
    for i in eachindex(disk)
        points = Vector{NSState{N,T}}(undef, length(disk[i]))
        for j in eachindex(points)
            points[j] = NSState(disk[i][j])
        end
        append!(newdisk, [points])
    end
    circles = [[paramise(newdisk[i], interp=interp)] for i in eachindex(newdisk)]
    NSVTwoDManifold(prob, circles, flawpoints)
end

function grow!(manifold::NSVTwoDManifold{F,S,N,T}; interp=LinearInterpolation) where {F,S,N,T}
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
    manifold
end

function growmanifold(prob::NSVTwoDManifoldProblem, disk, N; interp=LinearInterpolation)
    manifold = initialize(prob, disk, interp=interp)
    for i in 1:N
        grow!(manifold, interp=interp)
    end
    manifold
end

