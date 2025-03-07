"""
    TwoDManifoldProblem{F,T}

`TwoDManifoldProblem` is a struct to contain the main information for continuing the two-dimensional manifold of a nonlinear map.
# Fields
- `f` the nonlinear map, which should has the form `f(x,p)` and return a `SVector`;
- `para` the parameters of the nonlinear map;
- `amax` the maximum angle between points when continuing the manifold;
- `d` the maximum distance between points when continuing the manifold;
- `dcircle` the maximum distance between circles when continuing the manifold;
- `dsmin` the minimum arc length allowing; note that if in a continuation point, this value is achieved and the angle as well as the distance values are not achieved, then we will record this point as a [`FlawPoint`](@ref).

Convenient consturctors are `TwoDManifoldProblem(f)` and `TwoDManifoldProblem(f,para)`
"""
struct TwoDManifoldProblem{F,T}
    f::F
    para::Vector{T}
    amax::T
    d::T
    dcircle::T
    dsmin::T
end

function show(io::IO, m::MIME"text/plain", A::TwoDManifoldProblem)
    printstyled(io, "TwoDManifoldProblem:"; color=:cyan)
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
    print(io, "dcircle: ")
    show(io, m, A.dcircle)
    println(io)
    print(io, "dsmin: ")
    show(io, m, A.dsmin)
end

"""
    TwoDManifold{F,S,N,T}

`TwoDManifold` is a struct contains all the information of the two-dimensional numerical manifold of a nonlinear map.
# Fields
- `prob` the problem [`TwoDManifoldProblem`](@ref);
- `data` the numerical data that should be `Vector{S}`, where `S` is the interpolation curve (we use `DataInterpolation` in this package);
- `flawpoints` the flaw points generated during continuation.
"""
mutable struct TwoDManifold{F,S,N,T}
    prob::TwoDManifoldProblem{F,T}
    data::Vector{Vector{S}}
    flawpoints::Vector{FlawPoint{N,T}}
end

function TwoDManifoldProblem(f; amax=0.5, d=1e-3, dcircle=1e-2, dsmin=1e-6)
    TwoDManifoldProblem(f, Float64[], amax, d, dcircle, dsmin)
end


function TwoDManifoldProblem(f, para::AbstractVector{T};
    amax=T(0.5), d=T(1e-3), dcircle=T(1e-2), dsmin=T(1e-6)) where {T}
    TwoDManifoldProblem(f, para, amax, d, dcircle, dsmin)
end


function Base.show(io::IO, m::MIME"text/plain", A::TwoDManifold)
    m = 0
    n = length(A.data)
    for i in eachindex(A.data)
        n = n + length(A.data[i])
        for j in eachindex(A.data[i])
            m = m + length(A.data[i][j].t)
        end
    end
    k = length(A.flawpoints)
    printstyled(io, "Two-dimensional manifold:"; bold=true, color=:cyan)
    println(io)
    printstyled(io, "Circles number: "; color=:cyan)
    println(io, "$n")
    printstyled(io, "Points number: "; color=:cyan)
    println(io, "$m")
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
    InvariantManifolds.kd_distence

The function to measure the distance between two circles by using the package `NearestNeighbors.jl`.
"""
@inline function kd_distence(point::SVector{N,T}, someset::Vector{SVector{N,T}}) where {N,T}
    btree = KDTree(someset)
    nn(btree, point)[2]
end

@inline function kd_distence(someset1::Vector{SVector{N,T}}, someset2::Vector{SVector{N,T}}) where {N,T}
    btree = KDTree(someset1)
    dd = [nn(btree, x)[2] for x in someset2]
    maximum(dd)
end

"""
    InvariantManifolds.set_distence

The function to measure the distance between two sets by using the package `NearestNeighbors.jl`.
"""
@inline function set_distence(someset1::Vector{SVector{N,T}}, someset2::Vector{SVector{N,T}}) where {N,T}
    btree = KDTree(someset1)
    dd = [nn(btree, x)[2] for x in someset2]
    @show minimum(dd)
    minimum(dd)
end

function welldistributedpoints!(pcurve, points, para, d)
    olds = collect(para)
    n = length(points)
    i = 1
    while i + 1 <= n
        if norm(points[i] - points[i+1]) < d
            i = i + 1
        else
            s0 = olds[i]
            s1 = olds[i+1]
            paras = (s0 + s1) / 2
            addps = pcurve(paras)
            insert!(points, i + 1, addps)
            insert!(olds, i + 1, paras)
            n = n + 1
        end
    end
end
"""
    gen_disk(p, times)

`gen_disk` is a function to generate circles around the saddle, which represented as the local manifold.
# Parameters
- `p` the struct [`Saddle`](@ref) which should contains two unstable directions; the complex eigenvalues and eigenvectors are allowed.
# Keyword arguments
- `times` the iteration time, default to be `1`; for the computation of invariant manifolds of nonlinear map, this parameter is needed to adjust the torsion in different directions in the process of continuation.
- `n` the number of point in each circle, default to be `150`;
- `d` the max distance between points in a single circle, default to be `0.0002`;
- `r` the size of the disk, default to be `0.05`;
- `circles` the number of the circles, default to be `10`.
"""
function gen_disk(p::Saddle{N,T,S}; times=1, n=150, d=0.0002, r=0.01, circles=2) where {N,T,S}
    if S <: Real
        v1 = p.unstable_directions[1]
        v2 = p.unstable_directions[2]
        newv1 = normalize(v1)
        newv2 = normalize(v2)
        λ1 = p.unstable_eigen_values[1]
        λ2 = p.unstable_eigen_values[2]
        res = (λ1 / λ2)^(times)
        para = range(0, 2, n)
        result = [[p.saddle, p.saddle, p.saddle]]
        for i in 1:circles-1
            circle = Vector{SVector{N,T}}(undef, n)
            c1 = r * i / (circles - 1)
            c2 = r * i * res / (circles - 1)
            pcurve(s) = p.saddle + c1 * cospi(s) * newv1 +
                        c2 * sinpi(s) * newv2
            for j in eachindex(para)
                circle[j] = pcurve(para[j])
            end
            welldistributedpoints!(pcurve, circle, para, d)
            append!(result, [circle])
        end
        result
    else
        v1 = real(p.unstable_directions[1])
        v2 = imag(p.unstable_directions[1])
        A = hcat(v1, v2)
        Q = SMatrix(qr(A).Q)
        newv1 = Q[:, 1]
        newv2 = Q[:, 2]
        para = range(0, 2, n)
        result = [[p.saddle, p.saddle, p.saddle]]
        for i in 1:circles-1
            circle = Vector{SVector{N,T}}(undef, n)
            c1 = r * i / (circles - 1)
            c2 = r * i / (circles - 1)
            for j in eachindex(para)
                circle[j] = saddle + c1 * cospi(para[j]) * newv1 +
                            c2 * sinpi(para[j]) * newv2
            end
            append!(result, [circle])
        end
        result
    end
end

"""
    InvariantManifolds.addcircles!(f, para, d, circles, dsmin, αmax, dcircle, flawpoints; interp=LinearInterpolation)

Adds and refines circles in the two-dimensional manifold computation by iterating the map and ensuring proper point distribution.

# Arguments
- `f`: The nonlinear map function
- `para`: Vector of parameters for the map
- `d`: Maximum allowed distance between points in a circle
- `circles`: Vector of interpolated curves representing the current circles
- `dsmin`: Minimum allowed arc length between points
- `αmax`: Maximum allowed angle between consecutive points
- `dcircle`: Maximum allowed distance between consecutive circles
- `flawpoints`: Vector to store problematic points during computation
- `interp`: Interpolation method (default: QuadraticInterpolation)

# Returns
A vector of new interpolated curves representing the refined circles after one iteration of the map.

# Details
The function performs two main steps:
1. Iterates each circle forward under the map and refines point distribution within each circle
2. Adds intermediate circles where the distance between consecutive circles exceeds `dcircle`

Points are added to maintain proper spacing and curvature constraints specified by `d` and `αmax`.
"""
function addcircles!(f, para, d, circles, dsmin, αmax, dcircle, flawpoints; interp=QuadraticInterpolation)
    newdata = deepcopy(circles)
    k = length(newdata)
    Threads.@threads for i in 1:k
        states = similar(circles[i].u)
        ss = copy(circles[i].t)
        for j in eachindex(circles[i].u)
            states[j] = f(circles[i].u[j], para)
        end
        newdata[i] = interp(states, ss)
    end
    # first interpolate every circle
    Threads.@threads for j in 1:k
        ic = newdata[j].u
        olds = newdata[j].t
        oldcurve = circles[j]
        newpara = addpoints!(f, para, d, oldcurve, ic, olds, dsmin, αmax, flawpoints)
        olds .= newpara
    end
    # then interpolate between circles
    p = 1
    copydata = deepcopy(circles)
    @inbounds while p + 1 <= k
        dd = kd_distence(newdata[p].u, newdata[p+1].u)
        if dd > dcircle
            insert_circle = copy(copydata[p+1].u)
            pre_insert_circle = similar(copydata[p+1].u)
            kd_tree = KDTree(copydata[p].u)
            for m in eachindex(insert_circle)
                point = insert_circle[m]
                number = nn(kd_tree, point)[1]
                point2 = copydata[p].u[number]
                pre_insert_circle[m] = (point + point2) / 2
                insert_circle[m] = f((point + point2) / 2, para)
            end
            pre_curve = paramise(pre_insert_circle, interp=interp)
            insert_s = copy(pre_curve.t)
            newpara = addpoints!(f, para, d, pre_curve, insert_circle, insert_s, dsmin, αmax, flawpoints)
            insert_curve = interp(insert_circle, newpara)
            insert!(newdata, p + 1, insert_curve)
            insert!(copydata, p + 1, pre_curve)
            k = k + 1
        else
            p = p + 1
        end
    end
    newdata
end

function initialize(prob::TwoDManifoldProblem, disk::Vector{Vector{SVector{N,T}}}; interp=QuadraticInterpolation) where {N,T}
    para = prob.para
    f = prob.f
    αmax = prob.amax
    d = prob.d
    dsmin = prob.dsmin
    flawpoints = FlawPoint{N,T}[]
    dcircle = prob.dcircle
    circles = [paramise(disk[i], interp=interp) for i in eachindex(disk)]
    outercircle = last(disk)
    poutcircle = paramise(outercircle, interp=interp)
    newoutcircle = [f(outercircle[i], para) for i in eachindex(outercircle)]
    addpoints!(f, para, d, poutcircle, newoutcircle, copy(poutcircle.t), dsmin, αmax, flawpoints)
    pnewoutcircle = paramise(newoutcircle, interp=interp)
    result = [circles, [poutcircle, pnewoutcircle]]
    TwoDManifold(prob, result, flawpoints)
end

function grow!(manifold::TwoDManifold; interp=QuadraticInterpolation)
    prob = manifold.prob
    para = prob.para
    f = prob.f
    αmax = prob.amax
    d = prob.d
    dsmin = prob.dsmin
    flawpoints = manifold.flawpoints
    dsmin = prob.dsmin
    dcircle = prob.dcircle
    data = manifold.data
    circles = data[end]
    newcircles = addcircles!(f, para, d, circles, dsmin, αmax, dcircle, flawpoints; interp=interp)
    append!(data, [newcircles])
    manifold
end

function growmanifold(prob::TwoDManifoldProblem, disk, N; interp=QuadraticInterpolation)
    manifold = initialize(prob, disk, interp=interp)
    for i in 1:N
        grow!(manifold; interp=interp)
    end
    manifold
end
