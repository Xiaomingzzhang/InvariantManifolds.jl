"""
    TwoDManifoldProblem{F,T}

`TwoDManifoldProblem` is a struct to contain the main information for continuing the two-dimensional manifold of a nonlinear map.
# Fields
- `f` the nonlinear map, which should has the form `f(x,p)` and return a `SVector`;
- `para` the parameters of the nonlinear map;
- `amax` the maximum angle between points when continuing the manifold;
- `d` the maximum distance between points when continuing the manifold;
- `dcircle` the maximum distance between circles when continuing the manifold;
- `dsmin` the minimum arc length allowing; note that if in a continuation point, this value is achieved and the angle as well as the 
distance values are not achieved, then we will record this point as a [`FlawPoint`](@ref).

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

"""
    TwoDManifold{F,S,N,T}

`TwoDManifold` is a struct contains all the information of the two-dimensional numerical manifold.
# Fields
- `prob` the problem `TwoDManifoldProblem`;
- `data` the numerical data that should be `Vector{Vector{Vector{S}}}`, where `S` is the interpolation curve 
(we use `DataInterpolation` in this package);
- `flawpoints` the flaw points generated during continuation.
"""
mutable struct TwoDManifold{F,S,N,T}
    prob::TwoDManifoldProblem{F,T}
    data::Vector{Vector{S}}
    flawpoints::Vector{FlawPoint{N,T}}
end

function TwoDManifoldProblem(f; amax=0.5, d=0.001, dcircle=0.01, dsmin=1e-5)
    TwoDManifoldProblem(f, Float64[], amax, d, dcircle, dsmin)
end


function TwoDManifoldProblem(f, para::Vector{T};
    amax=T(0.5), d=T(0.001), dcircle=T(0.01), dsmin=T(1e-5)) where {T}
    TwoDManifoldProblem(f, para, amax, d, dcircle, dsmin)
end


function Base.show(io::IO, m::MIME"text/plain", A::TwoDManifold)
    m = 0
    n = 0
    for i in eachindex(A.data)
        n = n + length(A.data[i])
        for j in eachindex(A.data[i])
            m = m + length(A.data[i][j].t)
        end
    end
    println(io, "Two dimensional manifold with $n circles and $m points")
end

"""
    kd_distence

The function to measure the distance between two circles by using the package `NearestNeighbors`.
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
    disk(p, times)

`disk` is a function to generate circles around the saddle, which represented as the local manifold.
# Parameters
- `p` the struct `Saddle` which should contains two unstable directions; the complex eigenvalues and eigenvectors
are allowed.
- `times` the iteration times; this parameter is needed to adjust the torsion in different directions in the 
process of continuation.
# Keyword argument
- `n` the number of point in each circle, default to be `150`;
- `d` the max distance between points in a single circle, default to be `0.002`;
- `r` the size of the disk, default to be `0.05`;
- `circles` the number of the circles, default to be `10`.
"""
function disk(p::Saddle{N,T,S}, times; n=150, d=0.002, r=0.05, circles=10) where {N,T,S}
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

function addcircles!(f, para, d, circles, dsmin, αmax, dcircle, flawpoints; interp=LinearInterpolation)
    newdata = deepcopy(circles)
    k = length(newdata)
    for i in 1:k
        states = similar(circles[i].u)
        ss = copy(circles[i].t)
        for j in eachindex(circles[i].u)
            states[j] = f(circles[i].u[j], para)
        end
        newdata[i] = interp(states, ss)
    end
    # first interpolate every circle
    for j in 1:k
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

function inintialise(prob::TwoDManifoldProblem, disk::Vector{Vector{SVector{N,T}}}; interp=LinearInterpolation) where {N,T}
    para = prob.para
    f = prob.f
    αmax = prob.amax
    d = prob.d
    dsmin = prob.dsmin
    flawpoints = FlawPoint{N,T}[]
    dcircle = prob.dcircle
    circles = [paramise(disk[i], interp=interp) for i in eachindex(disk)]
    newcircles = addcircles!(f, para, d, circles, dsmin, αmax, dcircle, flawpoints; interp=interp)
    circle0 = first(disk)
    circlen = last(disk)
    j = 1
    dist = kd_distence(circle0, circlen)
    while kd_distence(newcircles[j].u, circle0) < dist
        j = j + 1
    end
    j = j + 1
    newcircles = newcircles[j:end]
    prepend!(newcircles, [paramise(circlen, interp=interp)])
    result = [circles, newcircles]
    TwoDManifold(prob, result, flawpoints)
end

function grow!(manifold::TwoDManifold; interp=LinearInterpolation)
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

function growmanifold(prob::TwoDManifoldProblem, disk, N; interp=LinearInterpolation)
    manifold = inintialise(prob, disk, interp=interp)
    for i in 1:N
        grow!(manifold; interp=interp)
    end
    manifold
end
