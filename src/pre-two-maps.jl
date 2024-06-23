module TwoDManifold
using StaticArrays, LinearAlgebra, NearestNeighbors, DataInterpolations

import DataInterpolations: LinearInterpolation

import Base: -, length, eltype, getindex, size, show

import NearestNeighbors: Metric, evaluate, euclidean

struct State{N,T<:Number} <: AbstractVector{T}
    state::SVector{N,T}
    s::T
end

function -(a, b)
    a.state - b.state
end

length(v::State{N,T}) where {N,T} = N

length(v::Type{State{N, T}}) where {N,T} = N

size(v::State{N,T}) where {N,T} = (N,)

eltype(v::State{N,T}) where {N,T} = T

getindex(v::State{N,T},k) where {N,T} = v.state[k]

function similar(v::State{N,T}) where {N,T}
    newv=similar(v.state)
    State(newv,T(0))
end

struct WEuclidean <: Metric end

function evaluate(d::WEuclidean, a,b)
    norm(a-b)
end

function State(v::SVector{N,T}, t::M) where {M<:Number,N,T<:Number}
    datetype = promote_type(M, T)
    newv = convert(SVector{N,datetype}, v)
    newt = convert(datetype, t)
    State(newv, newt)
end

struct IterationCurve{N,T<:Number,P}
    states::Vector{State{N,T}}
    pcurve::P
end

struct Annulus{N,T<:Number,P}
    circles::Vector{IterationCurve{N,T,P}}
end

function show(io::IO, v::IterationCurve)
    n = length(v.states)
    print(io, "IterationCurve with $n points")
end

function show(io::IO, v::Annulus)
    n = length(v.data)
    print(io, "$n IterationCurves")
end

function take_s(x)
    x.s
end

function take_state(x)
    x.state
end

function LinearInterpolation(v::Vector{State{N,T}}) where {N,T<:Number}
    LinearInterpolation(take_state.(v), take_s.(v))
end

@inline function distence(point, someset)
    btree = BallTree(someset, WEuclidean())
    nn(btree,point)[2]
end


# 为了节约计算时间, 不要内边界, 只要外边界; 下一次迭代时, 外边界就会成为下个环的外边界
# 为了使得离现在的内边界最近的圈距离内边界距离不是很远, 前一步的外边界要加到像上,距离判断完毕后再舍去即可.
function inintialise_mesh(f, p, saddle, v1, v2, λ1, λ2, n, r, k, N)
    newv1 = normalize(v1)
    newv2 = normalize(v2)
    rag = range(0, 1, length=n + 1)
    res = (λ2/λ1)^(N+1)
    datatype = typeof(saddle[1])
    circle1 = [State(saddle, t) for t in rag]
    data = [circle1]
    for _ in 1:k
        append!(data,[deepcopy(circle1)])
    end
    radius_range = range(0, 1, length=k + 1)
    for i in 1:k+1, j in 1:n+1
        data[i][j] = State(saddle + r * radius_range[i] * cospi(2 *rag[j]) * newv1 + 
        r *res* radius_range[i] * sinpi(2 * rag[j]) * newv2, rag[j])
    end
    newdata = deepcopy(data)
    for i in eachindex(data)
        for j in eachindex(data[i])
            newdata[i][j] = State(f(data[i][j].state, p), data[i][j].s)
        end
    end
    d = distence(State(saddle,datatype(0)), data[end])
    j = 1
    while distence(State(saddle,datatype(0)), newdata[j]) < d
        j = j + 1
    end
    result = newdata[j:end]
    append!(result, [data[end]])
    NN = length(result)
    final_result = [IterationCurve(result[1], LinearInterpolation(result[1]))]
    for i in 2:NN
        append!(final_result, [IterationCurve(result[i], LinearInterpolation(result[i]))])
    end
    para_data = [IterationCurve(circle1, LinearInterpolation(circle1))]
    for i in 2:k+1
        append!(para_data, [IterationCurve(data[i], LinearInterpolation(data[i]))])
    end
    [Annulus(para_data), Annulus(final_result)]
end

function grow_manifold!(annuluses, f, p, δ1, δ2)
    data = annuluses[end].circles
    newdata = deepcopy(data)
    k = length(newdata)
    for i in 1:k-1
        states = copy(data[i].states)
        ss = take_s.(data[i].states)
        for j in eachindex(data[i].states)
            states[j] = State(f(data[i].pcurve.u[j],p), ss[j])
        end
        newdata[i] = IterationCurve(states, LinearInterpolation(take_state.(states),ss))
    end
    prepend!(newdata, [data[end]])
    newdata
end

function fmap(X, p)
    a, b, c, α, β = p
    x, y, z = X
    SA[a*x, b*y, -a*(α^2)*x-(b^2)*(y^2)*β+c*(z+(α^2)*x+β*(y^2))]
end

para = SA[2.1, 15.3, 0.6, 0.2, -0.1]

v1, v2 = SA[1.0, 0, 0], SA[0.0, 1, 0]


result = inintialise_mesh(fmap, para, SA[0.0,0.0,0.0], v1,v2,1.0,1.0,100,0.05,10,3)

result[end].circles
grow_manifold!(result,fmap,para,0.1,0.1)

end

