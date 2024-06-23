module ExtendState
using StaticArrays, LinearAlgebra, NearestNeighbors

import Base: length, eltype, getindex, size

import NearestNeighbors: Metric, evaluate, euclidean

struct State{N,T<:Number} <: AbstractVector{T}
    state::SVector{N,T}
    s::T
end

length(v::Type{State{N, T}}) where {N,T} = N

size(v::State{N,T}) where {N,T} = (N,)

eltype(v::State{N,T}) where {N,T} = T

getindex(v::State{N,T},k) where {N,T} = v.state[k]

function similar(v::State{N,T}) where {N,T}
    newv=similar(v.state)
    State(newv,T(0))
end


function State(v::SVector{N,T}, t::M) where {M<:Number,N,T<:Number}
    datetype = promote_type(M, T)
    newv = convert(SVector{N,datetype}, v)
    newt = convert(datetype, t)
    State(newv, newt)
end

struct WEuclidean <: Metric end

# 这里的(a,b)是根据getindex方法得到的向量SVector, 而不是自己定义的类型 State.
function (d::WEuclidean)(a,b)
    Euclidean(0.0)(a,b)
end

data = [State(SA[1.0,2.0],0.0) for i in 1:12]

tree = BallTree(data, WEuclidean())

point = State(SA[2.0,2.0],0.0)

nn(tree,point)
    
end

BallTree




