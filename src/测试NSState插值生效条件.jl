module test
using DataInterpolations, StaticArrays

import Base: -,+,*,/, length, eltype, getindex, size, show

struct NSState{N,T<:Number} <: AbstractVector{T}
    state::SVector{N,T}
    event_t::Vector{T}
    event_state::Vector{SVector{N,T}}
    event_at::Vector{Int64}
end

function NSState{N,T}(v::NSState{N,T}) where {N,T}
    v
end

+(a,b)=a.state+b.state

-(a,b)=a.state-b.state

*(a::Number, b::NSState) = a * b.state

*(a::NSState, b::Number) = *(b,a)

/(a::NSState, b::Number) = a.state / b

length(v::NSState{N,T}) where {N,T} = N

eltype(v::NSState{N,T}) where {N,T} = T

getindex(v::NSState{N,T},k) where {N,T} = v.state[k]

size(v::NSState{N,T}) where {N,T} = (N,)



aa = [NSState(SA[1.0,i],Float64[],SVector{2,Float64}[],Int[]) for i in 1:20]

s0 = [(i-1)/19 for i in 1:20]

# 支持的插值类型: LinearInterpolation, CubicSpline, QuadraticInterpolation
result=CubicSpline(aa,s0)

typeof

result(0.2)

Vector{NSState{2,Float64}}(undef,3)
end

test.result.u