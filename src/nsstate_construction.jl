"""
    NSState{N,T<:Number} <: AbstractVector{T}

The struct `NSState` is to record the events data for a time-T-map.

# Fields
- `state` the final state of the time-T-map;
- `event_at` is a integer vector that contains the history of the events happened.

The construction of `NSState` allows to interpolate vectors consisting of `NSState`. Currently, `LinearInterpolation`, `CubicSpline`, and
`QuadraticInterpolation` in DataInterpolations.jl are supported.
"""
struct NSState{N,T<:Number} <: AbstractVector{T}
    state::SVector{N,T}
    event_at::Vector{Int64}
end

function NSState(v::SVector{N,T}) where {N,T}
    NSState(v, Int64[])
end


# Further definitions allow to interpolate vectors consisting of NSState.
# LinearInterpolation, CubicSpline, QuadraticInterpolation
function NSState{N,T}(v::NSState{N,T}) where {N,T}
    v
end

+(a::NSState, b::NSState) = a.state + b.state

-(a::NSState, b::NSState) = a.state - b.state

*(a::Number, b::NSState) = a * b.state

*(a::NSState, b::Number) = *(b, a)

/(a::NSState, b::Number) = a.state / b

Base.length(v::Type{NSState{N,T}}) where {N,T} = N

Base.eltype(v::NSState{N,T}) where {N,T} = T

Base.getindex(v::NSState{N,T}, k) where {N,T} = v.state[k]

Base.size(v::NSState{N,T}) where {N,T} = (N,)
