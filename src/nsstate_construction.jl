"""
    NSState{N,T<:Number} <: AbstractVector{T}

The struct `NSState` is to record the events data for a time-T-map.

# Fields
- `state` the final state of the time-T-map;
- `event_t` the times when events happen;
- `event_state` the solution's state when events happen;
- `event_at` is a vector that contains integers indicating which event happen;

The meaning of `event_at` in `PiecewiseV` and `BilliardV` systems is quite clear.
For a simple Fillippov system `SFilippovV`, we record three events:
- `1` represents the event that flow enters the sliding surface from `H<0` or `H>0`;
- `2` represents the event that flow cross the sliding surface;
- `3` represents the event that flow slides out the sliding surface from `H=0` to `H<0` or `H>0`.
"""
struct NSState{N,T<:Number} <: AbstractVector{T}
    state::SVector{N,T}
    event_t::Vector{T}
    event_state::Vector{SVector{N,T}}
    event_at::Vector{Int64}
end

function Base.show(io::IO, m::MIME"text/plain", A::NSState{N,T}) where {N,T}
    println(io, "NSState{$N, $T}: ")
    print(io, "state: ")
    show(io, m, A.state)
    println(io)
    print(io, "event_t: ")
    show(io, m, A.event_t)
    println(io)
    print(io, "event_state: ")
    show(io, m, A.event_state)
    println(io)
    print(io, "event_at: ")
    show(io, m, A.event_at)
end

# Further definitions allow to interpolate a NSState.
# LinearInterpolation, CubicSpline, QuadraticInterpolation

function NSState{N,T}(v::NSState{N,T}) where {N,T}
    v
end

+(a, b) = a.state + b.state

-(a, b) = a.state - b.state

*(a::Number, b::NSState) = a * b.state

*(a::NSState, b::Number) = *(b, a)

/(a::NSState, b::Number) = a.state / b

Base.length(v::NSState{N,T}) where {N,T} = N

Base.eltype(v::NSState{N,T}) where {N,T} = T

Base.getindex(v::NSState{N,T}, k) where {N,T} = v.state[k]

Base.size(v::NSState{N,T}) where {N,T} = (N,)