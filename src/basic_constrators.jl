"""
    State{N,T<:Number}

`State` is a struct contain point data in a numerical manifold.

# Fields
- `state` a `SVector` reprents the point in phase space;
- `s` the parameter of this point.
"""
struct State{N,T<:Number}
    state::SVector{N,T}
    s::T
end

function State(v::SVector{N,T}, t::M) where {M<:Number,N,T<:Number}
    datetype = promote_type(M, T)
    newv = convert(SVector{N,datetype}, v)
    newt = convert(datetype, t)
    State(newv, newt)
end

function show(io::IO, m::MIME"text/plain", A::State{N,T}) where {N,T}
    println(io, "State{$N, $T}: ")
    print(io, "state: ")
    show(io, m, A.state)
    println(io)
    print(io, "s: ")
    show(io, m, A.s)
end


"""
    IterationCurve{N,T<:Number}

`IterationCurve` is a struct that contains line data in a numerical manifold.

# Fields
- `states` a `Vector{State}` contains the points in phase space;
- `pcurve` the paramtric curve of these `states`, which is a `LinearInterpolation`.
"""
struct IterationCurve{N,T<:Number}
    states::Vector{State{N,T}}
    pcurve
end

function show(io::IO, v::IterationCurve{N,T}) where {N,T<:Number}
    n = length(v.states)
    print(io, "IterationCurve{$N,$T} with $n points")
end

function take_s(x)
    x.s
end

function take_state(x)
    x.state
end

function LinearInterpolation(v::Vector{State{N,T}}) where {N,T<:Number} # Vector{T} 中的T必须是具体的
    LinearInterpolation(take_state.(v), take_s.(v))
end

function -(a, b)
    a.state - b.state
end


"""
    NSState{N,T<:Number}

The struct `NSState` is to record the events data for a time-T-map.

# Fields
- `state` the final state of the time-T-map;
- `event_t` the times when events happen;
- `event_state` the solution's state when events happen;
- `event_at` is a vector that contains integers indicating which event happen;
- `s` the parameter of the `NSState`'s preimage.

The meaning of `event_at` in `PiecewiseV` and `BilliardV` systems is quite clear.
For a simple Fillippov system `SFilippovV`, we record three events:
- `1` reprents the event that flow enters the slding surface from `H<0` or `H>0`;
- `2` reprents the event that flow cross the slding surface;
- `3` reprents the event that flow slides out the slding surface from `H=0` to `H<0` or `H>0`.
"""
struct NSState{N,T<:Number}
    state::SVector{N,T}
    event_t::Vector{T}
    event_state::Vector{SVector{N,T}}
    event_at::Vector{Int64}
    s::T
end

function show(io::IO, m::MIME"text/plain", A::NSState{N,T}) where {N,T}
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
    println(io)
    print(io, "s: ")
    show(io, m, A.s)
end

"""
    NSSolution{N,T<:Number}

The `NSSolution` is a struct to contain all information of the solution of a nonsmooth ODE system.

# Fields
- `sol` `ODESolution` solved by `OrdinaryDiffEq`;
- `event_t` the times when events happen;
- `event_state` the solution's state when events happen;
- `event_at` is a vector that contains integers indicating which event happen.
"""
struct NSSolution{N,T<:Number}
    sol
    event_t::Vector{T}
    event_state::Vector{SVector{N,T}}
    event_at::Vector{Int64}
end

function show(io::IO, m::MIME"text/plain", A::NSSolution{N,T}) where {N,T}
    println(io, "NSSolution{$N, $T}: ")
    show(io, m, A.sol)
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

abstract type ContinuousVectorField end

abstract type JumpVectorField end

"""
    PiecewiseV

Piecewise smooth vector field. 

# Fields
- `fs` is a vector of smooth vector fields in different regions.
- `regions` is a vector of the region functions: `[r1,r2,...]`, where `r1(x,p,t)` should return a Bool value to indicate that `x` is in this region or not.
- `hypers` is a vector of the hypersurfaces separating the regions.

# Example

```julia
using StaticArrays, InvariantManifolds
f1(x,p,t)=SA[x[2],-2x[1]]
f2(x,p,t)=SA[x[2],-x[1]]
dom1(x,p,t)=x[1]>0
dom2(x,p,t)=x[2]<0
hyper(x,p,t)=x[1]
PiecewiseV([f1,f2],[dom1,dom2],[hyper])
```

The above codes generate a piecewise smooth vector field, which when `x[1]>0` is `f1`, and when `x[2]<0` is `f2`.
The hypersurface separating these smooth vector fields is `x[1]=0`.
"""
struct PiecewiseV <: ContinuousVectorField
    fs::Vector{Function}
    regions::Vector{Function}
    hypers::Vector{Function}
end

function show(io::IO, m::MIME"text/plain", A::PiecewiseV)
    println(io, "PiecewiseV: ")
    print(io, "fs: ")
    show(io, m, A.fs)
    println(io)
    print(io, "regions: ")
    show(io, m, A.regions)
    println(io)
    print(io, "hypers: ")
    show(io, m, A.hypers)
end

function (v::PiecewiseV)(x, p, t)
    n = Int(p[end])
    v.fs[n](x, p, t)
end


"""
    BilliardV

A vector field with multiple hypersurfaces such that the flow jump when hits these hypersurfaces.

# Fields
- `f` is the vector field, of type `f(x,p,t)`, and its output is a SVector;
- `hypers` is vector of hypersurfaces:`[h1,h2,...]`, `h1(x,p,t)`;
- `irules` is vector of rules on hypersurfaces:`[r1,r2,r3,...]`; the rules must be independent of time `t`, i.e., `r1(x,p)`.
"""
struct BilliardV <: JumpVectorField
    f
    hypers::Vector{Function}
    rules::Vector{Function}
end

function (v::BilliardV)(x, p, t)
    v.f(x, p, t)
end


"""
    SFilippovV

`SFilippovV` means simple Fillippov vector fields, which means that there only exists one hypersurface to separate the phase space.

# Fields
- `fs` vector fields in two sides of hypersurface. The slide vector field can be generated automatically.
- `hyper` hypersurface.
- `dhyper` grad of hypersurface. Warn!!! The grad must point to the second of `fs`.
- `exit` conditions to exit the hypersurface, which can also be generated automatically
"""
struct SFilippovV <: ContinuousVectorField
    fs::SVector{3,Function}
    hyper
    dhyper
    exit
end

function SFilippovV(fs, h, ∇h)
    f1 = fs[1]
    f2 = fs[2]
    function sv(x, p, t)
        α = dot(∇h(x, p, t), f1(x, p, t)) / dot(∇h(x, p, t), f1(x, p, t) - f2(x, p, t))
        (1 - α) * f1(x, p, t) + α * f2(x, p, t)
    end
    function exit(x, p, t)
        dot(∇h(x, p, t), f2(x, p, t)) * dot(∇h(x, p, t), f1(x, p, t))
    end
    SFilippovV(SA[f1, f2, sv], h, ∇h, exit)
end

function (v::SFilippovV)(x, p, t)
    n = Int(p[end])
    v.fs[n](x, p, t)
end


"""
    NSSetUp{T}

`NSSetUp` is a struct to contatin all the information needed in continuing the manifold of nonsmooth ODE.

# Fields
- `f::T` the Non-smooth vector field, like `PiecewiseV`;
- `timespan` the time span of time-T-map;
- `timetmap` the time-t-map of nonsmooth ODE, which maps a `State` and parameters of ODE to a `NSState`.
# Warn
In a `ContinuousVectorField`,this parameters must be long than the real parameter.
For example, if the parameter of your system is `[0.1,2.0]`, than you must set it to `[0.1,2.0,1.0]`.
The value of the last parameter is meaninglless, just for switching the vector fields.
"""
struct NSSetUp{T}
    f::T
    timespan
    timetmap
end

function show(io::IO, m::MIME"text/plain", A::NSSetUp{T}) where{T}
    println(io, "NSSetUp{$T}: ")
    print(io, "f: ")
    show(io, m, A.f)
    println(io)
    print(io, "timespan: ")
    show(io, m, A.timespan)
    println(io)
    print(io, "timetmap: ")
    show(io, m, A.timetmap)
end
