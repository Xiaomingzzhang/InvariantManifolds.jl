"""
    PiecewiseV

A callable struct to represent a piecewise smooth vector field. 

# Fields
- `fs` is a tuple of smooth vector fields in different regions.
- `regions` is a tuple of the region functions: `(r1,r2,...)`, where `r1(x,p,t)` should return a Bool value to indicate that `x` is in this region or not.
- `hypers` is a tuple of the hyper surfaces separating the regions.
- `n` is a integer to switch between vector fields. Default to be zero.

# Example

```julia
using StaticArrays, InvariantManifolds
f1(x,p,t)=SA[x[2],-2x[1]]
f2(x,p,t)=SA[x[2],-x[1]]
dom1(x,p,t)=x[1]>0
dom2(x,p,t)=x[2]<0
hyper(x,p,t)=x[1]
PiecewiseV((f1,f2),(dom1,dom2),(hyper,))
```

The above codes generate a piecewise smooth vector field, which when `x[1]>0` is `f1`, and when `x[2]<0` is `f2`.
The hyper surface separating these smooth vector fields is `x[1]=0`.
"""
mutable struct PiecewiseV{F1,F2,F3}
    fs::F1
    regions::F2
    hypers::F3
    n::Int
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
    println(io)
    print(io, "n: ")
    show(io, m, A.n)
end

id(x, p, t) = x

function PiecewiseV(f1, f2, f3)
    PiecewiseV(f1, f2, f3, 0)
end

function (v::PiecewiseV)(x, p, t)
    n = v.n
    v.fs[n](x, p, t)
end


"""
    BilliardV

A callable struct to represent a vector field with multiple hyper surfaces such that the flow jump when hits these hyper surfaces.

# Fields
- `f` is the vector field, of type `f(x,p,t)`, and its output is a SVector;
- `hypers` is tuple of hyper surfaces:`(h1,h2,...)`, `h1(x,p,t)`;
- `irules` is tuple of rules on hyper surfaces:`(r1,r2,r3,...)`.
"""
struct BilliardV{F1,F2,F3}
    f::F1
    hypers::F2
    rules::F3
end


function (v::BilliardV)(x, p, t)
    v.f(x, p, t)
end

"""
    PiecewiseImpactV

A callable struct to represent a vector field with both piecewise non-smoothness and impacts.

- `fs` is a tuple of smooth vector fields in different regions.
- `regions` is a tuple of the region functions: `(r1,r2,...)`, where `r1(x,p,t)` should return a Bool value to indicate that `x` is in this region or not.
- `hypers` is a tuple of the hyper surfaces separating the regions.
- `rules` is a tuple of rules on hyper surfaces:`(r1,r2,r3,...)`. Note that for hypersurfaces that only switch between two vector fields, we can set `r1=id`.
- `idxs` is a vector of integer to indicate hypersurfaces with impact effects.
- `n` is a integer to switch between vector fields. Default to be zero.
"""
mutable struct PiecewiseImpactV{F1,F2,F3,F4}
    fs::F1
    regions::F4
    hypers::F2
    rules::F3
    idxs::Vector{Int}
    n::Int
end

function (v::PiecewiseImpactV)(x, p, t)
    n = v.n
    v.fs[n](x, p, t)
end

function PiecewiseImpactV(fs, regions, hypers,  rules, idxs)
    PiecewiseImpactV(fs, regions, hypers, rules, idxs, 0)
end


"""
    NSSetUp{T}

`NSSetUp` is a struct to contain all the information needed in continuing the manifold of non-smooth ODE.

# Fields
- `f` the Non-smooth vector field, like [`PiecewiseV`](@ref);
- `timespan` the time span of time-T-map;
- `timetmap` the time-t-map of non-smooth ODE, which maps a [`NSState`](@ref) and parameters of ODE to a [`NSState`](@ref).
"""
struct NSSetUp{F1,T,F2}
    f::F1
    timespan::Tuple{T,T}
    timetmap::F2
end

function show(io::IO, m::MIME"text/plain", A::NSSetUp{T}) where {T}
    printstyled(io, "NSSetUp:"; color=:cyan)
    println(io)
    print(io, "f: ")
    show(io, A.f)
    println(io)
    print(io, "timespan: ")
    show(io, m, A.timespan)
end


"""
    NSSolution{N,T<:Number}

The `NSSolution` is a struct to contain all information of the solution of a non-smooth ODE system.

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

"""
    Saddle{N,T,S}

`Saddle` is a struct to contain the information of a saddle point needed in continuing the manifold of non-smooth ODE.
For an ODE's saddle, this struct can be constructed by the function [`findsaddle`](@ref).
# Fields
- `saddle` the location of the saddle point;
- `unstable_directions` the unstable directions;
- `unstable_eigen_values` eigenvalues of the linearized map in the saddle at the unstable eigenvectors.
"""
struct Saddle{N,T,S}
    saddle::SVector{N,T}
    unstable_directions::Vector{SVector{N,S}}
    unstable_eigen_values::Vector{S}
end

function show(io::IO, m::MIME"text/plain", A::Saddle{N,T,S}) where {N,T,S}
    println(io, "Saddle{$N,$T,$S}: ")
    print(io, "saddle: ")
    show(io, A.saddle)
    println(io)
    print(io, "unstable_directions: ")
    show(io, A.unstable_directions)
    println(io)
    print(io, "unstable_eigen_values: ")
    show(io, A.unstable_eigen_values)
end

