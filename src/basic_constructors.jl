abstract type ContinuousVectorField end

abstract type JumpVectorField end

"""
    PiecewiseV <: ContinuousVectorField

Piecewise smooth vector field. 

# Fields
- `fs` is a vector of smooth vector fields in different regions.
- `regions` is a vector of the region functions: `[r1,r2,...]`, where `r1(x,p,t)` should return a Bool value to indicate that `x` is in this region or not.
- `hypers` is a vector of the hyper surfaces separating the regions.
- `n` is a integer to switch between vector fields. It can be set to any integer when construct a `PiecewiseV`.

# Example

```julia
using StaticArrays, InvariantManifolds
f1(x,p,t)=SA[x[2],-2x[1]]
f2(x,p,t)=SA[x[2],-x[1]]
dom1(x,p,t)=x[1]>0
dom2(x,p,t)=x[2]<0
hyper(x,p,t)=x[1]
PiecewiseV((f1,f2),(dom1,dom2),(hyper,),1)
```

The above codes generate a piecewise smooth vector field, which when `x[1]>0` is `f1`, and when `x[2]<0` is `f2`.
The hyper surface separating these smooth vector fields is `x[1]=0`.
"""
mutable struct PiecewiseV{F1,F2,F3,F4,F5} <: ContinuousVectorField
    fs::F1
    regions::F2
    hypers::F3
    rules::F4
    mirrors::F5
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
    print(io, "rules: ")
    show(io, m, A.rules)
    println(io)
    print(io, "n: ")
    show(io, m, A.n)
end

id(x, p, t) = x

function PiecewiseV(f1, f2, f3)
    n = length(f3)
    ids = ntuple(_ -> id, n)
    PiecewiseV(f1, f2, f3, ids, ids, 0)
end

function (v::PiecewiseV)(x, p, t)
    n = v.n
    v.fs[n](x, p, t)
end


"""
    BilliardV <: JumpVectorField

A vector field with multiple hyper surfaces such that the flow jump when hits these hyper surfaces.

# Fields
- `f` is the vector field, of type `f(x,p,t)`, and its output is a SVector;
- `hypers` is tuple of hyper surfaces:`(h1,h2,...)`, `h1(x,p,t)`;
- `irules` is tuple of rules on hyper surfaces:`(r1,r2,r3,...)`.
"""
struct BilliardV{F1,F2,F3,F4} <: JumpVectorField
    f::F1
    hypers::F2
    rules::F3
    mirrors::F4
end

function BilliardV(f, hypers, rules)
    BilliardV(f, hypers, rules, nothing)
end

function (v::BilliardV)(x, p, t)
    v.f(x, p, t)
end


"""
    SFilippovV <: ContinuousVectorField

`SFilippovV` means simple Fillippov vector fields, which means that there only exists one hyper surface to separate the phase space.

# Fields
- `fs` vector fields in two sides of hyper surface. The slide vector field can be generated automatically.
- `hyper` hyper surface.
- `dhyper` grad of hyper surface. Warn!!! The grad must point to the second of `fs`.
- `exit` conditions to exit the hyper surface, which can also be generated automatically.
- `n` is a integer to switch between vector fields. It can be set to any integer when construct a `SFilippovV`.
"""
mutable struct SFilippovV{F,H,DH,E,G,Q} <: ContinuousVectorField
    fs::F
    hyper::H
    dhyper::DH
    exit::E
    rules::G
    mirrors::Q
    n::Int
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
    ids = (id, id, id)
    SFilippovV((f1, f2, sv), h, ∇h, exit, ids, ids, 0)
end

function (v::SFilippovV)(x, p, t)
    n = v.n
    v.fs[n](x, p, t)
end

# struct HyperPlaneMirror{T}
#     idx::Int64
#     value::T
# end

# function (v::SFilippovV)(x, p, t)
    
# end


"""
    NSSetUp{T}

`NSSetUp` is a struct to contain all the information needed in continuing the manifold of non-smooth ODE.

# Fields
- `f` the Non-smooth vector field, like `PiecewiseV`;
- `timespan` the time span of time-T-map;
- `timetmap` the time-t-map of non-smooth ODE, which maps a `State` and parameters of ODE to a `NSState`.
"""
struct NSSetUp{T}
    f::T
    timespan
    timetmap
end

function show(io::IO, m::MIME"text/plain", A::NSSetUp{T}) where {T}
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
