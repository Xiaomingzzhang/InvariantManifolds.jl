# monodromy.jl
#
# Purpose
# -------
# This file implements the derivative (monodromy matrix) of the time-T-map of a
# non-smooth ODE and uses it to locate saddle-type periodic orbits by Newton's
# method.  It backs the exported functions `ns_monodromy` and the non-smooth
# overloads of `findsaddle`.
#
# Main idea
# ---------
# A naive finite-difference derivative of the whole time-T-map is inaccurate near
# impacts or switching events.  Instead, we integrate an augmented system whose
# state is `(x, X)`, where `X` is the fundamental matrix of the variational
# equation `X' = DV(x,t) X`.  Whenever the flow hits a hypersurface transversally,
# the matrix part is multiplied from the left by the saltation (jump) matrix `S`.
#
# Two types of saltation matrices are implemented:
# - pure switching between two smooth vector fields: `S = I + (f_+ - f_-)∇Hᵀ/D`;
# - impact with reset map `g`: `S = Dg + (f_+ - Dg f_- - g_t)∇Hᵀ/D`,
# where `D = ∇H·f_- + H_t` is the transversal crossing speed.
#
# The returned monodromy matrix `X(T)` is then used by Newton's method
# `x -> x - (X(T)-I)^{-1}(F(x)-x)`.  Event times are stored in
# `NSMonodromyData` and can be inspected for validation.
#
# The gradients of `H`, the derivatives of the reset maps, and the derivatives of
# the smooth vector fields are obtained from user supplied Jacobians; for `H` and
# reset maps finite-difference fallbacks are used internally.

"""
    NSMonodromyData{N,T}

Internal result of integrating the augmented system consisting of the original
non-smooth ODE and the matrix variational equation.

# Fields
- `x1`: final state of the time-`T`-map;
- `D`: the derivative (monodromy matrix) of the time-`T`-map;
- `event_t`, `event_at`: event times and the indices of the hypersurfaces that
  were crossed.
"""
struct NSMonodromyData{N,T}
    x1::SVector{N,T}
    D::SMatrix{N,N,T}
    event_t::Vector{T}
    event_at::Vector{Int}
end

function Base.show(io::IO, m::MIME"text/plain", A::NSMonodromyData{N,T}) where {N,T}
    printstyled(io, "NSMonodromyData{$N,$T}:"; color=:cyan)
    println(io)
    print(io, "x1: ")
    show(io, m, A.x1)
    println(io)
    print(io, "D: ")
    show(io, m, A.D)
    println(io)
    print(io, "event_t: ")
    show(io, m, A.event_t)
    println(io)
    print(io, "event_at: ")
    show(io, m, A.event_at)
end

function _to_svector(u::AbstractVector, ::Val{N}) where {N}
    SVector{N}(ntuple(i -> u[i], N))
end

function _num_jacobian(f, x::AbstractVector)
    FiniteDiff.finite_difference_jacobian(u -> f(u), x)
end

function _num_gradient_scalar(h, x::AbstractVector, p, t)
    N = length(x)
    x0 = collect(x)
    δ = cbrt(eps(float(eltype(x0)))) * max(one(float(eltype(x0))), norm(x0))
    g = zeros(eltype(x0), N)
    @inbounds for i in 1:N
        xp = copy(x0)
        xm = copy(x0)
        xp[i] += δ
        xm[i] -= δ
        g[i] = (h(xp, p, t) - h(xm, p, t)) / (2δ)
    end
    g
end

function _num_time_derivative(f, x, p, t)
    δt = cbrt(eps(float(eltype(x)))) * max(one(float(eltype(x))), abs(t))
    (f(x, p, t + δt) - f(x, p, t - δt)) / (2δt)
end

"""
    InvariantManifolds._saltation_matrix(fpre, fpost, DG, gt, gradh, ht)

Return the saltation (jump) matrix for a transversal event at the state with
left/right vector fields `fpre` and `fpost`, reset derivative `DG` (identity for
pure switching), reset time derivative `gt` (zero for autonomous resets), and
event surface gradient/time derivative `gradh`/`ht`:

```math
S = DG + \\frac{\\left(f_+ - DG\\,f_- - g_t\\right)\\nabla H^{T}}
         {\\nabla H\\cdot f_- + H_t}.
```

This reduces to the formulas in the Floquet theory notes when `H` and `g` do not
depend explicitly on time.
"""
function _saltation_matrix(fpre, fpost, DG, gt, gradh, ht)
    denom = dot(gradh, fpre) + ht
    abs(denom) <= 10eps(float(denom)) * max(one(float(denom)), norm(gradh) * norm(fpre), abs(ht)) &&
        error("The orbit crosses the event surface tangentially; the saltation matrix is not defined")
    DG + (fpost - DG * fpre - gt) * gradh' ./ denom
end

function _monodromy_vectorfield_derivative(v::BilliardV, dv, x, p, t)
    dv(x, p, t)
end

function _monodromy_vectorfield_derivative(v::PiecewiseV, dvs::Tuple, x, p, t)
    dvs[v.n](x, p, t)
end

function _monodromy_vectorfield_derivative(v::PiecewiseImpactV, dvs::Tuple, x, p, t)
    dvs[v.n](x, p, t)
end

function _monodromy_augmented_dynamics!(dz, z, v, dvs, p, t, ::Val{N}) where {N}
    x = _to_svector(z[1:N], Val(N))
    X = reshape(view(z, N + 1:N + N * N), N, N)
    xdot = v(x, p, t)
    dX = _monodromy_vectorfield_derivative(v, dvs, x, p, t) * X
    dz[1:N] .= xdot
    dz[N+1:end] .= vec(dX)
    nothing
end

function _init_monodromy_region!(v::BilliardV, x, p, t0, region_detect=_region_detect)
    nothing
end

function _init_monodromy_region!(v::Union{PiecewiseV,PiecewiseImpactV}, x, p, t0,
    region_detect=_region_detect)
    v.n = region_detect(v.regions, x, p, t0)
    nothing
end

function _monodromy_saltation_at_event!(integrator, v::BilliardV, dvs, N, idx;
    event_t, event_at, repeat_nudge=1 // 100)
    t = integrator.t
    p = integrator.p
    u = integrator.u
    push!(event_t, t)
    push!(event_at, idx)
    xpre = _to_svector(u[1:N], Val(N))
    ypost = v.rules[idx](xpre, p, t)
    fpre = v(xpre, p, t)
    fpost = v(ypost, p, t)
    DG = _num_jacobian(x -> v.rules[idx](x, p, t), collect(xpre))
    gt = _num_time_derivative((u, pp, tt) -> v.rules[idx](u, pp, tt), xpre, p, t)
    gradh = _num_gradient_scalar(v.hypers[idx], xpre, p, t)
    ht = _num_time_derivative(v.hypers[idx], xpre, p, t)
    S = _saltation_matrix(fpre, fpost, DG, gt, gradh, ht)
    X = reshape(view(u, N + 1:N + N * N), N, N)
    u[1:N] .= ypost
    u[N+1:end] .= vec(S * X)
    nothing
end

function _monodromy_saltation_at_event!(integrator, v::PiecewiseV, dvs, N, idx;
    event_t, event_at, cross_time=0.01, region_detect=_region_detect,
    repeat_nudge=1 // 100)
    t = integrator.t
    p = integrator.p
    u = integrator.u
    push!(event_t, t)
    push!(event_at, idx)
    xpre = _to_svector(u[1:N], Val(N))
    t0 = t + cross_time
    z0 = integrator.sol(t0)
    x0 = _to_svector(z0[1:N], Val(N))
    npost = region_detect(v.regions, x0, p, t0)
    npre = v.n
    fpre = v(xpre, p, t)
    v.n = npost
    ypost = xpre
    fpost = v(ypost, p, t)
    DG = Matrix{eltype(xpre)}(I, N, N)
    gt = zeros(eltype(xpre), N)
    gradh = _num_gradient_scalar(v.hypers[idx], xpre, p, t)
    ht = _num_time_derivative(v.hypers[idx], xpre, p, t)
    S = _saltation_matrix(fpre, fpost, DG, gt, gradh, ht)
    X = reshape(view(u, N + 1:N + N * N), N, N)
    u[1:N] .= ypost
    u[N+1:end] .= vec(S * X)
    nothing
end

function _monodromy_saltation_at_event!(integrator, v::PiecewiseImpactV, dvs, N, idx;
    event_t, event_at, cross_time=0.01, region_detect=_region_detect,
    repeat_nudge=1 // 100)
    t = integrator.t
    p = integrator.p
    u = integrator.u
    push!(event_t, t)
    push!(event_at, idx)
    xpre = _to_svector(u[1:N], Val(N))
    npre = v.n
    if idx in v.idxs
        ypost = v.rules[idx](xpre, p, t)
        npost = npre
        fpre = v(xpre, p, t)
        v.n = npost
        fpost = v(ypost, p, t)
        DG = _num_jacobian(x -> v.rules[idx](x, p, t), collect(xpre))
        gt = _num_time_derivative((u, pp, tt) -> v.rules[idx](u, pp, tt), xpre, p, t)
    else
        t0 = t + cross_time
        z0 = integrator.sol(t0)
        x0 = _to_svector(z0[1:N], Val(N))
        npost = region_detect(v.regions, x0, p, t0)
        ypost = xpre
        fpre = v(xpre, p, t)
        v.n = npost
        fpost = v(ypost, p, t)
        DG = Matrix{eltype(xpre)}(I, N, N)
        gt = zeros(eltype(xpre), N)
    end
    gradh = _num_gradient_scalar(v.hypers[idx], xpre, p, t)
    ht = _num_time_derivative(v.hypers[idx], xpre, p, t)
    S = _saltation_matrix(fpre, fpost, DG, gt, gradh, ht)
    X = reshape(view(u, N + 1:N + N * N), N, N)
    u[1:N] .= ypost
    u[N+1:end] .= vec(S * X)
    nothing
end

function _monodromy_callback(v::BilliardV, dvs, N, event_t, event_at;
    repeat_nudge=1 // 100, extra...)
    nn = length(v.hypers)
    condition(out, u, t, integrator) = begin
        for i in eachindex(v.hypers)
            out[i] = v.hypers[i](u, integrator.p, t)
        end
    end
    affect!(integrator, idx) =
        _monodromy_saltation_at_event!(integrator, v, dvs, N, idx;
            event_t=event_t, event_at=event_at, repeat_nudge=repeat_nudge)
    VectorContinuousCallback(condition, affect!, nn; repeat_nudge=repeat_nudge)
end

function _monodromy_callback(v::Union{PiecewiseV,PiecewiseImpactV}, dvs, N,
    event_t, event_at;
    cross_time=0.01, region_detect=_region_detect, repeat_nudge=1 // 100, extra...)
    nn = length(v.hypers)
    condition(out, u, t, integrator) = begin
        for i in eachindex(v.hypers)
            out[i] = v.hypers[i](u, integrator.p, t)
        end
    end
    affect!(integrator, idx) = _monodromy_saltation_at_event!(
        integrator, v, dvs, N, idx;
        event_t=event_t, event_at=event_at, cross_time=cross_time,
        region_detect=region_detect, repeat_nudge=repeat_nudge)
    VectorContinuousCallback(condition, affect!, nn; repeat_nudge=repeat_nudge)
end

"""
    ns_monodromy(v, dvs, timespan, x, p; kwargs...)

Compute the derivative ``DF(x)`` of the time-`T`-map of a non-smooth ODE by
solving, together with the flow, the matrix variational equation and applying
the saltation (jump) matrix at every transversal event.

For a [`BilliardV`](@ref) vector field `dvs` is the derivative of its single
vector field, `dvs(x,p,t)`.  For [`PiecewiseV`](@ref) and
[`PiecewiseImpactV`](@ref), `dvs` must be a tuple containing the derivatives of
all vector fields in `fs`.

The return value is the monodromy matrix.  Event times and indices are recorded
inside the augmented integration used by [`findsaddle`](@ref).

The keyword arguments `alg`, `abstol`, `reltol`, `dtmax`, `cross_time`,
`region_detect` and `repeat_nudge` have the same meaning as in
[`setmap`](@ref)/[`ns_solver`](@ref).
"""
function ns_monodromy(v::BilliardV, dv, timespan::Tuple{T,T}, x::SVector{N,T}, p;
    alg=Tsit5(), abstol=1e-8, reltol=1e-8, repeat_nudge=1 // 100, extra...) where {N,T}
    data = _ns_monodromy_data(v, dv, timespan, x, p;
        alg=alg, abstol=abstol, reltol=reltol, repeat_nudge=repeat_nudge, extra...)
    data.D
end

function ns_monodromy(v::Union{PiecewiseV,PiecewiseImpactV}, dvs::Tuple,
    timespan::Tuple{T,T}, x::SVector{N,T}, p;
    alg=Tsit5(), abstol=1e-8, reltol=1e-8, cross_time=0.01,
    region_detect=_region_detect, repeat_nudge=1 // 100, extra...) where {N,T}
    data = _ns_monodromy_data(v, dvs, timespan, x, p;
        alg=alg, abstol=abstol, reltol=reltol, cross_time=cross_time,
        region_detect=region_detect, repeat_nudge=repeat_nudge, extra...)
    data.D
end

function _ns_monodromy_data(v, dvs, timespan::Tuple{T,T}, x::SVector{N,T}, p;
    alg=Tsit5(), abstol=1e-8, reltol=1e-8, cross_time=0.01,
    region_detect=_region_detect, repeat_nudge=1 // 100, extra...) where {N,T}
    _init_monodromy_region!(v, x, p, timespan[1], region_detect)
    event_t = T[]
    event_at = Int[]
    cb = _monodromy_callback(v, dvs, N, event_t, event_at;
        cross_time=cross_time, region_detect=region_detect, repeat_nudge=repeat_nudge)
    z0 = vcat(collect(x), vec(Matrix{T}(I, N, N)))
    prob = ODEProblem{true}((dz, z, p, t) -> _monodromy_augmented_dynamics!(
        dz, z, v, dvs, p, t, Val(N)), z0, timespan, p)
    sol = solve(prob, alg; callback=cb, abstol=abstol, reltol=reltol, extra...)
    z1 = sol.u[end]
    x1 = _to_svector(z1[1:N], Val(N))
    Dmat = reshape(z1[N+1:end], N, N)
    D = SMatrix{N,N,T}(Dmat)
    NSMonodromyData{N,T}(x1, D, event_t, event_at)
end

"""
    findsaddle(v::BilliardV, dv, timespan, x, p; kwargs...)

Locate a saddle-type periodic solution of the time-`T`-map of a non-smooth ODE
by Newton's method.  The Jacobian is obtained from the monodromy matrix, which
is computed with the variational equation and the saltation (jump) matrices at
the impact/switching events; see [`ns_monodromy`](@ref).

This is the recommended replacement for
`findsaddle(setup, x, p)` whenever the periodic orbit has impacts or switching
events, because the plain finite-difference version of the whole map is less
accurate and less efficient near such events.
"""
function findsaddle(v::BilliardV, dv, timespan, x::SVector{N,T}, p;
    n=100, abstol=1e-8, alg=Tsit5(), reltol=1e-8, repeat_nudge=1 // 100,
    extra...) where {N,T}
    ii = SMatrix{N,N,T}(I)
    function step(u)
        _ns_monodromy_data(v, dv, timespan, u, p;
            alg=alg, abstol=abstol, reltol=reltol,
            repeat_nudge=repeat_nudge, extra...)
    end
    nd = step(x)
    xn = x - inv(nd.D - ii) * (nd.x1 - x)
    itdata = [x, xn]
    i = 1
    while norm(itdata[2] - itdata[1]) > abstol && i <= n
        itdata[1] = itdata[2]
        nd = step(itdata[1])
        itdata[2] = itdata[1] - inv(nd.D - ii) * (nd.x1 - itdata[1])
        i += 1
    end
    nd = step(itdata[2])
    if norm(itdata[2] - itdata[1]) < abstol
        eigendata = eigen(nd.D)
        indics = findall(>(1), abs.(eigendata.values))
        directions = [SVector{N,T}(eigendata.vectors[:, j]) for j in indics]
        println("Fixed point found successfully!")
        println("Has contact with hypersurfaces: ", !isempty(nd.event_t))
        println("The fixed point: ")
        Saddle(itdata[2], directions, [eigendata.values[j] for j in indics])
    else
        println("Failed to find a fixed point after $n times iterations. The last point is:")
        Saddle(itdata[2], SVector{N,T}[], T[])
    end
end

function findsaddle(v::Union{PiecewiseV,PiecewiseImpactV}, dvs::Tuple,
    timespan, x::SVector{N,T}, p;
    n=100, abstol=1e-8, alg=Tsit5(), reltol=1e-8, cross_time=0.01,
    region_detect=_region_detect, repeat_nudge=1 // 100, extra...) where {N,T}
    ii = SMatrix{N,N,T}(I)
    function step(u)
        _ns_monodromy_data(v, dvs, timespan, u, p;
            alg=alg, abstol=abstol, reltol=reltol,
            cross_time=cross_time, region_detect=region_detect,
            repeat_nudge=repeat_nudge, extra...)
    end
    nd = step(x)
    xn = x - inv(nd.D - ii) * (nd.x1 - x)
    itdata = [x, xn]
    i = 1
    while norm(itdata[2] - itdata[1]) > abstol && i <= n
        itdata[1] = itdata[2]
        nd = step(itdata[1])
        itdata[2] = itdata[1] - inv(nd.D - ii) * (nd.x1 - itdata[1])
        i += 1
    end
    nd = step(itdata[2])
    if norm(itdata[2] - itdata[1]) < abstol
        eigendata = eigen(nd.D)
        indics = findall(>(1), abs.(eigendata.values))
        directions = [SVector{N,T}(eigendata.vectors[:, j]) for j in indics]
        println("Fixed point found successfully!")
        println("Has contact with hypersurfaces: ", !isempty(nd.event_t))
        println("The fixed point: ")
        Saddle(itdata[2], directions, [eigendata.values[j] for j in indics])
    else
        println("Failed to find a fixed point after $n times iterations. The last point is:")
        Saddle(itdata[2], SVector{N,T}[], T[])
    end
end
