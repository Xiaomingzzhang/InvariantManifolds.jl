struct Jacobi{S,F}
    sol::S
    dv::F
end

function (v::Jacobi)(x, p, t)
    v.dv(v.sol(t), p, t) * x
end

"""
    ODESolver{F1,F2,T}

A wrapper struct for solving ordinary differential equations (ODEs).

# Fields
- `f`: Vector field function of the ODE system in the form `f(x,p,t)`
- `timespan`: Time interval for solving the ODE, of type `Tuple{T,T}`
- `alg`: The numerical algorithm used for solving the ODE
- `abstol`: Absolute tolerance for the numerical solver
"""
struct ODESolver{F1,F2,T}
    f::F1
    timespan::Tuple{T,T}
    alg::F2
    abstol::T
end

function (v::ODESolver)(x, p)
    prob = ODEProblem{false}(v.f, x, v.timespan, p)
    solve(prob, v.alg, abstol=v.abstol)
end

"""
    findsaddle(v, dv, timespan, x, p)

`findsaddle` is a function to find the saddle of the time-T-map of smooth ODE systems, by using the Newton's method.
# Parameters
- `v` the vector field, which should be the form `f(x,p,t)` and return a `SVector`;
- `dv` the Jacobi matrix function of `v`, which should should be the form `dv(x,p,t)` and return a `SMatrix`;
- `timespan` the time span of the time-T-map;
- `x` the initial point to iterate.
# Keyword arguments
- `n` maximum iterate times, default to be 100;
- `abstol` absolute tolerance for the fixed point, default to be `1e-8`;
- `alg` the algorithm used to solve the ODE, default to be `Vern9()`.
We also provide the finite difference method to find the saddle of the time-T-map of nonsmooth ODE systems:
# Parameters
- `setup` the [`NSSetUp`](@ref) of the nonsmooth ODE system;
- `x` the initial point to iterate, of type `SVector{N,T}`.
- `p` the parameters of the nonsmooth ODE system.
# Keyword arguments
- `n` maximum iterate times, default to be 100;
- `abstol` absolute tolerance for the fixed point, default to be `1e-8`;
"""
function findsaddle(v, dv, timespan, x::SVector{N,T}, p; n=100, abstol=1e-8, alg=Vern9()) where {N,T}
    timemap = ODESolver(v, timespan, alg, abstol)
    ii = SMatrix{N,N,T}(I)
    function jac(u, p)
        sol = timemap(u, p)
        df = Jacobi(sol, dv)
        nprob = ODEProblem{false}(df, ii, timespan, p)
        solve(nprob, alg, abstol=abstol)[end]
    end
    xn = x - inv(jac(x, p) - ii) * (timemap(x, p)[end] - x)
    data = [x, xn]
    i = 1
    while norm(data[2] - data[1]) > abstol && i <= n
        data[1] = data[2]
        data[2] = data[1] - inv(jac(data[1], p) - ii) * (timemap(data[1], p)[end] - data[1])
        i += 1
    end
    if norm(data[2] - data[1]) < abstol
        eigendata = eigen(jac(data[2], p))
        indics = findall(>(1), abs.(eigendata.values))
        directions = [eigendata.vectors[:, i] for i in indics]
        println("Fixed point found successfully!")
        println("The fixed point: ")
        Saddle(data[2], directions, [eigendata.values[i] for i in indics])
    else
        println("Failed to find a fixed point after $n times iterations. The last point is:")
        Saddle(data[2], typeof(data[2])[], T[])
    end
end


function findsaddle(setup::NSSetUp, x::SVector{N,T}, p; n=100, abstol=1e-8) where{N,T}
    timemap = setup.timetmap
    ii = SMatrix{N,N,T}(I)
    function jac(u, p)
        ff(x) = timemap(NSState(x), p).state
        FiniteDiff.finite_difference_jacobian(ff, u)
    end
    xn = x - inv(jac(x, p) - ii) * (timemap(NSState(x), p).state - x)
    data = [x, xn]
    i = 1
    while norm(data[2] - data[1]) > abstol && i <= n
        data[1] = data[2]
        data[2] = data[1] - inv(jac(data[1], p) - ii) * (timemap(NSState(data[1]), p).state - data[1])
    end
    if norm(data[2] - data[1]) < abstol
        eigendata = eigen(jac(data[2], p))
        indics = findall(>(1), abs.(eigendata.values))
        directions = [eigendata.vectors[:, i] for i in indics]
        println("Fixed point found successfully!")
        println("Has contact with hypersurfaces: ", iscontact(setup, data[2], p))
        println("The fixed point: ")
        Saddle(data[2], directions, [eigendata.values[i] for i in indics])
    else
        println("Failed to find a fixed point after $n times iterations. The last point is:")
        Saddle(data[2], typeof(data[2])[], T[])
    end
end