using InvariantManifolds, LinearAlgebra, StaticArrays, OrdinaryDiffEq, FiniteDiff

@testset "thin_curve" begin
    # Dense, almost collinear points should be coarsened, keeping the ends.
    pts = [SA[0.02 * i - 1.0, 1e-5 * sin(4 * i)] for i in 0:100]
    curve = InvariantManifolds.paramise(pts)
    thinned = thin_curve(curve; amax=0.5, d=0.05)
    @test length(thinned.u) < length(curve.u)
    @test first(thinned.u) == first(curve.u)
    @test last(thinned.u) == last(curve.u)

    # A genuine sharp corner must not be destroyed by thinning.
    sharp = SVector{2,Float64}[SA[0.0, 0.0], SA[0.01, 0.0],
        SA[0.02, 0.0], SA[0.0201, 0.001], SA[0.0203, 0.003],
        SA[0.0206, 0.006], SA[0.021, 0.01], SA[0.03, 0.03]]
    sc = InvariantManifolds.paramise(sharp)
    sth = thin_curve(sc; amax=0.2, d=0.01)
    @test length(sth.u) >= 4
end

@testset "monodromy with jump matrix" begin
    f(x, p, t) = SA[x[2], sin(x[1]) + p[1] * x[2] + p[2] * sin(p[3] * t)]
    hyper1(x, p, t) = x[1] - p[4]
    hyper2(x, p, t) = x[1] + p[4]
    rule(x, p, t) = SA[x[1], -p[5] * x[2]]
    vf = BilliardV(f, (hyper1, hyper2), (rule, rule))
    para = [0.1, 0.1, 1.64, pi / 4, 0.9]
    timespan = (0.0, 2pi / para[3])
    setup = setmap(vf, timespan, Tsit5(), abstol=1e-8, reltol=1e-8)
    df(x, p, t) = SA[0 1; cos(x[1]) p[1]]
    xguess = SA[0.33105711548504757, -0.2795863847274121]

    DF_fd = FiniteDiff.finite_difference_jacobian(
        u -> setup.timetmap(NSState(u), para).state, xguess)
    data = InvariantManifolds._ns_monodromy_data(
        vf, df, timespan, xguess, para;
        alg=Tsit5(), abstol=1e-8, reltol=1e-8, dtmax=0.01)
    @test !isempty(data.event_t)
    @test data.D ≈ Matrix(DF_fd) rtol=1e-4 atol=1e-4

    saddle = findsaddle(vf, df, timespan, xguess, para;
        abstol=1e-8, reltol=1e-8, dtmax=0.01)
    @test norm(saddle.saddle - xguess) < 1e-5
    @test iscontact(setup, saddle, para)
end

@testset "piecewise monodromy with jump matrix" begin
    f1(x, p, t) = SA[x[2], -x[1]]
    f2(x, p, t) = SA[x[2], -2x[1]]
    hyper(x, p, t) = x[1]
    dom1(x, p, t) = x[1] < 0
    dom2(x, p, t) = x[1] > 0
    pw = PiecewiseV((f1, f2), (dom1, dom2), (hyper,))
    timespan = (0.0, 1.0)
    setup = setmap(pw, timespan, Tsit5(), abstol=1e-10, reltol=1e-10,
        cross_time=1e-5)
    dvs = (
        (x, p, t) -> SA[0 1; -1 0],
        (x, p, t) -> SA[0 1; -2 0],
    )
    x0 = SA[-0.5, 1.5]
    DF = ns_monodromy(pw, dvs, timespan, x0, Float64[];
        abstol=1e-10, reltol=1e-10, cross_time=1e-5)
    DF_fd = FiniteDiff.finite_difference_jacobian(
        u -> setup.timetmap(NSState(u), Float64[]).state, x0)
    @test Matrix(DF) ≈ Matrix(DF_fd) rtol=1e-5 atol=1e-5
end
