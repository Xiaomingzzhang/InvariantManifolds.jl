# Test piecewise vector field
@test_nowarn begin
    f1(x, p, t) = SA[x[2], p[1]*x[1]+p[4]*sin(2pi * t)]
    f2(x, p, t) = SA[x[2], -p[2]*x[1]+p[4]*sin(2pi * t)]
    f3(x, p, t) = SA[x[2], -p[3]*x[1]+p[4]*sin(2pi * t)]
    hyper1(x, p, t) = x[1] - p[5]
    hyper2(x, p, t) = x[1] + p[5]
    dom1(x, p, t) = -p[5] < x[1] < p[5]
    dom2(x, p, t) = x[1] > p[5]
    dom3(x, p, t) = x[1] < -p[5]
    vectorfield = PiecewiseV((f1, f2, f3), (dom1, dom2, dom3), (hyper1, hyper2))
    setup = setmap(vectorfield, (0.0, 1.0), Tsit5(), abstol=1e-8)
    function df(x, p, t)
        SA[0 1; p[1] 0]
    end
    para = [2, 5, 5, 0.7, 1.2]
    prob = NSOneDManifoldProblem(setup, para, ϵ=1e-3, amax=0.5, dsmin=1e-6, d=0.001)
    saddle = findsaddle(f1, df, (0.0,1.0), SA[0.0,0.0], para)
    seg = gen_segment(saddle)
    result = growmanifold(prob, seg, 5)
end
# Test vector field with impact
@test_nowarn begin
    f(x, p, t) = SA[x[2], sin(x[1])+p[1]*x[2]+p[2]*sin(2*pi*t)]
    hyper1(x, p, t) = x[1] - p[3]
    hyper2(x, p, t) = x[1] + p[3]
    rule(x, p, t) = SA[x[1], -p[4]*x[2]]
    vectorfield = BilliardV(f, (hyper1, hyper2), (rule,rule))
    setup = setmap(vectorfield, (0.0, 1.0), Tsit5(), abstol=1e-12)
    function df(x, p, t)
        SA[0 1; cos(x[1]) p[1]]
    end
    para = [0.1, 0.2, pi/4, 0.9]
    prob = NSOneDManifoldProblem(setup, para, ϵ=1e-5, d=0.0001, dsmin=1e-5)
    saddle = findsaddle(f, df, (0.0,1.0), SA[0.0,0.0], para)
    seg = gen_segment(saddle)
    result = growmanifold(prob, seg, 6)
end
# 