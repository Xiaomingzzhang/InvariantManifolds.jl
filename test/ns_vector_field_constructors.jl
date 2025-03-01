# Define all test functions locally
let
    f1(x, p, t) = SA[x[2], p[1]*x[1]]
    f2(x, p, t) = SA[x[2], -p[2]*x[1]]
    f3(x, p, t) = SA[sin(x[1]), cos(x[2])]
    f4(x, p, t) = SA[-x[2], x[1]^3]
    
    hyper(x, p, t) = x[1] - p[3]
    hyper2(x, p, t) = x[2]
    hyper3(x, p, t) = x[1]^2 + x[2]^2 - p[1]^2
    
    dom1(x, p, t) = p[3] < x[1]
    dom2(x, p, t) = x[1] > p[3]
    dom3(x, p, t) = x[1]^2 + x[2]^2 < p[1]^2
    dom4(x, p, t) = x[1]*x[2] > 0
    
    rule1(x, p, t) = SA[-x[1], x[2]]  # reflect x coordinate
    rule2(x, p, t) = SA[x[1], -x[2]]  # reflect y coordinate
    rule3(x, p, t) = -x  # reflect both coordinates

    # Basic PiecewiseV tests
    initialconditions = (SA[1.0, 5.0], SA[1.0, 1.0, 2.0], 0.0)
    pv = PiecewiseV((f1, f2), (dom1, dom2), (hyper,))
    
    @test pv.n == 0
    pv.n = 1
    @test pv(initialconditions...) == f1(initialconditions...)
    
    pv_1 = PiecewiseV((f1, f2), (dom1, dom2), (hyper,), 1)
    pv_2 = PiecewiseV((f1, f2), (dom1, dom2), (hyper,), 2)
    @test pv_1(initialconditions...) == f1(initialconditions...)
    @test pv_2(initialconditions...) == f2(initialconditions...)

    # Basic BilliardV tests
    bv = BilliardV(f1, (dom1, dom2), (hyper, hyper))
    @test bv(initialconditions...) == f1(initialconditions...)

    # Complex PiecewiseV tests
    pv_multi = PiecewiseV((f1, f2, f3), (dom1, dom2, dom3), (hyper, hyper2))
    @test pv_multi.n == 0
    @test length(pv_multi.fs) == 3
    @test length(pv_multi.hypers) == 2

    x1 = SA[2.0, 1.0]
    p1 = SA[1.0, 1.0, 0.0]
    t1 = 0.0
    pv_multi.n = 1
    @test pv_multi(x1, p1, t1) == f1(x1, p1, t1)
    
    pv_multi.n = 2
    @test pv_multi(x1, p1, t1) == f2(x1, p1, t1)

    # Complex BilliardV tests
    bv_multi = BilliardV(f1, (hyper, hyper2, hyper3), (rule1, rule2, rule3))
    @test length(bv_multi.hypers) == 3
    @test length(bv_multi.rules) == 3

    x2 = SA[1.0, -1.0]
    p2 = SA[2.0, 1.0, 0.0]
    t2 = 0.0
    @test bv_multi.rules[1](x2, p2, t2) == SA[-1.0, -1.0]
    @test bv_multi.rules[2](x2, p2, t2) == SA[1.0, 1.0]
    @test bv_multi.rules[3](x2, p2, t2) == SA[-1.0, 1.0]
    @test bv_multi(x2, p2, t2) == f1(x2, p2, t2)

    # PiecewiseImpactV tests
    # Test basic construction and field access
    piv = PiecewiseImpactV((f1, f2), (dom1, dom2), (hyper, hyper2), (rule1, rule2), [1])
    @test piv.n == 0
    @test length(piv.fs) == 2
    @test length(piv.regions) == 2
    @test length(piv.hypers) == 2
    @test length(piv.rules) == 2
    @test piv.idxs == [1]

    # Test vector field evaluation in different regions
    piv.n = 1
    @test piv(x1, p1, t1) == f1(x1, p1, t1)
    piv.n = 2
    @test piv(x1, p1, t1) == f2(x1, p1, t1)

    # Test impact rules
    @test piv.rules[1](x2, p2, t2) == SA[-1.0, -1.0]
    @test piv.rules[2](x2, p2, t2) == SA[1.0, 1.0]

    # Test complex configuration
    piv_complex = PiecewiseImpactV(
        (f1, f2, f3),  # multiple vector fields
        (dom1, dom2, dom3),  # multiple regions
        (hyper, hyper2, hyper3),  # multiple hypersurfaces
        (rule1, rule2, rule3),  # multiple rules
        [1, 3]  # multiple impact indices
    )
    @test length(piv_complex.fs) == 3
    @test length(piv_complex.regions) == 3
    @test length(piv_complex.hypers) == 3
    @test length(piv_complex.rules) == 3
    @test piv_complex.idxs == [1, 3]

    # Test error cases
    @test_throws ArgumentError PiecewiseV((f1,), (dom1, dom2), (hyper,))
    @test_throws ArgumentError BilliardV(f1, (hyper,), (rule1, rule2))
    @test_throws ArgumentError PiecewiseImpactV(
        (f1,),  # one field
        (dom1, dom2),  # two regions (mismatch)
        (hyper,),
        (rule1,),
        [1]
    )
end
