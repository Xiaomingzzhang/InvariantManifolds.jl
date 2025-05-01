@testset "OneDManifoldProblem" begin
    # Simple nonlinear map
    f(x, p) = SA[x[2], -x[1]+p[1]*sin(x[2])]

    # Test default constructor
    prob1 = OneDManifoldProblem(f)
    @test prob1.f === f
    @test isempty(prob1.para)
    @test prob1.amax ≈ 0.5
    @test prob1.d ≈ 1e-3
    @test prob1.dsmin ≈ 1e-6

    # Test constructor with parameters
    para = [0.1]
    prob2 = OneDManifoldProblem(f, para)
    @test prob2.para == para

    # Test constructor with custom settings
    prob3 = OneDManifoldProblem(f, para, amax=0.3, d=0.002, dsmin=1e-6)
    @test prob3.amax ≈ 0.3
    @test prob3.d ≈ 0.002
    @test prob3.dsmin ≈ 1e-6
end

@testset "NSOneDManifoldProblem" begin
    # Create a simple NSSetUp for testing
    f1(x, p, t) = SA[x[2], -x[1]]
    f2(x, p, t) = SA[x[2], x[1]]
    hyper(x, p, t) = x[1]
    region1(x, p, t) = x[1]>0
    region2(x, p, t) = x[1]<0
    vectorfield = PiecewiseV((f1,f2),(region1,region2),(hyper,))
    setup = NSSetUp(vectorfield, (0.0, 1.0), x -> x)

    # Test default constructor
    prob1 = NSOneDManifoldProblem(setup)
    @test prob1.f === setup
    @test isempty(prob1.para)
    @test prob1.amax ≈ 0.5
    @test prob1.d ≈ 1e-3
    @test prob1.ϵ ≈ 1e-5
    @test prob1.dsmin ≈ 1e-6

    # Test constructor with parameters
    para = [1.0]
    prob2 = NSOneDManifoldProblem(setup, para)
    @test prob2.para == para

    # Test constructor with custom settings
    prob3 = NSOneDManifoldProblem(setup, para, amax=0.3, d=0.002, ϵ=1e-6, dsmin=1e-6)
    @test prob3.amax ≈ 0.3
    @test prob3.d ≈ 0.002
    @test prob3.ϵ ≈ 1e-6
    @test prob3.dsmin ≈ 1e-6
end

@testset "VTwoDManifoldProblem" begin
    # Simple vector field
    f(x, p) = SA[x[2], -x[1]+p[1]*sin(x[2])]

    # Test default constructor
    prob1 = VTwoDManifoldProblem(f)
    @test prob1.f === f
    @test isempty(prob1.para)
    @test prob1.amax ≈ 0.5
    @test prob1.d ≈ 1e-3
    @test prob1.dsmin ≈ 1e-6

    # Test constructor with parameters
    para = [0.1]
    prob2 = VTwoDManifoldProblem(f, para)
    @test prob2.para == para

    # Test constructor with custom settings
    prob3 = VTwoDManifoldProblem(f, para, amax=0.3, d=0.002, dsmin=1e-4)
    @test prob3.amax ≈ 0.3
    @test prob3.d ≈ 0.002
    @test prob3.dsmin ≈ 1e-4
end

@testset "TwoDManifoldProblem" begin
    # Simple nonlinear map
    f(x, p) = SA[x[2], -x[1]+p[1]*sin(x[2])]

    # Test default constructor
    prob1 = TwoDManifoldProblem(f)
    @test prob1.f === f
    @test isempty(prob1.para)
    @test prob1.amax ≈ 0.5
    @test prob1.d ≈ 1e-3
    @test prob1.dcircle ≈ 0.01
    @test prob1.dsmin ≈ 1e-6

    # Test constructor with parameters
    para = [0.1]
    prob2 = TwoDManifoldProblem(f, para)
    @test prob2.para == para

    # Test constructor with custom settings
    prob3 = TwoDManifoldProblem(f, para, amax=0.3, d=0.002, dcircle=0.003, dsmin=1e-4)
    @test prob3.amax ≈ 0.3
    @test prob3.d ≈ 0.002
    @test prob3.dcircle ≈ 0.003
    @test prob3.dsmin ≈ 1e-4
end

@testset "NSVTwoDManifoldProblem" begin
    # Create a simple NSSetUp for testing
    f1(x, p, t) = SA[x[2], -x[1]]
    f2(x, p, t) = SA[x[2], x[1]]
    hyper(x, p, t) = x[1]
    region1(x, p, t) = x[1]>0
    region2(x, p, t) = x[1]<0
    vectorfield = PiecewiseV((f1,f2),(region1,region2),(hyper,))
    setup = NSSetUp(vectorfield, (0.0, 1.0), x -> x)

    # Test default constructor
    prob1 = NSVTwoDManifoldProblem(setup)
    @test prob1.f === setup
    @test isempty(prob1.para)
    @test prob1.amax ≈ 0.5
    @test prob1.d ≈ 1e-3
    @test prob1.ϵ ≈ 1e-5
    @test prob1.dsmin ≈ 1e-6

    # Test constructor with parameters
    para = [1.0]
    prob2 = NSVTwoDManifoldProblem(setup, para)
    @test prob2.para == para

    # Test constructor with custom settings
    prob3 = NSVTwoDManifoldProblem(setup, para, amax=0.3, d=0.002, ϵ=1e-6, dsmin=1e-5)
    @test prob3.amax ≈ 0.3
    @test prob3.d ≈ 0.002
    @test prob3.ϵ ≈ 1e-6
    @test prob3.dsmin ≈ 1e-5
end