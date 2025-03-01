@testset "gen_segment" begin
    # Test with Saddle object
    saddle_point = SA[0.0, 0.0]
    unstable_dir = SA[1.0, 0.0]
    saddle = Saddle(saddle_point, [unstable_dir], [1.0])
    
    # Test default parameters
    segment = gen_segment(saddle)
    @test length(segment) == 150  # default n=150
    @test segment[1] ≈ saddle_point
    @test norm(segment[end] - segment[1]) ≈ 0.01  # default d=0.01
    @test all(x -> length(x) == 2, segment)  # check dimension
    
    # Test with custom parameters
    n, d = 100, 0.02
    segment = gen_segment(saddle, n=n, d=d)
    @test length(segment) == n
    @test segment[1] ≈ saddle_point
    @test norm(segment[end] - segment[1]) ≈ d
    
    # Test that points are evenly distributed
    diffs = [norm(segment[i+1] - segment[i]) for i in 1:length(segment)-1]
    @test all(x -> isapprox(x, diffs[1], rtol=1e-10), diffs)
end

@testset "gen_disk" begin
    # Test with Saddle object in 2D
    saddle_point = SA[0.0, 0.0]
    unstable_dir1 = SA[1.0, 0.0]
    unstable_dir2 = SA[0.0, 1.0]
    eigen_values = [2.0, 1.5]  # λ₁ > λ₂
    saddle = Saddle(saddle_point, [unstable_dir1, unstable_dir2], eigen_values)
    
    # Test default parameters
    disk = gen_disk(saddle)
    @test length(disk) == 2  # default circles=2
    @test disk[1][1] ≈ saddle_point  # center point
    
    # Test with custom parameters
    n, d, r, circles = 100, 0.0001, 0.02, 3
    disk = gen_disk(saddle, n=n, d=d, r=r, circles=circles)
    @test length(disk) == circles
    @test disk[1][1] ≈ saddle_point
    
    # Test with different times parameter
    disk_evolved = gen_disk(saddle, times=2)
    @test length(disk_evolved) == 2  # default circles=2
    
end