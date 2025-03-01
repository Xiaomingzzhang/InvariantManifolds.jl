# constructor
@test NSState(SA[1.0,2.0]) == NSState(SA[1.0,2.0],Int[])
@test_throws MethodError NSState(SA[1.0,2.0],[1.2,2.0])

# more constructor tests
@test NSState(SA[1.0,2.0], Int[1,2]) == NSState(SA[1.0,2.0], [1,2])
@test NSState{2,Float64}(NSState(SA[1.0,2.0])) == NSState(SA[1.0,2.0])
@test_throws DimensionMismatch NSState(SA[1.0,2.0], Int[1,2]) + NSState(SA[1.0], Int[])

# arithmetic
@test NSState(SA[1.0,2.0])+NSState(SA[1.0,2.0]) == SA[2.0,4.0]
@test NSState(SA[1.0,2.0])+NSState(SA[1.0,2.0]) == SA[2.0,4.0]
@test pi*NSState(SA[1.0,2.0]) == pi*SA[1.0,2.0]
@test NSState(SA[1.0,2.0])*pi == SA[1.0,2.0]*pi
@test NSState(SA[1.0,2.0])/pi == SA[1.0,2.0]/pi

# additional arithmetic
@test NSState(SA[1.0,2.0]) - NSState(SA[0.5,1.0]) == SA[0.5,1.0]
@test -NSState(SA[1.0,2.0]) == SA[-1.0,-2.0]
@test 2 * NSState(SA[1.0,2.0]) == NSState(SA[2.0,4.0])

# broadcasting
@test NSState(SA[1.0,2.0]) .+ SA[1.0,2.0] == SA[2.0,4.0]
@test NSState(SA[1.0,2.0]) .* SA[2.0,3.0] == SA[2.0,6.0]

# behavior like Vectors

@test length(NSState(SA[1.0,2.0])) == 2
@test size(NSState(SA[1.0,2.0])) == (2,)
@test NSState(SA[1.0,2.0])[1] == 1.0

# additional vector behavior
@test eltype(NSState(SA[1.0,2.0])) == Float64
@test collect(NSState(SA[1.0,2.0])) == [1.0,2.0]
@test firstindex(NSState(SA[1.0,2.0])) == 1
@test lastindex(NSState(SA[1.0,2.0])) == 2

# type stability
let
    x = NSState(SA[1.0,2.0])
    @test typeof(x + x) == typeof(x.state)
    @test typeof(2 * x) == typeof(x.state)
end

# work with DataInterpolations
t = 0:0.01:1
nsstates = [NSState(SA[sin(2pi*t[i]),cos(2pi*t[i])]) for i in eachindex(t)]
linear_interp = LinearInterpolation(nsstates,t)
quad_interp= QuadraticInterpolation(nsstates,t)
qubic_interp= CubicSpline(nsstates,t)
@test linear_interp(0.2) ≈ SA[sin(2pi*0.2),cos(2pi*0.2)]
@test linear_interp(0.53) ≈ SA[sin(2pi*0.53),cos(2pi*0.53)]
@test quad_interp(0.2) ≈ SA[sin(2pi*0.2),cos(2pi*0.2)]
@test quad_interp(0.53) ≈ SA[sin(2pi*0.53),cos(2pi*0.53)]
@test qubic_interp(0.2) ≈ SA[sin(2pi*0.2),cos(2pi*0.2)]
@test qubic_interp(0.53) ≈ SA[sin(2pi*0.53),cos(2pi*0.53)]
