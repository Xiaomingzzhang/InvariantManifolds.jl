# constructor
@test NSState(SA[1.0,2.0]) == NSState(SA[1.0,2.0],Int[])
@test_throws MethodError NSState(SA[1.0,2.0],[1.2,2.0])


# arithmetic
@test NSState(SA[1.0,2.0])+NSState(SA[1.0,2.0]) == SA[2.0,4.0]
@test NSState(SA[1.0,2.0])+NSState(SA[1.0,2.0]) == SA[2.0,4.0]
@test pi*NSState(SA[1.0,2.0]) == pi*SA[1.0,2.0]
@test NSState(SA[1.0,2.0])*pi == SA[1.0,2.0]*pi
@test NSState(SA[1.0,2.0])/pi == SA[1.0,2.0]/pi

# behavior like Vectors

@test length(NSState(SA[1.0,2.0])) == 2
@test size(NSState(SA[1.0,2.0])) == (2,)
@test NSState(SA[1.0,2.0])[1] == 1.0

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

