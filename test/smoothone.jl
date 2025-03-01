@test_nowarn begin
    function henonmap(x, p)
        y1 = 1 - p[1] * x[1]^2 + x[2]
        y2 = p[2] * x[1]
        SA[y1, y2]
    end
    henonmap2(x, p) = henonmap(henonmap(x, p), p)
    henon_prob = OneDManifoldProblem(henonmap2, [1.4, 0.3])
    henon_seg = gen_segment(SA[-1.131354477, -0.339406343], SA[0.9957919906814664, 0.09164230079304274])
    henon_manifold = growmanifold(henon_prob, henon_seg, 2)
    henon_manifold = growmanifold(henon_prob, henon_seg, 2, interp=QuadraticInterpolation)
    henon_manifold = growmanifold(henon_prob, henon_seg, 2, interp=CubicSpline)
end