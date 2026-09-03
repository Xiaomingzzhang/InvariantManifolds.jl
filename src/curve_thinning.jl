"""
    InvariantManifolds._curve_angle(p0, p1, p2)

Return the local turning angle at `p1` formed by the three consecutive points
`p0`, `p1` and `p2`.  The definition is the same one as used in
[`addpoints!`](@ref) and [`ns_addpoints!`](@ref):

```math
\\alpha = 2\\arcsin\\frac{\\|\\bar p_0-p_0\\|}{2\\|p_1-p_0\\|},
\\qquad \\bar p_0=p_1+\\frac{p_1-p_2}{\\|p_1-p_2\\|}\\|p_1-p_0\\|.
```

Points that are equal give the angle `0`.
"""
@inline function _curve_angle(p0, p1, p2)
    d1 = norm(p0 - p1)
    d2 = norm(p1 - p2)
    d1 == 0 || d2 == 0 && return zero(float(d1 * d2))
    q = p1 + (p1 - p2) * (d1 / d2)
    x = norm(q - p0) / (2d1)
    2asin(clamp(x, -one(x), one(x)))
end

"""
    InvariantManifolds._angle_distance_tolerance(d, αmax)

Scale below which a large local turning angle over an extremely short segment is
treated as numerically/geometrically not resolvable, following the angle-times-length
control used in the adaptive methods of Goodman and Wróbel (2011).  A segment with
`α*d ≤ tolerance` has no meaningful contribution to the polygonal curve error even
if the raw angle exceeds `αmax`, so it is not recorded as a [`FlawPoint`](@ref).
"""
function _angle_distance_tolerance(d, αmax)
    max(d * αmax / 1000, 100eps(float(αmax * d)))
end

function _thin_curve_pass(points::Vector{P}, amax, d) where {P}
    n = length(points)
    keep = trues(n)
    tol = max(d * amax / 8, 100eps(float(amax * d)))
    @inbounds for i in 2:n-1
        if keep[i-1] && keep[i] && keep[i+1]
            a = points[i-1]
            b = points[i]
            c = points[i+1]
            ab = b - a
            ac = c - a
            L2 = dot(ac, ac)
            # sagitta-like flatness check: the removed point may not lie far from
            # the chord connecting its neighbours
            if L2 > 0
                t = clamp(dot(ab, ac) / L2, 0, 1)
                dist = norm(ab - t * ac)
                dist <= tol || continue
            end
            # points that resolve a genuinely sharp turn are retained
            _curve_angle(points[i-1], points[i], points[i+1]) <= amax || continue
            if i >= 3 && keep[i-2]
                _curve_angle(points[i-2], points[i-1], points[i+1]) <= amax || continue
            end
            if i + 2 <= n && keep[i+2]
                _curve_angle(points[i-1], points[i+1], points[i+2]) <= amax || continue
            end
            keep[i] = false
        end
    end
    points[keep]
end

"""
    thin_curve(curve; amax, d, interp=QuadraticInterpolation)

Remove unnecessary points of an already computed invariant curve after one
continuation step.  This follows the classical idea in Parker and Chua (1989)
that one should also delete unnecessary points while adding points by the
distance and turning-angle criteria.

A point is kept whenever deleting it would make any local turning angle larger
than `amax`, or when the point is farther from the chord connecting its
neighbours than the flatness tolerance `d*amax/8`.  Consequently, only redundant
points on nearly straight, over-resolved parts are removed, while points that
resolve genuine sharp turns are retained.

The parameterization of the returned curve is recomputed by arc length with
[`paramise`](@ref); no map evaluation is needed in this step.
"""
function thin_curve(curve; amax, d, interp=QuadraticInterpolation)
    points = curve.u
    n = length(points)
    n <= 4 && return curve
    newpoints = _thin_curve_pass(points, amax, d)
    if length(newpoints) == n
        curve
    else
        paramise(newpoints, interp=interp)
    end
end

thin_curve(curve::S, prob; interp=QuadraticInterpolation) where {S} =
    thin_curve(curve; amax=prob.amax, d=prob.d, interp=interp)
