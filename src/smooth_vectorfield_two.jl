"""
    AnnulusBoundaries

A struct contains data when generating the two dimensional manifold.

# Fields
- `inner` an `IterationCurve` represents inner boundary;
- `outer` an `IterationCurve` represents outer boundary;
"""
struct AnnulusBoundaries{T}
    inner::T
    outer::T
end

function inintialise_mesh(saddle, v1, v2, n, r; interp = LinearInterpolation)
    A = hcat(v1, v2)
    Q = SMatrix(qr(A).Q)
    newv1 = Q[:, 1]
    newv2 = Q[:, 2]
    rag = range(0, 1, length=n + 1)
    circle1 = [saddle for t in rag]
    circle2 = [saddle + r * cospi(2 * t / 1) * newv1 + r * sinpi(2 * t / 1) * newv2 for t in rag]
    ic1 = interp(circle1, rag)
    ic2 = interp(circle2, rag)
    [AnnulusBoundaries(ic1, ic2)]
end

function grow_surface!(f, p, annulus::Vector{AnnulusBoundaries{S}}, δ; interp = LinearInterpolation) where {S}
    oldu0 = copy(annulus[end].outer.u)
    olds0 = copy(annulus[end].outer.t)
    oldcurve = annulus[end].outer
    k = length(oldu0)
    newu0 = similar(oldu0)
    T = typeof(oldu0[1][1])
    newss = T[0]
    for i in eachindex(newu0)
        newu0[i] = f(oldu0[i], p)
    end
    j = 1
    while j + 1 <= k
        dist = norm(newu0[j] - newu0[j+1])
        if dist > δ
            m = ceil(Int, dist / δ)
            if m == 1
                m = 2
            end
            s1 = olds0[j+1]
            s0 = olds0[j]
            plengh = (s1 - s0) / m
            paras = [s0 + plengh * i for i in 1:m-1]
            preaddps = similar(oldu0, m - 1)
            addps = similar(oldu0, m - 1)
            for i in 1:m-1
                news = paras[i]
                preaddps[i] = oldcurve(news)
                addps[i] = f(oldcurve(news), p)
            end
            insert!(oldu0, j + 1, preaddps)
            insert!(newu0, j + 1, addps)
            insert!(olds0, j + 1, paras)
            k = k + m - 1
        else
            j = j + 1
            dd = newss[end]
            append!(newss, [dd + dist])
        end
    end
    ic1 = interp(newu0, newss)
    ic2 = interp(oldu0, olds0)
    newannulus = AnnulusBoundaries(ic2, ic1)
    append!(annulus, [newannulus])
    return nothing
end


"""
    generate_surface(f, p, saddle, v1, v2, N, r, δ)

Function to generate the two dimension manifold of a vector field.

# Parameters
- `f` the time-T-map a the vector field;
- `p` parameter of the system;
- `saddle` saddle fixed point of the ODE system;
- `v1` an instable direction of the `saddle`;
- `v2` another instable direction of the `saddle`;
- `N` iteration times;
- `r` the radius of origin disk to extend;
- `δ` the max distance between points when iterating.

# Keyword arguments
- `n=150` the number of points in the boundary of the origin disk to extend.
"""
function generate_surface(f, p, saddle, v1, v2, N, r, δ; n=150, interp = LinearInterpolation)
    @show myannulus = inintialise_mesh(saddle, v1, v2, n, r, interp=interp)
    for i in 1:N
        grow_surface!(f, p, myannulus,  δ; interp=interp)
    end
    myannulus
end