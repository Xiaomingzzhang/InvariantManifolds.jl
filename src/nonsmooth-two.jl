# The following algorithm only works with non-smooth vector fields's 
# time-T-map, and the vector fields must have nearly uniformly extension. The latter condition
# means that you have to change your vector field from x'=f(x,p,t) to x'=f(x,p,t)/sqrt(|f(x,p,t)|^2+0.1)

struct NSAnnulusBoundaries{N,T<:Number,P}
    inner::Vector{IterationCurve{N,T,P}}
    outer::Vector{IterationCurve{N,T,P}}
end

function ns_inintialise_mesh(saddle, v1, v2, n, r)
    A = hcat(v1, v2)
    Q = SMatrix(qr(A).Q)
    newv1 = Q[:, 1]
    newv2 = Q[:, 2]
    rag = range(0, 1, length=n + 1)
    circle1 = [State(saddle, t) for t in rag]
    circle2 = [State(saddle + r * cospi(2 * t / 1) * newv1 + r * sinpi(2 * t / 1) * newv2, t) for t in rag]
    ic1 = IterationCurve(circle1, LinearInterpolation(circle1))
    ic2 = IterationCurve(circle2, LinearInterpolation(circle2))
    [AnnulusBoundaries([ic1], [ic2])]
end

function grow_manifold!(annulus, f, p, δ)
    oldu0 = deepcopy(annulus[end].outer.states)
    oldcurve = annulus[end].outer.pcurve
    k = length(oldu0)
    newu0 = similar(oldu0)
    for i in eachindex(newu0)
        newu0[i] = State(f(oldu0[i].state, p), oldu0[i].s)
    end
    j = 1
    while j + 1 <= k
        dist = norm(newu0[j] - newu0[j+1])
        if dist > δ
            m = ceil(Int, dist / δ)
            if m == 1
                m = 2
            end
            s1 = oldu0[j+1].s
            s0 = oldu0[j].s
            plengh = (s1 - s0) / m
            paras = [s0 + plengh * i for i in 1:m-1]
            preaddps = similar(oldu0, m - 1)
            addps = similar(oldu0, m - 1)
            for i in 1:m-1
                news = paras[i]
                preaddps[i] = State(oldcurve(news), news)
                addps[i] = State(f(oldcurve(news), p), news)
            end
            insert!(oldu0, j + 1, preaddps)
            insert!(newu0, j + 1, addps)
            k = k + m - 1
        else
            j = j + 1
        end
    end
    finalnewu0 = similar(newu0)
    datatype = typeof(newu0[1].s)
    ss0 = convert(datatype, 0)
    finalnewu0[1] = State(newu0[1].state, ss0)
    n0 = length(finalnewu0)
    for i in 2:n0
        ss0 = ss0 + norm(newu0[i] - newu0[i-1])
        finalnewu0[i] = State(newu0[i].state, ss0)
    end
    ic1 = IterationCurve(oldu0, LinearInterpolation(oldu0))
    ic2 = IterationCurve(finalnewu0, LinearInterpolation(finalnewu0))
    newannulus = AnnulusBoundaries(ic1, ic2)
    append!(annulus, [newannulus])
end

function generate_surface(f, p, saddle, v1, v2, N, r, δ; n=150)
    myannulus = inintialise_mesh(saddle, v1, v2, n, r)
    for i in 1:N
        grow_manifold!(myannulus, f, p, δ)
    end
    myannulus
end


function generate_surface(v::NSSetUp{S}, p, saddle, v1, v2, N, r, δ; n=150) where {S<:ContinuousVectorField}
    myannulus = inintialise_mesh(saddle, v1, v2, n, r)
    for i in 1:N
        grow_manifold!(myannulus, f, p, δ)
    end
    myannulus
end

function generate_surface(v::NSSetUp{S}, p, saddle, v1, v2, N, r, δ; n=150) where {S<:JumpVectorField}
    myannulus = inintialise_mesh(saddle, v1, v2, n, r)
    for i in 1:N
        grow_manifold!(myannulus, f, p, δ)
    end
    myannulus
end