struct Annulus{T}
    circles::Vector{T}
end


@inline function kd_distence(point::SVector{N,T}, someset::Vector{SVector{N,T}}) where {N,T}
    btree = KDTree(someset)
    nn(btree, point)[2]
end

@inline function kd_distence(someset1::Vector{SVector{N,T}}, someset2::Vector{SVector{N,T}}) where {N,T}
    btree = KDTree(someset1)
    dd = [nn(btree, x)[2] for x in someset2]
    maximum(dd)
end

function inintialise_mesh(f, para, saddle::SVector{S,T}, v1, v2, λ1, λ2, n, r, N, δ; interp=LinearInterpolation) where{S,T}
    newv1 = normalize(v1)
    newv2 = normalize(v2)
    rag = range(0, 1, length=n)
    newv1 = normalize(v1)
    newv2 = normalize(v2)
    res = (λ1 / λ2)^(N)
    circle1 = [saddle]
    circle2 = Vector{SVector{S,T}}(undef, n)
    for i in 1:n
        circle2[i] = saddle + (r/λ1)  * cospi(2 * rag[i]) * newv2 +
        (r/λ2) * res * sinpi(2 * rag[i]) * newv1
    end
    oldcurve = paramise(circle2, interp = interp)
    olds = copy(oldcurve.t)
    circle3 = similar(circle2)
    for i in eachindex(circle3)
        circle3[i] = f(circle2[i], para)
    end
    addpoints!(f, para, δ, oldcurve, circle3, olds)
    result = [oldcurve, interp(circle3, olds)]
    pre_result = [interp(circle1, T[0]), oldcurve]
    [Annulus(pre_result), Annulus(result)]
end

function grow_surface!(f, para, annulus::Vector{Annulus{S}}, δ1, δ2; interp=LinearInterpolation) where {S}
    data = annulus[end].circles
    newdata = deepcopy(data)
    k = length(newdata)
    for i in 2:k
        states = similar(data[i].u)
        ss = copy(data[i].t)
        for j in eachindex(data[i].u)
            states[j] = f(data[i].u[j], para)
        end
        newdata[i] = interp(states, ss)
    end
    newdata[1] = deepcopy(data[end])
    # first interpolate every circle
    for j in 2:k
        ic = newdata[j].u
        olds = newdata[j].t
        oldcurve = data[j]
        newpara = addpoints!(f, para, δ1, oldcurve, ic, olds)
        olds .= newpara
    end
    # then interpolate between circles
    p = 1
    copydata = deepcopy(data)
    @inbounds while p + 1 <= k
        dd = kd_distence(newdata[p].u, newdata[p+1].u)
        if dd > δ2
            insert_circle = copy(copydata[p+1].u)
            pre_insert_circle = similar(copydata[p+1].u)
            kd_tree = KDTree(copydata[p].u)
            for m in eachindex(insert_circle)
                point = insert_circle[m]
                number = nn(kd_tree, point)[1]
                point2 = copydata[p].u[number]
                pre_insert_circle[m] = (point + point2) / 2
                insert_circle[m] = f((point + point2) / 2, para)
            end
            pre_curve = paramise(pre_insert_circle, interp=interp)
            insert_s = copy(pre_curve.t)
            newpara = addpoints!(f, para, δ1, pre_curve, insert_circle, insert_s)
            insert_curve = interp(insert_circle, newpara)
            insert!(newdata, p + 1, insert_curve)
            insert!(copydata, p + 1, pre_curve)
            k = k + 1
        else
            p = p + 1
        end
    end
    append!(annulus, [Annulus(newdata)])
    nothing
end

function generate_surface(f, p, saddle, v1, v2, λ1, λ2, N, r, δ1, δ2; n=150, interp=LinearInterpolation)
    myannulus = inintialise_mesh(f, p, saddle, v1, v2, λ1, λ2, n, r, N, δ1, interp=interp)
    for i in 1:N
        grow_surface!(f, p, myannulus, δ1, δ2; interp=interp)
    end
    myannulus
end
