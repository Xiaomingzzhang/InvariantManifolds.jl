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

# 第一次生成初始曲线的时候就要进行插值, 否则精度差太多
# function inintialise_mesh(f, para, saddle, v1, v2, λ1, λ2, n, r, k, N, δ1, δ2; interp=LinearInterpolation)
#     newv1 = normalize(v1)
#     newv2 = normalize(v2)
#     rag = range(0, 1, length=n)
#     res = (λ1 / λ2)^(N + 1)
#     circle1 = [copy(saddle) for i in 1:n]
#     data = [copy(circle1) for i in 1:k]
#     radius_range = range(0, 1, length=k)
#     for i in 1:k
#         for j in 1:n
#             data[i][j] = saddle + r * radius_range[i] * cospi(2 * rag[j]) * newv2 +
#                          r * res * radius_range[i] * sinpi(2 * rag[j]) * newv1
#         end
#     end
#     newdata = deepcopy(data)
#     for i in 2:k
#         for j in eachindex(data[i])
#             newdata[i][j] = f(data[i][j], para)
#         end
#     end
#     para_data = [paramise(circle1, interp=interp)]
#     for i in 2:k
#         append!(para_data, [paramise(data[i], interp=interp)])
#     end
#     para_newdata = [paramise(newdata[1], interp=interp)]
#     for i in 2:k
#         append!(para_newdata, [paramise(newdata[i], interp=interp)])
#     end
#     # first interpolate every circle in para_newdata
#     for j in 2:k
#         ic = para_newdata[j].u
#         olds = para_newdata[j].t
#         oldcurve = para_data[j]
#         addpoints!(f, para, δ1, oldcurve, ic, olds)
#     end
#     # then interpolate between circles
#     p = 1
#     q = 0
#     pre_curves = deepcopy(para_data)
#     @inbounds while p + 1 <= k && q < 50
#         dist = kd_distence(para_newdata[p].u, para_newdata[p+1].u)
#         @show dist
#         @show q
#         @show p
#         if dist > δ2
#             insert_circle = similar(pre_curves[p].u)
#             pre_insert_circle = similar(pre_curves[p].u)
#             insert_s = copy(pre_curves[p].t)
#             for m in eachindex(insert_circle)
#                 point = pre_curves[p].u[m]
#                 point2 = pre_curves[p+1].u[m]
#                 pre_insert_circle[m] = (point + point2) / 2
#                 insert_circle[m] = f(pre_insert_circle[m], para)
#             end
#             pre_insert_s = copy(insert_s)
#             pre_curve = interp(pre_insert_circle, pre_insert_s)
#             addpoints!(f, para, δ1, pre_curve, insert_circle, insert_s)
#             insert_curve = interp(insert_circle, insert_s)
#             insert!(para_newdata, p + 1, insert_curve)
#             insert!(pre_curves, p + 1, pre_curve)
#             k = k + 1
#             q = q + 1
#         else
#             p = p + 1
#             q = 0
#         end
#     end
#     d = kd_distence(saddle, para_data[end].u)
#     j = 1
#     while kd_distence(saddle, para_newdata[j].u) < d
#         j = j + 1
#     end
#     result = para_newdata[j:end]
#     prepend!(result, [para_data[end]])
#     [Annulus(para_data), Annulus(result)]
# end

# 越简单越好
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
        # @show length(states)
        # @show length(ss)
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
    # @show typeof(collect_s)
    # then interpolate between circles
    p = 1
    copydata = deepcopy(data)
    @inbounds while p + 1 <= k
        dd = kd_distence(newdata[p].u, newdata[p+1].u)
        # @show dd
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
