# The following algorithm only works with non-smooth vector fields's time-T-map


struct NSAnnulus{T}
    bkcircles::Vector{Vector{T}}
end


function ns_inintialise_mesh(f::NSSetUp, para, saddle::SVector{S,T},
    v1, v2, λ1, λ2, n, r, N, δ; dtimes=100, interp=LinearInterpolation) where {S,T}
    newv1 = normalize(v1)
    newv2 = normalize(v2)
    rag = range(0, 1, length=n)
    newv1 = normalize(v1)
    newv2 = normalize(v2)
    res = (λ1 / λ2)^(N)
    circle1 = [NSState(saddle)]
    circle2 = Vector{NSState{S,T}}(undef, n)
    for i in 1:n
        circle2[i] = NSState(saddle + (r / λ1) * cospi(2 * rag[i]) * newv2 +
                             (r / λ2) * res * sinpi(2 * rag[i]) * newv1)
    end
    oldcurve = paramise(circle2, interp=interp)
    olds = copy(oldcurve.t)
    circle3 = similar(circle2)
    tmap = f.timetmap
    tend = f.timespan[end]
    rules = f.f.rules
    for i in eachindex(circle3)
        circle3[i] = f.timetmap(circle2[i].state, para)
    end
    newpara = ns_addpoints!(tmap, para, δ, oldcurve, circle3, olds, dtimes, tend, rules)
    if ispartitioned(circle3) == false
        error("The initial radius has to be chosen more small")
    end
    result = [[oldcurve], [interp(circle3, newpara)]]
    pre_result = [[interp(circle1, T[0])], [oldcurve]]
    [NSAnnulus(pre_result), NSAnnulus(result)]
end

# 首先不论是外圈还是内圈都是已经分割的点列Vector{Vector{NSState{N,T}}}, 即每个Vector{NSState{N,T}}
# 中的事件都是一样的;
# 取内圈的一个片段, 寻找外圈中与其事件相同和相差为1的片段: 其中事件相差为1的片段要使用mirrors函数映射一下.
# 使用kd树判断内圈片段的每个点到外圈所取点的距离, 取最大值. 再对所有的内圈进行循环, 取整体距离的最大值即可.

# 插值时, 取内圈的一个片段, 取外圈与该内圈事件相同及相差1的片段,然后仍需镜像映射相差1的外圈片段, 在这两个片段之间进行插值.
# 要求1: 对于内圈与外圈中的两个点, 取中点得到的插值点必须与这两个点位于区域的同一侧, 否则我们需要作适当的处理使得其
# 位于同一侧.
# 要求2: 对于片段之间的端点, 需要使用特殊的插值使得插值点也在切换面附近

function take_event(v::Vector{Vector{NSState{N,T}}}) where {N,T}
    l = length(v)
    result = Vector{Vector{Int64}}(undef, l)
    for i in 1:l
        result[i] = v[i][1].event_at
    end
end

take_event(x::NSState) = x.event_at

@inline function kd_distence(f::NSSetUp{S}, para, innercircle::Vector{Vector{NSState{N,T}}},
    outercircle::Vector{Vector{NSState{N,T}}}) where {N,T,S<:JumpVectorField}
    ievent_data = take_event(innercircle)
    oevent_data = take_event(outercircle)
    i_maxl = maximum(length.(ievent_data))
    i_minl = minimum(length.(ievent_data))
    o_maxl = maximum(length.(oevent_data))
    o_minl = minimum(length.(oevent_data))
    if abs(i_maxl - o_maxl) >= 2 || abs(i_minl - o_minl) >= 2
        T(2)
    else
        mirrors = f.f.mirrors
        inverse_mirrors = f.f.inverse_mirrors
        dist_collect = Vector{T}(undef, length(innercircle))
        for i in eachindex(innercircle)
            kdtree = KDTree(innercircle[i])
            event = innercircle[i][1].event_at
            outpiece = Vector{NSState{N,T}}[]
            for j in eachindex(outercircle)
                piece = outercircle[j]
                if piece[1].event_at == event
                    append!(outpiece, [piece])
                elseif piece[1].event_at[1:end-1] == event
                    themap = inverse_mirrors[piece[1].event_at[end]]
                    data = similar(piece)
                    for k in eachindex(data)
                        data[k] = NSState(themap(piece[k].state, para), piece[k].event_at)
                    end
                    append!(outpiece, [data])
                elseif piece[1].event_at == event[1:end-1]
                    themap = mirrors[event[end]]
                    data = similar(piece)
                    for k in eachindex(data)
                        data[k] = NSState(themap(piece[k].state, para), piece[k].event_at)
                    end
                    append!(outpiece, [data])
                end
            end
            outpiece = vcat(outpiece...)
            dd = [nn(kdtree, x)[2] for x in outpiece]
            append!(dist_collect, [maximum(dd)])
        end
        maximum(dist_collect)
    end
end

@inline function kd_distence(f::NSSetUp{S}, innercircle::Vector{Vector{NSState{N,T}}},
    outercircle::Vector{Vector{NSState{N,T}}}) where {N,T,S<:ContinuousVectorField}
    newic = vcat(innercircle...)
    newoc = vcat(outercircle...)
    ievent_data = take_event(innercircle)
    oevent_data = take_event(outercircle)
    i_maxl = maximum(length.(ievent_data))
    i_minl = minimum(length.(ievent_data))
    o_maxl = maximum(length.(oevent_data))
    o_minl = minimum(length.(oevent_data))
    if abs(i_maxl - o_maxl) >= 2 || abs(i_minl - o_minl) >= 2
        T(2)
    else
        kdtree = KDTree(newic)
        dd = [nn(kdtree, x)[2] for x in newoc]
        maximum(dd)
    end
end

# from two bkcircles return an bkcircle between them
function insert_bkcircles(f::NSSetUp{S}, para, innercircle::Vector{Vector{NSState{N,T}}},
    outercircle::Vector{Vector{NSState{N,T}}}) where {N,T,S}
    mirrors = f.f.mirrors
    inverse_mirrors = f.f.inverse_mirrors
    dist_collect = Vector{T}(undef, length(innercircle))
    for i in eachindex(innercircle)
        kdtree = KDTree(innercircle[i])
        event = innercircle[i][1].event_at
        outpiece = Vector{NSState{N,T}}[]
        for j in eachindex(outercircle)
            piece = outercircle[j]
            if piece[1].event_at == event
                append!(outpiece, [piece])
            elseif piece[1].event_at[1:end-1] == event
                themap = inverse_mirrors[piece[1].event_at[end]]
                data = similar(piece)
                for k in eachindex(data)
                    data[k] = NSState(themap(piece[k].state, para), piece[k].event_at)
                end
                append!(outpiece, [data])
            elseif piece[1].event_at == event[1:end-1]
                themap = mirrors[event[end]]
                data = similar(piece)
                for k in eachindex(data)
                    data[k] = NSState(themap(piece[k].state, para), piece[k].event_at)
                end
                append!(outpiece, [data])
            end
        end
        outpiece = vcat(outpiece...)
        dd = [nn(kdtree, x)[2] for x in outpiece]
        append!(dist_collect, [maximum(dd)])
    end
    maximum(dist_collect)
end


function grow_surface!(f::NSSetUp, para, annulus::Vector{NSAnnulus{S}},
    δ1, δ2; dtimes=100, interp=LinearInterpolation) where {S}
    data = annulus[end].bkcircles
    newdata = deepcopy(data)
    k = length(newdata)
    tmap = f.timetmap
    tend = f.timespan[end]
    rules = f.f.rules
    irules = f.f.irules
    for i in 2:k
        for j in eachindex(data[i])
            states = similar(data[i][j].u)
            ss = copy(data[i][j].t)
            for m in eachindex(data[i][j].u)
                states[j] = tmap(data[i][j].u[m].state, para)
            end
            newdata[i][j] = interp(states, ss)
        end
    end
    newdata[1] = deepcopy(data[end])
    # first interpolate every piece of bkcircle but no partition yet
    for j in 2:k
        for m in newdata[j]
            ic = newdata[j][m].u
            olds = newdata[j][m].t
            oldcurve = data[j][m]
            newpara = ns_addpoints!(tmap, para, δ1, oldcurve, ic, olds, dtimes, tend, rules)
            olds .= newpara
        end
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
