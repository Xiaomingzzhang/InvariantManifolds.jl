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
    mirrors = f.f.mirrors
    for i in eachindex(circle3)
        circle3[i] = f.timetmap(circle2[i], para)
    end
    newpara = ns_addpoints!(tmap, para, δ, oldcurve, circle3, olds, dtimes, tend, mirrors)
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
    result
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
        dist_collect = T[]
        for i in eachindex(outercircle)
            event = outercircle[i][1].event_at
            innerpieces = Vector{NSState{N,T}}[]
            for j in eachindex(innercircle)
                piece = deepcopy(innercircle[j])
                if piece[1].event_at == event
                    append!(innerpieces, [piece])
                elseif piece[1].event_at[1:end-1] == event
                    idx = piece[1].event_at[end]
                    themap = inverse_mirrors[idx]
                    for k in eachindex(piece)
                        piece[k] = NSState(themap(piece[k].state, para), piece[k].event_at, false, idx, -1)
                    end
                    append!(innerpieces, [piece])
                elseif piece[1].event_at == event[1:end-1]
                    idx = event[end]
                    themap = mirrors[idx]
                    for k in eachindex(piece)
                        piece[k] = NSState(themap(piece[k].state, para), piece[k].event_at, false, idx, 1)
                    end
                    append!(innerpieces, [piece])
                end
            end
            finalinnerpieces = vcat(innerpieces...)
            kdtree = KDTree(finalinnerpieces)
            dd= Vector{T}(undef, length(outercircle[i]))
            for mn in eachindex(dd)
                dd[mn] = (nn(kdtree, outercircle[i][mn].state))[2]
            end
            append!(dist_collect, [maximum(dd)])
        end
        maximum(dist_collect)
    end
end

@inline function kd_distence(f::NSSetUp{S}, para, innercircle::Vector{Vector{NSState{N,T}}},
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

function partition(v::Vector{NSState{N,T}}; interp = LinearInterpolation) where {N,T}
    ctime1 = v[1].event_at
    n = length(v)
    a0 = NSState{N,T}[]
    result = Vector{NSState{N,T}}[]
    for j in 1:n
        if v[j].event_at == ctime1
            append!(a0, [v[j]])
        else
            append!(result, [a0])
            a0 = [v[j]]
            ctime1 = v[j].event_at
        end
    end
    append!(result, [a0])
    [paramise(result[i], interp=interp) for i in eachindex(result)]
end

function partitionbytoward(v::Vector{NSState{N,T}}; interp = LinearInterpolation) where {N,T}
    ctime1 = v[1].toward
    n = length(v)
    a0 = NSState{N,T}[]
    result = Vector{NSState{N,T}}[]
    for j in 1:n
        if v[j].toward == ctime1
            append!(a0, [v[j]])
        else
            append!(result, [a0])
            a0 = [v[j]]
            ctime1 = v[j].toward
        end
    end
    append!(result, [a0])
    [paramise(result[i], interp=interp) for i in eachindex(result)]
end

# from two bkcircles return an bkcircle between them
function insert_bkcircles(f::NSSetUp{S}, para, innercircle::Vector{Vector{NSState{N,T}}},
    outercircle::Vector{Vector{NSState{N,T}}}; interp=LinearInterpolation) where {N,T,S}
    mirrors = f.f.mirrors
    inverse_mirrors = f.f.inverse_mirrors
    pre_insercircle = NSState{N,T}[]
    hypers = f.f.hypers
    for i in eachindex(outercircle)
        outerpiece = outercircle[i]
        event = outerpiece[1].event_at
        innerpieces = Vector{NSState{N,T}}[]
        for j in eachindex(innercircle)
            piece = innercircle[j]
            if piece[1].event_at == event
                append!(innerpieces, [piece])
            elseif piece[1].event_at[1:end-1] == event
                idx = piece[1].event_at[end]
                themap = inverse_mirrors[idx]
                data = similar(piece)
                for k in eachindex(data)
                    data[k] = NSState(themap(piece[k].state, para), piece[k].event_at, true, idx, -1)
                end
                append!(innerpieces, [data])
            elseif piece[1].event_at == event[1:end-1]
                idx = event[end]
                themap = mirrors[idx]
                data = similar(piece)
                for k in eachindex(data)
                    data[k] = NSState(themap(piece[k].state, para), piece[k].event_at, true, idx, 1)
                end 
                append!(innerpieces, [data])
            end
        end
        newpiece = vcat(innerpieces...)
        kdtree = KDTree(newpiece)
        for kk in eachindex(outerpiece)
            p1 = outerpiece[kk]
            mm = nn(kdtree, p1)[1]
            p2 = newpiece[mm]
            pre_newp = (p1+p2)/2
            if p2.ismirrored == false
                sides = [ff(p1.state, para)>0 for ff in hypers]
                if sides == [ff(pre_newp, para)>0 for ff in hypers]
                    newp = NSState(pre_newp, event, false, 0, 0)
                    append!(pre_insercircle, [newp])
                else
                    error("cannot find the intersect point's domain")
                end
            elseif p2.ismirrored == true && p2.toward == 1
                idx = p2.idx
                side = hypers[idx](p1.state, para)>0
                if side == hypers[idx](pre_newp, para)>0
                    newp = NSState(pre_newp, event, false, 0, 1)
                    append!(pre_insercircle, [newp])
                else
                    newp = NSState(inverse_mirrors[idx](pre_newp, para), p2.event_at, false, 0, 2)
                    append!(pre_insercircle, [newp])
                end
            elseif p2.ismirrored == true && p2.toward == -1
                idx = p2.idx
                side = hypers[idx](p1.state, para)>0
                if side == hypers[idx](pre_newp, para)>0
                    newp = NSState(pre_newp, event, false, 0, -1)
                    append!(pre_insercircle, [newp])
                else
                    newp = NSState(mirrors[idx](pre_newp, para), p2.event_at, false, 0, -2)
                    append!(pre_insercircle, [newp])
                end
            end
        end
    end 
    partitionbytoward(pre_insercircle, interp=interp)
end


function testinsert_bkcircles(f::NSSetUp{S}, para, innercircle::Vector{Vector{NSState{N,T}}},
    outercircle::Vector{Vector{NSState{N,T}}}; interp=LinearInterpolation) where {N,T,S}
    mirrors = f.f.mirrors
    inverse_mirrors = f.f.inverse_mirrors
    pre_insercircle = NSState{N,T}[]
    hypers = f.f.hypers
    for i in eachindex(outercircle)
        outerpiece = outercircle[i]
        event = outerpiece[1].event_at
        innerpieces = Vector{NSState{N,T}}[]
        for j in eachindex(innercircle)
            piece = innercircle[j]
            if piece[1].event_at == event
                append!(innerpieces, [piece])
            elseif piece[1].event_at[1:end-1] == event
                idx = piece[1].event_at[end]
                themap = inverse_mirrors[idx]
                data = similar(piece)
                for k in eachindex(data)
                    data[k] = NSState(themap(piece[k].state, para), piece[k].event_at, true, idx, -1)
                end
                append!(innerpieces, [data])
            elseif piece[1].event_at == event[1:end-1]
                idx = event[end]
                themap = mirrors[idx]
                data = similar(piece)
                for k in eachindex(data)
                    data[k] = NSState(themap(piece[k].state, para), piece[k].event_at, true, idx, 1)
                end
                append!(innerpieces, [data])
            end
        end
        newpiece = vcat(innerpieces...)
        kdtree = KDTree(newpiece)
        for kk in eachindex(outerpiece)
            p1 = outerpiece[kk]
            mm = nn(kdtree, p1)[1]
            p2 = newpiece[mm]
            pre_newp = (p1+p2)/2
            if p2.ismirrored == false
                sides = [ff(p1.state, para)>0 for ff in hypers]
                if sides == [ff(pre_newp, para)>0 for ff in hypers]
                    newp = NSState(pre_newp, event, false, 0, 0)
                    append!(pre_insercircle, [newp])
                else
                    error("cannot find the intersect point's domain")
                end
            elseif p2.ismirrored == true && p2.toward == 1
                idx = p2.idx
                side = hypers[idx](p1.state, para)>0
                if side == hypers[idx](pre_newp, para)>0
                    newp = NSState(pre_newp, event, false, 0, 1)
                    append!(pre_insercircle, [newp])
                else
                    newp = NSState(inverse_mirrors[idx](pre_newp, para), p2.event_at, false, 0, 2)
                    append!(pre_insercircle, [newp])
                end
            elseif p2.ismirrored == true && p2.toward == -1
                idx = p2.idx
                side = hypers[idx](p1.state, para)>0
                if side == hypers[idx](pre_newp, para)>0
                    newp = NSState(pre_newp, event, false, 0, -1)
                    append!(pre_insercircle, [newp])
                else
                    newp = NSState(mirrors[idx](pre_newp, para), p2.event_at, false, 0, -2)
                    append!(pre_insercircle, [newp])
                end
            end
        end
    end 
    pre_insercircle
end

take_u(x) = x.u
take_t(x) = x.t

function grow_surface!(f::NSSetUp, para, annulus::Vector{NSAnnulus{S}},
    δ1, δ2; dtimes=100, interp=LinearInterpolation) where {S}
    data = annulus[end].bkcircles
    newdata = deepcopy(data)
    k = length(newdata)
    tmap = f.timetmap
    tend = f.timespan[end]
    mirrors = f.f.mirrors
    for i in 2:k
        for j in eachindex(data[i])
            states = copy(data[i][j].u)
            ss = copy(data[i][j].t)
            for m in eachindex(states)
                states[m] = tmap(states[m], para)
            end
            newdata[i][j] = interp(states, ss)
        end
    end
    newdata[1] = deepcopy(data[end])
    # first interpolate every piece of bkcircle
    for j in 2:k
        for m in eachindex(newdata[j])
            ic = newdata[j][m].u
            olds = newdata[j][m].t
            oldcurve = data[j][m]
            newpara = ns_addpoints!(tmap, para, δ1, oldcurve, ic, olds, dtimes, tend, mirrors)
            olds .= newpara
        end
    end
    # then partition every bkcircle of newpara
    for i in eachindex(newdata)
        prenewdata = S[]
        for ii in eachindex(newdata[i])
            append!(prenewdata, partition(newdata[i][ii].u, newdata[i][ii].t, interp=interp))
        end
        newdata[i] = prenewdata
    end
    # then interpolate between circles
    p = 1
    copydata = deepcopy(data)
    @inbounds while p + 1 <= k
        innercircle = take_u.(newdata[p])
        outercircle = take_u.(newdata[p+1])
        d = kd_distence(f, para, innercircle, outercircle)
        @show d
        if d > δ2
            pre_innercircle = take_u.(copydata[p])
            pre_outercircle = take_u.(copydata[p+1])
            pre_insert_circle = insert_bkcircles(f, para, pre_innercircle, pre_outercircle; interp=interp)
            d1 = kd_distence(f, para, pre_innercircle, take_u.(pre_insert_circle))
            d2 = kd_distence(f, para, take_u.(pre_insert_circle), pre_outercircle)
            if d1>δ2 || d2>δ2
                error("Incorect insert circle")
            end
            insert_circle = take_u.(deepcopy(pre_insert_circle))
            insert_circle_t = take_t.(deepcopy(pre_insert_circle))
            final_insert_circle = S[]
            for i in eachindex(insert_circle)
                for j in eachindex(insert_circle[i])
                    insert_circle[i][j] = tmap(insert_circle[i][j], para)
                end
            end
            for i in eachindex(insert_circle)
                newpara = ns_addpoints!(tmap, para, δ1, pre_insert_circle[i], insert_circle[i], insert_circle_t[i], dtimes, tend, mirrors)
                _result = partition(insert_circle[i], newpara, interp=interp)
                append!(final_insert_circle, _result)
            end
            insert!(newdata, p + 1, final_insert_circle)
            insert!(copydata, p + 1, pre_insert_circle)
            k = k + 1
        else
            p = p + 1
        end
    end
    append!(annulus, [NSAnnulus(newdata)])
    nothing
end


function generate_surface(f::NSSetUp, p, saddle, v1, v2, λ1, λ2, N, r, δ1, δ2;
    n=150, dtimes=100, interp=LinearInterpolation)
    myannulus = ns_inintialise_mesh(f, p, saddle, v1, v2, λ1, λ2, n, r, N, δ1, dtimes=dtimes, interp=interp)
    for i in 1:N
        grow_surface!(f, p, myannulus, δ1, δ2; interp=interp, dtimes=dtimes)
    end
    myannulus
end