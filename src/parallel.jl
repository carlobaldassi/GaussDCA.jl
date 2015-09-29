module AuxFunctions

include("common.jl")

export use_threading, compute_weights, compute_DI, compute_FN, remove_duplicate_seqs

@compat use_threading(x::Bool) = blas_set_num_threads(x ? Int(get(ENV, "OMP_NUM_THREADS", CPU_CORES)) : 1)

@compat typealias TriuInd Tuple{Tuple{Int,Int},Tuple{Int,Int},Int}

function ptriu(sz::Int, RT::Type, func::Function, args...)

    tot_inds = div(sz * (sz-1), 2)
    nw = nworkers()

    if tot_inds >= nw
        inds_dist = diff(round(Int, linspace(1, tot_inds+1, nw+1)))
    else
        inds_dist = [ones(Int, tot_inds), zeros(Int, nw-tot_inds)]
    end

    inds = Array(TriuInd, nw)
    i0, j0 = 1, 2
    for p = 1:nw
        l = inds_dist[p]
        I1 = (i0, j0)
        I3 = l
        i1, j1 = i0, j0
        while l > 0
            if sz-j1+1 < l
                l -= sz-j1+1
                i1 += 1
                j1 = i1 + 1
            else
                j1 += l - 1
                l = 0
            end
        end
        I2 = (i1, j1)
        inds[p] = (I1, I2, I3)
        i0, j0 = i1, j1
        j0 += 1
        if j0 > sz
            i0 += 1
            j0 = i0 + 1
        end
    end

    ret = Array(RT, nw)
    wrk = workers()

    use_threading(false)

    @sync begin
        for p = 1:nw
            @async begin
                ret[p] = remotecall_fetch(wrk[p], func, inds[p], args...)
            end
        end
    end

    use_threading(true)

    return ret, inds
end

function ptriu_compose{T}(src::Vector{Vector{T}}, sz::Int, inds::Vector{TriuInd})

    dest = Array(T, sz, sz)

    @assert length(src) == length(inds)

    for p = 1:length(src)
        i0, j0 = inds[p][1]
        i1, j1 = inds[p][2]
        srcp = src[p]
        l = 1
        for i = i0 : i1
            jj0 = i==i0 ? j0 : i+1
            jj1 = i==i1 ? j1 : sz
            for j = jj0 : jj1
                x = srcp[l]
                dest[i,j] = x
                dest[j,i] = x
                l += 1
            end
        end
    end
    for i = 1:sz
        dest[i,i] = 0
    end

    return dest
end

function compute_theta_chunk(inds::TriuInd, cZ::Vector{Vector{UInt64}}, N::Int, M::Int)

    const cl = clength(N)
    const cr = 5 * (packfactor - crest(N)) + packrest

    const kmax = div(cl - 1, 31)
    const rmax = (cl - 1) % 31

    meanfracid = 0.0
    i0, j0 = inds[1]
    i1, j1 = inds[2]
    for i = i0 : i1
        cZi = unsafe(cZ[i])
        nids::UInt64 = 0
        jj0 = i==i0 ? j0 : i+1
        jj1 = i==i1 ? j1 : M
        for j = jj0 : jj1
            cZj = unsafe(cZ[j])

            czi = start(cZi)
            czj = start(cZj)

            z::UInt64 = 0

            for k = 1:kmax
                z = 0
                for r = 1:31
                    zi, czi = next(cZi, czi)
                    zj, czj = next(cZj, czj)

                    ny = (~zi) $ zj
                    z += nz_aux(ny)
                end
                t = @collapse(z)
                nids += t
            end
            z = 0
            for r = 1:rmax
                zi, czi = next(cZi, czi)
                zj, czj = next(cZj, czj)
                ny = (~zi) $ zj
                z += nz_aux(ny)
            end
            zi, czi = next(cZi, czi)
            zj, czj = next(cZj, czj)
            ny = (~zi) $ zj
            z += nz_aux2(ny, cr)
            t = @collapse(z)
            nids += t
        end
        fracid = nids / N
        meanfracid += fracid
    end
    return meanfracid
end

@compat function compute_theta{T<:Union{Int8,UInt64}}(cZ::Vector{Vector{T}}, N::Int, M::Int)

    chunk_means, _ = ptriu(M, Float64, compute_theta_chunk, cZ, N, M)
    meanfracid = sum(chunk_means) / (0.5 * M * (M-1))
    theta = min(0.5, 0.38 * 0.32 / meanfracid)

    return theta
end

# slow fallback
function compute_theta_chunk(inds::TriuInd, ZZ::Vector{Vector{Int8}}, N::Int, M::Int)
    meanfracid = 0.0
    i0, j0 = inds[1]
    i1, j1 = inds[2]
    for i = i0 : i1
        Zi = ZZ[i]
        jj0 = i==i0 ? j0 : i+1
        jj1 = i==i1 ? j1 : M
        for j = jj0 : jj1
            Zj = ZZ[j]
            nids = 0
            for k = 1:N
                nids += Zi[k] == Zj[k]
            end
            fracid = nids / N
            meanfracid += fracid
        end
    end
    return meanfracid
end

function compute_weights_chunk(inds::TriuInd, cZ::Vector{Vector{UInt64}}, thresh::Real, N::Int, M::Int)

    const cl = clength(N)
    const kmax = div(cl - 1, 31)
    const rmax = (cl - 1) % 31 + 1

    W = zeros(M)
    i0, j0 = inds[1]
    i1, j1 = inds[2]

    for i = i0:i1
        cZi = unsafe(cZ[i])
        jj0 = i==i0 ? j0 : i+1
        jj1 = i==i1 ? j1 : M
        for j = jj0 : jj1
            cZj = unsafe(cZ[j])

            czi = start(cZi)
            czj = start(cZj)
            dist::UInt64 = 0
            z::UInt64 = 0
            for k = 1:kmax
                z = 0
                for r = 1:31
                    zi, czi = next(cZi, czi)
                    zj, czj = next(cZj, czj)

                    y = zi $ zj
                    z += nnz_aux(y)
                end
                t = @collapse(z)
                dist += t
                dist >= thresh && break
            end
            if dist < thresh
                z = 0
                for r = 1:rmax
                    zi, czi = next(cZi, czi)
                    zj, czj = next(cZj, czj)

                    y = zi $ zj
                    z += nnz_aux(y)
                end
                t = @collapse(z)
                dist += t
            end
            if dist < thresh
                W[i] += 1
                W[j] += 1
            end
        end
    end
    return W
end

@compat function compute_weights{T<:Union{Int8,UInt64}}(cZ::Vector{Vector{T}}, theta::Real, N::Int, M::Int)

    theta = Float64(theta)

    Meff = 0.0

    thresh = floor(theta * N)
    println("theta = $theta threshold = $thresh")

    if theta == 0
        println("M = $M N = $N Meff = $M")
        return W, Float64(M)
    end

    Ws, _ = ptriu(M, Vector{Float64}, compute_weights_chunk, cZ, thresh, N, M)
    W = sum(Ws)

    for i = 1:M
        W[i] = 1 / (1 + W[i])
        Meff += W[i]
    end
    println("M = $M N = $N Meff = $Meff")
    return W, Meff
end

# slow fallback
function compute_weights_chunk(inds::TriuInd, ZZ::Vector{Vector{Int8}}, thresh::Real, N::Int, M::Int)
    W = zeros(M)
    i0, j0 = inds[1]
    i1, j1 = inds[2]
    for i = i0 : i1
        Zi = ZZ[i]
        jj0 = i==i0 ? j0 : i+1
        jj1 = i==i1 ? j1 : M
        for j = jj0 : jj1
            Zj = ZZ[j]
            dist = 0
            k = 1
            while dist < thresh && k <= N
                dist += Zi[k] != Zj[k]
                k += 1
            end
            if dist < thresh
                W[i] += 1.
                W[j] += 1.
                #ret[j-i] = true
            end
        end
    end
    return W
end

@compat function compute_dists_chunk(inds::TriuInd, cZ::Vector{Vector{UInt64}}, N::Int, M::Int)

    const cl = clength(N)
    const kmax = div(cl - 1, 31)
    const rmax = (cl - 1) % 31 + 1

    D = zeros(Float16, inds[3])

    i0, j0 = inds[1]
    i1, j1 = inds[2]

    l = 1
    for i = i0:i1
        cZi = unsafe(cZ[i])
        jj0 = i==i0 ? j0 : i+1
        jj1 = i==i1 ? j1 : M
        for j = jj0 : jj1
            cZj = unsafe(cZ[j])

            czi = start(cZi)
            czj = start(cZj)
            dist::UInt64 = 0
            z::UInt64 = 0
            for k = 1:kmax
                z = 0
                for r = 1:31
                    zi, czi = next(cZi, czi)
                    zj, czj = next(cZj, czj)

                    y = zi $ zj
                    z += nnz_aux(y)
                end
                t = @collapse(z)
                dist += t
            end
            z = 0
            for r = 1:rmax
                zi, czi = next(cZi, czi)
                zj, czj = next(cZj, czj)

                y = zi $ zj
                z += nnz_aux(y)
            end
            t = @collapse(z)
            dist += t
            D[l] = dist / N
            l += 1
        end
    end

    return D
end

@compat function compute_dists(cZ::Vector{Vector{UInt64}}, N::Int, M::Int)

    Ds, inds = ptriu(M, Vector{Float16}, compute_dists_chunk, cZ, N, M)
    D = ptriu_compose(Ds, M, inds)

    return D
end

function compute_DI_chunk(inds::TriuInd, N::Int, s::Integer, iKs::Vector{Matrix{Float64}}, mJ::Matrix{Float64}, z::Float64)

    DI = Array(Float64, inds[3])
    i0, j0 = inds[1]
    i1, j1 = inds[2]
    l = 1
    for i = i0 : i1
        row = (i-1)*s + 1 : i*s

        invsqrtKi = iKs[i]

        jj0 = i==i0 ? j0 : i+1
        jj1 = i==i1 ? j1 : N
        for j = jj0 : jj1
            col = (j-1)*s + 1 : j*s

            invsqrtKj = iKs[j]

            mJij = mJ[row, col]
            if sum(mJij) != 0
                MM = invsqrtKi * mJij * invsqrtKj
                #V = Is + 4 * (MM * MM')
                V = MM * MM'
                #X = Is + sqrtm(Symmetric(V))
                eigV = eig(V)[1]
                eigX = sqrt(1 .+ 4 * eigV)
                #DI[l] = z + 0.5 * log(det(X))
                DI[l] = z + 0.5 * sum(log(1 .+ eigX))
            else
                DI[l] = 0
            end
            l += 1
        end
    end
    return DI
end

if VERSION ≥ v"0.4-"
    typealias KT Symmetric{Float64,Matrix{Float64}}
else
    typealias KT Matrix{Float64}
end

function compute_DI(mJ::Matrix{Float64}, C::Matrix{Float64}, N::Int, q::Integer)

    s = q - 1

    iKs = Array(KT, N)
    rowi = 0
    for i = 1:N
        row = rowi + (1:s)
        rowi += s

        iKs[i] = sqrtm(Symmetric(C[row,row]))
    end

    z = 0.5 * s * log(0.5)

    DIs, inds = ptriu(N, Vector{Float64}, compute_DI_chunk, N, s, iKs, mJ, z)
    DI = ptriu_compose(DIs, N, inds)

    return DI
end

end
