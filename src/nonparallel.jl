module AuxFunctions

include("common.jl")

export use_threading, compute_weights, compute_DI, compute_FN, remove_duplicate_seqs

use_threading(x) = nothing

function compute_theta(cZ::Vector{Vector{UInt64}}, N::Int, M::Int)

    const cl = clength(N)
    const cr = 5 * (packfactor - crest(N)) + packrest

    const kmax = div(cl - 1, 31)
    const rmax = (cl - 1) % 31

    meanfracid = 0.0

    # count seqs below theta dist
    for i = 1:M-1
        cZi = unsafe(cZ[i])
        nids::UInt64 = 0
        for j = i+1:M
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
    meanfracid /= 0.5 * M * (M-1)
    theta = min(0.5, 0.38 * 0.32 / meanfracid)

    return theta
end

function compute_theta(ZZ::Vector{Vector{Int8}}, N::Int, M::Int)
    meanfracid = 0.0
    for i = 1:M-1
        Zi = ZZ[i]
        for j = i+1:M
            Zj = ZZ[j]
            nids = 0
            for k = 1:N
                nids += Zi[k] == Zj[k]
            end
            fracid = nids / N
            meanfracid += fracid
        end
    end
    meanfracid /= 0.5 * M * (M-1)
    theta = min(0.5, 0.38 * 0.32 / meanfracid)
    return theta
end

@compat function compute_weights(cZ::Vector{Vector{UInt64}}, theta::Real, N::Int, M::Int)

    theta = Float64(theta)

    const cl = clength(N)
    const kmax = div(cl - 1, 31)
    const rmax = (cl - 1) % 31 + 1

    Meff = 0.0

    W = ones(M)

    thresh = floor(theta * N)
    println("theta = $theta threshold = $thresh")

    if theta == 0
        println("M = $M N = $N Meff = $M")
        return W, Float64(M)
    end

    for i = 1:M-1
        cZi = unsafe(cZ[i])
        for j = i+1:M
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
    for i = 1:M
        W[i] = 1 / W[i]
    end
    Meff = sum(W)
    println("M = $M N = $N Meff = $Meff")
    return W, Meff
end

function compute_weights(ZZ::Vector{Vector{Int8}}, theta::Float64, N::Int, M::Int)

    Meff = 0.

    W = ones(M)

    thresh = floor(theta * N)
    println("theta = $theta threshold = $thresh")

    if theta == 0
        println("M = $M N = $N Meff = $M")
        return W, float64(M)
    end

    for i = 1:M-1
        Zi = ZZ[i]
        for j = i+1:M
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
            end
        end
    end
    for i = 1:M
        W[i] = 1. / W[i]
        Meff += W[i]
    end
    println("M = $M N = $N Meff = $Meff")
    return W, Meff
end

function compute_dists(cZ::Vector{Vector{UInt64}}, N::Int, M::Int)

    const cl = clength(N)
    const kmax = div(cl - 1, 31)
    const rmax = (cl - 1) % 31 + 1

    D = zeros(Float16, M, M)

    for i = 1:M-1
        cZi = unsafe(cZ[i])
        for j = i+1:M
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
            D[i,j] = dist / N
            D[j,i] = D[i,j]
        end
    end
    return D
end

if VERSION ≥ v"0.4-"
    typealias KT Symmetric{Float64,Matrix{Float64}}
else
    typealias KT Matrix{Float64}
end

function compute_DI(mJ::Matrix{Float64}, C::Matrix{Float64}, N::Int, q::Integer)

    DI = zeros(N, N)
    s = q - 1
    #Is = eye(s)

    iKs = Array(KT, N)
    rowi = 0
    for i = 1:N
        row = rowi + (1:s)
        rowi += s

        iKs[i] = sqrtm(Symmetric(C[row,row]))
    end

    z = 0.5 * s * log(0.5)

    rowi = 0
    for i = 1:N-1
        row = rowi + (1:s)
        rowi += s

        invsqrtKi = iKs[i]

        coli = rowi
        for j = i + 1 : N
            col = coli + (1:s)
            coli += s

            invsqrtKj = iKs[j]

            mJij = mJ[row, col]
            if sum(mJij) != 0
                MM = invsqrtKi * mJij * invsqrtKj
                #V = Is + 4 * (MM * MM')
                V = MM * MM'
                #X = Is + sqrtm(Symmetric(V))
                eigV = eig(V)[1]
                eigX = sqrt(1 .+ 4 * eigV)
                #DI[i,j] = z + 0.5 * log(det(X))
                DI[i,j] = z + 0.5 * sum(log(1 .+ eigX))
                DI[j,i] = DI[i,j]
            end
        end
    end

    return DI
end

end
