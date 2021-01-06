module AuxFunctions

using LinearAlgebra
using Base.Threads

include("common.jl")

export use_threading, compute_weights, compute_DI, compute_FN, remove_duplicate_seqs

use_threading(x) = nothing

const chunk_size = 64


function compute_theta(cZ::Vector{Vector{UInt64}}, N::Int, M::Int)
    cl = clength(N)
    cr = 5 * (packfactor - crest(N)) + packrest

    kmax = (cl - 1) ÷ 31
    rmax = (cl - 1) % 31

    meanfracid = 0.0

    perthread_fracids = zeros(nthreads())

    @inbounds for i = 1:M-1
        cZi = cZ[i]
        # nids::UInt64 = 0
        Threads.@threads for j1 = i+1:chunk_size:M
            for j = j1:min(j1+chunk_size-1, M)
                cZj = cZ[j]
                czi, czj = 1, 1
                z::UInt64 = 0
                for k = 1:kmax
                    z = 0
                    for r = 1:31
                        zi, zj = cZi[czi], cZj[czj]
                        czi += 1
                        czj += 1

                        ny = (~zi) ⊻ zj
                        z += nz_aux(ny)
                    end
                    t = collapse(z)
                    # nids += t
                    perthread_fracids[Threads.threadid()] += t / N
                end
                z = 0
                for r = 1:rmax
                    zi, zj = cZi[czi], cZj[czj]
                    czi += 1
                    czj += 1
                    ny = (~zi) ⊻ zj
                    z += nz_aux(ny)
                end
                zi, zj = cZi[czi], cZj[czj]
                czi += 1
                czj += 1
                ny = (~zi) ⊻ zj
                z += nz_aux2(ny, cr)
                t = collapse(z)
                # nids += t
                perthread_fracids[Threads.threadid()] += t / N
            end
        end
        # fracid = nids / N
        # meanfracid += fracid
    end
    meanfracid = sum(perthread_fracids) / (M * (M - 1) ÷ 2)
    # meanfracid /= 0.5 * M * (M-1)
    θ = min(0.5, 0.38 * 0.32 / meanfracid)

    return θ
end

# slow fallback
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
    θ = min(0.5, 0.38 * 0.32 / meanfracid)
    return θ
end

function compute_weights(cZ::Vector{Vector{UInt64}}, θ::Real, N::Int, M::Int)
    θ = Float64(θ)

    cl = clength(N)
    kmax = (cl - 1) ÷ 31
    rmax = (cl - 1) % 31 + 1

    Meff = 0.0

    W = ones(M)

    thresh = floor(θ * N)
    println("θ = $θ threshold = $thresh")

    if θ == 0
        println("M = $M N = $N Meff = $M")
        return W, Float64(M)
    end

    @inbounds for i = 1:M-1
        cZi = cZ[i]
        for j = i+1:M
            cZj = cZ[j]
            czi, czj = 1, 1
            dist::UInt64 = 0
            z::UInt64 = 0
            for k = 1:kmax
                z = 0
                for r = 1:31
                    zi, zj = cZi[czi], cZj[czj]
                    czi += 1
                    czj += 1

                    y = zi ⊻ zj
                    z += nnz_aux(y)
                end
                t = collapse(z)
                dist += t
                dist >= thresh && break
            end
            if dist < thresh
                z = 0
                for r = 1:rmax
                    zi, zj = cZi[czi], cZj[czj]
                    czi += 1
                    czj += 1

                    y = zi ⊻ zj
                    z += nnz_aux(y)
                end
                t = collapse(z)
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

# slow fallback
function compute_weights(ZZ::Vector{Vector{Int8}}, θ::Float64, N::Int, M::Int)
    Meff = 0.0

    W = ones(M)

    thresh = floor(θ * N)
    println("θ = $θ threshold = $thresh")

    if θ == 0
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
    cl = clength(N)
    kmax = (cl - 1) ÷ 31
    rmax = (cl - 1) % 31 + 1

    D = zeros(Float16, M, M)

    @inbounds for i = 1:M-1
        cZi = cZ[i]
        for j = i+1:M
            cZj = cZ[j]
            czi, czj = 1, 1
            dist::UInt64 = 0
            z::UInt64 = 0
            for k = 1:kmax
                z = 0
                for r = 1:31
                    zi, zj = cZi[czi], cZj[czj]
                    czi += 1
                    czj += 1

                    y = zi ⊻ zj
                    z += nnz_aux(y)
                end
                t = collapse(z)
                dist += t
            end
            z = 0
            for r = 1:rmax
                zi, zj = cZi[czi], cZj[czj]
                czi += 1
                czj += 1

                y = zi ⊻ zj
                z += nnz_aux(y)
            end
            t = collapse(z)
            dist += t
            D[i,j] = dist / N
            D[j,i] = D[i,j]
        end
    end
    return D
end

const KT = Symmetric{Float64,Matrix{Float64}}

function compute_DI(mJ::Matrix{Float64}, C::Matrix{Float64}, N::Int, q::Integer)
    DI = zeros(N, N)
    s = q - 1
    iKs = Array{KT}(undef, N)
    rowi = 0
    for i = 1:N
        row = rowi .+ (1:s)
        rowi += s
        iKs[i] = √(Symmetric(C[row,row]))
    end

    z = 0.5 * s * log(0.5)

    rowi = 0
    for i = 1:N-1
        row = rowi .+ (1:s)
        rowi += s
        invsqrtKi = iKs[i]
        coli = rowi
        for j = i + 1 : N
            col = coli .+ (1:s)
            coli += s
            invsqrtKj = iKs[j]
            mJij = mJ[row, col]
            if sum(mJij) ≠ 0
                MM = invsqrtKi * mJij * invsqrtKj
                V = MM * MM'
                eigV = eigvals(V)
                DI[i,j] = z + 0.5 * sum(@. log1p(√(1 + 4 * eigV)))
                DI[j,i] = DI[i,j]
            end
        end
    end

    return DI
end

end # module
