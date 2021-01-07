module AuxFunctions

using LinearAlgebra
using Base.Threads
# using Distributed

export compute_weights, compute_DI, compute_FN, remove_duplicate_seqs

## function borrowed from stdlib/Distributed test file
function get_num_blas_threads()::Int
    blas = LinearAlgebra.BLAS.vendor()
    # Wrap in a try to catch unsupported blas versions
    try
        if blas == :openblas
            return ccall((:openblas_get_num_threads, Base.libblas_name), Cint, ())
        elseif blas == :openblas64
            return ccall((:openblas_get_num_threads64_, Base.libblas_name), Cint, ())
        elseif blas == :mkl
            return ccall((:MKL_Get_Max_Num_Threads, Base.libblas_name), Cint, ())
        end

        # OSX BLAS looks at an environment variable
        if Sys.isapple()
            return tryparse(Cint, get(ENV, "VECLIB_MAXIMUM_THREADS", "1"))
        end
    catch
    end

    return Int(get(ENV, "OMP_NUM_THREADS", Sys.CPU_THREADS))
end

# create a single-threaded-BLAS scope with do..end
function disable_blas_threads(f)
    old_num_threads = get_num_blas_threads()
    try
        LinearAlgebra.BLAS.set_num_threads(1)
        return f()
    finally
        LinearAlgebra.BLAS.set_num_threads(old_num_threads)
    end
end

const PACKBITS = 64
const _u = 0x0084210842108421
const _alt = 0x007c1f07c1f07c1f

const packfactor = PACKBITS ÷ 5
const packrest = PACKBITS % 5

const _msk = ~(UInt64(0));

clength(l::Int) = (l-1) ÷ packfactor + 1
crest(l::Int)   = (l-1) % packfactor + 1

nnz_aux(x::UInt64) = ((x | (x >>> 1) | (x >>> 2) | (x >>> 3) | (x >>> 4)) & _u)

nz_aux(nx::UInt64) = (nx & (nx >>> 1) & (nx >>> 2) & (nx >>> 3) & (nx >>> 4) & _u)

nz_aux2(nx::UInt64, s) = nz_aux(nx) & (_msk >>> s)

@inline function collapse(z)
    t = (z & _alt) + ((z >> 5) & _alt)
    t += t >> 10
    t += t >> 20
    t += t >> 40
    t &= 0x3ff
    return t
end

function compress_Z(Z::Matrix{Int8})
    N, M = size(Z)
    ZZ = Vector{Int8}[Z[:,i] for i = 1:M]

    cl = clength(N)

    cZ = [zeros(UInt64, cl) for i=1:M]

    @inbounds for i = 1:M
        cZi = cZ[i]
        ZZi = ZZ[i]
        for k = 1:N
            k0 = (k-1) ÷ packfactor + 1
            k1 = (k-1) % packfactor
            cZi[k0] |= (UInt64(ZZi[k]) << (5 * k1))
        end
    end

    return cZ
end

macro hash(x)
    tst = get(ENV, "GDCA_TESTING", "false") == "true"
    return Expr(:call, tst ? :sum : :hash, esc(x))
end

function remove_duplicate_seqs(Z::Matrix{Int8})
    N, M = size(Z)
    hZ = Array{UInt}(undef, M)
    @inbounds for i = 1:M
        hZ[i] = @hash(@view Z[:,i])
    end
    print("removing duplicate sequences... ")

    ref_seq_ind = Array{Int}(undef, M)
    ref_seq = Dict{UInt,Int}()
    @inbounds for i = 1:M
        ref_seq_ind[i] = get!(ref_seq, hZ[i], i)
    end
    uniqueseqs = collect(values(ref_seq))

    # Check for collisions
    old_collided = trues(M)
    collided = falses(M)
    while true
        fill!(collided, false)
        @inbounds for i = 1:M
            k = ref_seq_ind[i]
            (!old_collided[i] || k == i) && continue
            collided[i] = @views Z[:,i] ≠ Z[:,k]
        end
        any(collided) || break

        # Collect index of first row for each collided hash
        empty!(ref_seq)
        @inbounds for i = 1:M
            collided[i] || continue
            ref_seq_ind[i] = get!(ref_seq, hZ[i], i)
        end
        for v in values(ref_seq)
            push!(uniqueseqs, v)
        end
        old_collided, collided = collided, old_collided
    end
    sort!(uniqueseqs)

    newM = length(uniqueseqs)
    newZ = Z[:,uniqueseqs]
    println("done: $M -> $newM")
    return newZ, uniqueseqs
end

compute_weights(Z::Matrix{Int8}, θ) = compute_weights(Z, maximum(Z), θ)

compute_weights(Z::Matrix{Int8}, q, θ) = throw(ArgumentError("θ must be either :auto or a single real value"))

function Z_to_cZ(Z, q)
    M = size(Z, 2)
    fast = q ≤ 32 && get(ENV, "GDCA_FORCE_FALLBACK", "false") ≠ "true"
    fast || println("GaussDCA: using slower fallbacks")
    return fast ? compress_Z(Z) : [Z[:,i]::Vector{Int8} for i = 1:M]
end

function compute_weights(Z::Matrix{Int8}, q, θ::Symbol)
    θ ≠ :auto && throw(ArgumentError("θ must be either :auto or a single real value"))

    N, M = size(Z)
    cZ = Z_to_cZ(Z, q)
    θ = compute_theta(cZ, N, M)
    return compute_weights(cZ, θ, N, M)
end

function compute_weights(Z::Matrix{Int8}, q, θ::Real)
    N, M = size(Z)
    cZ = Z_to_cZ(Z, q)
    return compute_weights(cZ, θ, N, M)
end

function compute_FN(mJ::Matrix{Float64}, N::Int, q::Integer)
    q = Int(q)
    s = q - 1

    mJij = Array{Float64}(undef, s, s)
    amJi = Array{Float64}(undef, s)
    amJj = Array{Float64}(undef, s)
    fs = Float64(s)
    fs2 = Float64(s^2)

    FN = zeros(N, N)

    for i = 1:N-1
        row0 = (i-1) * s
        for j = i+1:N
            col0 = (j-1) * s

            ## devectorize for speed and memory efficicency
            # mJij = mJ[row,col]
            # mK = mJij .- mean(mJij, 1) .- mean(mJij, 2) .+ mean(mJij)

            amJ = 0.0
            fill!(amJi, 0.0)
            fill!(amJj, 0.0)
            for b = 1:s, a = 1:s
                x = mJ[row0 + a, col0 + b]
                mJij[a,b] = x
                amJi[b] += x / fs
                amJj[a] += x / fs
                amJ += x / fs2
            end
            fn = 0.0
            for b = 1:s, a = 1:s
                fn += (mJij[a,b] - amJi[b] - amJj[a] + amJ)^2
            end

            FN[i,j] = √fn
            FN[j,i] = FN[i,j]
        end
    end
    return FN
end

## parallel code

const TriuInd = Tuple{Tuple{Int,Int},Tuple{Int,Int},Int}

function ptriu(sz::Int, RT::Type, func::Function, args...)
    tot_inds = (sz * (sz-1)) ÷ 2
    # nw = nworkers() # use this for multi-processing
    nw = nthreads()

    if tot_inds ≥ nw
        inds_dist = diff([round(Int,x) for x in range(1, tot_inds+1, length=nw+1)])
    else
        inds_dist = [ones(Int, tot_inds); zeros(Int, nw-tot_inds)]
    end

    inds = Array{TriuInd}(undef, nw)
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

    ret = Array{RT}(undef, nw)

    disable_blas_threads() do
        ## use this for multi-processing
        # wrk = workers()
        # @sync for p = 1:nw
        #     @async ret[p] = remotecall_fetch(func, wrk[p], inds[p], args...)
        # end

        Threads.@threads for p = 1:nw
            ret[p] = func(inds[p], args...)
        end
    end

    return ret, inds
end

function ptriu_compose(src::Vector{Vector{T}}, sz::Int, inds::Vector{TriuInd}) where {T}
    dest = Array{T}(undef, sz, sz)

    @assert length(src) == length(inds)

    for p = 1:length(src)
        i0, j0 = inds[p][1]
        i1, j1 = inds[p][2]
        srcp = src[p]
        l = 1
        for i = i0:i1
            jj0 = i==i0 ? j0 : i+1
            jj1 = i==i1 ? j1 : sz
            for j = jj0:jj1
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
    cl = clength(N)
    cr = 5 * (packfactor - crest(N)) + packrest

    kmax = (cl - 1) ÷ 31
    rmax = (cl - 1) % 31

    meanfracid = 0.0
    i0, j0 = inds[1]
    i1, j1 = inds[2]
    @inbounds for i = i0:i1
        cZi = cZ[i]
        nids::UInt64 = 0
        jj0 = i==i0 ? j0 : i+1
        jj1 = i==i1 ? j1 : M
        for j = jj0:jj1
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
                nids += t
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
            nids += t
        end
        fracid = nids / N
        meanfracid += fracid
    end
    return meanfracid
end

function compute_theta(cZ::Vector{Vector{T}}, N::Int, M::Int) where {T<:Union{Int8,UInt64}}
    chunk_means, _ = ptriu(M, Float64, compute_theta_chunk, cZ, N, M)
    meanfracid = sum(chunk_means) / (0.5 * M * (M-1))
    θ = min(0.5, 0.38 * 0.32 / meanfracid)
    return θ
end

# slow fallback
function compute_theta_chunk(inds::TriuInd, ZZ::Vector{Vector{Int8}}, N::Int, M::Int)
    meanfracid = 0.0
    i0, j0 = inds[1]
    i1, j1 = inds[2]
    @inbounds for i = i0:i1
        Zi = ZZ[i]
        jj0 = i==i0 ? j0 : i+1
        jj1 = i==i1 ? j1 : M
        for j = jj0:jj1
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
    cl = clength(N)
    kmax = (cl - 1) ÷ 31
    rmax = (cl - 1) % 31 + 1

    W = zeros(M)
    i0, j0 = inds[1]
    i1, j1 = inds[2]

    @inbounds for i = i0:i1
        cZi = cZ[i]
        jj0 = i==i0 ? j0 : i+1
        jj1 = i==i1 ? j1 : M
        for j = jj0:jj1
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
    return W
end

function compute_weights(cZ::Vector{Vector{T}}, θ::Real, N::Int, M::Int) where {T<:Union{Int8,UInt64}}
    θ = Float64(θ)

    Meff = 0.0

    thresh = floor(θ * N)
    println("θ = $θ threshold = $thresh")

    if θ == 0
        println("M = $M N = $N Meff = $M")
        W = ones(M)
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
    @inbounds for i = i0:i1
        Zi = ZZ[i]
        jj0 = i==i0 ? j0 : i+1
        jj1 = i==i1 ? j1 : M
        for j = jj0:jj1
            Zj = ZZ[j]
            dist = 0
            k = 1
            while dist < thresh && k <= N
                dist += Zi[k] ≠ Zj[k]
                k += 1
            end
            if dist < thresh
                W[i] += 1
                W[j] += 1
            end
        end
    end
    return W
end

function compute_dists_chunk(inds::TriuInd, cZ::Vector{Vector{UInt64}}, N::Int, M::Int)
    cl = clength(N)
    kmax = (cl - 1) ÷ 31
    rmax = (cl - 1) % 31 + 1

    D = zeros(Float16, inds[3])

    i0, j0 = inds[1]
    i1, j1 = inds[2]

    l = 1
    for i = i0:i1
        # cZi = unsafe(cZ[i])
        cZi = cZ[i]
        jj0 = i==i0 ? j0 : i+1
        jj1 = i==i1 ? j1 : M
        for j = jj0:jj1
            # cZj = unsafe(cZ[j])
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
            D[l] = dist / N
            l += 1
        end
    end

    return D
end

function compute_dists(cZ::Vector{Vector{UInt64}}, N::Int, M::Int)
    Ds, inds = ptriu(M, Vector{Float16}, compute_dists_chunk, cZ, N, M)
    D = ptriu_compose(Ds, M, inds)
    return D
end

const KT = Symmetric{Float64,Matrix{Float64}}

function compute_DI_chunk(inds::TriuInd, N::Int, s::Integer, iKs::Vector{KT}, mJ::Matrix{Float64}, z::Float64)
    DI = Array{Float64}(undef, inds[3])
    i0, j0 = inds[1]
    i1, j1 = inds[2]
    l = 1
    @inbounds for i = i0:i1
        row = ((i-1)*s+1):i*s
        invsqrtKi = iKs[i]
        jj0 = i==i0 ? j0 : i+1
        jj1 = i==i1 ? j1 : N
        for j = jj0:jj1
            col = ((j-1)*s+1):j*s
            invsqrtKj = iKs[j]
            mJij = mJ[row, col]
            if sum(mJij) ≠ 0
                MM = invsqrtKi * mJij * invsqrtKj
                V = MM * MM'
                eigV = eigvals(V)
                DI[l] = z + 0.5 * sum(@. log1p(√(1 + 4 * eigV)))
            else
                DI[l] = 0
            end
            l += 1
        end
    end
    return DI
end

function compute_DI(mJ::Matrix{Float64}, C::Matrix{Float64}, N::Int, q::Integer)
    s = q - 1
    iKs = Array{KT}(undef, N)
    rowi = 0
    for i = 1:N
        row = rowi .+ (1:s)
        rowi += s
        iKs[i] = (√(Symmetric(@view C[row,row])))::KT
    end

    z = 0.5 * s * log(0.5)

    DIs, inds = ptriu(N, Vector{Float64}, compute_DI_chunk, N, s, iKs, mJ, z)
    DI = ptriu_compose(DIs, N, inds)

    return DI
end

end # module
