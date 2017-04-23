include("unsafe_array.jl")
using .Unsafe
using Compat

const PACKBITS = 64
const _u = 0x0084210842108421
const _alt = 0x007c1f07c1f07c1f

const packfactor = div(PACKBITS, 5)
const packrest = PACKBITS % 5

const _msk = ~(UInt64(0));

clength(l::Int) = div(l-1, packfactor) + 1
crest(l::Int)   = (l-1) % packfactor + 1

nnz_aux(x::UInt64) = ((x | (x >>> 1) | (x >>> 2) | (x >>> 3) | (x >>> 4)) & _u)

nz_aux(nx::UInt64) = (nx & (nx >>> 1) & (nx >>> 2) & (nx >>> 3) & (nx >>> 4) & _u)

nz_aux2(nx::UInt64, s) = nz_aux(nx) & (_msk >>> s)

macro collapse(z)
    quote
        t = ($(esc(z)) & _alt) + (($(esc(z)) >> 5) & _alt)
        t += t >> 10
        t += t >> 20
        t += t >> 40
        t &= 0x3ff
    end
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
            k0 = div(k-1, packfactor)+1
            k1 = (k-1) % packfactor
            cZi[k0] |= (UInt64(ZZi[k]) << 5*k1)
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
    hZ = Array{UInt}(M)
    @inbounds for i = 1:M
        hZ[i] = @hash(Z[:,i])
    end
    print("removing duplicate sequences... ")

    ref_seq_ind = Array{Int}(M)
    ref_seq = Dict{UInt,Int}()
    @inbounds for i = 1:M
        ref_seq_ind[i] = get!(ref_seq, hZ[i], i)
    end
    uniqueseqs = collect(values(ref_seq))

    # Check for collisions
    collided = falses(M)
    @inbounds for i = 1:M
        k = ref_seq_ind[i]
        k == i && continue
        for j = 1:N
            Z[j,i] != Z[j,k] && (collided[i] = true; break)
        end
    end

    if any(collided)
        nowcollided = BitArray(M)
        while any(collided)
            # Collect index of first row for each collided hash
            empty!(ref_seq)
            @inbounds for i = 1:M
                collided[i] || continue
                ref_seq_ind[i] = get!(ref_seq, hZ[i], i)
            end
            for v in values(ref_seq)
                push!(uniqueseqs, v)
            end

            # Check for collisions
            fill!(nowcollided, false)
            @inbounds for i = 1:M
                k = ref_seq_ind[i]
                (!collided[i] || k == i) && continue
                for j = 1:N
                    Z[j,i] != Z[j,k] && (nowcollided[i] = true; break)
                end
            end
            collided, nowcollided = nowcollided, collided
        end
    end

    sort!(uniqueseqs)

    newM = length(uniqueseqs)
    newZ = Z[:,uniqueseqs]
    println("done: $M -> $newM")
    return newZ, uniqueseqs
end

compute_weights(Z::Matrix{Int8}, theta) = compute_weights(Z, maximum(Z), theta)

compute_weights(Z::Matrix{Int8}, q, theta) = error("theta must be either :auto or a single real value")

function Z_to_cZ(Z, q)
    M = size(Z, 2)
    fast = q ≤ 32 && get(ENV, "GDCA_FORCE_FALLBACK", "false") ≠ "true"
    fast || println("GaussDCA: using slower fallbacks")
    return fast ? compress_Z(Z) : [Z[:,i]::Vector{Int8} for i = 1:M]
end

function compute_weights(Z::Matrix{Int8}, q, theta::Symbol)
    theta != :auto && return invoke(compute_weights, (Matrix{Int8}, Any, Any), Z, theta)

    N, M = size(Z)
    cZ = Z_to_cZ(Z, q)
    theta = compute_theta(cZ, N, M)
    return compute_weights(cZ, theta, N, M)
end

function compute_weights(Z::Matrix{Int8}, q, theta::Real)
    N, M = size(Z)
    cZ = Z_to_cZ(Z, q)
    return compute_weights(cZ, theta, N, M)
end

#function matsqrt(K::Matrix{Float64})
    #vL, U = eig(K)
    #L = diagm(sqrt(vL))
    #U * L * U'
#end

function compute_FN(mJ::Matrix{Float64}, N::Int, q::Integer)

    q = Int(q)
    s = q - 1

    mJij = Array{Float64}(s, s)
    amJi = Array{Float64}(s)
    amJj = Array{Float64}(s)
    fs = Float64(s)
    fs2 = Float64(s^2)

    FN = zeros(N, N)

    for i = 1:N-1
        #row = (i-1)*s + 1 : i*s
        row0 = (i-1)*s

        for j = i+1:N
            #col = (j-1)*s + 1 : j*s
            col0 = (j-1)*s

            # devectorize for speed and memory efficicency
            #mJij = mJ[row,col]
            #mK = mJij .- mean(mJij, 1) .- mean(mJij, 2) .+ mean(mJij)

            amJ = 0.0
            for a = 1:s
                amJi[a] = 0.0
                amJj[a] = 0.0
            end
            for b = 1:s, a = 1:s
                x = mJ[row0 + a, col0 + b]
                mJij[a,b] = x
                amJi[b] += x / fs
                amJj[a] += x / fs
                amJ += x / fs2
            end
            fn = 0.0
            for b = 1:s, a = 1:s
                fn += (mJij[a,b] - amJi[b] - amJj[a] + amJ) ^ 2
            end

            FN[i,j] = sqrt(fn)
            #FN[i,j] = normfro(mK)
            FN[j,i] = FN[i,j]
        end
    end
    return FN
end
