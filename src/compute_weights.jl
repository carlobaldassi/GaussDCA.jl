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

compute_weights(Z::Matrix{Int8}, θ) = compute_weights(Z, maximum(Z), θ)

compute_weights(Z::Matrix{Int8}, q, θ) = throw(ArgumentError("θ must be either :auto or a single real value"))

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

