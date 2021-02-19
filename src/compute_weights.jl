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

"""
    compute_weights(Z::Matrix{Int8}, [q::Integer,] θ; verbose::Bool = true) -> (Vector{Float64}, Float64)

This function computes the reweighting vector. It retuns the vector and its sum `Meff` (the latter
represents the number of "effective" sequences).

`Z` is an \$N × M\$ multiple sequence aLignment (see [`read_fasta_alignment`](@ref)). \$N\$
is the length of each sequence and \$M\$ the number of sequences.

`q` is the maximum value in the alphabet, if omitted it's computed from `maximum(Z)`.

`θ` is the distance threshold: for any sequence, the number \$n\$ of sequences (including itself)
that are at normalized distance smaller than \$⌊θN⌋\$ is counted, and the weight of that sequence
is then \$1/n\$.

`θ` can be a real value between 0 and 1, or the symbol `:auto`, in which case the
[`compute_theta`](@ref) function is used.

This function can use multiple threads if available.
"""
function compute_weights end

function compute_weights(cZ::Vector{<:Vector{<:Union{Int8,UInt64}}}, θ::Real, N::Int, M::Int; verbose::Bool = true)
    θ = Float64(θ)

    Meff = 0.0

    thresh = floor(θ * N)
    verbose && println("θ = $θ threshold = $thresh")

    if θ == 0
        verbose && println("M = $M N = $N Meff = $M")
        W = ones(M)
        return W, Float64(M)
    end

    Ws, _ = ptriu(M, Vector{Float64}, compute_weights_chunk, cZ, thresh, N, M)
    W = sum(Ws)

    for i = 1:M
        W[i] = 1 / (1 + W[i])
        Meff += W[i]
    end
    verbose && println("M = $M N = $N Meff = $Meff")
    return W, Meff
end

compute_weights(Z::Matrix{Int8}, θ; kw...) = compute_weights(Z, maximum(Z), θ; kw...)

compute_weights(Z::Matrix{Int8}, q, θ; kw...) = throw(ArgumentError("θ must be either :auto or a single real value"))

function compute_weights(Z::Matrix{Int8}, q::Integer, θ::Symbol; verbose::Bool = true)
    θ ≠ :auto && throw(ArgumentError("θ must be either :auto or a single real value"))

    N, M = size(Z)
    cZ = Z_to_cZ(Z, q)
    θ = compute_theta(cZ, N, M)
    return compute_weights(cZ, θ, N, M; verbose)
end

function compute_weights(Z::Matrix{Int8}, q::Integer, θ::Real; verbose::Bool = true)
    N, M = size(Z)
    cZ = Z_to_cZ(Z, q)
    return compute_weights(cZ, θ, N, M; verbose)
end

