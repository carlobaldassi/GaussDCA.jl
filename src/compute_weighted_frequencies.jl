using LinearAlgebra

function compute_Pij_chunk(inds::TriuInd, N::Int, s::Integer, ZZ::Vector{Vector{Int8}}, M::Int, W::Vector{Float64}, Meff::Float64)
    q = s + 1
    Pij = Array{Matrix{Float64}}(undef, inds[3])
    i0, j0 = inds[1]
    i1, j1 = inds[2]
    l = 1
    @inbounds for i = i0:i1
        Zi = ZZ[i]
        jj0 = i==i0 ? j0 : i+1
        jj1 = i==i1 ? j1 : N
        for j = jj0:jj1
            Zj = ZZ[j]
            block = zeros(s, s)
            for k = 1:M
                a, b = Zi[k], Zj[k]
                (a == q || b == q) && continue
                block[a, b] += W[k]
            end
            block ./= Meff
            Pij[l] = block
            l += 1
        end
    end
    return Pij
end

"""
    compute_weighted_frequencies(Z::Matrix{Int8}, W::Vector{Float64} [, q]) -> (Vector{Float64}, Matrix{Float64})

Given a multiple sequence alignment matrix `Z` (see [`read_fasta_alignment`](@ref)), and a vector
of weights (see [`compute_weights`](@ref)), returns the empirical one- and two-point frequencies
\$P_i\$ and \$P_{ij}\$.

`q` is the size of the alphabet, if not given it's computed as `maximum(Z)`.

If `Z` has size \$N × M\$ (i.e. \$M\$ sequences of length \$N\$), the resulting vector \$P_i\$ has
length \$N (q-1)\$ and contains \$N\$ blocks (one for each residue position), each block containing
the frequencies of the amino-acids, weighted according to `W`. The frequency of the last symbol,
which usually represents the gap, is omitted and can be recovered by normalization. The resulting
matrix \$P_{ij}\$ has size \$N (q-1) × N (q-1)\$ and it also has a block structure, with \$N × N\$
blocks, one for each pair of residues (the last row and column of each block are omitted and can be
recovered by normalization).

    compute_weighted_frequencies(Z::Matrix{Int8}, [q,] θ) -> (Vector{Float64}, Matrix{Float64}, Float64, Vector{Float64})

This form of the function just calls [`compute_weights`](@ref) with the given values of `θ` and `q`
and then uses the result to call the version desrcibed above.

Besides returning the one- and two-point frequencies, it also returns the result of
`compute_weights`: the `Meff` and the reweighting vector.
"""
function compute_weighted_frequencies end

compute_weighted_frequencies(Z::Matrix{Int8}, θ::Union{Real,Symbol}) = compute_weighted_frequencies(Z, maximum(Z), θ)

function compute_weighted_frequencies(Z::Matrix{Int8}, q::Integer, θ::Union{Real,Symbol})
    W, Meff = compute_weights(Z, q, θ)
    Pi_true, Pij_true = compute_weighted_frequencies(Z, W, Meff, q)
    return Pi_true, Pij_true, Meff, W
end

compute_weighted_frequencies(Z::Matrix{Int8}, W::Vector{Float64}) = compute_weighted_frequencies(Z, W, maximum(Z))
compute_weighted_frequencies(Z::Matrix{Int8}, W::Vector{Float64}, q::Integer) = compute_weighted_frequencies(Z, W, sum(W), q)
compute_weighted_frequencies(Z::Matrix{Int8}, W::Vector{Float64}, Meff::Float64) = compute_weighted_frequencies(Z, W, Meff, maximum(Z))

function compute_weighted_frequencies(Z::Matrix{Int8}, W::Vector{Float64}, Meff::Float64, q::Integer)
    N, M = size(Z)
    s = q - 1

    Ns = N * s

    Pij = zeros(Ns, Ns)
    Pi = zeros(Ns)

    ZZ = Vector{Int8}[Z[i,:] for i = 1:N]

    i0 = 0
    for i = 1:N
        Zi = ZZ[i]
        for k = 1:M
            a = Zi[k]
            a == q && continue
            Pi[i0 + a] += W[k]
        end
        i0 += s
    end
    Pi ./= Meff

    Pijs, inds = ptriu(N, Vector{Matrix{Float64}}, compute_Pij_chunk, N, s, ZZ, M, W, Meff)
    Pij = ptriu_compose_blocks(Pijs, N, s, inds)

    Pij[diagind(Pij)] .= Pi

    return Pi, Pij
end
