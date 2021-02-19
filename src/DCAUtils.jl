module DCAUtils

export read_fasta_alignment,
       remove_duplicate_sequences,
       compute_theta,
       compute_weights,
       compute_dists,
       compute_weighted_frequencies,
       add_pseudocount,
       compute_DI_gauss,
       compute_FN

include("read_fasta_alignment.jl")
using .ReadFastaAlignment

include("compress_sequencies.jl")
using .CompressSequencies

include("parallel_triu.jl")
using .ParallelTriU

macro hash(x)
    tst = get(ENV, "DCAUTILS_TESTING", "false") == "true"
    return Expr(:call, tst ? :sum : :hash, esc(x))
end

"""
    remove_duplicate_sequences(Z::Matrix{Int8}; verbose::Bool = true) -> (Matrix{Int8}, Vector{Int})

Takes a matrix representing a mutiple sequence alignment (see [`read_fasta_alignment`](@ref))
and returns a new matrix with all duplicated sequences removed. It also returns a vector of column
indices with the positions of the unique sequences in the input matrix.
"""
function remove_duplicate_sequences(Z::Matrix{Int8}; verbose::Bool = true)
    N, M = size(Z)
    hZ = Array{UInt}(undef, M)
    @inbounds for i = 1:M
        hZ[i] = @hash(@view Z[:,i])
    end
    verbose && print("removing duplicate sequences... ")

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
    verbose && println("done: $M -> $newM")
    return newZ, uniqueseqs
end

"""
    add_pseudocount(Pi::Vector{Float64}, Pij::Matrix{Float64}, pc::Float64, q::Integer = 21) -> (Vector{Float64}, Matrix{Float64})

This function takes one- and two-points frequencies (see [`compute_weighted_frequencies`](@ref))
and returns the corresponding frequencies with a pseudocount `pc` added.

The resulting frequencies are the same that would be obtained by mixing the original data with
weight `(1-pc)` with a uniform distribution with weight `pc`. So `pc` must be between 0 (returns a
copy of the original data) and 1 (returns the frequencies for the uniform distribution).

The integer `q` is used to specify the block size of the matrices (each block has size
\$(q-1)×(q-1)\$).
"""
function add_pseudocount(Pi_true::Vector{Float64}, Pij_true::Matrix{Float64}, pc::Float64, q::Integer = 21)
    Nq = length(Pi_true)
    s = q - 1
    Nq % s == 0 || throw(ArgumentError("incompatible length of Pi_true with q"))
    size(Pij_true) == (Nq, Nq) || throw(ArgumentError("incompatible sizes of Pi_true and Pij_true"))
    N = Nq ÷ s

    pcq = pc / q

    Pij = @. (1 - pc) * Pij_true + pcq / q
    Pi = @. (1 - pc) * Pi_true + pcq

    i0 = 0
    for i = 1:N
        xr = i0 .+ (1:s)
        Pij[xr, xr] .= @. (1 - pc) * @view Pij_true[xr, xr]
        for α = 1:s
            x = i0 + α
            Pij[x, x] += pcq
        end
        i0 += s
    end

    return Pi, Pij
end

include("compute_theta.jl")
include("compute_weights.jl")
include("compute_dists.jl")
include("compute_weighted_frequencies.jl")
include("compute_DI_gauss.jl")
include("compute_FN.jl")

end # module
