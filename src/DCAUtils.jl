module DCAUtils

export read_fasta_alignment,
       remove_duplicate_sequences,
       compute_theta,
       compute_weights,
       compute_dists,
       compute_DI,
       compute_FN,
       compute_weighted_frequencies,
       add_pseudocount

include("read_fasta_alignment.jl")
using .ReadFastaAlignment

include("compress_sequencies.jl")
using .CompressSequencies

include("parallel_triu.jl")
using .ParallelTriU

macro hash(x)
    tst = get(ENV, "GDCA_TESTING", "false") == "true"
    return Expr(:call, tst ? :sum : :hash, esc(x))
end

function remove_duplicate_sequences(Z::Matrix{Int8})
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

function compute_weighted_frequencies(Z::Matrix{Int8}, q, θ)
    W, Meff = compute_weights(Z, q, θ)
    Pi_true, Pij_true = compute_weighted_frequencies(Z, W, Meff)
    return Pi_true, Pij_true, Meff, W
end

function compute_weighted_frequencies(Z::Matrix{Int8}, W::Vector{Float64}, Meff::Float64)
    N, M = size(Z)
    q = maximum(Z)
    s = q - 1

    Ns = N * s

    Pij = zeros(Ns, Ns)
    Pi = zeros(Ns)

    ZZ = Vector{Int8}[vec(Z[i,:]) for i = 1:N]

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
    Pi /= Meff

    i0 = 0
    for i = 1:N
        Zi = ZZ[i]
        j0 = i0
        for j = i:N
            Zj = ZZ[j]
            for k = 1:M
                a = Zi[k]
                b = Zj[k]
                (a == q || b == q) && continue
                Pij[i0+a, j0+b] += W[k]
            end
            j0 += s
        end
        i0 += s
    end
    for i = 1:Ns
        Pij[i,i] /= Meff
        for j = i+1:Ns
            Pij[i,j] /= Meff
            Pij[j,i] = Pij[i,j]
        end
    end

    return Pi, Pij
end

function add_pseudocount(Pi_true::Vector{Float64}, Pij_true::Matrix{Float64}, pc::Float64, N::Int, q::Int)
    pcq = pc / q

    Pij = (1 - pc) * Pij_true .+ pcq / q
    Pi = (1 - pc) * Pi_true .+ pcq

    s = q - 1

    i0 = 0
    for i = 1:N
        xr = i0 .+ (1:s)
        Pij[xr, xr] = (1 - pc) * Pij_true[xr, xr]
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
include("compute_DI.jl")
include("compute_FN.jl")

end # module
