module GaussDCA

export gDCA, printrank

using LinearAlgebra, Printf

include("DCAUtils.jl")
using .DCAUtils

function gDCA(
        filename::AbstractString;
        pseudocount::Real = 0.8,
        θ = :auto,
        max_gap_fraction::Real = 0.9,
        score::Symbol = :frob,
        min_separation::Integer = 5,
        remove_dups::Bool = false
    )

    check_arguments(filename, pseudocount, θ, max_gap_fraction, score, min_separation)

    Z = read_fasta_alignment(filename, max_gap_fraction)
    if remove_dups
        Z, _ = remove_duplicate_sequences(Z)
    end
    N, M = size(Z)
    q = Int(maximum(Z))
    q > 32 && error("parameter q=$q is too big (max 32 is allowed)")

    Pi_true, Pij_true, Meff, _ = compute_weighted_frequencies(Z, q, θ)

    Pi, Pij = add_pseudocount(Pi_true, Pij_true, Float64(pseudocount), N, q)

    C = compute_C(Pi, Pij)

    mJ = inv(cholesky(C))

    if score == :DI
        S = compute_DI(mJ, C, N, q)
    else
        S = compute_FN(mJ, N, q)
    end

    S = correct_APC(S)

    R = compute_ranking(S, min_separation)

    return R
end

function check_arguments(filename, pseudocount, θ, max_gap_fraction, score, min_separation)
    aerror(s) = throw(ArgumentError(s))
    0 <= pseudocount <= 1 ||
        aerror("invalid pseudocount value: $pseudocount (must be between 0 and 1)")
    θ == :auto || (θ isa Real && 0 <= θ <= 1) ||
        aerror("invalid θ value: $θ (must be either :auto, or a number between 0 and 1)")
    0 <= max_gap_fraction <= 1 ||
        aerror("invalid max_gap_fraction value: $max_gap_fraction (must be between 0 and 1)")
    score in [:DI, :frob] ||
        aerror("invalid score value: $score (must be either :DI or :frob)")
    min_separation >= 1 ||
        aerror("invalid min_separation value: $min_separation (must be >= 1)")
    isfile(filename) ||
        aerror("cannot open file $filename")

    return true
end

function printrank(io::IO, R::Vector{Tuple{Int,Int,Float64}})
    for I in R
        @printf(io, "%i %i %e\n", I[1], I[2], I[3])
    end
end
printrank(R::Vector{Tuple{Int,Int,Float64}}) = printrank(STDOUT, R)

printrank(outfile::AbstractString, R::Vector{Tuple{Int,Int,Float64}}) = open(f->printrank(f, R), outfile, "w")

compute_C(Pi::Vector{Float64}, Pij::Matrix{Float64}) = Pij - Pi * Pi'

function correct_APC(S::Matrix)
    N = size(S, 1)
    Si = sum(S, dims=1)
    Sj = sum(S, dims=2)
    Sa = sum(S) * (1 - 1/N)

    S -= (Sj * Si) / Sa
    return S
end

function compute_ranking(S::Matrix{Float64}, min_separation::Int = 5)
    N = size(S, 1)
    R = Array{Tuple{Int,Int,Float64}}(undef, ((N-min_separation)*(N-min_separation+1)) ÷ 2)
    counter = 0
    for i = 1:N-min_separation, j = i+min_separation:N
        counter += 1
        R[counter] = (i, j, S[j,i])
    end

    sort!(R, by=x->x[3], rev=true)
    return R
end

end # module
