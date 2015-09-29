module GaussDCA

export gDCA, printrank

include("read_fasta_alignment.jl")

using Compat
using .ReadFastaAlignment

if nprocs() > 2 && parse(get(ENV, "PARALLEL_GDCA", "true"))
    include("parallel.jl")
else
    include("nonparallel.jl")
end
using .AuxFunctions

if VERSION < v"0.3-"
    Base.sqrtm{T<:Real}(A::Symmetric{T}) = sqrtm(A, false)
end

@compat function gDCA(filename::AbstractString;
                      pseudocount::Real = 0.8,
                      theta = :auto,
                      max_gap_fraction::Real = 0.9,
                      score::Symbol = :frob,
                      min_separation::Integer = 5,
                      remove_dups::Bool = false)


    check_arguments(filename, pseudocount, theta, max_gap_fraction, score, min_separation)

    use_threading(true)

    Z = read_fasta_alignment(filename, max_gap_fraction)
    if remove_dups
        Z, _ = remove_duplicate_seqs(Z)
    end
    N, M = size(Z)
    q = Int(maximum(Z))
    q > 32 && error("parameter q=$q is too big (max 32 is allowed)")

    Pi_true, Pij_true, Meff, _ = compute_new_frequencies(Z, q, theta)

    Pi, Pij = add_pseudocount(Pi_true, Pij_true, Float64(pseudocount), N, q)

    C = compute_C(Pi, Pij)

    mJ = inv(cholfact(C))

    if score == :DI
        S = compute_DI(mJ, C, N, q)
    else
        S = compute_FN(mJ, N, q)
    end

    S = correct_APC(S)

    R = compute_ranking(S, min_separation)

    use_threading(false)

    return R
end

function check_arguments(filename, pseudocount, theta, max_gap_fraction, score, min_separation)

    aerror(s) = throw(ArgumentError(s))
    0 <= pseudocount <= 1 ||
        aerror("invalid pseudocount value: $pseudocount (must be between 0 and 1)")
    theta == :auto || (isa(theta, Real) && 0 <= theta <= 1) ||
        aerror("invalid theta value: $theta (must be either :auto, or a number between 0 and 1)")
    0 <= max_gap_fraction <= 1 ||
        aerror("invalid max_gap_fraction value: $max_gap_fraction (must be between 0 and 1)")
    score in [:DI, :frob] ||
        aerror("invalid score value: $score (must be either :DI or :frob)")
    min_separation >= 1 ||
        aerror("invalid min_separation value: $min_separation (must be >= 1)")
    isreadable(filename) ||
        aerror("cannot open file $filename")

    return true
end

@compat function printrank(io::IO, R::Vector{Tuple{Int,Int,Float64}})
    for I in R
        @printf(io, "%i %i %e\n", I[1], I[2], I[3])
    end
end
@compat printrank(R::Vector{Tuple{Int,Int,Float64}}) = printrank(STDOUT, R)

@compat printrank(outfile::AbstractString, R::Vector{Tuple{Int,Int,Float64}}) = open(f->printrank(f, R), outfile, "w")

function compute_new_frequencies(Z::Matrix{Int8}, q, theta)

    W, Meff = compute_weights(Z, q, theta)
    Pi_true, Pij_true = compute_freqs(Z, W, Meff)

    return Pi_true, Pij_true, Meff, W
end


function compute_freqs(Z::Matrix{Int8}, W::Vector{Float64}, Meff::Float64)
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
	xr = i0 + (1:s)
	Pij[xr, xr] = (1 - pc) * Pij_true[xr, xr]
        for alpha = 1:s
            x = i0 + alpha
            Pij[x, x] += pcq
	end
        i0 += s
    end

    return Pi, Pij
end

compute_C(Pi::Vector{Float64}, Pij::Matrix{Float64}) = Pij - Pi * Pi'

function correct_APC(S::Matrix)
    N = size(S, 1)
    Si = sum(S, 1)
    Sj = sum(S, 2)
    Sa = sum(S) * (1 - 1/N)

    S -= (Sj * Si) / Sa
    return S
end

@compat function compute_ranking(S::Matrix{Float64}, min_separation::Int = 5)

    N = size(S, 1)
    R = Array(Tuple{Int,Int,Float64}, div((N-min_separation)*(N-min_separation+1), 2))
    counter = 0
    for i = 1:N-min_separation, j = i+min_separation:N
        counter += 1
        R[counter] = (i, j, S[j,i])
    end

    sort!(R, by=x->x[3], rev=true)
    return R

end

function _warmup()
    myid() == 1 || return
    while any(map(w->!(remotecall_fetch(w, isdefined, :_GaussDCAloaded)), workers()))
        sleep(0.1)
    end
    smallfile = joinpath(dirname(Base.source_path()), "..", "data", "small.fasta.gz")
    ostdout = STDOUT
    try
        (rd, wr) = redirect_stdout()
        gDCA(smallfile)
        gDCA(smallfile, score=:DI)
    finally
        redirect_stdout(ostdout)
    end
end

end

module _GaussDCAloaded
end

VERSION < v"0.4-" && GaussDCA._warmup() # TODO: doesn't work properly on 0.4, should investigate (or rather use precompilation)
