"""
    compute_FN(J::Matrix{Float64}, q::Integer = 21) -> Matrix{Float64}

Compute the matrix of Frobenius norms.

`J` is the inverse of the covariance matrix \$C_{ij} = P_{ij} - P_i P_j\$.

The integer `q` is used to specify the block size of the matrices (each block has size
\$(q-1)×(q-1)\$). The result has one entry per block.

This function can use multiple threads if available.
"""
function compute_FN(J::Matrix{Float64}, q::Integer = 21)
    Nq = size(J, 1)
    size(J, 2) == Nq || throw(ArgumentError("J matrix must be square"))
    s = q - 1
    Nq % s == 0 || throw(ArgumentError("incompatible size of J with q"))

    N = Nq ÷ s
    Jij = Array{Float64}(undef, s, s)
    aJi = Array{Float64}(undef, s)
    aJj = Array{Float64}(undef, s)
    fs = Float64(s)
    fs2 = Float64(s^2)

    FN = zeros(N, N)

    for i = 1:N-1
        row0 = (i-1) * s
        for j = i+1:N
            col0 = (j-1) * s

            ## devectorize for speed and memory efficicency
            # Jij = J[row,col]
            # mK = Jij .- mean(Jij, 1) .- mean(Jij, 2) .+ mean(Jij)

            aJ = 0.0
            fill!(aJi, 0.0)
            fill!(aJj, 0.0)
            for b = 1:s, a = 1:s
                x = J[row0 + a, col0 + b]
                Jij[a,b] = x
                aJi[b] += x / fs
                aJj[a] += x / fs
                aJ += x / fs2
            end
            fn = 0.0
            for b = 1:s, a = 1:s
                fn += (Jij[a,b] - aJi[b] - aJj[a] + aJ)^2
            end

            FN[i,j] = √fn
            FN[j,i] = FN[i,j]
        end
    end
    return FN
end
