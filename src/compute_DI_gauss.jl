using LinearAlgebra

const KT = Symmetric{Float64,Matrix{Float64}}

function compute_DI_gauss_chunk(inds::TriuInd, N::Int, s::Integer, iKs::Vector{KT}, J::Matrix{Float64}, z::Float64)
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
            Jij = J[row, col]
            if sum(Jij) ≠ 0
                MM = invsqrtKi * Jij * invsqrtKj
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

"""
    compute_DI_gauss(J::Matrix{Float64}, C::Matrix{Float64}, q::Integer = 21) -> Matrix{Float64}

Compute the Direct Information matrix, assuming a Gaussian model.

`C` is the covariance matrix \$C_{ij} = P_{ij} - P_i P_j\$, and `J` its inverse.

The integer `q` is used to specify the block size of the matrices (each block has size
\$(q-1)×(q-1)\$). The result has one entry per block.

This function can use multiple threads if available.
"""
function compute_DI_gauss(J::Matrix{Float64}, C::Matrix{Float64}, q::Integer = 21)
    Nq = size(J, 1)
    size(J, 2) == Nq || throw(ArgumentError("J matrix must be square"))
    size(C) == size(J) || throw(ArgumentError("incompatible J and C matrices"))
    s = q - 1
    Nq % s == 0 || throw(ArgumentError("incompatible size of J with q"))

    N = Nq ÷ s
    iKs = Array{KT}(undef, N)
    rowi = 0
    for i = 1:N
        row = rowi .+ (1:s)
        rowi += s
        iKs[i] = (√(Symmetric(@view C[row,row])))::KT
    end

    z = 0.5 * s * log(0.5)

    DIs, inds = ptriu(N, Vector{Float64}, compute_DI_gauss_chunk, N, s, iKs, J, z)
    DI = ptriu_compose(DIs, N, inds)

    return DI
end
