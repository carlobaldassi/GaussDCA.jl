using LinearAlgebra

const KT = Symmetric{Float64,Matrix{Float64}}

function compute_DI_chunk(inds::TriuInd, N::Int, s::Integer, iKs::Vector{KT}, mJ::Matrix{Float64}, z::Float64)
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
            mJij = mJ[row, col]
            if sum(mJij) ≠ 0
                MM = invsqrtKi * mJij * invsqrtKj
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

function compute_DI(mJ::Matrix{Float64}, C::Matrix{Float64}, N::Int, q::Integer)
    s = q - 1
    iKs = Array{KT}(undef, N)
    rowi = 0
    for i = 1:N
        row = rowi .+ (1:s)
        rowi += s
        iKs[i] = (√(Symmetric(@view C[row,row])))::KT
    end

    z = 0.5 * s * log(0.5)

    DIs, inds = ptriu(N, Vector{Float64}, compute_DI_chunk, N, s, iKs, mJ, z)
    DI = ptriu_compose(DIs, N, inds)

    return DI
end

