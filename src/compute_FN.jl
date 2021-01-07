function compute_FN(mJ::Matrix{Float64}, N::Int, q::Integer)
    q = Int(q)
    s = q - 1

    mJij = Array{Float64}(undef, s, s)
    amJi = Array{Float64}(undef, s)
    amJj = Array{Float64}(undef, s)
    fs = Float64(s)
    fs2 = Float64(s^2)

    FN = zeros(N, N)

    for i = 1:N-1
        row0 = (i-1) * s
        for j = i+1:N
            col0 = (j-1) * s

            ## devectorize for speed and memory efficicency
            # mJij = mJ[row,col]
            # mK = mJij .- mean(mJij, 1) .- mean(mJij, 2) .+ mean(mJij)

            amJ = 0.0
            fill!(amJi, 0.0)
            fill!(amJj, 0.0)
            for b = 1:s, a = 1:s
                x = mJ[row0 + a, col0 + b]
                mJij[a,b] = x
                amJi[b] += x / fs
                amJj[a] += x / fs
                amJ += x / fs2
            end
            fn = 0.0
            for b = 1:s, a = 1:s
                fn += (mJij[a,b] - amJi[b] - amJj[a] + amJ)^2
            end

            FN[i,j] = √fn
            FN[j,i] = FN[i,j]
        end
    end
    return FN
end
