function compute_dists_chunk(inds::TriuInd, cZ::Vector{Vector{UInt64}}, N::Int, M::Int)
    cl = clength(N)
    kmax = (cl - 1) ÷ 31
    rmax = (cl - 1) % 31 + 1

    D = zeros(Float16, inds[3])

    i0, j0 = inds[1]
    i1, j1 = inds[2]

    l = 1
    for i = i0:i1
        # cZi = unsafe(cZ[i])
        cZi = cZ[i]
        jj0 = i==i0 ? j0 : i+1
        jj1 = i==i1 ? j1 : M
        for j = jj0:jj1
            # cZj = unsafe(cZ[j])
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
            end
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
            D[l] = dist / N
            l += 1
        end
    end

    return D
end

function compute_dists(cZ::Vector{Vector{UInt64}}, N::Int, M::Int)
    Ds, inds = ptriu(M, Vector{Float16}, compute_dists_chunk, cZ, N, M)
    D = ptriu_compose(Ds, M, inds)
    return D
end

