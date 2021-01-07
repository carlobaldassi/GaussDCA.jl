function compute_theta_chunk(inds::TriuInd, cZ::Vector{Vector{UInt64}}, N::Int, M::Int)
    cl = clength(N)
    cr = 5 * (packfactor - crest(N)) + packrest

    kmax = (cl - 1) ÷ 31
    rmax = (cl - 1) % 31

    meanfracid = 0.0
    i0, j0 = inds[1]
    i1, j1 = inds[2]
    @inbounds for i = i0:i1
        cZi = cZ[i]
        nids::UInt64 = 0
        jj0 = i==i0 ? j0 : i+1
        jj1 = i==i1 ? j1 : M
        for j = jj0:jj1
            cZj = cZ[j]
            czi, czj = 1, 1
            z::UInt64 = 0
            for k = 1:kmax
                z = 0
                for r = 1:31
                    zi, zj = cZi[czi], cZj[czj]
                    czi += 1
                    czj += 1

                    ny = (~zi) ⊻ zj
                    z += nz_aux(ny)
                end
                t = collapse(z)
                nids += t
            end
            z = 0
            for r = 1:rmax
                zi, zj = cZi[czi], cZj[czj]
                czi += 1
                czj += 1
                ny = (~zi) ⊻ zj
                z += nz_aux(ny)
            end
            zi, zj = cZi[czi], cZj[czj]
            czi += 1
            czj += 1
            ny = (~zi) ⊻ zj
            z += nz_aux2(ny, cr)
            t = collapse(z)
            nids += t
        end
        fracid = nids / N
        meanfracid += fracid
    end
    return meanfracid
end

# slow fallback
function compute_theta_chunk(inds::TriuInd, ZZ::Vector{Vector{Int8}}, N::Int, M::Int)
    meanfracid = 0.0
    i0, j0 = inds[1]
    i1, j1 = inds[2]
    @inbounds for i = i0:i1
        Zi = ZZ[i]
        jj0 = i==i0 ? j0 : i+1
        jj1 = i==i1 ? j1 : M
        for j = jj0:jj1
            Zj = ZZ[j]
            nids = 0
            for k = 1:N
                nids += Zi[k] == Zj[k]
            end
            fracid = nids / N
            meanfracid += fracid
        end
    end
    return meanfracid
end

function compute_theta(cZ::Vector{Vector{T}}, N::Int, M::Int) where {T<:Union{Int8,UInt64}}
    chunk_means, _ = ptriu(M, Float64, compute_theta_chunk, cZ, N, M)
    meanfracid = sum(chunk_means) / (0.5 * M * (M-1))
    θ = min(0.5, 0.38 * 0.32 / meanfracid)
    return θ
end

