module CompressSequencies

export clength, crest,
       packfactor, packrest,
       nnz_aux, nz_aux, nz_aux2,
       collapse,
       compress_Z, Z_to_cZ

const PACKBITS = 64
const _u = 0x0084210842108421
const _alt = 0x007c1f07c1f07c1f

const packfactor = PACKBITS ÷ 5
const packrest = PACKBITS % 5

const _msk = ~(UInt64(0));

clength(l::Int) = (l-1) ÷ packfactor + 1
crest(l::Int)   = (l-1) % packfactor + 1

nnz_aux(x::UInt64) = ((x | (x >>> 1) | (x >>> 2) | (x >>> 3) | (x >>> 4)) & _u)

nz_aux(nx::UInt64) = (nx & (nx >>> 1) & (nx >>> 2) & (nx >>> 3) & (nx >>> 4) & _u)

nz_aux2(nx::UInt64, s) = nz_aux(nx) & (_msk >>> s)

@inline function collapse(z)
    t = (z & _alt) + ((z >> 5) & _alt)
    t += t >> 10
    t += t >> 20
    t += t >> 40
    t &= 0x3ff
    return t
end

function compress_Z(Z::Matrix{Int8})
    N, M = size(Z)
    ZZ = Vector{Int8}[Z[:,i] for i = 1:M]

    cl = clength(N)

    cZ = [zeros(UInt64, cl) for i=1:M]

    @inbounds for i = 1:M
        cZi = cZ[i]
        ZZi = ZZ[i]
        for k = 1:N
            k0 = (k-1) ÷ packfactor + 1
            k1 = (k-1) % packfactor
            cZi[k0] |= (UInt64(ZZi[k]) << (5 * k1))
        end
    end

    return cZ
end

function Z_to_cZ(Z, q)
    M = size(Z, 2)
    fast = q ≤ 32 && get(ENV, "DCAUTILS_FORCE_FALLBACK", "false") ≠ "true"
    fast || println("GaussDCA: using slower fallbacks")
    return fast ? compress_Z(Z) : [Z[:,i]::Vector{Int8} for i = 1:M]
end

end # module
