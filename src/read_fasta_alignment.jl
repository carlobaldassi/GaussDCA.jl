module ReadFastaAlignment

using FastaIO

# this is silly but it's required for correct cross-linking in docstrings, apparently
using ...DCAUtils

export read_fasta_alignment

"""
    read_fasta_alignment(filename::AbstractString, max_gap_fraction::Real) -> Matrix{Int8}

Parses a FASTA file containing a multiple sequence alignment, and returns a matrix of integers
that represents one sequence per column.

The mapping between the aminoacid symbols and the integers uses this table:

```
  A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
```

Any unrecognized capital letter and the gap symbol `-` are mapped to the value 21; any other
symbol or lowercase letter is ignored.

If a sequence contains a fraction of gaps that exceeds `max_gap_fraction`, it is discarded. Set this
value to 1 to keep all the sequences.

The input file can be plaintext (ASCII) or gzip-compressed plaintext (with the extension ".gz")

See also [`remove_duplicate_sequences`](@ref).
"""
function read_fasta_alignment(filename::AbstractString, max_gap_fraction::Real)
    max_gap_fraction = Float64(max_gap_fraction)
    return FastaReader(filename) do f
        # pass 1
        seqs = Int[]
        inds = Int[]
        fseqlen = 0
        for (name, seq) in f
            ngaps = 0
            if f.num_parsed == 1
                ls = length(seq)
                resize!(inds, ls)
                for i = 1:ls
                    c = seq[i]
                    if c ≠ '.' && !islowercase(c)
                        fseqlen += 1
                        inds[fseqlen] = i
                        c == '-' && (ngaps += 1)
                    end
                end
            else
                ls = length(seq)
                ls == length(inds) || error("inputs are not aligned")
                tstfseqlen = 0
                for i = 1:ls
                    c = seq[i]
                    if c ≠ '.' && !islowercase(c)
                        tstfseqlen += 1
                        inds[tstfseqlen] == i || error("inconsistent inputs")
                        c == '-' && (ngaps += 1)
                    end
                end
                tstfseqlen == fseqlen || error("inconsistent inputs")
            end
            ngaps / fseqlen <= max_gap_fraction && push!(seqs, f.num_parsed)
        end

        length(seqs) > 0 || error("Out of $(f.num_parsed) sequences, none passed the filter (max_gap_fraction=$max_gap_fraction)")

        # pass 2
        Z = Array{Int8}(undef, fseqlen, length(seqs))
        seqid = 1
        for (name, seq) in f
            seqs[end] < f.num_parsed && break
            seqs[seqid] == f.num_parsed || continue
            for i = 1:fseqlen
                c = seq[inds[i]]
                Z[i, seqid] = letter2num(c)
            end
            seqid += 1
        end
        @assert seqid == length(seqs) + 1

        return Z
    end
end

let alphabet = [ 1,21, 2, 3, 4, 5, 6, 7, 8,21, 9,10,11,12,21,13,14,15,16,17,21,18,19,21,20]
               # A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y
    global function letter2num(c::Union{Char,UInt8})
        i = UInt8(c) - 0x40
        1 <= i <= 25 && return alphabet[i]
        return 21
    end
end

end # module
