module ParallelTriU

export TriuInd, ptriu, ptriu_compose, ptriu_compose_blocks

using LinearAlgebra
using Base.Threads
# using Distributed

## function borrowed from stdlib/Distributed test file
function get_num_blas_threads()::Int
    blas = LinearAlgebra.BLAS.vendor()
    # Wrap in a try to catch unsupported blas versions
    try
        if blas == :openblas
            return ccall((:openblas_get_num_threads, Base.libblas_name), Cint, ())
        elseif blas == :openblas64
            return ccall((:openblas_get_num_threads64_, Base.libblas_name), Cint, ())
        elseif blas == :mkl
            return ccall((:MKL_Get_Max_Num_Threads, Base.libblas_name), Cint, ())
        end

        # OSX BLAS looks at an environment variable
        if Sys.isapple()
            return tryparse(Cint, get(ENV, "VECLIB_MAXIMUM_THREADS", "1"))
        end
    catch
    end

    return Int(get(ENV, "OMP_NUM_THREADS", Sys.CPU_THREADS))
end

# create a single-threaded-BLAS scope with do..end
function disable_blas_threads(f)
    old_num_threads = get_num_blas_threads()
    try
        LinearAlgebra.BLAS.set_num_threads(1)
        return f()
    finally
        LinearAlgebra.BLAS.set_num_threads(old_num_threads)
    end
end


const TriuInd = Tuple{Tuple{Int,Int},Tuple{Int,Int},Int}

function ptriu(sz::Int, RT::Type, func::Function, args...)
    tot_inds = (sz * (sz-1)) ÷ 2
    # nw = nworkers() # use this for multi-processing
    nw = nthreads()

    if tot_inds ≥ nw
        inds_dist = diff([round(Int,x) for x in range(1, tot_inds+1, length=nw+1)])
    else
        inds_dist = [ones(Int, tot_inds); zeros(Int, nw-tot_inds)]
    end

    inds = Array{TriuInd}(undef, nw)
    i0, j0 = 1, 2
    for p = 1:nw
        l = inds_dist[p]
        I1 = (i0, j0)
        I3 = l
        i1, j1 = i0, j0
        while l > 0
            if sz-j1+1 < l
                l -= sz-j1+1
                i1 += 1
                j1 = i1 + 1
            else
                j1 += l - 1
                l = 0
            end
        end
        I2 = (i1, j1)
        inds[p] = (I1, I2, I3)
        i0, j0 = i1, j1
        j0 += 1
        if j0 > sz
            i0 += 1
            j0 = i0 + 1
        end
    end

    ret = Array{RT}(undef, nw)

    disable_blas_threads() do
        ## use this for multi-processing
        # wrk = workers()
        # @sync for p = 1:nw
        #     @async ret[p] = remotecall_fetch(func, wrk[p], inds[p], args...)
        # end

        Threads.@threads for p = 1:nw
            ret[p] = func(inds[p], args...)
        end
    end

    return ret, inds
end

function ptriu_compose(src::Vector{Vector{T}}, sz::Int, inds::Vector{TriuInd}) where {T}
    dest = Array{T}(undef, sz, sz)

    @assert length(src) == length(inds)

    for p = 1:length(src)
        i0, j0 = inds[p][1]
        i1, j1 = inds[p][2]
        srcp = src[p]
        l = 1
        for i = i0:i1
            jj0 = i==i0 ? j0 : i+1
            jj1 = i==i1 ? j1 : sz
            for j = jj0:jj1
                x = srcp[l]
                dest[i,j] = x
                dest[j,i] = x
                l += 1
            end
        end
    end
    for i = 1:sz
        dest[i,i] = 0
    end

    return dest
end

function ptriu_compose_blocks(src::Vector{Vector{Matrix{T}}}, sz::Int, bsz::Int, inds::Vector{TriuInd}) where {T}
    dest = zeros(T, sz * bsz, sz * bsz)

    @assert length(src) == length(inds)

    for p = 1:length(src)
        i0, j0 = inds[p][1]
        i1, j1 = inds[p][2]
        srcp = src[p]
        l = 1
        for i = i0:i1
            jj0 = i==i0 ? j0 : i+1
            jj1 = i==i1 ? j1 : sz
            for j = jj0:jj1
                block = srcp[l]
                row = (i-1)*bsz .+ (1:bsz)
                col = (j-1)*bsz .+ (1:bsz)
                dest[row, col] .= block
                dest[col, row] .= block'
                l += 1
            end
        end
    end

    return dest
end

end # module
