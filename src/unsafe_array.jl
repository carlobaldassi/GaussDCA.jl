module Unsafe

export unsafe
import Base: getindex, length, size, iterate

# this is unsafe in 2 ways:
#  1) no boundary checking
#  2) calling functions must keep a reference to the
#     original Array to avoid garbage collection
#     (e.g. A=UnsafeArray(A) is dangerous)
# it must be used with care!
struct UnsafeArray{T} <: AbstractVector{T}
    p::Ptr{T}
    l::Int
    UnsafeArray(A::Array{T}) where {T} = new{T}(pointer(A), length(A))
end

getindex(uA::UnsafeArray, i::Int) = unsafe_load(uA.p, i)
length(uA::UnsafeArray) = uA.l
size(uA::UnsafeArray) = (uA.l,)

function iterate(uA::UnsafeArray, ind=1)
    length(uA) < ind && return nothing
    return uA[ind], ind+1
end

unsafe(x) = x
unsafe(A::Array) = UnsafeArray(A)

end
