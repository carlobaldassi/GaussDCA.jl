module Unsafe

export unsafe

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

Base.getindex(uA::UnsafeArray, i::Int) = unsafe_load(uA.p, i)
Base.length(uA::UnsafeArray) = uA.l
Base.size(uA::UnsafeArray) = (uA.l,)

@static if VERSION < v"0.7.0-DEV.5126"
    Base.start(uA::UnsafeArray) = 1
    Base.next(uA::UnsafeArray, ind::Int) = (uA[ind], ind+1)
    Base.done(uA::UnsafeArray, ind::Int) = ind > length(uA)
else
    function Base.iterate(uA::UnsafeArray, ind=1)
        length(uA) < ind && return nothing
        return uA[ind], ind+1
    end
end

unsafe(x) = x
unsafe(A::Array) = UnsafeArray(A)

end
