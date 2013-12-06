module Unsafe

export unsafe

# this is unsafe in 2 ways:
#  1) no boundary checking
#  2) calling functions must keep a reference to the
#     original Array to avoid garbage collection
#     (e.g. A=UnsafeArray(A) is dangerous)
# it must be used with care!
immutable UnsafeArray{T} <: AbstractArray
    p::Ptr{T}
    l::Int
    function UnsafeArray(A::Array{T})
        new(pointer(A), length(A))
    end
end

UnsafeArray{T}(A::Array{T}) = UnsafeArray{T}(A)

Base.getindex(uA::UnsafeArray, i::Int) = unsafe_load(uA.p, i)
Base.length(uA::UnsafeArray) = uA.l

Base.start(uA::UnsafeArray) = 1
Base.next(uA::UnsafeArray, ind::Int) = (uA[ind], ind+1)
Base.done(uA::UnsafeArray, ind::Int) = ind > length(uA)

unsafe(x) = x
unsafe(A::Array) = UnsafeArray(A)

end
