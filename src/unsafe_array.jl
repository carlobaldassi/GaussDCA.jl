module Unsafe

export unsafe

annotate_new!(T, ex) = ex
function annotate_new!(T, ex::Expr)
    if ex.head == :call && ex.args[1] == :new
        ex.args[1] = Expr(:curly, :new, T)
    end
    map!(x->annotate_new!(T, x), ex.args, ex.args)
    return ex
end

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

Base.start(uA::UnsafeArray) = 1
Base.next(uA::UnsafeArray, ind::Int) = (uA[ind], ind+1)
Base.done(uA::UnsafeArray, ind::Int) = ind > length(uA)

unsafe(x) = x
unsafe(A::Array) = UnsafeArray(A)

end
