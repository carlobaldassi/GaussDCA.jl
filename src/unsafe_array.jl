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

# horrible macro to keep compatibility with both julia 0.5 and 0.6,
# while avoiding some even more horrible syntax
macro inner(T, ex)
    VERSION < v"0.6-" && return esc(ex)
    @assert Base.Meta.isexpr(ex, [:(=), :function])
    @assert length(ex.args) == 2
    @assert isa(ex.args[1], Expr) && ex.args[1].head == :call
    @assert isa(ex.args[1].args[1], Symbol)
    fn = ex.args[1].args[1]
    fargs = ex.args[1].args[2:end]
    body = ex.args[2]
    annotate_new!(T, body)

    return esc(Expr(ex.head, Expr(:where, Expr(:call, Expr(:curly, fn, T), fargs...), T), body))
end

# this is unsafe in 2 ways:
#  1) no boundary checking
#  2) calling functions must keep a reference to the
#     original Array to avoid garbage collection
#     (e.g. A=UnsafeArray(A) is dangerous)
# it must be used with care!
immutable UnsafeArray{T} <: AbstractVector{T}
    p::Ptr{T}
    l::Int
    @inner T UnsafeArray(A::Array{T}) = new(pointer(A), length(A))
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
