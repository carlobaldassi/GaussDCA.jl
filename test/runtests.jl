module GaussDCATests

using GaussDCA
using Compat
using Compat.Test
import Compat: @info

datadir = joinpath(dirname(@__FILE__), "data")
isdir(datadir) || (datadir = joinpath(Pkg.dir(), "GaussDCA.jl", "test", "data"))
isdir(datadir) || error("data directory not found")

const fastafile_s = joinpath(datadir, "small.fasta.gz")
const FNRfile_s = joinpath(datadir, "small.FNRout.txt")
const DIRfile_s = joinpath(datadir, "small.DIRout.txt")
const DIRfile_s2 = joinpath(datadir, "small.DIRout2.txt")

const fastafile_l = joinpath(datadir, "large.fasta.gz")
const DIRfile_l = joinpath(datadir, "large.DIRout.txt")

macro tostring(ex)
    @assert ex.head == :call
    f = esc(ex.args[1])
    a = map(esc, ex.args[2:end])
    newcall = Expr(:call, f, :io, a...)
    quote
        io = IOBuffer()
        $newcall
        String(take!(io))
    end
end

function todict(r)
    d = Dict{NTuple{2,Int},Float64}()
    for l in Compat.split(r, ['\r', '\n'], keepempty = false)
        sl = split(l)
        @test length(sl) == 3
        i, j, x = parse(Int, sl[1]), parse(Int, sl[2]), parse(Float64, sl[3])
        @test !haskey(d, (i,j))
        d[(i,j)] = x
    end
    return d
end

function compare_results(r1, r2)
    d1 = todict(r1)
    d2 = todict(r2)
    allk = sort!(collect(keys(d1)))
    @test allk == sort!(collect(keys(d2)))

    for k in allk
        @test d1[k] ≈ d2[k]
    end
    return true
end


function test1()
    FNR = gDCA(fastafile_s)
    results_FNR = @tostring printrank(FNR)
    expected_results_FNR = read(FNRfile_s, String)
    compare_results(results_FNR, expected_results_FNR)

    DIR = gDCA(fastafile_s, pseudocount = 0.2, score = :DI, remove_dups = true);
    results_DIR = @tostring printrank(DIR)
    expected_results_DIR = read(DIRfile_s, String)
    compare_results(results_DIR, expected_results_DIR)

    DIR2 = gDCA(fastafile_s, pseudocount = 0.2, score = :DI, theta = 0.0, max_gap_fraction = 0.8, min_separation = 4)
    results_DIR2 = @tostring printrank(DIR2)
    expected_results_DIR2 = read(DIRfile_s2, String)
    compare_results(results_DIR2, expected_results_DIR2)
    return true
end

function test2()
    DIR = gDCA(fastafile_l, pseudocount = 0.2, score = :DI, remove_dups = true);
    results_DIR = @tostring printrank(DIR)
    expected_results_DIR = read(DIRfile_l, String)
    compare_results(results_DIR, expected_results_DIR)
    return true
end

function test3()
    ENV["GDCA_FORCE_FALLBACK"] = "true"
    DIR = gDCA(fastafile_s, pseudocount = 0.2, score = :DI, remove_dups = true);
    results_DIR = @tostring printrank(DIR)
    expected_results_DIR = read(DIRfile_s, String)
    compare_results(results_DIR, expected_results_DIR)
    delete!(ENV, "GDCA_FORCE_FALLBACK")
    return true
end

for t in [test1, test2, test3]
    @info("Running $t")
    t()
end

end
