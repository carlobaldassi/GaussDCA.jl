Gaussian Direct Coupling Analysis
=================================

Overview
--------

This is the code which accompanies the paper "Fast and accurate multivariate
Gaussian modeling of protein families: Predicting residue contacts and
protein-interaction partners" by Carlo Baldassi, Marco Zamparo, Christoph
Feinauer, Andrea Procaccini, Riccardo Zecchina, Martin Weight and Andrea
Pagnani, submitted to PLOS ONE (2013).

This code is released under the GPL version 3 (or later) license; see the
`LICENCE.md` file for details.

The code is written in [Julia](www.julialang.org) and requires julia version
0.2 or later; it provides a function which reads a multiple sequence alignment
(in FASTA format) and returns a ranking of all pairs of residue positions in
the aligned amino-acid sequences.

Installation
------------

The recommended way to install this software is by giving this command in the
julia command line:

  ```
  julia> Pkg.clone("https://github.com/carlobaldassi/GaussDCA.jl")
  ```

in this way, all dependencies will be satisfied automatically, and the code
will be upgraded every time the `Pkg.update()` command is used.

Alternatively, you can manually copy the whole directory structure to your
julia package directory (use `Pkg.dir()` to locate it), and then run
`Pkg.update()` to download the dependencies.

Usage
-----

To load the code, just type `using GaussDCA`.  

This software provides one main function, `gDCA(filname::String, ...)`. This
function takes the name of a (possibly gzipped) FASTA file, and returns a
predicted contact ranking, in the form of a Vector of triples, each triple
containing two indices `i` and `j` (with `i < j + 4`) and a score. The indices
start counting from 1, and indicate pair of residue positions in the given
alignment. The triples are sorted by score in descending order, such that
predicted contacts should come up on top.  

For convenience, a utility function is also provided, `printrank(output, R)`,
which prints the result of `gDCA` either in a file or to a stream, given as
first argument.  If the first argument `output` is omitted, the standard
terminal output will be used.

The `gDCA` function takes some additional, optional keyword arguments:

 * `pseudocount`: the value of the pseudo-count parameter, between `0` and `1`.
                  the default is `0.8`, which gives good results when the
                  Frobenius norm score is used (see below); a good value for the
                  Direct Information score is `0.2`.
 * `theta`: the value of the similarity threshold. By default it is `:auto`,
            which means it will be automatically computed (this takes additional
            time); otherwise, a real value between `0` and `1` can be given.
 * `max_gap_fraction`: maximum fraction of gap symbols in a sequence; sequences
                       which exceed this threshold are discarded. The default
                       value is `0.9`.
 * `score`: the scoring function to use. There are two possibilities, `:DI` for
            the Direct Information, and `:frob` for the Frobenius norm. The
            default is `:frob`. (Note the leading colon: this argument is passed
            as a symbol).

The code will be parallelized if more than one julia worker (as obtained by the
`nworkers()` function) is available (by either launching julia with the `-p`
option from the command line, or by using the `addprocs` function). See the
"Additional thechnical notes" section at the end of this document.  

Examples
--------

Here is a basic usage example, assuming an alignment in FASTA format is found
in the file "alignment.fasta.gz":

  ```
  julia> using GaussDCA

  julia> FNR = gDCA("alignment.fasta.gz");

  julia> printrank("results_FN.txt", FNR)
  ```

The above uses the Frobenius norm ranking with default parameters.
This is how to get the Direct Information ranking instead:

  ```
  julia> DIR = gDCA("alignment.fasta.gz", pseudocount = 0.2, score = :DI);

  julia> printrank("results_DI.txt", DIR)
  ```

Additional technical details
----------------------------

The parallelization can be forcefully disabled even in presence of extra
workers, by setting the environment variable `PARALLEL_GDCA` to `false`
before loading the `GaussDCA` module.

In julia 0.2, when using workers, and using either OpenBLAS - which is the
default - or MKL as the BLAS backend, the default behaviour is to disable
threading in BLAS libraries. In this case, i.e. when many workers are found and
parallelization is not manually disabled, the `gDCA` function overrides the
default julia behaviour and sets the number of threads to match the number of
workers (except when running the parallel portions of the code). It then resets
the number of threads to 1 when finished. The number of cores used in the
non-parallel portions of the code can be explicitly controlled by the user via
the `OMP_NUM_THREADS` environment variable.
