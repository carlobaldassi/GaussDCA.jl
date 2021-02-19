Gaussian Direct Coupling Analysis for protein contacts predicion
================================================================

![CI][CI-url] [![][codecov-img]][codecov-url] [![CODECOV][codecov-img]][codecov-url]

Overview
--------

This is the code which accompanies the paper ["Fast and accurate multivariate
Gaussian modeling of protein families: Predicting residue contacts and
protein-interaction partners"][paper]
by Carlo Baldassi, Marco Zamparo, Christoph Feinauer, Andrea Procaccini,
Riccardo Zecchina, Martin Weigt and Andrea Pagnani, (2014)
PLoS ONE 9(3): e92721. doi:10.1371/journal.pone.0092721

See also [this Wikipedia article][wikiDCA] for a general overview of the Direct
Coupling Analysis technique.

This code is released under the GPL version 3 (or later) license; see the
`LICENSE.md` file for details.

The code is written in [Julia][julia] and requires julia version
1.5 or later; it provides a function which reads
a multiple sequence alignment (in FASTA format) and returns a ranking of all
pairs of residue positions in the aligned amino-acid sequences.

Since version 2, most of the internal functions used to parse and manipulate
the data have been factored out into the package [DCAUtils.jl][DCAUtils].
The code in this module is essentially a wrapper around those utilities.

[paper]: http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0092721
[julia]: https://www.julialang.org
[wikiDCA]: https://en.wikipedia.org/wiki/Direct_coupling_analysis

[CI-url]: https://github.com/carlobaldassi/GaussDCA.jl/workflows/CI/badge.svg

[codecov-img]: https://codecov.io/gh/carlobaldassi/GaussDCA.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/carlobaldassi/GaussDCA.jl

[DCAUtils]: https://github.com/carlobaldassi/DCAUtils.jl

Installation
------------

To install the package, enter in Pkg mode by pressing the <kbd>]</kbd> key,
then in the pkg prompt enter

```
(@v1.5) pkg> add https://github.com/carlobaldassi/GaussDCA.jl"
```

Usage
-----

To load the code, just type `using GaussDCA`.

This software provides one main function, `gDCA(filname::String, ...)`. This
function takes the name of a (possibly gzipped) FASTA file, and returns a
predicted contact ranking, in the form of a Vector of triples, each triple
containing two indices `i` and `j` (with `i` &lt; `j`) and a score. The indices
start counting from 1, and denote pair of residue positions in the given
alignment; pairs which are separated by less than a given number of residues
(by default 5) are filtered out. The triples are sorted by score in descending
order, such that predicted contacts should come up on top.

For convenience, a utility function is also provided, `printrank(output, R)`,
which prints the result of `gDCA` either in a file or to a stream, given as
first argument.  If the first argument `output` is omitted, the standard
terminal output will be used.

The `gDCA` function takes some additional, optional keyword arguments:

 * `pseudocount`: the value of the pseudo-count parameter, between `0` and `1`.
                  the default is `0.8`, which gives good results when the
                  Frobenius norm score is used (see below); a good value for the
                  Direct Information score is `0.2`.
 * `θ`: the value of the similarity threshold. By default it is `:auto`,
      which means it will be automatically computed (this takes additional
      time); otherwise, a real value between `0` and `1` can be given.
 * `max_gap_fraction`: maximum fraction of gap symbols in a sequence; sequences
                       that exceed this threshold are discarded. The default
                       value is `0.9`.
 * `score`: the scoring function to use. There are two possibilities, `:DI` for
            the Direct Information, and `:frob` for the Frobenius norm. The
            default is `:frob`. (Note the leading colon: this argument is passed
            as a symbol).
 * `min_separation`: the minimum separation between residues in the output
                     ranking. Must be ≥ `1`. The default
                     is `5`.

The code is multi-threaded: if you start julia with the `-t` option, for example
as `julia -t 8`, the computations will run in parallel on the given number of
threads.

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
