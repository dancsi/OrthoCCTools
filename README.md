# OrthoCCTools - tools for constructing orthogonal sets of coiled coils

This repo contains several interconnected tools for constructing orthogonal (sub)sets of coiled coils.
 - `fastscore`: a tool that takes a list of peptide sequences in a .fasta file, and scores their pairwise interactions
 - `pyccsscore`, `jsccscore`: Python and Javascript bindings for the scoring functions used by `fastscore`
 - `setbuilder`: attempts to build a large orthogonal set from a smaller one
 - `solver`: takes an interaction matrix produced by `fastscore` and outputs the largest orthogonal set
 
 # Build instructions
 
The code requires a functional C++17 compiler. It has no external dependencies, and can be compiled using CMake as follows:
```shell
mkdir build
cd build
cmake ..
cmake --build .
```
The only exception to this is `jsccscore`, which is built separately using the [WASI SDK](https://github.com/WebAssembly/wasi-sdk). See its [README.md](src/jsccscore/README.md) for more details.
