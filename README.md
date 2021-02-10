# OrthoCCTools - tools for constructing orthogonal sets of coiled coils

This repo contains several interconnected tools for constructing orthogonal (sub)sets of coiled coils.
 - `fastscore`: a tool that takes a list of peptide sequences in a .fasta file, and scores their pairwise interactions
 - `pyccsscore`, `jsccscore`: Python and Javascript bindings for the scoring functions used by `fastscore`
 - `setbuilder`: attempts to build a large orthogonal set from a smaller one
 - `solver`: takes an interaction matrix produced by `fastscore` and outputs the largest orthogonal set
 
# Build instructions
 
The code requires a functional C++20 compiler, Python 3 development libraries and OpenMP. It can be build with CMake in the following way:
```shell
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX="$PWD" ..
cmake --build . --config Release
cmake --install . --config Release
```
The only exception to this is `jsccscore`, which is built separately using the [WASI SDK](https://github.com/WebAssembly/wasi-sdk). See its [README.md](src/jsccscore/README.md) for more details.

The build was tested on Windows 10 running Visual Studio 16.8.5 and Python 3.8.5 (Anaconda) as well as Ubuntu 20.10 running GCC 10.2.0 and Python 3.8.6.

# Usage instructions

First, we generate an initial set of 4-heptad peptides.
```shell
mkdir data
python src/util/01_generate.py >data/full4096.fasta
```

Next, we compute the pairwise scores using `fastscore`
```shell
./build/fastscore data/full4096.fasta
```
You will see 3 additional files in your `data/` directory: `full4096.bin`, `full4096.align.bin` and `full4096.orientation.bin`. These files contain the matrices of interaction scores, as well as the "best" alignment and orientation, respectively.
You can view the interaction scores with the `02_plot_score_matrices.py` utility:
```shell
python src/util/02_plot_score_matrices.py data/full4096.bin
```

In order to construct the interaction graph, you need to compute the interaction score thresholds. They can be computed using `03_recommend_cutoff.py`:
```shell
python src/util/03_recommend_cutoff.py data/full4096.bin
```
Now, you can use these cutoffs to construct the interaction graph and find an orthogonal set in it:
```shell
./build/solver data/full4096.bin --binding-cutoff=-8.5 --nonbinding-cutoff=-7
```
The solver will output the orthogonal set to a `.pairs` file with the same name as the score matrix. In our case, the largest orthogonal set will be saved in `data/full4096.pairs`. 
