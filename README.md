# CoiledCoilOrtho
Program to compute an orthogonal set of peptides.

Compile with

```shell
g++ -std=c++11 src/solver.cpp -o solver
```

```shell
USAGE: solver FILE [OPTIONS...]
Available options:
        --homo-only
        --hetero-only
        --binding-cutoff=PARAM
        --nonbinding-cutoff=PARAM
```