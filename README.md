# SimpleMultigrid

| **Build Status** | **Coverage** |
|------------------|--------------|
| [![Build Status](https://travis-ci.org/PieterjanRobbe/SimpleMultigrid.jl.png)](https://travis-ci.org/PieterjanRobbe/SimpleMultigrid.jl) [![Build status](https://ci.appveyor.com/api/projects/status/e2ui5mkwewwgtny5?svg=true)](https://ci.appveyor.com/project/PieterjanRobbe/simplemultigrid-jl) | [![Coverage](https://codecov.io/gh/PieterjanRobbe/SimpleMultigrid.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PieterjanRobbe/SimpleMultigrid.jl) [![Coverage Status](https://coveralls.io/repos/github/PieterjanRobbe/SimpleMultigrid.jl/badge.svg)](https://coveralls.io/github/PieterjanRobbe/SimpleMultigrid.jl) |

Simple implementation of geometric multigrid methods. 

## Installation

From the Julia REPL, type ] to enter Pkg mode and run

```julia
pkg> add https://github.com/PieterjanRobbe/SimpleMultigrid.jl
```

Load the package in Julia by

```julia
julia> using SimpleMultigrid
```

## Details

Basic implementation of a Geometric Multigrid V-cycle, W-cycle and FMG-cycle.

Implemented smoothers use the `jacobi!` and `gauss_seidel!` methods from [IterativeSolvers.jl](https://github.com/JuliaMath/IterativeSolvers.jl).
