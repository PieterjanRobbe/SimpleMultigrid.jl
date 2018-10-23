# SimpleMultigrid
[![Build Status](https://travis-ci.org/PieterjanRobbe/SimpleMultigrid.jl.png)](https://travis-ci.org/PieterjanRobbe/SimpleMultigrid.jl)

Simple implementation of geometric multigrid methods. 

## Installation

```julia
] add https://github.com/PieterjanRobbe/SimpleMultigrid.jl
```

Load the package in Julia by

```julia
using SimpleMultigrid
```

## Details

Basic implementation of a Geometric Multigrid V-cycle, W-cycle and FMG-cycle.

Implemented smoothers use the `jacobi!` and `gauss_seidel!` methods from the [IterativeSolvers](https://github.com/JuliaMath/IterativeSolvers.jl) package.
