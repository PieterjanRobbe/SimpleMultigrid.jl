# SimpleMultigrid
[![Build Status](https://travis-ci.org/PieterjanRobbe/SimpleMultigrid.jl.png)](https://travis-ci.org/PieterjanRobbe/SimpleMultigrid.jl)

Simple implementation of geometric multigrid methods. 

## Installation

```julia
Pkg.clone("https://github.com/PieterjanRobbe/SimpleMultigrid.jl")
```

Load the package in Julia by

```julia
using SimpleMultigrid.jl
```

## Details

The package provides a basic implementation of a Multigrid V-cycle, W-cycle and FMG-cycle. As input, these methods need an array of `Grid`s.

Implemented smoothers use the `jacobi` and `gauss_seidel` methods from the [IterativeSolvers](https://github.com/JuliaMath/IterativeSolvers.jl) package.

There is also a method `FMG_solve(grids, tol)`, that solves the grid hierarchy `grids` up to a certain tolerance `tol` on the residu.

## Usage

We provide a sample implementation of the poisson problem and a general elliptic PDE in 2d using finite differences. See [test](https://github.com/PieterjanRobbe/SimpleMultigrid.jl/tree/master/test) for more usage examples.