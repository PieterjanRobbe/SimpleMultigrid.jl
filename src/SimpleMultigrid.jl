module SimpleMultigrid

# dependencies
using IterativeSolvers, LinearAlgebra, StaticArrays, SparseArrays, Printf

# export statements
export laplace1d, laplace2d, laplace3d, elliptic1d, elliptic2d, elliptic3d
export Injection, FullWeighting, Cubic
export GaussSeidel, Jacobi, RedBlackGaussSeidel
export coarsen
export MultigridMethod, V, W, F, V_cycle, W_cycle, F_cycle, \, size

# import statements
import Base: expand, show, iterate, \, size

# include source files
include("compose_matrix.jl")
include("grid_transfer_operators.jl")
include("iterative_methods.jl")
include("grids.jl")
include("multigrid.jl")

end # module
