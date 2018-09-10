__precompile__()
module SimpleMultigrid

# dependencies
using IterativeSolvers, StaticArrays

# export statements
export laplace1d, laplace2d, laplace3d, elliptic1d, elliptic2d, elliptic3d
export Injection, FullWeighting, Cubic
export GaussSeidel, Jacobi
export coarsen
export MultigridMethod, V_cycle, W_cycle, F_cycle, \

# import statements
import Base: expand, show, start, next, done, \

# include source files
include("compose_matrix.jl")
include("grid_transfer_operators.jl")
include("iterative_methods.jl")
include("grids.jl")
include("multigrid.jl")

end # module
