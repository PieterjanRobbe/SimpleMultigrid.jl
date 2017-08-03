using SimpleMultigrid
using Base.Test

verbose = true

# problems
include("poisson.jl")
include("elliptic.jl")

# test discretization
include("test_elliptic.jl")

# test multigrid
include("test_V_cycle.jl")
include("test_W_cycle.jl")
include("test_FMG.jl")
include("test_non_constant_rhs.jl")
include("test_sin.jl")
include("test_asymmetric_grid.jl")
include("test_FMG_solve.jl")
