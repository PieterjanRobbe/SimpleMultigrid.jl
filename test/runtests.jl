using SimpleMultigrid, Base.Test, Requires

include("test_utils.jl")
include("test_problems.jl")

include("test_convergence.jl")
include("test_multigrid.jl")
include("test_direct_discretization.jl")
