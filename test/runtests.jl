using SimpleMultigrid, Base.Test, Requires

include("test_utils.jl")

# plot
include("test_convergence.jl")
include("test_multigrid.jl")

