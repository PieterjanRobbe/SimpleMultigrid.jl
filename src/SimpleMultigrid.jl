module SimpleMultigrid

# dependencies
using Interpolations
using IterativeSolvers
using KrylovMethods

# export statements
export discr_norm # from utils.jl
export V_cycle!, W_cycle! # from V_cycle.jl
export FMG!, FMG_solve # from FMG.jl
export Grid # from grids.jl

# import statements
import Interpolations.interpolate

# include source files
include("utils.jl")
include("grids.jl")
include("solvers.jl")
include("V_cycle.jl")
include("FMG.jl")

end # module
