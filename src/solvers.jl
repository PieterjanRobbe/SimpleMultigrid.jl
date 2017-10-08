## solver.jl : collection of relaxation methods for multigrid package

# from IterativeSolvers
jacobi(grid::Grid,nsweeps) = jacobi!(grid.v, grid.A, grid.f, maxiter=nsweeps)

gauss_seidel(grid::Grid,nsweeps) = gauss_seidel!(grid.v, grid.A, grid.f, maxiter=nsweeps)

# from KrylovMethods
# has better sparse properties for stationary smoothers
function ssor(grid::Grid,nsweeps)
    grid.v,flag,err,iter,resvec = KrylovMethods.ssor(grid.A, grid.f, x=grid.v, maxIter=nsweeps, out=-1, omega=2/3)
end

# default
relax(grid::Grid,nsweeps::Integer) = ssor(grid,nsweeps)
