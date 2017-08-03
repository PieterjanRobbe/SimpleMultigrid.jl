#
# SOLVERS.JL
#
# Collection of relaxation methods for multigrid package.

# from IterativeSolvers
jacobi(grid::Grid,nsweeps) = jacobi!(grid.v, grid.A, grid.f, maxiter=nsweeps)

gauss_seidel(grid::Grid,nsweeps) = gauss_seidel!(grid.v, grid.A, grid.f, maxiter=nsweeps)

# from KrylovMethods
# has better sparse porperties for stationary smoothers
function ssor(grid::Grid,nsweeps)
	grid.v,flag,err,iter,resvec = KrylovMethods.ssor(grid.A, grid.f, x=grid.v, maxIter=nsweeps, out=-1)
end

# defaults to gauss-seidel
#relax(grid::Grid,nsweeps::Integer) = gauss_seidel(grid,nsweeps)
relax(grid::Grid,nsweeps::Integer) = ssor(grid,nsweeps)
