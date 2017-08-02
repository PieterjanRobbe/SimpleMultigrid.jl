#
# SOLVERS.JL
#
# Collection of relaxation methods for multigrid package.



jacobi(grid::Grid,nsweeps) = jacobi!(grid.v, grid.A, grid.f, maxiter=nsweeps)

gauss_seidel(grid::Grid,nsweeps) = gauss_seidel!(grid.v, grid.A, grid.f, maxiter=nsweeps)

# defaults to gauss-seidel
relax(grid::Grid,nsweeps::Integer) = gauss_seidel(grid,nsweeps)
