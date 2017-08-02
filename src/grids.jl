#
# GRIDS.JL 
# 
# Functions and definitions for dealing with grid-related operators
# such as interpolate, restrict

mutable struct Grid
	A # system matrix
	f # right-hand side
	v # current approximation
	x # points where the solution is defined
end

"""
	prolong(coarse_grid, fine_grid)

	Prolongate (interpolate) the solution from the coarse_grid \Omega^2h to the fine_grid \Omega^h
"""
prolong(coarse_grid::Grid, fine_grid::Grid) = interpolate(coarse_grid.x, coarse_grid.v, fine_grid.x)

"""
	restrict_residu(fine_grid, coarse_grid)

	Restrict the residu from the fine_grid \Omega^h to the coarse_grid \Omega^2h 
"""
restrict_residu(fine_grid::Grid, coarse_grid::Grid) = interpolate(fine_grid.x, fine_grid.f-fine_grid.A*fine_grid.v, coarse_grid.x)


"""
	restrict_rhs(fine_grid, coarse_grid)

	Restrict the right-hand side from the fine_grid \Omega^h to the coarse_grid \Omega^2h 
"""
restrict_rhs(fine_grid::Grid, coarse_grid::Grid) = interpolate(fine_grid.x, fine_grid.f, coarse_grid.x)
