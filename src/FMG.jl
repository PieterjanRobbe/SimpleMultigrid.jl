#
# FMG.JL Full Multigrid
#
# Implementation of a full multigrid cycle.

function FMG!(grids::Array{Grid}, nu0::Integer; nu1::Integer=2, nu2::Integer=2)
	if length(grids) != 1
		grids[2].f = restrict_rhs(grids[1],grids[2])
		FMG!(grids[2:end],nu0)
		grids[1].v = prolong(grids[2],grids[1])
	else
		grids[1].v = zero(grids[1].v)
	end
	for i = 1:nu0
		V_cycle!(grids, nu1, nu2)
	end
end

"""
    FMG_solve(grids, tol)

	Solve the grid hierarchy in grids using full multigrid and repeated V-cycles until the norm of the residu is smaller then tol
"""
function FMG_solve(grids, tol)
	FMG!(grids,1)
	residu = norm(grids[1].f-grids[1].A*grids[1].v)
	while residu > tol
		V_cycle!(grids,2,1)
		residu = norm(grids[1].f-grids[1].A*grids[1].v)
	end
end

