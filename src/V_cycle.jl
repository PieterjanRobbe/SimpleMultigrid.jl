#
# V_CYCLE.JL
#
# Implementation of a multigrid V-cycle.

function V_cycle!(grids::Array{Grid}, nu1::Integer, nu2::Integer)
	relax(grids[1],nu1)
	if length(grids) != 1
		grids[2].f = restrict_residu(grids[1],grids[2])
		grids[2].v = zero(grids[2].v)
		V_cycle!(grids[2:end], nu1, nu2)
		grids[1].v += prolong(grids[2],grids[1])
	else # direct solve
		grids[1].v = grids[1].A\grids[1].f
	end
	relax(grids[1],nu2)
end

function W_cycle!(grids::Array{Grid}, mu::Integer; nu1::Integer=2, nu2::Integer=1)
	relax(grids[1],nu1)
	if length(grids) != 1
		grids[2].f = restrict_residu(grids[1],grids[2])
		grids[2].v = zero(grids[2].v)
		for i = 1:mu
			W_cycle!(grids[2:end], mu, nu1=nu1, nu2=nu2)
		end
		grids[1].v += prolong(grids[2],grids[1])
	else # direct solve
		grids[1].v = grids[1].A\grids[1].f
	end
	relax(grids[1],nu2)
end
