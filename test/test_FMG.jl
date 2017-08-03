#
# TEST_FMG.JL
#
# Test FMG implementation for simple poisson problem.

## 1D POISSON PROBLEM
verbose && println("## testing FMG, 1D poisson problem")

# grid specification
ngrids = 10
grids = Grid[]
for i = ngrids:-1:1
	m = 2^i
	grid = Grid(poisson1d(m),ones(m-1),zeros(m-1),[1/m:1/m:1-1/m])
	push!(grids,grid)
end

# FMG
FMG!(grids,1)

# extract solution
u_hat = grids[1].v

# test value at mid point
verbose && println(@sprintf("comparing exact value %0.6f with computed value %0.6f, error is %3.2e",1/8,u_hat[2^(ngrids-1)],discr_norm(1/8-u_hat[2^(ngrids-1)],grids[1].x[1][1],2)))
@test 1/8 ≈u_hat[2^(ngrids-1)] atol=sqrt(sqrt(eps()))

## 2D POISSON PROBLEM
verbose && println("## testing FMG, 2D poisson problem")

# grid specification
ngrids = 5
grids = Grid[]
for i = ngrids:-1:2
	m = 2^i
	grid = Grid(poisson2d(m,m),ones((m-1)*(m-1)),zeros((m-1)*(m-1)),[1/m:1/m:1-1/m,1/m:1/m:1-1/m])
	push!(grids,grid)
end

# FMG
FMG!(grids,1)

# extract solution
u_hat = reshape(grids[1].v,Tuple(length.(grids[1].x)))
Q_hat = u_hat[2^(ngrids-1),2^(ngrids-1)]

# "exact" discrete solution
u_discr = reshape(grids[1].A\grids[1].f,Tuple(length.(grids[1].x)))
Q_discr = u_discr[2^(ngrids-1),2^(ngrids-1)]

# test value at mid point
verbose && println(@sprintf("comparing exact value %0.6f with computed value %0.6f, error is %3.2e",Q_discr,Q_hat,discr_norm(Q_discr-Q_hat,grids[1].x[1][1],2)))
@test Q_hat ≈Q_discr atol=sqrt(sqrt(eps()))

