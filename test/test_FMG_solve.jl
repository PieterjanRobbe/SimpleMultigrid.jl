#
# TEST_FMG_SOLVE.JL
#
# Test FMG for general elliptic PDEs

## 1D
verbose && println("## testing FMG solver for 1D poisson problem")

# grid specification
ngrids = 10
grids = Grid[]
for i = ngrids:-1:1
	m = 2^i
	grid = Grid(poisson1d(m),ones(m-1),zeros(m-1),[1/m:1/m:1-1/m])
	push!(grids,grid)
end

# FMG solve
tol = 1e-6
FMG_solve(grids,tol)
u_hat = grids[1].v[2^(ngrids-1)]

@test u_hat ≈1/8 atol=tol

## 2D
verbose && println("## testing FMG solver for 2D elliptic problem")

# grid specification
ngrids = 5
grids = Grid[]
for i = ngrids:-1:1
	n = 2^i
	m = 2^i
	kx = 1+sin.(5π*(1/2/m:1/m:1-1/2/m))*sin.(5π*(1/n:1/n:1-1/n))'
	ky = 1+sin.(5π*(1/m:1/m:1-1/m))*sin.(5π*(1/2/n:1/n:1-1/2/n))'
	A = elliptic2d(kx,ky)
	f = ones((n-1)*(m-1))
	grid = Grid(A,f,zeros((m-1)*(n-1)),[1/m:1/m:1-1/m,1/n:1/n:1-1/n])
	push!(grids,grid)
end

# "exact" discrete solution
exact = reshape(grids[1].A\grids[1].f,Tuple(length.(grids[1].x)))
Q_exact = exact[2^(ngrids-1),2^(ngrids-1)]
residu = norm(grids[1].f-grids[1].A*exact[:])

# FMG solver
FMG_solve(grids,residu)
u_hat = reshape(grids[1].v,Tuple(length.(grids[1].x)))
Q_hat = u_hat[2^(ngrids-1),2^(ngrids-1)]

# test value at mid point
verbose && println(@sprintf("comparing exact value %0.6f with computed value %0.6f, error is %3.2e",Q_exact,Q_hat,discr_norm(Q_exact-Q_hat,grids[1].x[1][1],1)))
@test Q_hat ≈Q_exact atol=sqrt(sqrt(eps()))

## TIMING
verbose && println("## timing test")

function make_grids(ngrids)
	grids = Grid[]
	for i = ngrids:-1:1
		n = 2^(i)
		m = 2^(i)
		kx = 1+sin.(5π*(1/2/m:1/m:1-1/2/m))*sin.(5π*(1/n:1/n:1-1/n))'
		ky = 1+sin.(5π*(1/m:1/m:1-1/m))*sin.(5π*(1/2/n:1/n:1-1/2/n))'
		A = elliptic2d(kx,ky)
		f = ones((n-1)*(m-1))
		grid = Grid(A,f,zeros((m-1)*(n-1)),[1/m:1/m:1-1/m,1/n:1/n:1-1/n])
		push!(grids,grid)
	end
	return grids
end

function time_direct(ngrids)
	grids = make_grids(ngrids)
	return @elapsed grids[1].A\grids[1].f
end

function time_FMG(ngrids)
	grids = make_grids(ngrids)
	return @elapsed FMG_solve(grids,1e-3)
end

relax(grid,nsweeps) = jacobi(grid,nsweeps)

time_direct(3)
t_direct = time_direct(8)
verbose && println("direct solve takes $(t_direct) s")
time_FMG(3)
t_FMG = time_FMG(8)
verbose && println("FMG solve takes $(t_FMG) s")
