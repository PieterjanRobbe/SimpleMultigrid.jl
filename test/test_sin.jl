#
# TEST_SIN.JL
#
# Test multigrid for general elliptic PDEs.

## 1D
verbose && println("## testing V-cyle, 1D elliptic problem")

# grid specification
ngrids = 10
grids = Grid[]
for i = ngrids:-1:1
	m = 2^i
	x = 1/2/m:1/m:1-1/2/m
	grid = Grid(elliptic1d(1+1/10*sin.(5π*x)),ones(m-1),zeros(m-1),(1/m:1/m:1-1/m,))
	push!(grids,grid)
end

# "exact" discrete solution
exact = grids[1].A\grids[1].f
Q_exact = exact[2^(ngrids-1)]

# repeat V-cycles
ncycles = 10
residu = zeros(ncycles)
error = zeros(ncycles)
for i = 1:ncycles
	V_cycle!(grids,2,1)
	residu[i] = discr_norm(grids[1].f-grids[1].A*grids[1].v,grids[1].x[1][1],1)
	error[i] = discr_norm(exact-grids[1].v,grids[1].x[1][1],1)
end
u_hat = grids[1].v
Q_hat = u_hat[2^(ngrids-1)]

# compute ratios
residu_ratio = residu[2:end]./residu[1:end-1]
prepend!(residu_ratio,0)
error_ratio = error[2:end]./error[1:end-1]
prepend!(error_ratio,0)

# print convergence history
verbose && begin
	println("V-cycle   |r|_h      ratio    |e|_h      ratio")
	println("----------------------------------------------")
	for i = 1:ncycles
		println(@sprintf("%3i",i)*"       "*@sprintf("%3.2e",residu[i])*"   "*
			@sprintf("%0.2f",residu_ratio[i])*"     "*@sprintf("%3.2e",error[i])*
			"   "*@sprintf("%0.2f",error_ratio[i]))
	end
end

# test value at mid point
@test Q_exact ≈Q_hat atol=sqrt(sqrt(eps()))

# now test FMG for same problem
verbose && println("## testing FMG, 1D elliptic problem")

# grid specification
grids = Grid[]
for i = ngrids:-1:1
	m = 2^i
	x = 1/2/m:1/m:1-1/2/m
	grid = Grid(elliptic1d(1+1/10*sin.(5π*x)),ones(m-1),zeros(m-1),(1/m:1/m:1-1/m,))
	push!(grids,grid)
end

# FMG
FMG!(grids,1)

# extract solution
u_hat = grids[1].v
Q_hat = u_hat[2^(ngrids-1)]

# test value at mid point
verbose && println(@sprintf("comparing exact value %0.6f with computed value %0.6f, error is %3.2e",Q_exact,Q_hat,discr_norm(Q_exact-Q_hat,grids[1].x[1][1],1)))
@test Q_hat ≈Q_exact atol=sqrt(sqrt(eps()))

## 2D
verbose && println("## testing V-cyle, 2D elliptic problem")

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
	grid = Grid(A,f,zeros((m-1)*(n-1)),(1/m:1/m:1-1/m,1/n:1/n:1-1/n))
	push!(grids,grid)
end

# "exact" discrete solution
exact = reshape(grids[1].A\grids[1].f,length.(grids[1].x))
Q_exact = exact[2^(ngrids-1),2^(ngrids-1)]

# repeat V-cycles
ncycles = 10
residu = zeros(ncycles)
error = zeros(ncycles)
for i = 1:ncycles
	V_cycle!(grids,2,1)
	residu[i] = discr_norm(grids[1].f-grids[1].A*grids[1].v,grids[1].x[1][1],2)
	error[i] = discr_norm(exact[:]-grids[1].v,grids[1].x[1][1],2)
end
u_hat = reshape(grids[1].v,length.(grids[1].x))
Q_hat = u_hat[2^(ngrids-1),2^(ngrids-1)]

# compute ratios
residu_ratio = residu[2:end]./residu[1:end-1]
prepend!(residu_ratio,0)
error_ratio = error[2:end]./error[1:end-1]
prepend!(error_ratio,0)

# print convergence history
verbose && begin
	println("V-cycle   |r|_h      ratio    |e|_h      ratio")
	println("----------------------------------------------")
	for i = 1:ncycles
		println(@sprintf("%3i",i)*"       "*@sprintf("%3.2e",residu[i])*"   "*
			@sprintf("%0.2f",residu_ratio[i])*"     "*@sprintf("%3.2e",error[i])*
			"   "*@sprintf("%0.2f",error_ratio[i]))
	end
end

# test value at mid point
@test Q_exact ≈Q_hat atol=sqrt(sqrt(eps()))

# now test FMG for same problem
verbose && println("## testing FMG, 2D elliptic problem")

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
	grid = Grid(A,f,zeros((m-1)*(n-1)),(1/m:1/m:1-1/m,1/n:1/n:1-1/n))
	push!(grids,grid)
end

# FMG
FMG!(grids,2)

# extract solution
u_hat = reshape(grids[1].v,length.(grids[1].x))
Q_hat = u_hat[2^(ngrids-1),2^(ngrids-1)]

# test value at mid point
verbose && println(@sprintf("comparing exact value %0.6f with computed value %0.6f, error is %3.2e",Q_exact,Q_hat,discr_norm(Q_exact-Q_hat,grids[1].x[1][1],1)))
@test Q_hat ≈Q_exact atol=sqrt(sqrt(eps()))

