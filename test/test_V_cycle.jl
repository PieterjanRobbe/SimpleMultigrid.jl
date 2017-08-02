#
# TEST_V_CYCLE.JL
#
# Test V-cycle implementation for simple poisson problem.

## 1D POISSON PROBLEM
verbose && println("## testing V-cyle, 1D poisson problem")

# grid specification
ngrids = 10
grids = Grid[]
for i = ngrids:-1:1
	m = 2^i
	grid = Grid(poisson1d(m),ones(m-1),zeros(m-1),(1/m:1/m:1-1/m,))
	push!(grids,grid)
end

# "exact" discrete solution
exact = grids[1].A\grids[1].f
Q_exact = 1/8

# repeat V-cycles
ncycles = 15
residu = zeros(ncycles)
error = zeros(ncycles)
for i = 1:ncycles
	V_cycle!(grids,2,1)
	residu[i] = discr_norm(grids[1].f-grids[1].A*grids[1].v,grids[1].x[1][1],1)
	error[i] = discr_norm(1/8-grids[1].v,grids[1].x[1][1],1)
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
@test Q_exact ≈Q_hat atol=sqrt(eps())

## 2D POISSON PROBLEM
verbose && println("## testing V-cyle, 2D poisson problem")

# grid specification
ngrids = 5
grids = Grid[]
for i = ngrids:-1:1
	m = 2^i
	grid = Grid(poisson2d(m,m),ones((m-1)*(m-1)),zeros((m-1)*(m-1)),(1/m:1/m:1-1/m,1/m:1/m:1-1/m))
	push!(grids,grid)
end

# "exact" discrete solution
exact = reshape(grids[1].A\grids[1].f,length.(grids[1].x))
Q_exact = exact[2^(ngrids-1),2^(ngrids-1)]

# repeat V-cycles
ncycles = 15
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
@test Q_exact ≈Q_hat atol=sqrt(eps())

