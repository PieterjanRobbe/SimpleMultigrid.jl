#
# TEST_NON_CONSTANT_RHS.JL
#
# Test example from "A Multigrid Tutorial", page 64.

## V-CYCLES
verbose && println("## testing non_constant_rhs, V-cycles")

# specify grid size
ell = 5
n = 2^ell

# compute grids
grids = Grid[]
for i = ell:-1:1
	m = 2^i
	x = 1/m:1/m:1-1/m
	y = x'
	f = @. 2((1-6x^2)*y^2*(1-y^2)+(1-6y^2)*x^2*(1-x^2))
	grid = Grid(poisson2d(m,m),f[:],zeros((m-1)*(m-1)),[1/m:1/m:1-1/m,1/m:1/m:1-1/m])
	push!(grids,grid)
end

# exact solution
x = grids[1].x[1]
y = x'
exact = @. (x^2-x^4)*(y^4-y^2)

# repeat V-cycles
ncycles = 15
residu = zeros(ncycles)
error = zeros(ncycles)
for i = 1:ncycles
	V_cycle!(grids,2,1)
	residu[i] = discr_norm(grids[1].f-grids[1].A*grids[1].v,grids[1].x[1][1],1)
	error[i] = discr_norm(exact[:]-grids[1].v,grids[1].x[1][1],1)
end

# compute ratios
residu_ratio = residu[2:end]./residu[1:end-1]
prepend!(residu_ratio,0)
error_ratio = error[2:end]./error[1:end-1]
prepend!(error_ratio,0)

# print
verbose && begin println("V-cycle   |r|_h      ratio    |e|_h      ratio")
	println("----------------------------------------------")
	for i = 1:ncycles
		println(@sprintf("%3i",i)*"       "*@sprintf("%3.2e",residu[i])*"   "*
			@sprintf("%0.2f",residu_ratio[i])*"     "*@sprintf("%3.2e",error[i])*
			"   "*@sprintf("%0.2f",error_ratio[i]))
	end
end

# extract solution computed by V-cycles
u_v_cycles = reshape(grids[1].v,Tuple(length.(grids[1].x)))

# test norm of residu
@test residu[end] < sqrt(eps())

## FMG
verbose && println("## testing non_constant_rhs, FMG")

# compute grids
grids = Grid[]
for i = ell:-1:1
	m = 2^i
	x = 1/m:1/m:1-1/m
	y = x'
	f = @. 2((1-6x^2)*y^2*(1-y^2)+(1-6y^2)*x^2*(1-x^2))
	grid = Grid(poisson2d(m,m),f[:],zeros((m-1)*(m-1)),[1/m:1/m:1-1/m,1/m:1/m:1-1/m])
	push!(grids,grid)
end

# FMG
FMG!(grids,1)

# extract solution computed by FMG
u_fmg = reshape(grids[1].v,Tuple(length.(grids[1].x)))

# compare solution of FMG to exact and V-cycle solution
verbose && println(@sprintf("difference with exact solution is %3.2e",discr_norm(u_fmg-exact,grids[1].x[1][1],2)))
verbose && println(@sprintf("difference with V-cycle solution is %3.2e",discr_norm(u_v_cycles-u_fmg,grids[1].x[1][1],2)))

@test discr_norm(u_v_cycles-u_fmg,grids[1].x[1][1],2) < sqrt(sqrt(eps()))
