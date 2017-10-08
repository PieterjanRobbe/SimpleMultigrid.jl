## test_elliptic.jl : test finite difference discretization of the general elliptic PDE in 1D and 2D

## 1D

# test if poisson equation equals general elliptic with k=1
verbose && println("## test if 1D poisson equation equals general elliptic with k=1")
n = 1000
P = poisson1d(n)
E = elliptic1d(ones(n))
@test all(P .≈E)

# test solution of 1d elliptic PDE with k = 1+sin(πx)
verbose && println("## test 1D elliptic PDE with k=1+sin(\pi x)")
n = 200
x = 1/2/n:1/n:1-1/2/n
A = elliptic1d(1+sin.(π*x))
f = ones(n-1)
u = A\f

# test grid independence of the solution
verbose && println("## test grid independence of the solution")
ns = 15
us = zeros(ns)
for i in 1:ns
	n = 2^i
	x = 1/2/n:1/n:1-1/2/n
	A = elliptic1d(1+sin.(π*x))
	f = ones(n-1)
	u = A\f
	us[i]=u[2^(i-1)]
end

err = abs.(us[1:ns-1]-us[ns])
rat = zeros(ns-1)
rat[2:end] = err[2:end]./err[1:end-1]
verbose && begin
	println("    n     error      ratio")
	println("--------------------------")
	for i = 1:ns-1
		println(@sprintf("%6i",i)*"    "*@sprintf("%3.2e",err[i])*"   "*
			@sprintf("%0.2f",rat[i]))
	end
end
@test rat[10] ≈0.25 atol=0.001

## 2D

# test if poisson equation equals general elliptic with k=1
verbose && println("## test if 2D poisson equation equals general elliptic with k=1")
n = 10
m = 10
P = poisson2d(n,m)
E = elliptic2d(ones(m,n-1),ones(m-1,n))

@test all(P .≈E)

# test solution of 2d elliptic PDE with k = 1+sin(πx)
verbose && println("## test 2D elliptic PDE with k=1+sin(\pi x)*sin(\pi y)")
n = 2
m = 2
kx = 1+sin.(π*(1/2/m:1/m:1-1/2/m))*sin.(π*(1/n:1/n:1-1/n))'
ky = 1+sin.(π*(1/m:1/m:1-1/m))*sin.(π*(1/2/n:1/n:1-1/2/n))'

A = elliptic2d(kx,ky)
f = ones((n-1)*(m-1))
u = A\f

# test grid independence of the solution
verbose && println("## test grid independence of the solution")
ns = 10
us = zeros(ns)
for i in 1:ns
	n = 2^i
	m = 2^i
	kx = 1+sin.(π*(1/2/m:1/m:1-1/2/m))*sin.(π*(1/n:1/n:1-1/n))'
	ky = 1+sin.(π*(1/m:1/m:1-1/m))*sin.(π*(1/2/n:1/n:1-1/2/n))'
	A = elliptic2d(kx,ky)
	f = ones((n-1)*(m-1))
	u = reshape(A\f,(n-1,m-1))
	us[i]=u[2^(i-1),2^(i-1)]
end

err = abs.(us[1:ns-1]-us[ns])
rat = zeros(ns-1)
rat[2:end] = err[2:end]./err[1:end-1]
verbose && begin
	println("    n     error      ratio")
	println("--------------------------")
	for i = 1:ns-1
		println(@sprintf("%6i",i)*"    "*@sprintf("%3.2e",err[i])*"   "*
			@sprintf("%0.2f",rat[i]))
	end
end
@test rat[6] ≈0.25 atol=0.001

