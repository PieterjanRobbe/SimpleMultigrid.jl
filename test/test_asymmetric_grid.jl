## test_asymmetric_grid.jl : test implementation for domains of n-by-m points, with m != n

## POISSON PROBLEM
verbose && println("## testing asymmetric grid, 2D poisson problem")

# grid specification
ngrids = 3
val = zeros(ngrids,ngrids)

for i = 1:ngrids, j = 1:ngrids
	m = 8*2^i
	n = 8*2^j
	A = poisson2d(m,n)
	f = ones((m-1)*(n-1))
	u = reshape(A\f,m-1,n-1)
	val[i,j] = u[8*2^(i-1),8*2^(j-1)]
end

# "exact value"
m = 8*2^(ngrids+1)
n = 8*2^(ngrids+1)
A = poisson2d(m,n)
f = ones((m-1)*(n-1))
u = reshape(A\f,m-1,n-1)
exact = u[8*2^ngrids,8*2^ngrids]

err = abs.(val-exact)./abs(exact)

verbose && begin
	println("grid     rel. error")
	println("-------------------")
	for i = 1:ngrids, j=1:ngrids
		println(@sprintf("[%2i,%2i]",i,j)*"  "*@sprintf("%3.2e",err[i,j]))
	end
end

for i = 1:ngrids, j = 1:ngrids
	@test err[i,j] ≈err[j,i] atol=sqrt(eps())
end

## GENERAL ELLIPTIC PROBLEM
verbose && println("## testing asymmetric grid, general elliptic problem")

# grid specification
ngrids = 3
val = zeros(ngrids,ngrids)

for i = 1:ngrids, j = 1:ngrids
	m = 8*2^i
	n = 8*2^j
	kx = 1+sin.(π*(1/2/m:1/m:1-1/2/m))*sin.(π*(1/n:1/n:1-1/n))'
	ky = 1+sin.(π*(1/m:1/m:1-1/m))*sin.(π*(1/2/n:1/n:1-1/2/n))'
	A = elliptic2d(kx,ky)
	f = ones((m-1)*(n-1))
	u = reshape(A\f,(m-1,n-1))
	val[i,j] = u[8*2^(i-1),8*2^(j-1)]
end

# "exact value"
m = 8*2^(ngrids+1)
n = 8*2^(ngrids+1)
kx = 1+sin.(π*(1/2/m:1/m:1-1/2/m))*sin.(π*(1/n:1/n:1-1/n))'
ky = 1+sin.(π*(1/m:1/m:1-1/m))*sin.(π*(1/2/n:1/n:1-1/2/n))'
A = elliptic2d(kx,ky)
f = ones((m-1)*(n-1))
u = reshape(A\f,(m-1,n-1))
exact = u[8*2^ngrids,8*2^ngrids]

err = abs.(val-exact)./abs(exact)

verbose && begin
	println("grid     rel. error")
	println("-------------------")
	for i = 1:ngrids, j=1:ngrids
		println(@sprintf("[%2i,%2i]",i,j)*"  "*@sprintf("%3.2e",err[i,j]))
	end
end

for i = 1:ngrids, j = 1:ngrids
	@test err[i,j] ≈err[j,i] atol=sqrt(eps())
end

