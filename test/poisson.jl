#
# POISSON.JL
#
# Generate the system matrix for the poisson equation as model problem


"""
    laplacian(n)

	Return 1d laplacian operator in n points.
"""
function laplacian(n::Integer)
	d1 = 2*ones(n-1)
	d2 = -ones(n-2)
	return spdiagm((d2,d1,d2),(-1,0,1))
end

"""
    poisson1d(n)

	Return system matrix A for the 1d poisson problem on an n-point grid. 
"""
function poisson1d(n::Integer)
	A = n^2*laplacian(n)	
end

"""
    poisson2d(n,m)

	Return system matrix A for the 2d poisson problem on an n-by-m grid. 
"""
function poisson2d(n::Integer,m::Integer)
	A = n^2*kron(laplacian(n),speye(m-1)) + m^2*kron(speye(n-1),laplacian(m))	
end
