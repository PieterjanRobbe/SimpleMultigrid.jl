# test_problems.jl : collection of test problems

# example rhs from A Multigrid Tutorial, scnd. ed.
h(x,y) = 2((1-6x^2)*y^2*(1-y^2)+(1-6y^2)*x^2*(1-x^2))

# default test examples 
g(x) = Float64[1 + sin(π*xᵢ) for xᵢ in x]
g(x,y) = Float64[1 + sin(π*xᵢ)*sin(π*yᵢ) for xᵢ in x, yᵢ in y]
g(x,y,z) = Float64[1 + sin(π*xᵢ)*sin(π*yᵢ)*sin(π*zᵢ) for xᵢ in x, yᵢ in y, zᵢ in z]

f(n) = g(range(1/n, stop=1-1/n, length=n-1))
f(n,m) = g(range(1/n, stop=1-1/n, length=n-1),range(1/m, stop=1-1/m, length=m-1))
f(n,m,l) = g(range(1/n, stop=1-1/n, length=n-1),range(1/m, stop=1-1/m, length=m-1),range(1/l, stop=1-1/l, length=l-1))

p₁₁(n) = laplace1d(n)
p₁₂(n) = laplace2d(n,n)
p₁₃(n) = laplace3d(n,n,n)

p₂₁(n) = elliptic1d(f(n))
p₂₂(n) = elliptic2d(f(n,n))
p₂₃(n) = elliptic3d(f(n,n,n))
