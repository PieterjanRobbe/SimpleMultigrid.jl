# test_direct_discretization.jl

# direct discretization for elliptic problem in 2D/3D
function get_problem_2d(n)
	A = MultigridMethod((n,m) -> p₂₂(n), (n, n), V(3, 2))
	b = fill(1,(n-1)*(n-1))
	return (A,b)
end

@testset "Direct discretization of elliptic PDE, 2D" begin
	for n in [16 32 64 128 256 512 1024]
		(A, b) = get_problem_2d(n)
		A.grids[1].b .= b # copy rhs
		push!(A.resnorm,SimpleMultigrid.norm_of_residu(A.grids[1])) # log convergence history
		for i in 1:15
			next(A, i)
		end
		@test A.resnorm[end] < 1/n^2
		log(A,2)
	end
end

function get_problem_3d(n)
	A = MultigridMethod((n,m,k) -> p₂₃(n), (n, n, n), V(3, 2))
	b = fill(1,(n-1)*(n-1)*(n-1))
	return (A,b)
end

@testset "Direct discretization of elliptic PDE, 2D" begin
	for n in [4 8 16 32 64 128]
		(A, b) = get_problem_3d(n)
		A.grids[1].b .= b # copy rhs
		push!(A.resnorm,SimpleMultigrid.norm_of_residu(A.grids[1])) # log convergence history
		for i in 1:15
			next(A, i)
		end
		@test A.resnorm[end] < 1/n^3
		log(A,3)
	end
end
