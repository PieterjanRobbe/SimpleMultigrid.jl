# test_convergence.jl

# example problem from A Multigrid Tutorial, scnd. ed.

f(x,y) = 2((1-6x^2)*y^2*(1-y^2)+(1-6y^2)*x^2*(1-x^2))

function get_problem(ns)
    n = ns
    m = n
    L = laplace2d(n,m)
    A = V_cycle(L, (n,m))
    pts = 1/n:1/n:1-1/n
    b = Float64[f(x,y) for x in pts, y in pts][:]
    return (A,b)
end

@testset "Convergence for example from Briggs et. al." begin
    for n in [16 32 64 128 256 512 1024]
        (A,b) = get_problem(n)
        A.grids[1].b .= b # copy rhs
        tol = 1/prod(A.grids[1].sz)*SimpleMultigrid.norm_of_residu(A.grids[1]) # target norm of residu is O(h²)
        push!(A.resnorm,SimpleMultigrid.norm_of_residu(A.grids[1])) # log convergence history
        for i in 1:15
            SimpleMultigrid.cycle!(A)
            push!(A.resnorm,SimpleMultigrid.norm_of_residu(A.grids[1])) # log convergence history
        end
        @test A.resnorm[end] < 1/n^2
        log(A)
    end
end
