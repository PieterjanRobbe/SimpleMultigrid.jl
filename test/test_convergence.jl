# test_convergence.jl

function get_problem(ns)
    n = m = ns
    L = laplace2d(n, m)
    A = V_cycle(L, (n, m))
    pts = 1/n:1/n:1-1/n
    b = Float64[h(x,y) for x in pts, y in pts][:]
    return (A,b)
end

@testset "Convergence for example from Briggs et. al." begin
    for n in [16 32 64 128 256 512 1024]
        (A, b) = get_problem(n)
        A.grids[1].b .= b # copy rhs
        for item in Base.Iterators.take(A,15) end # iterate
        @test A.resnorm[end] < 1/n^2
        log(A,2)
    end
end
