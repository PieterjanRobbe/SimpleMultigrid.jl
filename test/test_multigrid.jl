# test_multigrid.jl : test multigrid implementation

g(x,y) = 1 + sin.(π*x)*sin.(π*y)'
f(n,m) = g(linspace(0,1,n+1),linspace(0,1,m+1))
p₁(n,m) = laplace2d(n,m)
p₂(n,m) = elliptic2d(f(n,m))

@testset "Default Multigrid methods for 2D laplace and elliptic" begin
    for (problem,problem_name) in zip([p₁,p₂],["LAPLACE 2D ","ELLIPTIC 2D"])
        println("********************")
        println("*   $(problem_name)    *")
        println("********************")
        for (method,method_name) in zip([V_cycle, W_cycle, F_cycle],["V-CYCLE", "W-CYCLE", "  FMG  "])
            println("+------------------+")
            println("|      $(method_name)     |")
            println("+------------------+")
            for n in 2.^(2:10)
                print("n = $(n)...")
                A = method(problem(n,n),(n,n))
                b = ones((n-1)*(n-1))
                x = A\b
                @test A.resnorm[end] < 1/n^2
                println("done")
            end
        end
    end
end
