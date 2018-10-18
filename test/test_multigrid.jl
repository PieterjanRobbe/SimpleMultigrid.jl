# test_multigrid.jl : test multigrid implementation

for (d,d_name) in zip(1:3,["₁","₂","₃"])
    @testset "Default Multigrid methods for $(d)D laplace and elliptic problem" begin
        for (problem,problem_name) in zip(["p₁","p₂"],["LAPLACE $(d)D ","ELLIPTIC $(d)D"])
            println("********************")
            println("*   $(problem_name)    *")
            println("********************")
            for (method,method_name) in zip([V_cycle, W_cycle, F_cycle],["V-CYCLE", "W-CYCLE", "  FMG  "])
                println("+------------------+")
                println("|      $(method_name)     |")
                println("+------------------+")
                for n in 2 .^(1:(14-3d))
                    print("n = $(n)...")
                    A = method(eval(Symbol(problem,d_name))(n),ntuple(i->n,d))
                    b = ones(prod([n-1 for i in 1:d]))
                    x = A\b
                    @test A.resnorm[end] < 1/n^d
                    println("done ($(length(A.resnorm)-1) iterations)")
                end
            end
        end
    end
end
