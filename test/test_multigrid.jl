# test_multigrid.jl : test multigrid implementation

g(x) = Float64[1 + sin(π*xᵢ) for xᵢ in x]
g(x,y) = Float64[1 + sin(π*xᵢ)*sin(π*yᵢ) for xᵢ in x, yᵢ in y]
g(x,y,z) = Float64[1 + sin(π*xᵢ)*sin(π*yᵢ)*sin(π*zᵢ) for xᵢ in x, yᵢ in y, zᵢ in z]

f(n) = g(linspace(0,1,n))
f(n,m) = g(linspace(0,1,n),linspace(0,1,m))
f(n,m,l) = g(linspace(0,1,n),linspace(0,1,m),linspace(0,1,l))

p₁₁(n) = laplace1d(n)
p₁₂(n) = laplace2d(n,n)
p₁₃(n) = laplace3d(n,n,n)

p₂₁(n) = elliptic1d(f(n+1))
p₂₂(n) = elliptic2d(f(n+1,n+1))
p₂₃(n) = elliptic3d(f(n+1,n+1,n+1))

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
                for n in 2.^(2:6)
                    print("n = $(n)...")
                    A = method(eval(Symbol(problem,d_name))(n),ntuple(i->n,d))
                    b = ones(prod([n-1 for i in 1:d]))
                    x = A\b
                    @test A.resnorm[end] < 1/n^d
                    println("done")
                end
            end
        end
    end
end
