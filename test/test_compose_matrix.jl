using Base.Test
include(joinpath(Pkg.dir("SimpleMultigrid"),"test","compose_matrix.jl"))

# print util
function print_err(n,err)
    println("------------------------")
    println(" N     | RELATIVE ERROR ")
    println("------------------------")
    for i in 1:min(length(n),length(err))
        println(@sprintf(" %6i | %6.3e",n[i],err[i]))
    end
end

@testset "Test Laplace 1d              " begin
    val = Float64[]
    ns = [5,10,20,50,100,200,1000]
    for n in ns
        println(@sprintf("Solving 1d Laplacian with n = %i...",n))
        L = laplace1d(n)
        y = ones(n-1)
        x = L\y
        push!(val,maximum(x))
    end
    err = abs.(val[1:end-1]-val[end])/val[end]
    print_err(ns,err)
end

@testset "Test Laplace 2d              " begin
    val = Float64[]
    ns = [5,10,20,50,100,200,500]
    for n in ns
        println(@sprintf("Solving 2d Laplacian with n = %i...",n))
        L = laplace2d(n,n)
        y = ones((n-1)*(n-1))
        x = L\y
        push!(val,maximum(x))
    end
    err = abs.(val[1:end-1]-val[end])/val[end]
    print_err(ns,err)
end

f(x,y) = 1+sin.(π*x)*sin.(π*y)'
g(n,m)=f(linspace(0,1,n+1),linspace(0,1,m+1))

@testset "Test Elliptic 2d              " begin
    val = Float64[]
    ns = [5,10,20,50,100,200,500]
    for n in ns
        println(@sprintf("Solving 2d Elliptic with n = %i...",n))
        k = g(n,n)
        L = elliptic2d(k)
        y = ones((n-1)*(n-1))
        x = L\y
        push!(val,maximum(x))
    end
    err = abs.(val[1:end-1]-val[end])/val[end]
    print_err(ns,err)
end
