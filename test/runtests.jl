using SimpleMultigrid, Base.Test, Requires

using PyPlot

n = 10
x = linspace(0,1,n)
y = sin.(π*x)

x₂ = linspace(0,1,2n+1)
@show P = SimpleMultigrid.P(Cubic(),n+1)
@show full(P)
@show size(P)
y₂ = P*y
@show size(y)
@show size(y₂)

plot(x,y)
plot(x₂,y₂)



#=
include("test_utils.jl")

#include("test_convergence.jl")
include("test_multigrid.jl")
=#
