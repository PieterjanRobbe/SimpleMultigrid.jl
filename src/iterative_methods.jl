# iterative_methods.jl : convenience wrappers for smoothers in the IterativeSolvers package

abstract type Smoother end

struct Jacobi <: Smoother end

struct GaussSeidel <: Smoother end

struct RedBlackGaussSeidel <: Smoother end

smooth!(x::AbstractVector, A::AbstractMatrix, b::AbstractVector, n::Int, ::Jacobi) = jacobi!(x,A,b,maxiter=n)

smooth!(x::AbstractVector, A::AbstractMatrix, b::AbstractVector, n::Int, ::GaussSeidel) = gauss_seidel!(x,A,b,maxiter=n)
