# iterative_methods.jl : convenience wrappers for smoothers in the IterativeSolvers package

abstract type Smoother end

struct Jacobi <: Smoother end

struct GaussSeidel <: Smoother end

struct LineSmoother{S<:Smoother} <: Smoother
    smoother::S
end

LineSmoother(smoother::Smoother) = smoother isa LineSmoother ? throw(ArgumentError("cannot create a LineSmoother for another LineSmoother!")) : LineSmoother(smoother)

smooth!(x::AbstractVector, A::AbstractSparseMatrix, b::AbstractVector, n::Int, ::Jacobi) = jacobi!(x,A,b,maxiter=n)

smooth!(x::AbstractVector, A::AbstractSparseMatrix, b::AbstractVector, n::Int, ::GaussSeidel) = gauss_seidel!(x,A,b,maxiter=n)
