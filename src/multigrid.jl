# multigrid.jl : define multigrid cycle operations

# Multigrid cycles
abstract type MultigridCycle end

struct V <: MultigridCycle 
    ν₁::Int
    ν₂::Int
end
V() = V(2,1)

struct W <: MultigridCycle
    ν₁::Int
    ν₂::Int
end
W() = W(2,1)

struct F <: MultigridCycle
    ν₀::Int
    ν₁::Int
    ν₂::Int
end
F() = F(2,2,1)

# Multigrid struct
struct MultigridIterable{C<:MultigridCycle,Gs,S,V}
    grids::Gs
    max_iter::Int
    cycle_type::C
    smoother::S
    resnorm::V
end
show(io::IO, mg::MultigridIterable) = print(io, string("$(length(mg.grids))-level Multigrid method"))

"""
    MultigridMethod(A, sz, cycle_type)
    MultigridMethod(A, sz, cycle_type; kwargs...)

Geometric Multigrid method of type `cycle_type` for matrix `A`, that results from a discretization of a PDE on [0,1]^d using an `sz`-point mesh.

Inputs
======
A          : SparseMatrixCSC, sparse matrix with discretized PDE
sz         : NTuple, PDE grid size, e.g., `(n,m)`
cycle_type : Multigrid cycle type, can be `V(ν₁,ν₂)` (V-cycle), `W(ν₁,ν₂)` (W-cycle) or `F(ν₀,ν₁,ν₂)` (FMG)

Options
=======
* max_iter : maximum number of iterations. Default value is `20` For FMG using the `F()`-cycle, this is set to 1.
* R_op     : restriction operator type, can be `Injection()` or `FullWeighting()` (default) 
* P_op     : interpolation operator type, can be `Injection()`, `FullWeighting()` (default), or `Cubic()` 
* ngrids   : total number of grids to use, default is `min.(⌊log₂(sz)⌋)`
* smoother : smoother, can be `GaussSeidel()` of `Jacobi()`
"""
MultigridMethod(A::AbstractMatrix, sz::NTuple, cycle_type::MultigridCycle; max_iter::Int=20, R_op::TransferKind=FullWeighting(), P_op::TransferKind=FullWeighting(), ngrids::Int=minimum(floor.(Int,log2.(sz))), smoother::Smoother=GaussSeidel()) = MultigridIterable(coarsen(A,sz,R_op,P_op,ngrids),max_iter,cycle_type,smoother,Float64[])

"""
    V_cycle(A, sz)
    V_cycle(A, sz; kwargs...)

Geometric Multigrid V(2,1)-cycle for matrix `A`, that results from a discretization of a PDE on [0,1]^d using an `sz`-point mesh.

Inputs
======
A          : SparseMatrixCSC, sparse matrix with discretized PDE
sz         : NTuple, PDE grid size, e.g., `(n,m)`

For other options, see `MultigridMethod`.
"""
V_cycle(A::AbstractMatrix, sz::NTuple; kwargs...) = MultigridMethod(A,sz,V(); kwargs...)

"""
    W_cycle(A, sz)
    W_cycle(A, sz; kwargs...)

Geometric Multigrid W(2,1)-cycle for matrix `A`, that results from a discretization of a PDE on [0,1]^d using an `sz`-point mesh.

Inputs
======
A          : SparseMatrixCSC, sparse matrix with discretized PDE
sz         : NTuple, PDE grid size, e.g., `(n,m)`

For other options, see `MultigridMethod`.
"""
W_cycle(A::AbstractMatrix, sz::NTuple; kwargs...) = MultigridMethod(A,sz,W(); kwargs...)

"""
F_cycle(A, sz)
F_cycle(A, sz; kwargs...)

Geometric Full Multigrid for matrix `A`, that results from a discretization of a PDE on [0,1]^d using an `sz`-point mesh.

Inputs
======
A          : SparseMatrixCSC, sparse matrix with discretized PDE
sz         : NTuple, PDE grid size, e.g., `(n,m)`

For other options, see `MultigridMethod`.
"""
F_cycle(A::AbstractMatrix, sz::NTuple; kwargs...) = MultigridMethod(A,sz,F(); max_iter=1, kwargs...)

# solver
"""
    \\(M,b)

Solve the matrix problem `Mx=b` using Multigrid.

Inputs
======
M : MultigridIterable
b : right-hand side
"""
function \(mg::MultigridIterable, b::AbstractVector)
    length(b) == length(mg.grids[1].b) || throw(DimensionMismatch("Right-hand side b has length $(length(b)), but needs $(length(mg.grids[1].b))")) # check dimensions
    mg.grids[1].b .= b # copy rhs
    push!(mg.resnorm,norm_of_residu(mg.grids[1])) # log convergence history

    for item in Base.Iterators.take(mg,mg.max_iter) end # iterate

    mg.resnorm[end] >= 1/prod(mg.grids[1].sz) && warn(2sprintf("maximum number of iterations reached, norm of residual is %6.3e > %6.3e",mg.resnorm[end],1/prod(mg.grids[1].sz))) # check convergence with max_iter igterations

    return mg.grids[1].x
end

# iterator commands
start(mg::MultigridIterable) = 1

done(mg::MultigridIterable, it::Int) = mg.resnorm[end] < 1/prod(mg.grids[1].sz) # target norm of residu is O(h²)

function next(mg::MultigridIterable, it::Int)
    cycle!(mg)
    push!(mg.resnorm,norm_of_residu(mg.grids[1])) # log convergence history

    nothing, it+1
end

# multigrid cycles
cycle!(mg::MultigridIterable{C} where {C<:V}) = μ_cycle!(mg.grids,1,mg.cycle_type.ν₁,mg.cycle_type.ν₂,1,mg.smoother)
cycle!(mg::MultigridIterable{C} where {C<:W}) = μ_cycle!(mg.grids,2,mg.cycle_type.ν₁,mg.cycle_type.ν₂,1,mg.smoother)
cycle!(mg::MultigridIterable{C} where {C<:F}) = F_cycle!(mg.grids,mg.cycle_type.ν₀,mg.cycle_type.ν₁,mg.cycle_type.ν₂,1,mg.smoother)

function μ_cycle!(grids::Vector{G} where {G<:Grid}, μ::Int, ν₁::Int, ν₂::Int, grid_ptr::Int, smoother::Smoother)
    smooth!(grids[grid_ptr],ν₁,smoother)
    if grid_ptr == length(grids)
        grids[grid_ptr].x .= grids[grid_ptr].A\grids[grid_ptr].b # exact solve
    else
        grids[grid_ptr+1].b .= grids[grid_ptr].R*residu(grids[grid_ptr])
        grids[grid_ptr+1].x .= zeros(grids[grid_ptr+1].x)
        μ_cycle!(grids,μ,ν₁,ν₂,grid_ptr+1,smoother)
        grids[grid_ptr].x .+= grids[grid_ptr+1].P*grids[grid_ptr+1].x
    end
    smooth!(grids[grid_ptr],ν₂,smoother)
end

function F_cycle!(grids::Vector{G} where {G<:Grid}, ν₀::Int, ν₁::Int, ν₂::Int, grid_ptr::Int, smoother::Smoother)
    if grid_ptr == length(grids)
        grids[grid_ptr].x .= grids[grid_ptr].A\grids[grid_ptr].b # exact solve
    else
        grids[grid_ptr+1].b .= grids[grid_ptr].R*grids[grid_ptr].b
        F_cycle!(grids,ν₀,ν₁,ν₂,grid_ptr+1,smoother)
        grids[grid_ptr].x .= P(Cubic(),grids[grid_ptr+1].sz...)*grids[grid_ptr+1].x # FMG with cubic interpolation
    end
    for i in 1:ν₀
        μ_cycle!(grids,1,ν₁,ν₂,grid_ptr,smoother)
    end
end
