# transfer_operators.jl : tranform quantities from fine to coarse grid and back

import Base.expand

# TransferKind
abstract type TransferKind end

struct Injection <: TransferKind end

struct FullWeighting <: TransferKind end

struct HalfWeighting <: TransferKind end

# TransferOperator
abstract type TransferOperator end

struct Restriction <: TransferOperator end

struct Interpolation <: TransferOperator end

# GridTransferOperator
struct GridTransferOperator{d,K<:TransferKind,O<:TransferOperator}
    kind::K
    op::O
end

GridTransferOperator(n::N where {N<:Integer},kind::K,op::O) where {K<:TransferKind,O<:TransferOperator} = GridTransferOperator{n,K,O}(kind,op)

# stencil definitions
stencil(op::GridTransferOperator{1,Injection,Restriction}) = [0,1,0] 

stencil(op::GridTransferOperator{1,FullWeighting,Restriction}) = 1/4*[1,2,1] 

stencil(op::GridTransferOperator{1,FullWeighting,Interpolation}) = 1/2*[1,2,1] 

stencil(op::GridTransferOperator{2,Injection,Restriction}) = [0 0 0; 0 1 0; 0 0 0] 

stencil(op::GridTransferOperator{2,FullWeighting,Restriction}) = 1/16*[1 2 1; 2 4 2; 1 2 1] 

stencil(op::GridTransferOperator{2,FullWeighting,Interpolation}) = 1/4*[1 2 1; 2 4 2; 1 2 1] 

stencil(op::GridTransferOperator{2,HalfWeighting,Restriction}) = 1/8*[0 1 0; 1 4 1; 0 1 0] 

# restriction
restrict(u::Vector{T},op::TransferKind,n::N...) where {T<:AbstractFloat,N<:Integer} = expand(GridTransferOperator(length(n),op,Restriction()),n...)*u

# interpolation
interpolate(u::Vector{T},op::TransferKind,n::N...) where {T<:AbstractFloat,N<:Integer} = expand(GridTransferOperator(length(n),op,Interpolation()),2.*n...)*u

# prolongation = interpolation
prolongate(u,op) = interpolate(u,op)

# expand
function expand(op::GridTransferOperator{1,K,O} where {K,O},n)
    Is = Int64[]; Js = Int64[]; Vs = Float64[]
    for (idx,i) in enumerate(2:2:n-1)
        (is,js,vs) = stencil2mat(stencil(op),n,i)
        push!(Is,idx*ones(length(is))...)
        push!(Js,js...)
        push!(Vs,vs...)
    end
    return sparse_matrix_from_op(op,Is,Js,Vs) 
end

sparse_matrix_from_op(op::GridTransferOperator{d,K,Interpolation} where {d,K},Is,Js,Vs) = sparse(Js,Is,Vs)
sparse_matrix_from_op(op::GridTransferOperator{d,K,Restriction} where {d,K},Is,Js,Vs) = sparse(Is,Js,Vs)

function expand(op::GridTransferOperator{2,K,O} where {K,O},n,m)
    Is = Int64[]; Js = Int64[]; Vs = Float64[]
    for (idx,(j,i)) in enumerate(Base.Iterators.product(2:2:m-1,2:2:n-1))
        (is,js,vs) = stencil2mat(stencil(op),n,m,i,j)
        push!(Is,idx*ones(length(is))...)
        push!(Js,js...)
        push!(Vs,vs...)
    end
    return sparse_matrix_from_op(op,Is,Js,Vs) 
end
