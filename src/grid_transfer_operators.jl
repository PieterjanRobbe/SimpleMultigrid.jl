# transfer_operators.jl : tranform quantities from fine to coarse grid and back

# TransferKind
abstract type TransferKind end

struct Injection <: TransferKind end

struct FullWeighting <: TransferKind end

struct HalfWeighting <: TransferKind end

struct Cubic <: TransferKind end

# TransferOperator
abstract type TransferOperator end

struct Restriction <: TransferOperator end

struct Interpolation <: TransferOperator end

# GridTransferOperator
struct GridTransferOperator{K<:TransferKind,O<:TransferOperator}
    kind::K
    op::O
end

# stencil definitions
stencil(op::GridTransferOperator{Injection,Restriction}) = @SVector [0,1.,0] 

stencil(op::GridTransferOperator{FullWeighting,Restriction}) = @SVector [1/4,1/2,1/4] 

stencil(op::GridTransferOperator{FullWeighting,Interpolation}) = @SVector [1/2,1,1/2]

stencil(op::GridTransferOperator{Cubic,Interpolation}) = @SVector [-1/16,0,9/16,1,9/16,0,-1/16]

# TODO: half weighting needs a 2d stencil...

# restriction
R₁(op::TransferKind,n::Int) = expand(GridTransferOperator(op,Restriction()),n) # 1d restriction operator
R(op::TransferKind,n::Int) = R₁(op,n)
R(op::TransferKind,n::Int...) = kron(R₁.(op,reverse(n))...)
restrict(u::Vector{T},op::TransferKind,n::Int...) where {T<:AbstractFloat} = R(op,n...)*u

# interpolation
P₁(op::TransferKind,n::Int) = expand(GridTransferOperator(op,Interpolation()),n) # 1d interpolation operator
P(op::TransferKind,n::Int) = P₁(op,2n)
P(op::TransferKind,n::Int...) = kron(P₁.(op,2.*reverse(n))...)
interpolate(u::Vector{T},op::TransferKind,n::Int...) where {T<:AbstractFloat} = P(op,n...)*u

# prolongation = interpolation
prolongate(u,op) = interpolate(u,op)

# expand
function expand(op::GridTransferOperator{K,O} where {K,O},n)
    Is = Int64[]; Js = Int64[]; Vs = Float64[]
    for (idx,i) in enumerate(2:2:n-1)
        (is,js,vs) = stencil2mat(stencil(op),n,i)
        push!(Is,idx*ones(length(is))...)
        push!(Js,js...)
        push!(Vs,vs...)
    end
    S = sparse_matrix_from_op(op,Is,Js,Vs)
    return apply_boundary!(op,S) 
end

sparse_matrix_from_op(op::GridTransferOperator{K,Interpolation} where {K},Is,Js,Vs) = sparse(Js,Is,Vs)
sparse_matrix_from_op(op::GridTransferOperator{K,Restriction} where {K},Is,Js,Vs) = sparse(Is,Js,Vs)

# correction for boundary in cubic interpolation operator
function apply_boundary!(op::GridTransferOperator{Cubic,Interpolation},S::SparseMatrixCSC)
    S[1,1] .+= 1/16 
    S[end,end] .+= 1/16 
    return S
end
apply_boundary!(op::GridTransferOperator,S::SparseMatrixCSC) = S