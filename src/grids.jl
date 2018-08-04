# grids.jl : introduces the Grid type for easy use of multigrid algorithms

# Ax=b
struct Grid{matT<:SparseMatrixCSC,vecT<:AbstractVector,opeT,tupT<:NTuple}
    A::matT # matrix A
    x::vecT # solution vector
    b::vecT # rhs vector
    R::opeT # restriction matrix
    P::opeT # prolongation matrix
    sz::tupT # PDE grid size
end
show(io::IO, grid::Grid) = print(io, string(join(grid.sz," x ")," grid"))

# coarsen the problem A into a sequence of coarser grids
"""
    coarsen(A, sz, R_op, P_op, ngrids)

Coarsen the matrix `A` using `R_op` as restriction operator and `P_op` as interpolation operator. The matrix `A` is the discrete version of a PDE defined on an `sz`-point mesh, and `ngrids` is the number of levels. 
"""
function coarsen(A::SparseMatrixCSC, sz::NTuple, R_op::TransferKind, P_op::TransferKind, ngrids::Int)
    grids = Vector{Grid{typeof(A),Vector{eltype(A)},typeof(A),typeof(sz)}}(ngrids)
    grids[1] = Grid(A,zero_x(A),zero_x(A),R(R_op,sz...),spzeros(0,0),sz)
    for i in 2:ngrids
        sz_c = grids[i-1].sz.>>1
        R_mat = R(R_op,sz_c...)
        P_mat = P(P_op,sz_c...)
        A_c = grids[i-1].R*grids[i-1].A*P_mat
        grids[i] = Grid(A_c,zero_x(A_c),zero_x(A_c),R_mat,P_mat,sz_c)
    end
    grids
end

zero_x(A::SparseMatrixCSC) = zeros(eltype(A),size(A,1)) 

δ(i,n) = Int[j==i for j=1:n]

# compute the residu on this level
residu(grid::Grid) = grid.b - grid.A*grid.x

# compute the h-norm of the residu on this level
h_norm(v::AbstractVector,sz) = 1./prod(sz) * sqrt(sum(v.*v))
norm_of_residu(grid::Grid) = h_norm(residu(grid),grid.sz)

# apply smoother
smooth!(grid::Grid,ν::Int,smoother::Smoother) = smooth!(grid.x,grid.A,grid.b,ν,smoother)
