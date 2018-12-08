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

ndims(grid::Grid) = length(size(grid))

size(grid::Grid) = grid.sz
size(grid::Grid, n::Integer) = grid.sz[n]

show(io::IO, grid::Grid) = print(io, string(join(grid.sz," x ")," grid"))

# coarsen the problem A into a sequence of coarser grids
"""
    coarsen(A, sz, R_op, P_op, ngrids)

Coarsen the matrix `A` using `R_op` as restriction operator and `P_op` as interpolation operator. The matrix `A` is the discrete version of a PDE defined on an `sz`-point mesh, and `ngrids` is the number of levels.  A Galerkin approach is used to compose the coarse matrices, unless a function `f` is provided for direct discretization.
"""
function coarsen(A::SparseMatrixCSC, sz::NTuple, R_op::TransferKind, P_op::TransferKind, ngrids::Int)
    grids = Vector{Grid{typeof(A),Vector{eltype(A)},typeof(A),typeof(sz)}}(undef, ngrids)
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

function coarsen(f::Function, sz::NTuple, R_op::TransferKind, P_op::TransferKind, ngrids::Int)
	A0 = f(sz...)
    grids = Vector{Grid{typeof(A0),Vector{eltype(A0)},typeof(A0),typeof(sz)}}(undef, ngrids)
    grids[1] = Grid(A0,zero_x(A0),zero_x(A0),R(R_op,sz...),spzeros(0,0),sz)
    for i in 2:ngrids
        sz_c = grids[i-1].sz.>>1
        R_mat = R(R_op,sz_c...)
        P_mat = P(P_op,sz_c...)
		A_c = f(sz_c...)
        grids[i] = Grid(A_c,zero_x(A_c),zero_x(A_c),R_mat,P_mat,sz_c)
    end
    grids
end

# automatically determine number of grids
# TODO: because of the way the grid transfer operators are implemented
# we can only do prolongation for grids with an even number of points...
function factor_twos(n::Int)
    i = 0
    while iseven(n) && n > 1
        n >>= 1
        i+=1
    end
    i
end

zero_x(A::SparseMatrixCSC) = zeros(eltype(A),size(A,1)) 

# compute the residu on this level
residu(grid::Grid) = grid.b - grid.A*grid.x

# compute the h-norm of the residu on this level
h_norm(v::AbstractVector,sz) = 1.0/prod(sz) * sqrt(sum(v.*v))
norm_of_residu(grid::Grid) = h_norm(residu(grid),grid.sz)

# apply smoother
smooth!(grid::Grid, ν::Int, smoother::Smoother) = smooth!(grid.x, grid.A, grid.b, ν, smoother)

# apply line smoother
function smooth!(grid::Grid, ν::Int, smoother::LineSmoother)
    R = CartesianIndices(size(grid))
    L = LinearIndices(R)
    for dir in 1:ndims(grid)
        for νᵢ in 1:ν
            for i in 1:size(grid, dir)
                idcs = selectdim(R, dir, i)
                smooth!(view(grid.x, L[idcs]), view(grid.A, L[idcs], L[idcs]), view(grid.b, L[idcs]), 1, smoother.smoother)
            end
        end
    end
end


