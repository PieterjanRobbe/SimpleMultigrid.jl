## compose_matrix.jl : methods to compose example matrices

"""
    laplace1d(n)

1d Laplacian on [0,1] using n points
"""
function laplace1d(n)
    stencil = @SVector [-n^2,2n^2,-n^2]
    Is = Int64[]
    Js = Int64[]
    Vs = Float64[]
    for i = 1:n-1
        (is,js,vs) = stencil2mat(stencil,n,i)
        push!(Is,is...)
        push!(Js,js...)
        push!(Vs,vs...)
    end
    return sparse(Is,Js,Vs)
end

"""
    laplace2d(n,m)

2d Laplacian on [0,1]^2 using n x m points 
"""
function laplace2d(n,m)
    stencil = @SMatrix [0 -m^2 0; -n^2 2n^2+2m^2 -n^2; 0 -m^2 0]
    Is = Int64[]
    Js = Int64[]
    Vs = Float64[]
    for i in 1:n-1, j in 1:m-1
        (is,js,vs) = stencil2mat(stencil,n,m,i,j)
        push!(Is,is...)
        push!(Js,js...)
        push!(Vs,vs...)
    end
    return sparse(Is,Js,Vs)
end

"""
    elliptic1d(k)

Generate system matrix for general 1d elliptic PDE defined on [0,1].
The diffusion coefficient k must given on an equidistant grid with nodes {0,1,...,n}, 
boundaries included. The step size is thus h = 1/n.
"""
function elliptic1d(k)
    n = length(k)-1
    Is = Int64[]
    Js = Int64[]
    Vs = Float64[]
    for i = 1:n-1
        stencil = @SVector [-k[i]*n^2,2*k[i+1]*n^2,-k[i+2]*n^2]
        (is,js,vs) = stencil2mat(stencil,n,i)
        push!(Is,is...)
        push!(Js,js...)
        push!(Vs,vs...)
    end
    return sparse(Is,Js,Vs)
end

"""
    elliptic2d(k)

Generate system matrix for general 2d elliptic PDE defined on [0,1]^2.
The diffusion coefficient k must given on an equidistant grid with nodes {0,1,...,n}x{0,1,..,m}, 
boundaries included. The step sizes are thus h_x = 1/n and h_y = 1/m.
"""
function elliptic2d(k)
    (n,m) = size(k).-1
    Is = Int64[]
    Js = Int64[]
    Vs = Float64[]
    for i in 1:n-1, j in 1:m-1
        stencil = @SMatrix [0 -m^2*k[i+1,j+2] 0; -n^2*k[i,j+1] n^2*(k[i,j+1]+k[i+2,j+1])+m^2*(k[i+1,j+2]+k[i+1,j]) -n^2*k[i+2,j+1]; 0 -m^2*k[i+1,j] 0]
        (is,js,vs) = stencil2mat(stencil,n,m,i,j)
        push!(Is,is...)
        push!(Js,js...)
        push!(Vs,vs...)
    end
    return sparse(Is,Js,Vs)
end

# convert stencil to sparse matrix - 1d - length 3
function stencil2mat(stencil::SVector{3},n::Int,i::Int)
    els = [i-1,i,i+1]
    idcs = lexicographic_idcs_3(i,n)
    is = [i for k in 1:length(idcs)]
    js = els[idcs]
    vs = stencil[idcs]
    return (is,js,vs)
end

# convert stencil to sparse matrix - 1d - length 7
function stencil2mat(stencil::SVector{7},n::Int,i::Int)
    els = [i-3,i-2,i-1,i,i+1,i+2,i+3]
    idcs = lexicographic_idcs_7(i,n)
    is = [i for k in 1:length(idcs)]
    js = els[idcs]
    vs = stencil[idcs]
    return (is,js,vs)
end

# convert stencil to sparse matrix - 2d
function stencil2mat(stencil::SMatrix{3,3},n::Int,m::Int,i::Int,j::Int)
    el = sub2ind((n-1,m-1),i,j)
    els = [el+n-2 el+n-1 el+n; el-1 el el+1; el-n el-n+1 el-n+2]
    idcs = lexicographic_idcs_3(el,n,m)
    is = [el for k in 1:prod(length.(idcs))]
    js = els[idcs...]
    vs = stencil[idcs...]
    return (is,js,vs)
end

# helper function for lexicographic ordening of the indices
function lexicographic_idcs_3(el,n)
    if n == 2 # hack for problem with n = 2
        return 2:2
    elseif el == 1
        return 2:3
    elseif el == n-1
        return 1:2
    else
        return 1:3
    end
end

function lexicographic_idcs_7(el,n)
    if n == 4
        return 3:5
    elseif el == 1 || el == 2
        return 3:7
    elseif el == n-1 || el == n-2
        return 1:5
    else
        return 1:7
    end
end

function lexicographic_idcs_3(el,n,m)
    r = n == 2 ? 2 : 1 # catch matrices of size 2xN or Mx2 or 2x2
    s = m == 2 ? 2 : 1
    t = n == 2 ? 2 : 3
    if el == 1 # SW corner
        return (s:2,2:t)
    elseif el == n-1 # SE corner
        return (s:2,1:2)
    elseif el == (n-1)*(m-2)+1 # NW corner
        return (2:3,2:t)
    elseif el == (n-1)*(m-1) # NE corner
        return (2:3,1:2)
    elseif el < n-1 # S border
        return (s:2,1:3)
    elseif mod(el,n-1) == 1 # W border
        return (1:3,2:3)
    elseif mod(el,n-1) == 0 # E border
        return (1:3,r:2)
    elseif el > (n-1)*(m-2)+1 # N border
        return (2:3,1:3)
    else # mid
        return (1:3,1:3)
    end
end
