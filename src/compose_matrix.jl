## compose_matrix.jl : methods to compose example matrices

"""
    laplace1d(n)

1d Laplacian on [0,1] using n points
"""
laplace1d(n) = stencil2mat(SArray{Tuple{3}}([-n^2 2n^2 -n^2]), n)

"""
    laplace2d(n, m)

2d Laplacian on [0,1]^2 using n x m points 
"""
laplace2d(n, m) = stencil2mat(SArray{Tuple{3,3}}([0 -m^2 0 -n^2 2(n^2+m^2) -n^2 0 -m^2 0]), n, m)

"""
    laplace3d(n, m, l)

3d Laplacian on [0,1]^3 using n x m x l points 
"""
laplace3d(n, m, l) = stencil2mat(SArray{Tuple{3,3,3}}([0 0 0 0 -l^2 0 0 0 0 0 -m^2 0 -n^2 2(m^2+n^2+l^2) -n^2 0 -m^2 0 0 0 0 0 -l^2 0 0 0 0]), n, m, l)
    
"""
    elliptic1d(k)

Generate system matrix for general 1d elliptic PDE defined on [0,1].
The diffusion coefficient k must given on an equidistant grid with nodes {0,1,...,n}, 
boundaries included. The step size is thus h = 1/n.
"""
function elliptic1d(k)
    n = length(k)-1
    stencil2mat(SArray{Tuple{3}}([-n^2 2n^2 -n^2]), k)
end

"""
    elliptic2d(k)

Generate system matrix for general 2d elliptic PDE defined on [0,1]^2.
The diffusion coefficient k must given on an equidistant grid with nodes {0,1,...,n}x{0,1,..,m}, 
boundaries included. The step sizes are thus h_x = 1/n and h_y = 1/m.
"""
function elliptic2d(k)
    n,m = size(k).-1
    stencil2mat(SArray{Tuple{3,3}}([0 -m^2 0 -n^2 2(n^2+m^2) -n^2 0 -m^2 0]), k)
end

"""
    elliptic3d(k)

Generate system matrix for general 3d elliptic PDE defined on [0,1]^3.
The diffsion coefficient k must given on an equidistant grid with nodes {0,1,...,n}x{0,1,..,m}x{0,1,...,l}, 
boundaries included. The step sizes are thus h_x = 1/n, h_y = 1/m and h_z = 1/l.
"""
function elliptic3d(k)
    n,m,l = size(k).-1
    stencil2mat(SArray{Tuple{3,3,3}}([0 0 0 0 -l^2 0 0 0 0 0 -m^2 0 -n^2 2(n^2+m^2+l^2) -n^2 0 -m^2 0 0 0 0 0 -l^2 0 0 0 0]), k)
end

# construction for constant stencil problems
function stencil2mat(stencil::SArray,n::Int...)
    R = CartesianRange(n.-1)
    I1, Iend = first(R), last(R)
    Is = Int64[]
    Js = Int64[]
    Vs = Float64[]
    @inbounds for I in R
        is,js,vs = _stencil2mat(stencil,I,I1,Iend)
        push!(Is,is...)
        push!(Js,js...)
        push!(Vs,vs...)
    end
    sparse(Is,Js,Vs)
end

# main driver code - constant stencil
@noinline function _stencil2mat(stencil::SArray,I,I1,Iend)
    R = CartesianRange(max(I1, I-I1), min(Iend, I+I1))
    Is = fill(0,length(R))
    Js = fill(0,length(R))
    Vs = fill(0.,length(R))
    @inbounds for (i,J) in enumerate(R)
        idx = J-I+2I1
        Is[i] = sub2ind(Iend.I,I.I...)
        Js[i] = sub2ind(Iend.I,J.I...)
        Vs[i] = stencil[idx]
    end
    Is,Js,Vs
end

# main driver code for length 7 constant 1d stencils (used in `Cubic()` interpolation)
function _stencil2mat(stencil::SVector{7},I,I1,Iend)
    R = CartesianRange(max(I1, I-3I1), min(Iend, I+3I1))
    Is = fill(0,length(R))
    Js = fill(0,length(R))
    Vs = fill(0.,length(R))
    for (i,J) in enumerate(R)
        idx = J-I+4I1
        Is[i] = sub2ind(Iend.I,I.I...)
        Js[i] = sub2ind(Iend.I,J.I...)
        Vs[i] = stencil[idx]
    end
    if I.I[1] == 2 # correction for boundary
        Vs[1] .+= 1/16
    elseif I.I[1] == Iend.I[1]-1
        Vs[end] .+= 1/16
    end
    Is,Js,Vs
end

# construction for non-constant stencil problems
function stencil2mat(stencil::SArray,k::AbstractArray)
    n = size(k).-2
    R = CartesianRange(n)
    I1, Iend = first(R), last(R)
    Is = Int64[]
    Js = Int64[]
    Vs = Float64[]
    for I in R
        is,js,vs = _stencil2mat(k,stencil,I,I1,Iend)
        push!(Is,is...)
        push!(Js,js...)
        push!(Vs,vs...)
    end
    sparse(Is,Js,Vs)
end

# main driver code - non-constant stencil
function _stencil2mat(k::AbstractArray,stencil::SArray,I,I1,Iend)
    R = CartesianRange(max(I1, I-I1), min(Iend, I+I1))
    Is = fill(0,length(R))
    Js = fill(0,length(R))
    Vs = fill(0.,length(R))
    for (i,J) in enumerate(R)
        idx = J-I+2I1
        Is[i] = sub2ind(Iend.I,I.I...)
        Js[i] = sub2ind(Iend.I,J.I...)
        Vs[i] = k[J+I1]*stencil[idx]
    end
    Is,Js,Vs
end
