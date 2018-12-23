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
The diffusion coefficient k must given on an equidistant grid with nodes {1,...,n-1}.
The step size is thus h = 1/n.
"""
function elliptic1d(k)
    n = length(k) + 1
    stencil2mat(SArray{Tuple{3}}([-n^2 2n^2 -n^2]), k)
end

"""
    flowcell(k)

Returns the discretization for an elliptic PDE with flow cell geometry.
The diffusion coefficient k must given on an equidistant grid with nodes {1,...,n-1} (1d)
or {1,...,n-1} x {0,..,m} (2d). The step sizes are thus h_x = 1/n and h_y = 1/m.

"""
function flowcell(k::AbstractVector)
    n = size(k) + 1
    A = elliptic1d(k)
    b = fill(0., size(A, 1))
    b[1] = k[1]*n^2

    A, b
end

function flowcell(k::AbstractMatrix)
    n = size(k, 1) + 1
    m = size(k, 2) + 2
    A = _elliptic2d(k, n, m)
    b = fill(0., size(A, 1))
    b[1:n-1] .= k[:, 1]*m^2

    R = CartesianIndices(size(k))
    L = LinearIndices(R)
    for i in 1:size(k, 2)
        A[L[1, i], L[2, i]] -= k[1, i]*n^2
        A[L[end, i], L[end-1, i]] -= k[end, i]*n^2
    end

    A, b
end

"""
    elliptic2d(k)

Generate system matrix for general 2d elliptic PDE defined on [0,1]^2.
The diffusion coefficient k must given on an equidistant grid with nodes {1,...,n-1} x {1,..,m-1}, 
boundaries included. The step sizes are thus h_x = 1/n and h_y = 1/m.
"""
elliptic2d(k) = _elliptic2d(k, size(k).+1...)

_elliptic2d(k, n, m) = stencil2mat(SArray{Tuple{3,3}}([0 -m^2 0 -n^2 2(n^2+m^2) -n^2 0 -m^2 0]), k)

"""
    elliptic3d(k)

Generate system matrix for general 3d elliptic PDE defined on [0,1]^3.
The diffsion coefficient k must given on an equidistant grid with nodes {0,1,...,n}x{0,1,..,m}x{0,1,...,l}, 
boundaries included. The step sizes are thus h_x = 1/n, h_y = 1/m and h_z = 1/l.
"""
function elliptic3d(k)
    n,m,l = size(k).+1
    stencil2mat(SArray{Tuple{3,3,3}}([0 0 0 0 -l^2 0 0 0 0 0 -m^2 0 -n^2 2(n^2+m^2+l^2) -n^2 0 -m^2 0 0 0 0 0 -l^2 0 0 0 0]), k)
end

# construction for constant stencil problems
function stencil2mat(stencil::SArray, n::Int...)
    R = CartesianIndices(n.-1)
    I1, Iend = extrema(R)
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
    R = CartesianIndices(UnitRange.(max(I1, I-I1).I, min(Iend, I+I1).I)) # TODO: #29440 replaces UnitRange with : (colon), see https://discourse.julialang.org/t/psa-replacement-of-ind2sub-sub2ind-in-julia-0-7/14666/9 [repeats 3 times]
    Is = fill(0,length(R))
    Js = fill(0,length(R))
    Vs = fill(0.,length(R))
    sz = LinearIndices(Iend.I)
    @inbounds for (i,J) in enumerate(R)
        idx = J-I+2I1
        Is[i] = sz[I]
        Js[i] = sz[J]
        Vs[i] = stencil[idx]
    end
    Is,Js,Vs
end

# main driver code for length 7 constant 1d stencils (used in `Cubic()` interpolation)
function _stencil2mat(stencil::SVector{7},I,I1,Iend)
    R = CartesianIndices(UnitRange.(max(I1, I-3I1).I, min(Iend, I+3I1).I))
    Is = fill(0,length(R))
    Js = fill(0,length(R))
    Vs = fill(0.,length(R))
    sz = LinearIndices(Iend.I)
    @inbounds for (i,J) in enumerate(R)
        idx = J-I+4I1
        Is[i] = sz[I]
        Js[i] = sz[J]
        Vs[i] = stencil[idx]
    end
    if I.I[1] == 2 # correction for boundary
        Vs[1] += 1/16
    elseif I.I[1] == Iend.I[1]-1
        Vs[end] += 1/16
    end
    Is,Js,Vs
end

# construction for non-constant stencil problems
function stencil2mat(stencil::SArray, k::AbstractArray)
    n = size(k)
    R = CartesianIndices(n)
    I1, Iend = extrema(R)
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
    R = CartesianIndices(UnitRange.(max(I1, I-I1).I, min(Iend, I+I1).I))
    #R = max(I1, I-I1):min(Iend, I+I1)
    Is = fill(0, length(R))
    Js = fill(0, length(R))
    Vs = fill(0., length(R))
    sz = LinearIndices(Iend.I)
    for (i,J) in enumerate(R)
        idx = J-I+2I1
        Is[i] = sz[I]
        Js[i] = sz[J]
        Vs[i] = k[J]*stencil[idx]
    end
    Is,Js,Vs
end
