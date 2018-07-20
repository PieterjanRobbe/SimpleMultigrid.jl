using PyPlot

include(joinpath(Pkg.dir("SimpleMultigrid"),"test","compose_matrix.jl"))
include(joinpath(Pkg.dir("SimpleMultigrid"),"test","grid_transfer_operators.jl"))

function mplot(x...)
    n = maximum(length.(x))+1
    figure(figsize=(5,5))
    for i in 1:length(x)
        factor = round(Int64,n/(length(x[i])+1))
        plot(factor:factor:n-1,x[i],"-o")
    end
end

function msurf(s...)
    n = maximum(size.(s,1))+1
    m = maximum(size.(s,2))+1
    figure(figsize=(5,5))
    alphas = linspace(0.2,1,length(s))
    for i in 1:length(s)
        nfactor = round(Int64,n/(size(s[i],1)+1))
        mfactor = round(Int64,m/(size(s[i],2)+1))
        x = nfactor:nfactor:n-1
        y = mfactor:mfactor:m-1
        xgrid = repmat(x,1,length(y))
        ygrid = repmat(y',length(x),1)
        surf(xgrid,ygrid,s[i],alpha=alphas[i])
    end
end

##################
#      1D        #
##################
n = 16
k(n) = 1+sin.(π*linspace(0,1,n+1))
L = elliptic1d(k(n))
b = ones(n-1)
x = L\b

rx1 = restrict(x,Injection(),n)
rx2 = restrict(rx1,Injection(),n>>1)
rx3 = restrict(rx2,Injection(),n>>2)
mplot(x,rx1,rx2,rx3)

rx1 = restrict(x,FullWeighting(),n)
rx2 = restrict(rx1,FullWeighting(),n>>1)
rx3 = restrict(rx2,FullWeighting(),n>>2)
mplot(x,rx1,rx2,rx3)

ix1 = interpolate(x,FullWeighting(),n)
ix2 = interpolate(ix1,FullWeighting(),2n)
mplot(ix2,ix1,x)

##################
#      2D        #
##################
f(x,y) = 1+sin.(π*x)*sin.(π*y)'
g(n,m) = f(linspace(0,1,n+1),linspace(0,1,m+1))
n = 16
m = 16
L = elliptic2d(g(n,m))
b = ones((n-1)*(m-1))
x = L\b

rx1 = restrict(x,Injection(),n,m)
rx2 = restrict(rx1,Injection(),n>>1,m>>1)
rx3 = restrict(rx2,Injection(),n>>2,m>>2)
msurf(reshape(x,(n-1,m-1)),reshape(rx1,(n>>1-1,m>>1-1)),reshape(rx2,(n>>2-1,m>>2-1)),reshape(rx3,(n>>3-1,m>>3-1)))

rx1 = restrict(x,FullWeighting(),n,m)
rx2 = restrict(rx1,FullWeighting(),n>>1,m>>1)
rx3 = restrict(rx2,FullWeighting(),n>>2,m>>2)
msurf(reshape(x,(n-1,m-1)),reshape(rx1,(n>>1-1,m>>1-1)),reshape(rx2,(n>>2-1,m>>2-1)),reshape(rx3,(n>>3-1,m>>3-1)))

rx1 = restrict(x,HalfWeighting(),n,m)
rx2 = restrict(rx1,HalfWeighting(),n>>1,m>>1)
rx3 = restrict(rx2,HalfWeighting(),n>>2,m>>2)
msurf(reshape(x,(n-1,m-1)),reshape(rx1,(n>>1-1,m>>1-1)),reshape(rx2,(n>>2-1,m>>2-1)),reshape(rx3,(n>>3-1,m>>3-1)))

ix1 = interpolate(x,FullWeighting(),n,m)
ix2 = interpolate(ix1,FullWeighting(),2n,2m)
msurf(reshape(ix2,(4n-1,4m-1)),reshape(ix1,(2n-1,2m-1)),reshape(x,(n-1,m-1)))
